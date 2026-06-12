"""
WGBS 分析流程编排
================
SampleManifest → Pipeline → 完整分析流程。
"""

import os
import sys
import yaml
from dataclasses import dataclass, field
from typing import Optional, List

from .config import AnalysisConfig, DEFAULT_REF_DIR, DEFAULT_OUTPUT_DIR, GENOME_FNA
from .modules.genome import build_bismark_index, compute_and_cache_c_stats
from .modules.trimming import run_fastp
from .modules.mapping import run_bismark_alignment
from .modules.dedup import deduplicate
from .modules.extractor import run_methylation_extractor
from .modules.caller import run_methylation_calling
from .modules.summary import generate_summary


@dataclass
class SampleManifest:
    """单个样本的描述。"""
    name: str
    read1: str
    read2: Optional[str] = None       # None = 单端
    skip_trimming: bool = False       # _clean 文件跳过修剪
    metadata: dict = field(default_factory=dict)

    @property
    def paired_end(self) -> bool:
        return bool(self.read2)


def load_samples(yaml_path: str) -> List[SampleManifest]:
    """
    从 YAML 文件加载样本清单。

    YAML 格式:
        samples:
          - name: NC2M_1
            read1: /path/to/NC2M_1_clean_1.fq.gz
            read2: /path/to/NC2M_1_clean_2.fq.gz
            skip_trimming: true
          - name: SRR1566116
            read1: /path/to/SRR1566116_1.fastq.gz
            # read2 省略 = 单端
    """
    with open(yaml_path) as fh:
        data = yaml.safe_load(fh)

    samples = []
    for entry in data.get("samples", []):
        # 自动检测 _clean 文件
        skip = entry.get("skip_trimming", False)
        if not skip and "_clean" in entry["read1"]:
            skip = True

        samples.append(SampleManifest(
            name=entry["name"],
            read1=entry["read1"],
            read2=entry.get("read2"),
            skip_trimming=skip,
            metadata=entry.get("metadata", {}),
        ))

    return samples


class Pipeline:
    """WGBS 分析管线。"""

    def __init__(
        self,
        ref_dir: str = DEFAULT_REF_DIR,
        output_dir: str = DEFAULT_OUTPUT_DIR,
        config: Optional[AnalysisConfig] = None,
        dry_run: bool = False,
    ):
        self.ref_dir = ref_dir
        self.output_dir = output_dir
        self.config = config or AnalysisConfig()
        self.dry_run = dry_run

        # 预计算
        self.fasta_path = os.path.join(ref_dir, GENOME_FNA)

        if not os.path.exists(self.fasta_path):
            raise FileNotFoundError(
                f"参考基因组不存在: {self.fasta_path}\n"
                f"请确认 --ref_dir 参数指向包含 {GENOME_FNA} 的目录"
            )

    def setup_genome(self) -> dict:
        """一次性: 构建 Bismark 索引 + 计算 C 统计。"""
        print("=" * 50)
        print("[pipeline] 基因组准备...")

        # 构建 Bismark 索引
        build_bismark_index(self.ref_dir, dry_run=self.dry_run)

        # 计算 C 统计
        c_stats = compute_and_cache_c_stats(self.ref_dir, GENOME_FNA)
        print(f"[pipeline] 基因组 C 统计: "
              f"CG={c_stats['totals']['C_CG']:,} "
              f"CHG={c_stats['totals']['C_CHG']:,} "
              f"CHH={c_stats['totals']['C_CHH']:,} "
              f"Total={c_stats['totals']['total_C']:,}")
        return c_stats

    def run_sample(self, sample: SampleManifest, c_stats: dict) -> dict:
        """
        运行单个样本的完整分析流程。

        返回:
            {"mc_level_file": path, "summary_file": path}
        """
        sample_dir = os.path.join(self.output_dir, sample.name)
        os.makedirs(sample_dir, exist_ok=True)

        print(f"\n{'=' * 50}")
        print(f"[pipeline] 处理样本: {sample.name} "
              f"({'PE' if sample.paired_end else 'SE'})")
        print(f"{'=' * 50}")

        cfg = self.config

        # --- Step 1: Trimming ---
        if sample.skip_trimming:
            print(f"[pipeline:{sample.name}] 跳过修剪 (_clean 文件)")
            trimmed_r1 = sample.read1
            trimmed_r2 = sample.read2 or ""
        else:
            trim_dir = os.path.join(sample_dir, "trimmed")
            print(f"[pipeline:{sample.name}] fastp 修剪...")
            trimmed_r1, trimmed_r2 = run_fastp(
                read1=sample.read1,
                read2=sample.read2 or "",
                output_dir=trim_dir,
                sample_name=sample.name,
                threads=min(cfg.total_cores, 32),
                fastp_params=cfg.fastp_params,
                dry_run=self.dry_run,
            )

        # --- Step 2: Bismark Alignment ---
        bam_dir = os.path.join(sample_dir, "bam")
        print(f"[pipeline:{sample.name}] Bismark 比对...")
        bam_file = run_bismark_alignment(
            read1=trimmed_r1,
            read2=trimmed_r2 if sample.paired_end else "",
            ref_dir=self.ref_dir,
            output_dir=bam_dir,
            sample_name=sample.name,
            parallel_instances=cfg.bismark_parallel_instances,
            bowtie2_threads=cfg.bismark_bowtie2_threads,
            score_min=cfg.bowtie2_score_min,
            max_insert=cfg.max_insert,
            dovetail=cfg.dovetail,
            dry_run=self.dry_run,
        )

        # --- Step 3: Deduplication ---
        dedup_dir = os.path.join(sample_dir, "bam")
        print(f"[pipeline:{sample.name}] 去重...")
        dedup_bam = deduplicate(
            bam_file=bam_file,
            output_dir=dedup_dir,
            sample_name=sample.name,
            paired_end=sample.paired_end,
            dry_run=self.dry_run,
        )

        # --- Step 4: Methylation Extraction ---
        extract_dir = os.path.join(sample_dir, "meth_extract")
        print(f"[pipeline:{sample.name}] 甲基化提取...")
        run_methylation_extractor(
            dedup_bam=dedup_bam,
            output_dir=extract_dir,
            sample_name=sample.name,
            paired_end=sample.paired_end,
            multicore=cfg.methyl_extractor_multicore,
            ignore_r1=cfg.ignore_r1,
            ignore_r2=cfg.ignore_r2,
            dry_run=self.dry_run,
        )

        # --- Step 5: Methylation Calling ---
        result_dir = os.path.join(self.output_dir, sample.name)
        print(f"[pipeline:{sample.name}] 甲基化调用 (Python)...")
        mc_file = run_methylation_calling(
            extract_dir=extract_dir,
            output_dir=result_dir,
            sample_name=sample.name,
            fasta_path=self.fasta_path,
            non_conversion_rate=cfg.bisulfite_non_conversion_rate,
            bh_grouping=cfg.bh_correction_grouping,
            cores=cfg.methyl_extractor_multicore,
        )

        # --- Step 6: Summary ---
        print(f"[pipeline:{sample.name}] 生成摘要...")
        summary_file = generate_summary(
            mc_level_file=mc_file,
            genome_c_stats=c_stats,
            output_dir=result_dir,
            sample_name=sample.name,
        )

        print(f"\n[pipeline:{sample.name}] ✓ 完成!")
        print(f"  mC_level: {mc_file}")
        print(f"  summary:  {summary_file}")

        return {"mc_level_file": mc_file, "summary_file": summary_file}
