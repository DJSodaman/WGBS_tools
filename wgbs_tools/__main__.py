"""
WGBS_tools CLI 入口
===================
用法:
    python -m wgbs_tools --ref_dir /path/to/ref --output_dir /path/to/output --samples samples.yaml

环境变量:
    WGBS_TOOLS_DRY_RUN=1  → 仅打印命令，不执行
"""

import argparse
import os
import sys
import time

from .config import AnalysisConfig, DEFAULT_REF_DIR, DEFAULT_OUTPUT_DIR
from .pipeline import Pipeline, load_samples


def main():
    parser = argparse.ArgumentParser(
        description="WGBS_tools - WGBS 甲基化分析管线",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 分析用户自己的数据
  python -m wgbs_tools \\
      --ref_dir /path/to/ref \\
      --output_dir /path/to/output \\
      --samples my_samples.yaml

  # dry-run 模式 (仅打印命令)
  WGBS_TOOLS_DRY_RUN=1 python -m wgbs_tools --samples my_samples.yaml
        """,
    )

    parser.add_argument(
        "--ref_dir", type=str, default=DEFAULT_REF_DIR,
        help=f"参考基因组目录 (含 GCF_*_genomic.fna). 默认: {DEFAULT_REF_DIR}"
    )
    parser.add_argument(
        "--output_dir", type=str, default=DEFAULT_OUTPUT_DIR,
        help=f"输出目录. 默认: {DEFAULT_OUTPUT_DIR}"
    )
    parser.add_argument(
        "--samples", type=str, required=True,
        help="样本清单 YAML 文件路径"
    )
    parser.add_argument(
        "--cores", type=int, default=112,
        help="可用 CPU 核心总数 (默认: 112)"
    )

    args = parser.parse_args()

    # Dry-run 模式
    dry_run = os.environ.get("WGBS_TOOLS_DRY_RUN", "") == "1"

    if dry_run:
        print("⚠ DRY-RUN 模式 — 仅打印命令，不执行")

    # 验证 ref_dir
    ref_dir = os.path.abspath(args.ref_dir)
    if not os.path.isdir(ref_dir):
        print(f"[错误] ref_dir 不存在: {ref_dir}")
        sys.exit(1)

    # 验证 samples YAML
    if not os.path.exists(args.samples):
        print(f"[错误] 样本清单文件不存在: {args.samples}")
        sys.exit(1)

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # 加载样本
    print(f"[main] 加载样本清单: {args.samples}")
    samples = load_samples(args.samples)
    print(f"[main] 共 {len(samples)} 个样本")

    for i, s in enumerate(samples):
        mode = "PE" if s.paired_end else "SE"
        trim = "skip" if s.skip_trimming else "trim"
        print(f"  [{i+1}] {s.name} ({mode}, {trim})")
        if s.paired_end:
            print(f"      R1: {s.read1}")
            print(f"      R2: {s.read2}")
        else:
            print(f"      Read: {s.read1}")

    # 构建配置
    config = AnalysisConfig(
        total_cores=args.cores,
    )

    # 自动调整并行参数
    if args.cores >= 112:
        config.bismark_parallel_instances = 8
        config.bismark_bowtie2_threads = 13
        config.methyl_extractor_multicore = 16
    elif args.cores >= 56:
        config.bismark_parallel_instances = 6
        config.bismark_bowtie2_threads = 8
        config.methyl_extractor_multicore = 12
    elif args.cores >= 28:
        config.bismark_parallel_instances = 4
        config.bismark_bowtie2_threads = 6
        config.methyl_extractor_multicore = 8
    else:
        config.bismark_parallel_instances = 2
        config.bismark_bowtie2_threads = 4
        config.methyl_extractor_multicore = 4

    # 创建管线
    pipeline = Pipeline(
        ref_dir=ref_dir,
        output_dir=output_dir,
        config=config,
        dry_run=dry_run,
    )

    # 基因组准备 (一次性)
    t0_genome = time.time()
    c_stats = pipeline.setup_genome()
    print(f"[main] 基因组准备完成 ({time.time() - t0_genome:.0f}s)")

    # 运行每个样本
    print(f"\n{'=' * 60}")
    print(f"[main] 开始处理 {len(samples)} 个样本...")
    print(f"{'=' * 60}")

    t0_all = time.time()
    results = {}
    for i, sample in enumerate(samples):
        print(f"\n{'#' * 50}")
        print(f"### 样本 [{i+1}/{len(samples)}]: {sample.name}")
        print(f"{'#' * 50}")

        t0_sample = time.time()
        try:
            result = pipeline.run_sample(sample, c_stats)
            results[sample.name] = {"status": "success", **result}
            print(f"  耗时: {time.time() - t0_sample:.0f}s")
        except Exception as e:
            print(f"[错误] 样本 {sample.name} 失败: {e}")
            results[sample.name] = {"status": "error", "error": str(e)}

    total_elapsed = time.time() - t0_all

    # 汇总
    print(f"\n{'=' * 60}")
    print(f"[main] 分析完成")
    print(f"{'=' * 60}")
    success = sum(1 for r in results.values() if r.get("status") == "success")
    failed = len(results) - success
    print(f"  成功: {success}, 失败: {failed}, 总耗时: {total_elapsed:.0f}s")
    print(f"  输出目录: {output_dir}")

    for name, res in results.items():
        status = "✓" if res.get("status") == "success" else "✗"
        print(f"  {status} {name}")
        if res.get("mc_level_file"):
            print(f"    mC_level: {res['mc_level_file']}")
        if res.get("summary_file"):
            print(f"    summary:  {res['summary_file']}")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
