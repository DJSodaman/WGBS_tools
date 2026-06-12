"""
甲基化摘要统计模块
生成 *.methylation_summary.stat.txt (公司格式)
"""

import os
from collections import defaultdict


def generate_summary(
    mc_level_file: str,
    genome_c_stats: dict,  # {"per_chrom": ..., "totals": ...}
    output_dir: str,
    sample_name: str,
) -> str:
    """
    从 mC_level_Identification_stat.txt 生成摘要统计文件。

    参数:
        mc_level_file: caller.py 输出的位点级文件
        genome_c_stats: genome.py 计算的基因组 C 统计
        output_dir: 输出目录
        sample_name: 样本名

    返回:
        输出的 methylation_summary.stat.txt 路径
    """
    os.makedirs(output_dir, exist_ok=True)

    # 基因组 C 统计
    totals = genome_c_stats["totals"]
    total_C_genome = totals["total_C"]
    C_CG_genome = totals["C_CG"]
    C_CHG_genome = totals["C_CHG"]
    C_CHH_genome = totals["C_CHH"]

    # 从 mC_level 文件统计
    mC_total = 0
    mC_CG = 0
    mC_CHG = 0
    mC_CHH = 0

    with open(mc_level_file) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            context = parts[6]
            # 统计位点数 (sites)，不是 read counts
            mC_total += 1
            if context == "CG":
                mC_CG += 1
            elif context == "CHG":
                mC_CHG += 1
            elif context == "CHH":
                mC_CHH += 1

    # 计算百分比
    pct_total = (mC_total / total_C_genome * 100) if total_C_genome > 0 else 0
    pct_CG = (mC_CG / C_CG_genome * 100) if C_CG_genome > 0 else 0
    pct_CHG = (mC_CHG / C_CHG_genome * 100) if C_CHG_genome > 0 else 0
    pct_CHH = (mC_CHH / C_CHH_genome * 100) if C_CHH_genome > 0 else 0

    # 组成百分比 (各上下文占总 methylated C 的比例)
    comp_CG = (mC_CG / mC_total * 100) if mC_total > 0 else 0
    comp_CHG = (mC_CHG / mC_total * 100) if mC_total > 0 else 0
    comp_CHH = (mC_CHH / mC_total * 100) if mC_total > 0 else 0

    # 写文件 (精确匹配公司格式)
    output_path = os.path.join(
        output_dir, f"{sample_name}.methylation_summary.stat.txt"
    )

    with open(output_path, "w") as fh:
        fh.write(f"Cytosines in the genome:\t{total_C_genome}\n")
        fh.write(f"Methylated C:\t{mC_total}\n")
        fh.write(f"Methylated C percent:\t{pct_total:.2f}%\n")
        fh.write(f"C in CG:\t{C_CG_genome}\n")
        fh.write(f"mC in CG:\t{mC_CG}\n")
        fh.write(f"Methylated CG percent:\t{pct_CG:.2f}%\n")
        fh.write(f"C in CHG:\t{C_CHG_genome}\n")
        fh.write(f"mC in CHG:\t{mC_CHG}\n")
        fh.write(f"Methylated CHG percent:\t{pct_CHG:.2f}%\n")
        fh.write(f"C in CHH:\t{C_CHH_genome}\n")
        fh.write(f"mC in CHH:\t{mC_CHH}\n")
        fh.write(f"Methylated CHH percent:\t{pct_CHH:.2f}%\n")
        fh.write(
            f"mC:\t{mC_total}\t100%\t"
            f"mCG:\t{mC_CG}\t{comp_CG:.2f}%\t"
            f"mCHG:\t{mC_CHG}\t{comp_CHG:.2f}%\t"
            f"mCHH:\t{mC_CHH}\t{comp_CHH:.2f}%\n"
        )

    return output_path
