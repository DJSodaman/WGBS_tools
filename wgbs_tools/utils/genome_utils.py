"""
参考基因组处理工具
- FASTA 读取
- 双链 C 位点上下文判定 (CG/CHG/CHH)
- 基因组 C 统计计算
"""

from typing import Dict, Tuple, List
import os


def load_fasta(fasta_path: str) -> Dict[str, str]:
    """
    加载 FASTA 文件，返回 {chromosome_id: sequence} dict。
    所有序列大写。
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA 文件不存在: {fasta_path}")

    seqs: Dict[str, str] = {}
    current_chr = ""
    current_seq: List[str] = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_chr:
                    seqs[current_chr] = "".join(current_seq)
                current_chr = line[1:].split()[0]  # 仅取 accession
                current_seq = []
            else:
                current_seq.append(line.upper())

    if current_chr:
        seqs[current_chr] = "".join(current_seq)

    return seqs


def get_chromosome_order(fasta_path: str) -> List[str]:
    """返回 FASTA 中染色体出现顺序。"""
    order: List[str] = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                order.append(line[1:].split()[0])
    return order


def determine_context(seq: str, pos: int, strand: str) -> str:
    """
    判定给定位置 C 的三核苷酸上下文。

    参数:
        seq: 正向链序列 (大写)
        pos: 0-based 位置 (正向链)
        strand: '+' 或 '-'

    返回: 'CG', 'CHG', 'CHH'

    正向链 (+):
        位置 pos 处碱基为 C:
        - pos+1 为 G → CG
        - pos+1 非 G, pos+2 为 G → CHG
        - pos+2 非 G → CHH

    反向链 (-):
        位置 pos 处的 C 位于互补链:
        - 相当于正向链 pos 处必须为 G (互补 C)
        - context 由正向链上游碱基决定:
          - pos-1 为 C (即 CpG 的反向互补) → CG
          - pos-2 为 C, pos-1 非 C → CHG
          - pos-2 非 C → CHH
    """
    if strand == "+":
        # 检查正向链下游
        if pos + 1 < len(seq) and seq[pos + 1] == "G":
            return "CG"
        elif pos + 2 < len(seq) and seq[pos + 2] == "G":
            return "CHG"
        else:
            return "CHH"
    else:
        # 反向链: 反向互补 C 位于 pos
        # CG 互补是 CpG → 正向链 pos 处为 G, pos-1 处为 C
        if pos - 1 >= 0 and seq[pos - 1] == "C":
            return "CG"
        # CHG: C?G → 正向链 pos 为 G, pos-1 非 C, pos-2 为 C
        elif pos - 2 >= 0 and seq[pos - 2] == "C":
            return "CHG"
        else:
            return "CHH"


def compute_genome_c_stats(fasta_path: str) -> Dict[str, Dict[str, int]]:
    """
    计算基因组中每条染色体的 C 统计 (双链)。

    返回:
        {
            chr: {
                "total_C": int,     # 双链 C 总数
                "C_CG": int,        # CG 上下文 C
                "C_CHG": int,       # CHG 上下文 C
                "C_CHH": int,       # CHH 上下文 C
            },
            ...
        }

    逻辑:
        对于正向链每个 C 位点:
        - 判定其上下文 → 计为正向链 C
        - 其对侧(反向链)在该位置的互补 C 也计入对应上下文
    """
    seqs = load_fasta(fasta_path)
    chrom_order = get_chromosome_order(fasta_path)

    stats: Dict[str, Dict[str, int]] = {}

    for chr_name in chrom_order:
        seq = seqs.get(chr_name, "")
        if not seq:
            continue

        chr_stats = {"total_C": 0, "C_CG": 0, "C_CHG": 0, "C_CHH": 0}

        for pos in range(len(seq)):
            base = seq[pos]

            # 正向链 C
            if base == "C":
                ctx = determine_context(seq, pos, "+")
                chr_stats["total_C"] += 1
                if ctx == "CG":
                    chr_stats["C_CG"] += 1
                elif ctx == "CHG":
                    chr_stats["C_CHG"] += 1
                else:
                    chr_stats["C_CHH"] += 1

            # 反向链 C (正向链 G → 反向链 C)
            if base == "G":
                # 该 G 在反向链上对应一个 C
                # context 由 determine_context(seq, pos, '-') 判定
                ctx = determine_context(seq, pos, "-")
                chr_stats["total_C"] += 1
                if ctx == "CG":
                    chr_stats["C_CG"] += 1
                elif ctx == "CHG":
                    chr_stats["C_CHG"] += 1
                else:
                    chr_stats["C_CHH"] += 1

        stats[chr_name] = chr_stats

    return stats


def get_genome_totals(stats: Dict[str, Dict[str, int]]) -> Dict[str, int]:
    """汇总所有染色体的 C 统计。"""
    total = {"total_C": 0, "C_CG": 0, "C_CHG": 0, "C_CHH": 0}
    for chr_stats in stats.values():
        for key in total:
            total[key] += chr_stats[key]
    return total
