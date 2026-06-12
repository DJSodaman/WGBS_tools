"""
甲基化调用核心模块 (v2 - strand-merged key)
"""
import gzip, os, time
from collections import defaultdict
from multiprocessing import Pool
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import binom as scipy_binom

PREFIX_CTX = {
    "CpG_context":"CG","CpG_OT":"CG","CpG_CTOB":"CG","CpG_CTOT":"CG","CpG_OB":"CG",
    "CHG_context":"CHG","CHG_OT":"CHG","CHG_CTOB":"CHG","CHG_CTOT":"CHG","CHG_OB":"CHG",
    "CHH_context":"CHH","CHH_OT":"CHH","CHH_CTOB":"CHH","CHH_CTOT":"CHH","CHH_OB":"CHH",
}
METH = {"Z","X","H"}; UNMETH = {"z","x","h"}

def find_extract_files(extract_dir: str) -> List[Tuple[str, str]]:
    files = []
    for fname in os.listdir(extract_dir):
        for prefix, ctx in PREFIX_CTX.items():
            if fname.startswith(prefix):
                files.append((os.path.join(extract_dir, fname), ctx)); break
    return files

def _parse_one(args):
    fp, ctx = args[0], args[1]
    local = defaultdict(lambda: [0,0])
    opener = gzip.open if fp.endswith(".gz") else open
    with opener(fp, "rt") as fh:
        for line in fh:
            if not line.strip(): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5: continue
            try: pos = int(p[3])
            except: continue
            call = p[4].strip()
            key = (p[2], pos, ctx)  # strand-merged
            if call in METH: local[key][0] += 1
            elif call in UNMETH: local[key][1] += 1
    return dict(local)

def run_methylation_calling(
    extract_dir, output_dir, sample_name, fasta_path,
    non_conversion_rate=0.01, bh_grouping="by_context",
    min_coverage=0, cores=16, min_mc=2, min_depth=5,
):
    t0 = time.time()
    os.makedirs(output_dir, exist_ok=True)

    # Load reference
    ref_seqs = {}
    with open(fasta_path) as f:
        cur, cur_seq = "", []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur: ref_seqs[cur] = "".join(cur_seq)
                cur = line[1:].split()[0]; cur_seq = []
            else: cur_seq.append(line.upper())
        if cur: ref_seqs[cur] = "".join(cur_seq)
    chrom_order = list(ref_seqs.keys())
    chrom_rank = {ch: i for i, ch in enumerate(chrom_order)}

    def get_genomic_strand(chrom, pos):
        seq = ref_seqs.get(chrom, "")
        if pos < 1 or pos > len(seq): return "+"
        return "+" if seq[pos-1] == "C" else "-"

    # Find & parse files
    extract_files = find_extract_files(extract_dir)
    if not extract_files:
        raise FileNotFoundError(f"在 {extract_dir} 中未找到 Bismark 输出文件")
    print(f"[caller:{sample_name}] 找到 {len(extract_files)} 个 extract 文件")

    print(f"[caller:{sample_name}] 并行解析 (cores={cores})...")
    tasks = [(fp, ctx) for fp, ctx in extract_files]
    with Pool(min(cores, len(tasks))) as pool:
        partial = pool.map(_parse_one, tasks)

    merged = defaultdict(lambda: [0,0])
    for d in partial:
        for k, (mc, umc) in d.items():
            merged[k][0] += mc; merged[k][1] += umc
    counts = {k: (v[0], v[1]) for k, v in merged.items()}
    print(f"[caller:{sample_name}] 总位点数: {len(counts):,}")

    # Sort & compute
    sorted_items = sorted(counts.items(), key=lambda x: (chrom_rank.get(x[0][0], 999999), x[0][1]))

    n_sites = len(sorted_items)
    k_arr = np.empty(n_sites, dtype=np.int32)
    n_arr = np.empty(n_sites, dtype=np.int32)
    ctx_arr = np.empty(n_sites, dtype='<U3')
    lev_arr = np.empty(n_sites, dtype=np.float64)
    records = []

    for idx, ((chrom, pos, ctx), (mc, umc)) in enumerate(sorted_items):
        total = mc + umc
        k_arr[idx] = mc; n_arr[idx] = max(total, 1)
        ctx_arr[idx] = ctx
        lev_arr[idx] = mc/total if total>0 else 0.0
        records.append({"chrom":chrom,"pos":pos,"strand":get_genomic_strand(chrom,pos),
                        "mc":mc,"umc":umc,"level":lev_arr[idx],"context":ctx})

    print(f"[caller:{sample_name}] 向量化二项检验 + BH 校正...")
    pvals = scipy_binom.sf(k_arr-1, n_arr, non_conversion_rate)

    corrected = np.full(n_sites, np.nan)
    for ctx in ("CG","CHG","CHH"):
        mask = ctx_arr == ctx
        if not mask.any(): continue
        cp = pvals[mask].copy()
        order = np.argsort(cp)
        ranks = np.arange(1, mask.sum()+1)
        corr = cp[order] * mask.sum() / ranks
        for i in range(len(corr)-2, -1, -1):
            if corr[i] > corr[i+1]: corr[i] = corr[i+1]
        corr = np.minimum(corr, 1.0)
        result = np.empty_like(corr); result[order] = corr
        corrected[mask] = result

    # Filter & write
    print(f"[caller:{sample_name}] 过滤 (mC>={min_mc}, depth>={min_depth})...")
    out_path = os.path.join(output_dir, f"{sample_name}.mC_level_Identification_stat.txt")
    written = 0
    with open(out_path, "w") as fh:
        for i, rec in enumerate(records):
            if rec["mc"] < min_mc: continue
            if rec["mc"] + rec["umc"] < min_depth: continue
            lv = f"{rec['level']:.6f}".rstrip("0")
            if lv.endswith("."): lv += "0"
            fh.write(f"{rec['chrom']}\t{rec['pos']}\t{rec['strand']}\t{rec['mc']}\t{rec['umc']}\t{lv}\t{rec['context']}\t{pvals[i]:.15e}\t{corrected[i]:.15e}\n")
            written += 1

    elapsed = time.time()-t0
    print(f"[caller:{sample_name}] ✓ 输出 {written:,} 位点 ({elapsed:.0f}s)")
    return out_path
