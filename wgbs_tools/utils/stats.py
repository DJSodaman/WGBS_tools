"""
统计工具: 二项检验 + Benjamini-Hochberg 校正
"""

import numpy as np
from scipy.stats import binom
from typing import List, Tuple


def binomial_test_right_tailed(k: int, n: int, p_null: float = 0.01) -> float:
    """
    右侧二项检验。
    H0: 甲基化概率 = p_null (非转化率)
    H1: 甲基化概率 > p_null
    返回 p-value = P(X >= k | n, p_null)
    """
    if n == 0:
        return 1.0
    if k == 0:
        return 1.0
    # scipy binom.sf(k-1, n, p) = P(X >= k)
    return binom.sf(k - 1, n, p_null)


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR 校正。
    输入: 原始 p-value 数组
    输出: BH 校正后的 p-value 数组 (保持原顺序)
    处理 NaN 值。
    """
    n = len(pvalues)
    if n == 0:
        return pvalues

    # 找出非 NaN 的索引
    valid_mask = ~np.isnan(pvalues)
    valid_idx = np.where(valid_mask)[0]
    n_valid = len(valid_idx)

    if n_valid == 0:
        return pvalues

    result = np.full_like(pvalues, np.nan, dtype=float)
    valid_pvals = pvalues[valid_idx]

    # 按 p-value 排序
    sort_idx = np.argsort(valid_pvals)
    sorted_pvals = valid_pvals[sort_idx]

    # BH 校正: padj[i] = p[i] * n / rank[i], 然后限制为 1.0
    ranks = np.arange(1, n_valid + 1)
    corrected = sorted_pvals * n_valid / ranks

    # 确保单调性 (从后往前)
    for i in range(n_valid - 2, -1, -1):
        if corrected[i] > corrected[i + 1]:
            corrected[i] = corrected[i + 1]

    # 限制为 1.0
    corrected = np.minimum(corrected, 1.0)

    # 恢复原顺序
    result[valid_idx[sort_idx]] = corrected
    return result


def compute_effective_size(counts: List[int]) -> float:
    """
    计算有效样本量 (harmonic mean)。
    用于稳健方差估计。
    """
    if not counts or sum(counts) == 0:
        return 0.0
    arr = np.array(counts)
    arr = arr[arr > 0]
    if len(arr) == 0:
        return 0.0
    return len(arr) / np.sum(1.0 / arr)
