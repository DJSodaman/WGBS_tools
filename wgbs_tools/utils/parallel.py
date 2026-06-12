"""
CPU 资源管理与并行调度
"""

import os
import math
from typing import List, Callable, Any
from concurrent.futures import ThreadPoolExecutor, as_completed


def detect_cores() -> int:
    """自动检测可用 CPU 核心数。"""
    return len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count() or 1


def allocate_bismark_resources(total_cores: int) -> dict:
    """
    分配 Bismark 并行资源。
    返回 {"parallel": int, "bowtie2_threads": int}

    Bismark 总核数 ≈ parallel * (bowtie2_threads + 1 overhead)
    """
    best = {"parallel": 1, "bowtie2_threads": 4}
    best_usage = 0

    for parallel in range(1, 21):
        for bt2 in range(2, 33):
            used = parallel * (bt2 + 1)
            if used <= total_cores and used > best_usage:
                best = {"parallel": parallel, "bowtie2_threads": bt2}
                best_usage = used

    return best


def run_parallel(func: Callable, items: List[Any], max_workers: int = 4,
                 desc: str = "") -> List[Any]:
    """
    用 ThreadPoolExecutor 并行执行函数 (适合 I/O 密集操作)。
    返回结果列表，失败项为 None。
    """
    results = [None] * len(items)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_idx = {
            executor.submit(func, item): i
            for i, item in enumerate(items)
        }
        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            try:
                results[idx] = future.result()
            except Exception as e:
                print(f"  [并行任务 {desc}] 项 {idx} 失败: {e}")
                results[idx] = None
    return results
