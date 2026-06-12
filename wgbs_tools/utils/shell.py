"""
子进程执行与日志工具
"""

import subprocess
import sys
import time
import os
from typing import Optional, List


class ShellError(Exception):
    """子进程执行失败异常。"""
    def __init__(self, cmd: str, returncode: int, stderr: str):
        self.cmd = cmd
        self.returncode = returncode
        self.stderr = stderr
        super().__init__(f"命令失败 (rc={returncode}): {cmd}\n{stderr}")


def run_cmd(
    cmd: str,
    log_prefix: str = "",
    cwd: Optional[str] = None,
    env: Optional[dict] = None,
    timeout: Optional[int] = None,
    dry_run: bool = False,
) -> str:
    """
    执行 shell 命令并实时输出。
    返回 stdout (如有)。
    失败时抛出 ShellError。
    """
    header = f"[{log_prefix}] " if log_prefix else ""
    print(f"{header}执行: {cmd}", flush=True)

    if dry_run:
        return ""

    t0 = time.time()
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            cwd=cwd,
            env=env or os.environ.copy(),
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        raise ShellError(cmd, -1, f"超时 ({elapsed:.0f}s)")

    elapsed = time.time() - t0

    if result.stdout.strip():
        # 仅输出最后几行，避免刷屏
        lines = result.stdout.strip().split("\n")
        for line in lines[-10:]:
            print(f"{header}  {line}", flush=True)

    if result.returncode != 0:
        raise ShellError(cmd, result.returncode, result.stderr)

    print(f"{header}✓ 完成 ({elapsed:.0f}s)", flush=True)
    return result.stdout


def run_cmd_silent(
    cmd: str,
    cwd: Optional[str] = None,
    env: Optional[dict] = None,
    timeout: Optional[int] = None,
) -> subprocess.CompletedProcess:
    """静默执行命令，返回 CompletedProcess。"""
    return subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
        cwd=cwd,
        env=env or os.environ.copy(),
        timeout=timeout,
    )


def check_tool(tool_name: str) -> Optional[str]:
    """检查工具是否在 PATH 中，返回路径或 None。"""
    result = subprocess.run(
        f"which {tool_name}",
        shell=True, capture_output=True, text=True
    )
    if result.returncode == 0:
        return result.stdout.strip()
    return None
