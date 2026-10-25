# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import os
import subprocess

VERBOSE = False
LOG_PATH: Path | None = None

def available_cpu_threads() -> int:
    """Return logical CPU threads available to this workflow process."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)

def automatic_thread_count(cpu_threads: int) -> int:
    """Keep two logical CPU threads free when possible."""
    if cpu_threads <= 2:
        return 1
    return max(1, cpu_threads - 2)

def log_line(message: str) -> None:
    """Append one line to the active run log when logging is enabled."""
    if LOG_PATH is None:
        return
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")

def section(title: str) -> None:
    print(f"\n{title}")
    print("─" * len(title))

def executable_version(executable: str) -> str:
    """Return a best-effort external tool version without failing a run."""
    for command in ([executable, "--version"], [executable, "-version"], [executable, "version"]):
        try:
            result = subprocess.run(command, check=False, capture_output=True, text=True, timeout=10)
        except Exception:
            continue
        output = (result.stdout or result.stderr or "").strip()
        if output:
            return output.splitlines()[0][:300]
    return "unknown"
