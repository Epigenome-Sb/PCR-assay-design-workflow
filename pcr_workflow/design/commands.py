# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import shutil

from ..core.commands import run_logged_command

def resolve_executable(value: str) -> str:
    """Resolve an executable path or fail with an actionable message."""
    resolved = shutil.which(value)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found: {value}")
    return resolved

def check_required_tools(tools: list[tuple[str, str]]) -> None:
    missing = [label for label, executable in tools if shutil.which(executable) is None]
    if missing:
        raise RuntimeError("Required program(s) not found: " + ", ".join(missing) + ".")

def run_command(
    command: list[str],
    *,
    stdout_file: Path | None = None,
    env: dict[str, str] | None = None,
) -> None:
    run_logged_command(command, stdout_file=stdout_file, env=env)
