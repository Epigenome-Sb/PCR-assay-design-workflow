# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path

from ..core.commands import run_logged_command

def run_command(command: list[str]) -> None:
    run_logged_command(command)

def run_command_to_file(command: list[str], output_file: Path) -> None:
    run_logged_command(command, stdout_file=output_file, remove_output_on_error=True)
