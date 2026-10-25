# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import subprocess

from . import runtime

def run_logged_command(
    command: list[str],
    *,
    stdout_file: Path | None = None,
    env: dict[str, str] | None = None,
    remove_output_on_error: bool = False,
) -> None:
    """Run a command with the workflows' shared logging and error behavior."""
    command_text = " ".join(command)
    if stdout_file is not None:
        command_text += f" > {stdout_file}"
    runtime.log_line("COMMAND: " + command_text)
    if runtime.VERBOSE:
        print("\nCommand:", command_text)
    try:
        if stdout_file is None:
            subprocess.run(command, check=True, env=env)
        else:
            stdout_file.parent.mkdir(parents=True, exist_ok=True)
            with stdout_file.open("w", encoding="utf-8") as handle:
                subprocess.run(command, check=True, stdout=handle, env=env)
    except subprocess.CalledProcessError as error:
        if stdout_file is not None and remove_output_on_error:
            stdout_file.unlink(missing_ok=True)
        raise RuntimeError(
            f"Command failed with exit code {error.returncode}: " + command_text
        ) from error
