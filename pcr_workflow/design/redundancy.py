# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path

from .commands import run_command

def cd_hit_word_size(identity: float) -> int:
    """Return a CD-HIT-EST word size compatible with the identity threshold."""

    if 0.95 <= identity <= 1.0:
        return 10
    if 0.90 <= identity < 0.95:
        return 8
    if 0.88 <= identity < 0.90:
        return 7
    if 0.85 <= identity < 0.88:
        return 6
    if 0.80 <= identity < 0.85:
        return 5
    raise ValueError("CD-HIT-EST identity must be between 0.80 and 1.0.")
def run_seqkit_deduplication(
    input_file: Path,
    output_file: Path,
    seqkit_executable: str,
) -> None:
    """Remove exact same-strand duplicate sequences with SeqKit."""

    run_command(
        [
            seqkit_executable,
            "rmdup",
            "-s",
            "-P",
            "-o",
            str(output_file),
            str(input_file),
        ]
    )
def run_cdhit_clustering(
    input_file: Path,
    output_file: Path,
    *,
    identity: float,
    mode: str,
    strand: str,
    threads: int,
    cdhit_executable: str,
) -> int:
    """Cluster nucleotide sequences with CD-HIT-EST."""

    word_size = cd_hit_word_size(identity)
    g_value = "1" if mode == "accurate" else "0"
    r_value = "1" if strand == "both" else "0"

    run_command(
        [
            cdhit_executable,
            "-i",
            str(input_file),
            "-o",
            str(output_file),
            "-c",
            str(identity),
            "-n",
            str(word_size),
            "-g",
            g_value,
            "-r",
            r_value,
            "-T",
            str(threads),
            "-d",
            "0",
        ]
    )
    return word_size
