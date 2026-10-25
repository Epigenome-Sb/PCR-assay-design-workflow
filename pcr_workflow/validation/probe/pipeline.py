# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path

from .blast import run_probe_blast
from .summary import summarize_probe_blast

def analyse_probe_with_blast(
    scheme: str,
    probe: str,
    valid_amplicons_fasta: Path,
    valid_rows: list[dict[str, object]],
    scheme_dir: Path,
    total_target_sequences: int,
    threads: int,
) -> dict[str, object]:
    """Run and summarize BLAST analysis for the original VarVAMP probe."""
    raw_file, query_file, variant_count = run_probe_blast(
        scheme,
        probe,
        valid_amplicons_fasta,
        valid_rows,
        scheme_dir,
        threads,
    )

    return summarize_probe_blast(
        scheme,
        probe,
        raw_file,
        query_file,
        variant_count,
        valid_rows,
        scheme_dir,
        total_target_sequences,
    )
def empty_probe_analysis() -> dict[str, object]:
    """Return a schema-compatible empty probe result."""
    return {
        "amplicon_total": 0,
        "amplicon_exact_0": 0,
        "amplicon_full_length": 0,
        "amplicon_1_mismatch": 0,
        "amplicon_2_mismatch": 0,
        "amplicon_ge3_mismatch": 0,
        "exact_0": 0,
        "mismatch_1": 0,
        "mismatch_2": 0,
        "mismatch_ge3": 0,
        "partial_only": 0,
        "no_hit": 0,
        "full_length_hitids": 0,
        "pcr_positive_hitids": 0,
        "hitid_summary": "",
        "blast_summary": "",
        "alignment_hitid_file": "",
        "mismatch_positions_file": "",
        "variant_statistics_file": "",
        "variant_hitid_matrix_file": "",
    }
