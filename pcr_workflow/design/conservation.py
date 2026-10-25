# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from collections import Counter
from math import log2
from pathlib import Path

import pandas as pd

from .alignment import read_alignment
from .constants import GAP_CHARACTERS, VALID_BASES

def calculate_real_lengths(alignment) -> list[int]:
    """Calculate sequence lengths after removing alignment gaps."""

    lengths: list[int] = []
    for record in alignment:
        sequence = str(record.seq).upper()
        ungapped = "".join(
            char for char in sequence if char not in GAP_CHARACTERS
        )
        lengths.append(len(ungapped))
    return lengths
def calculate_shannon_entropy(base_counts: Counter[str]) -> float:
    """Calculate Shannon entropy for A, C, G and T counts."""

    total = sum(base_counts.values())
    if total == 0:
        return 0.0

    entropy = 0.0
    for count in base_counts.values():
        if count:
            frequency = count / total
            entropy -= frequency * log2(frequency)
    return entropy
def analyse_alignment_positions(
    alignment,
    min_occupancy: float,
    min_major_frequency: float,
) -> tuple[pd.DataFrame, str]:
    """Analyse each alignment column and build a conservative consensus."""

    number_of_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()
    results: list[dict[str, object]] = []
    consensus_bases: list[str] = []

    for position_index in range(alignment_length):
        column = [
            str(record.seq[position_index]).upper()
            for record in alignment
        ]

        gap_count = sum(char in GAP_CHARACTERS for char in column)
        non_gap_characters = [
            char for char in column if char not in GAP_CHARACTERS
        ]
        occupied_count = len(non_gap_characters)
        occupancy = occupied_count / number_of_sequences

        valid_bases = [char for char in non_gap_characters if char in VALID_BASES]
        ambiguous_count = sum(
            char not in VALID_BASES for char in non_gap_characters
        )
        base_counts = Counter(valid_bases)

        if base_counts:
            major_base, major_base_count = base_counts.most_common(1)[0]
        else:
            major_base = "N"
            major_base_count = 0

        # Preserve the behavior of the previous workflow: ambiguous characters
        # are part of occupied_count and therefore reduce major-base frequency.
        major_base_frequency = (
            major_base_count / occupied_count if occupied_count else 0.0
        )

        sufficiently_represented = occupancy >= min_occupancy
        strictly_conserved = (
            gap_count == 0
            and ambiguous_count == 0
            and major_base_count == number_of_sequences
        )
        highly_conserved = (
            sufficiently_represented
            and major_base_frequency >= min_major_frequency
        )

        consensus_base = major_base if highly_conserved else "N"
        consensus_bases.append(consensus_base)

        results.append(
            {
                "position_alignment": position_index + 1,
                "A_count": base_counts.get("A", 0),
                "C_count": base_counts.get("C", 0),
                "G_count": base_counts.get("G", 0),
                "T_count": base_counts.get("T", 0),
                "valid_base_count": len(valid_bases),
                "gap_count": gap_count,
                "ambiguous_count": ambiguous_count,
                "occupied_count": occupied_count,
                "occupancy": occupancy,
                "major_base": major_base,
                "major_base_count": major_base_count,
                "major_base_frequency": major_base_frequency,
                "shannon_entropy": calculate_shannon_entropy(base_counts),
                "sufficiently_represented": sufficiently_represented,
                "strictly_conserved": strictly_conserved,
                "highly_conserved": highly_conserved,
                "consensus_base": consensus_base,
            }
        )

    return pd.DataFrame(results), "".join(consensus_bases)
def save_consensus_fasta(
    consensus_sequence: str,
    output_file: Path,
    project_name: str,
    min_occupancy: float,
    min_major_frequency: float,
) -> None:
    """Save the conservative consensus sequence."""

    header = (
        f">{project_name}_consensus_"
        f"occupancy{min_occupancy:.3f}_major{min_major_frequency:.3f}"
    )
    with output_file.open("w", encoding="utf-8") as handle:
        handle.write(header + "\n")
        for start in range(0, len(consensus_sequence), 80):
            handle.write(consensus_sequence[start : start + 80] + "\n")
def perform_conservation_analysis(
    alignment_file: Path,
    results_dir: Path,
    project_name: str,
    min_occupancy: float,
    min_major_frequency: float,
) -> None:
    """Generate conservation table, conservative consensus and text summary."""

    alignment = read_alignment(alignment_file)
    number_of_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()

    real_lengths = calculate_real_lengths(alignment)
    minimum_length = min(real_lengths)
    maximum_length = max(real_lengths)
    average_length = sum(real_lengths) / len(real_lengths)

    results_df, consensus_sequence = analyse_alignment_positions(
        alignment,
        min_occupancy,
        min_major_frequency,
    )

    total_alignment_characters = number_of_sequences * alignment_length
    total_gap_characters = int(results_df["gap_count"].sum())
    total_ambiguous_characters = int(results_df["ambiguous_count"].sum())

    gap_percentage = 100 * total_gap_characters / total_alignment_characters
    ambiguous_percentage = (
        100 * total_ambiguous_characters / total_alignment_characters
    )

    represented_count = int(results_df["sufficiently_represented"].sum())
    strict_count = int(results_df["strictly_conserved"].sum())
    conserved_count = int(results_df["highly_conserved"].sum())

    represented_percentage = 100 * represented_count / alignment_length
    strict_percentage = 100 * strict_count / alignment_length
    conserved_total_percentage = 100 * conserved_count / alignment_length
    conserved_represented_percentage = (
        100 * conserved_count / represented_count if represented_count else 0.0
    )

    results_dir.mkdir(parents=True, exist_ok=True)
    position_file = results_dir / f"{project_name}_conservation_by_position.csv"
    consensus_file = results_dir / f"{project_name}_consensus.fasta"
    summary_file = results_dir / f"{project_name}_conservation_summary.txt"

    results_df.to_csv(position_file, index=False)
    save_consensus_fasta(
        consensus_sequence,
        consensus_file,
        project_name,
        min_occupancy,
        min_major_frequency,
    )

    summary_lines = [
        "ALIGNMENT CONSERVATION ANALYSIS",
        "=" * 55,
        "",
        "ANALYSED DATA",
        f"Sequences analysed: {number_of_sequences}",
        f"Alignment length: {alignment_length} positions",
        f"Ungapped sequence length: {minimum_length}-{maximum_length} nt",
        f"Mean ungapped sequence length: {average_length:.2f} nt",
        "",
        "ALIGNMENT QUALITY",
        (
            f"Gaps: {total_gap_characters} / {total_alignment_characters} "
            f"({gap_percentage:.2f} %)"
        ),
        (
            f"Ambiguous characters: {total_ambiguous_characters} / "
            f"{total_alignment_characters} ({ambiguous_percentage:.4f} %)"
        ),
        "",
        "SUFFICIENTLY REPRESENTED POSITIONS",
        (
            f"Occupancy >= {min_occupancy:.0%}: {represented_count} / "
            f"{alignment_length} ({represented_percentage:.2f} %)"
        ),
        "",
        "CONSERVATION",
        (
            f"Strictly conserved positions: {strict_count} / {alignment_length} "
            f"({strict_percentage:.2f} %)"
        ),
        (
            f"Positions with occupancy >= {min_occupancy:.0%} and major-base "
            f"frequency >= {min_major_frequency:.0%}: {conserved_count} / "
            f"{alignment_length} ({conserved_total_percentage:.2f} %)"
        ),
        (
            "Conservation among sufficiently represented positions: "
            f"{conserved_count} / {represented_count} "
            f"({conserved_represented_percentage:.2f} %)"
        ),
        "",
        "COORDINATE NOTE",
        "Positions are multiple-sequence-alignment coordinates and are not",
        "automatically equivalent to coordinates in a reference genome.",
    ]

    summary_text = "\n".join(summary_lines)
    summary_file.write_text(summary_text + "\n", encoding="utf-8")
    print("\n" + summary_text)
