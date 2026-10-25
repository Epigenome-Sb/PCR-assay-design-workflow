# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from collections import Counter
from pathlib import Path

from Bio import AlignIO, SeqIO

from .constants import ALLOWED_NUCLEOTIDES, GAP_CHARACTERS, VALID_BASES

def validate_fraction(
    value: float,
    name: str,
    *,
    allow_zero: bool = False,
) -> None:
    """Validate a fraction in [0,1] or (0,1]."""

    minimum_ok = value >= 0 if allow_zero else value > 0
    if not minimum_ok or value > 1:
        interval = "[0, 1]" if allow_zero else "(0, 1]"
        raise ValueError(f"{name} must be in {interval}.")
def validate_nucleotide_fasta(file_path: Path) -> None:
    """Check nucleotide content and require unique FASTA record IDs."""

    if not file_path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {file_path}")

    sequence_count = 0
    sequence_characters = 0
    invalid_characters: set[str] = set()
    ids: list[str] = []

    with file_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                sequence_count += 1
                ids.append(stripped[1:].split()[0] if stripped[1:] else "")
                continue

            upper_line = stripped.upper().replace(" ", "")
            sequence_characters += len(upper_line)
            invalid_characters.update(
                char for char in upper_line if char not in ALLOWED_NUCLEOTIDES
            )

    if sequence_count == 0:
        raise ValueError("The file does not contain any FASTA headers.")
    if sequence_characters == 0:
        raise ValueError("The FASTA file does not contain any sequence characters.")
    if any(not identifier for identifier in ids):
        raise ValueError("At least one FASTA header has no sequence identifier.")

    duplicate_ids = [identifier for identifier, count in Counter(ids).items() if count > 1]
    if duplicate_ids:
        preview = ", ".join(duplicate_ids[:10])
        raise ValueError(
            "FASTA record IDs must be unique so sequences can be tracked across "
            f"orientation/rotation steps. Duplicate IDs include: {preview}"
        )

    if invalid_characters:
        invalid_display = ", ".join(sorted(invalid_characters))
        raise ValueError(
            "The file contains non-nucleotide characters: "
            f"{invalid_display}. Allowed characters are standard DNA/RNA and "
            "IUPAC ambiguity codes plus '-' and '.'."
        )
def count_fasta_sequences(file_path: Path) -> int:
    """Count FASTA records."""

    if not file_path.is_file():
        raise FileNotFoundError(f"File not found: {file_path}")
    return sum(1 for _ in SeqIO.parse(str(file_path), "fasta"))
def calculate_fasta_statistics(file_path: Path) -> dict[str, float | int]:
    """Calculate lightweight statistics for an unaligned/preprocessed FASTA."""

    records = list(SeqIO.parse(str(file_path), "fasta"))
    if not records:
        raise ValueError(f"No sequences found in FASTA file: {file_path}")

    ungapped_sequences: list[str] = []
    lengths: list[int] = []
    total_gaps = 0
    ambiguous_count = 0
    canonical_count = 0
    gc_count = 0

    for record in records:
        raw = str(record.seq).upper()
        total_gaps += sum(char in GAP_CHARACTERS for char in raw)
        sequence = "".join(char for char in raw if char not in GAP_CHARACTERS)
        ungapped_sequences.append(sequence)
        lengths.append(len(sequence))

        for char in sequence:
            if char in VALID_BASES:
                canonical_count += 1
                if char in {"G", "C"}:
                    gc_count += 1
            else:
                ambiguous_count += 1

    sorted_lengths = sorted(lengths)
    middle = len(sorted_lengths) // 2
    if len(sorted_lengths) % 2:
        median_length = float(sorted_lengths[middle])
    else:
        median_length = (sorted_lengths[middle - 1] + sorted_lengths[middle]) / 2

    unique_count = len(set(ungapped_sequences))
    exact_duplicate_records = len(records) - unique_count
    total_ungapped = sum(lengths)

    return {
        "sequence_count": len(records),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": total_ungapped / len(lengths),
        "median_length": median_length,
        "total_ungapped_nt": total_ungapped,
        "gap_characters": total_gaps,
        "ambiguous_characters": ambiguous_count,
        "ambiguous_percentage": (
            100.0 * ambiguous_count / total_ungapped if total_ungapped else 0.0
        ),
        "gc_percentage": (100.0 * gc_count / canonical_count if canonical_count else 0.0),
        "unique_exact_sequences": unique_count,
        "exact_duplicate_records": exact_duplicate_records,
    }
def print_fasta_statistics(
    stats: dict[str, float | int],
    *,
    title: str = "Dataset statistics",
) -> None:
    """Display concise FASTA statistics before/after preprocessing steps."""

    print(f"\n{title}")
    print("=" * len(title))
    print(f"Number of sequences:          {int(stats['sequence_count'])}")
    print(
        "Ungapped length range:       "
        f"{int(stats['min_length'])}-{int(stats['max_length'])} nt"
    )
    print(f"Mean ungapped length:         {float(stats['mean_length']):.2f} nt")
    print(f"Median ungapped length:       {float(stats['median_length']):.2f} nt")
    print(f"Mean GC (A/C/G/T only):       {float(stats['gc_percentage']):.2f} %")
    print(
        f"Ambiguous characters:         {int(stats['ambiguous_characters'])} "
        f"({float(stats['ambiguous_percentage']):.4f} %)"
    )
    print(f"Gap characters present:       {int(stats['gap_characters'])}")
    print(f"Unique exact sequences:       {int(stats['unique_exact_sequences'])}")
    print(f"Exact duplicate records:      {int(stats['exact_duplicate_records'])}")
