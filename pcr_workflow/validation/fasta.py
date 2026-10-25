# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path

from .constants import IUPAC_DNA

def read_fasta(path: Path) -> list[tuple[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    records: list[tuple[str, str]] = []
    current: str | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal current, seq_parts
        if current is None:
            return

        sequence = "".join(seq_parts).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"Empty sequence for FASTA record: {current}")

        invalid = sorted(set(sequence) - IUPAC_DNA)
        if invalid:
            raise ValueError(
                f"Invalid nucleotide character(s) in {current}: "
                + ", ".join(invalid)
            )

        records.append((current, sequence))

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">"):
                flush()
                current = line[1:].split()[0]
                seq_parts = []
            else:
                if current is None:
                    raise ValueError(
                        "Sequence encountered before the first FASTA header."
                    )
                seq_parts.append(line)

    flush()

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return records
def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
def read_fasta_with_full_headers(path: Path) -> list[tuple[str, str]]:
    """
    Read FASTA records while preserving the complete header text.

    This is used for database filtering/formatting checks so descriptive
    headers are not silently shortened.
    """
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal header, seq_parts
        if header is None:
            return

        sequence = "".join(seq_parts).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"Empty sequence for FASTA record: {header}")

        invalid = sorted(set(sequence) - IUPAC_DNA)
        if invalid:
            raise ValueError(
                f"Invalid nucleotide character(s) in {header}: "
                + ", ".join(invalid)
            )

        records.append((header, sequence))

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                seq_parts = []
            else:
                if header is None:
                    raise ValueError(
                        "Sequence encountered before the first FASTA header."
                    )
                seq_parts.append(line)

    flush()

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return records
def write_fasta_with_full_headers(
    path: Path,
    records: list[tuple[str, str]],
) -> None:
    """Write FASTA records while preserving complete header descriptions."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
def fasta_statistics(path: Path) -> dict[str, float | int]:
    records = read_fasta_with_full_headers(path)
    lengths = [len(seq) for _, seq in records]
    sorted_lengths = sorted(lengths)
    n = len(sorted_lengths)

    if n % 2:
        median_length = float(sorted_lengths[n // 2])
    else:
        median_length = (
            sorted_lengths[n // 2 - 1] + sorted_lengths[n // 2]
        ) / 2.0

    ambiguous = sum(
        sum(base not in {"A", "C", "G", "T"} for base in seq)
        for _, seq in records
    )

    return {
        "sequences": len(records),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": sum(lengths) / len(lengths),
        "median_length": median_length,
        "ambiguous": ambiguous,
    }
def print_fasta_statistics(path: Path, title: str) -> None:
    stats = fasta_statistics(path)
    print(f"\n{title}")
    print("=" * len(title))
    print(f"FASTA:                 {path}")
    print(f"Sequences:             {stats['sequences']}")
    print(
        f"Sequence length range: {stats['min_length']}-"
        f"{stats['max_length']} nt"
    )
    print(f"Mean sequence length:  {stats['mean_length']:.2f} nt")
    print(f"Median sequence length:{stats['median_length']:.2f} nt")
    print(f"Ambiguous bases:       {stats['ambiguous']}")
