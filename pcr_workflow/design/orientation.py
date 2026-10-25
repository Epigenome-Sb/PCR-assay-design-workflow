# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..core import runtime
from .commands import run_command
from .constants import GAP_CHARACTERS
from .fasta import count_fasta_sequences, validate_nucleotide_fasta

def calculate_cyclic_rotation_statistics(
    before_file: Path,
    after_file: Path,
) -> dict[str, float | int]:
    """Compare MARS input/output and summarize detectable cyclic shifts."""

    before = {
        record.id: "".join(
            char for char in str(record.seq).upper() if char not in GAP_CHARACTERS
        )
        for record in SeqIO.parse(str(before_file), "fasta")
    }
    after = {
        record.id: "".join(
            char for char in str(record.seq).upper() if char not in GAP_CHARACTERS
        )
        for record in SeqIO.parse(str(after_file), "fasta")
    }

    shifts: list[int] = []
    unchanged = 0
    unverified = 0

    for record_id, sequence_after in after.items():
        sequence_before = before.get(record_id)
        if (
            sequence_before is None
            or len(sequence_before) != len(sequence_after)
            or not sequence_before
        ):
            unverified += 1
            continue

        if sequence_before == sequence_after:
            unchanged += 1
            shifts.append(0)
            continue

        position = (sequence_before + sequence_before).find(sequence_after)
        if 0 < position < len(sequence_before):
            shifts.append(position)
        else:
            unverified += 1

    rotated = sum(shift > 0 for shift in shifts)
    positive_shifts = [shift for shift in shifts if shift > 0]

    return {
        "sequence_count": len(after),
        "rotated_sequences": rotated,
        "unchanged_sequences": unchanged,
        "unverified_sequences": unverified,
        "mean_shift": (
            sum(positive_shifts) / len(positive_shifts) if positive_shifts else 0.0
        ),
        "min_shift": min(positive_shifts) if positive_shifts else 0,
        "max_shift": max(positive_shifts) if positive_shifts else 0,
    }
def print_rotation_statistics(stats: dict[str, float | int]) -> None:
    """Display the effect of MARS cyclic start-position normalization."""

    count = int(stats["sequence_count"])
    rotated = int(stats["rotated_sequences"])
    percentage = 100.0 * rotated / count if count else 0.0

    print("\nMARS rotation statistics")
    print("========================")
    print(f"Sequences evaluated:          {count}")
    print(f"Cyclic shifts detected:       {rotated} ({percentage:.2f} %)")
    print(f"Start position unchanged:     {int(stats['unchanged_sequences'])}")
    print(f"Shift could not be verified:  {int(stats['unverified_sequences'])}")
    if rotated:
        print(
            "Detected shift range:         "
            f"{int(stats['min_shift'])}-{int(stats['max_shift'])} nt"
        )
        print(f"Mean detected shift:          {float(stats['mean_shift']):.2f} nt")
    print(
        "Note: MARS changes cyclic start positions, not sequence length or "
        "nucleotide composition."
    )
def copy_fasta_without_gaps(input_file: Path, output_file: Path) -> int:
    """Write an ungapped FASTA copy and return the number of removed gaps."""

    output_file.parent.mkdir(parents=True, exist_ok=True)
    records: list[SeqRecord] = []
    removed = 0

    for record in SeqIO.parse(str(input_file), "fasta"):
        sequence = str(record.seq)
        cleaned = "".join(char for char in sequence if char not in GAP_CHARACTERS)
        removed += len(sequence) - len(cleaned)
        records.append(
            SeqRecord(
                Seq(cleaned),
                id=record.id,
                name=record.name,
                description=record.description,
            )
        )

    if not records:
        raise ValueError(f"No sequences found in FASTA file: {input_file}")

    SeqIO.write(records, str(output_file), "fasta")
    return removed
def normalize_orientation_with_mafft(
    input_file: Path,
    workdir: Path,
    project_name: str,
    strategy: str,
    threads: int,
    mafft_executable: str,
) -> tuple[Path, int]:
    """Orient sequences with MAFFT and return output plus reversal count."""

    temporary_alignment = workdir / f"{project_name}_orientation_alignment.fasta"
    oriented_file = workdir / f"{project_name}_oriented.fasta"
    report_file = workdir / f"{project_name}_orientation_report.txt"

    option = "--adjustdirection" if strategy == "adjust" else "--adjustdirectionaccurately"
    run_command(
        [
            mafft_executable,
            option,
            "--thread",
            str(threads),
            str(input_file),
        ],
        stdout_file=temporary_alignment,
    )

    original_descriptions = {
        record.id: record.description
        for record in SeqIO.parse(str(input_file), "fasta")
    }

    output_records: list[SeqRecord] = []
    reversed_ids: list[str] = []

    for record in SeqIO.parse(str(temporary_alignment), "fasta"):
        mafft_id = record.id
        reversed_sequence = mafft_id.startswith("_R_")
        original_id = mafft_id[3:] if reversed_sequence else mafft_id

        if reversed_sequence:
            reversed_ids.append(original_id)

        description = original_descriptions.get(original_id, original_id)
        ungapped = "".join(
            char for char in str(record.seq) if char not in GAP_CHARACTERS
        )
        output_records.append(
            SeqRecord(
                Seq(ungapped),
                id=original_id,
                name=original_id,
                description=description,
            )
        )

    if len(output_records) != count_fasta_sequences(input_file):
        raise RuntimeError(
            "MAFFT orientation output does not contain the same number of "
            "sequences as the input."
        )

    SeqIO.write(output_records, str(oriented_file), "fasta")

    report_lines = [
        "SEQUENCE ORIENTATION REPORT",
        "=" * 55,
        f"Method: MAFFT {option}",
        f"Input sequences: {len(output_records)}",
        f"Reverse-complemented sequences: {len(reversed_ids)}",
        "",
    ]
    if reversed_ids:
        report_lines.append("Reverse-complemented record IDs:")
        report_lines.extend(reversed_ids)
    else:
        report_lines.append("No sequence required reverse-complementation.")

    report_file.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    if runtime.VERBOSE:
        print(f"Oriented FASTA: {oriented_file.resolve()}")
    return oriented_file, len(reversed_ids)
def run_mars_rotation(
    input_file: Path,
    workdir: Path,
    project_name: str,
    mars_executable: str,
    threads: int,
) -> Path:
    """Run MARS on ungapped nucleotide sequences."""

    mars_input = workdir / f"{project_name}_mars_input_ungapped.fasta"
    rotated_file = workdir / f"{project_name}_rotated.fasta"

    removed_gaps = copy_fasta_without_gaps(input_file, mars_input)
    if removed_gaps:
        print(
            f"Removed {removed_gaps} gap characters before MARS. "
            "Only temporary/alignment gaps should be present at this stage."
        )

    run_command(
        [
            mars_executable,
            "-a",
            "DNA",
            "-i",
            str(mars_input),
            "-o",
            str(rotated_file),
            "-T",
            str(threads),
        ]
    )

    validate_nucleotide_fasta(rotated_file)
    if count_fasta_sequences(rotated_file) != count_fasta_sequences(input_file):
        raise RuntimeError(
            "MARS output contains a different number of sequences than its input."
        )

    if runtime.VERBOSE:
        print(f"MARS-rotated FASTA: {rotated_file.resolve()}")
    return rotated_file
