# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""Final MAFFT alignment, quality control and optional trimAl trimming."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pandas as pd
from Bio import AlignIO

from . import run_management as _run_management
from .fasta import GAP_CHARACTERS
from .questions import ask_required_choice, ask_yes_no, validate_fraction
from .run_management import log_line, run_command


MAFFT_STRATEGIES: dict[str, dict[str, object]] = {
    "auto": {
        "description": (
            "Auto — MAFFT automatically selects an appropriate strategy "
            "according to dataset size."
        ),
        "arguments": ["--auto"],
    },
    "fft-ns-1": {
        "description": (
            "FFT-NS-1 — very fast progressive alignment; useful for very "
            "large datasets."
        ),
        "arguments": ["--retree", "1", "--maxiterate", "0"],
    },
    "fft-ns-2": {
        "description": (
            "FFT-NS-2 — fast progressive alignment with two guide-tree "
            "calculations."
        ),
        "arguments": ["--retree", "2", "--maxiterate", "0"],
    },
    "fft-ns-i-2": {
        "description": "FFT-NS-i (2 cycles) — fast iterative refinement.",
        "arguments": ["--retree", "2", "--maxiterate", "2"],
    },
    "fft-ns-i-1000": {
        "description": (
            "FFT-NS-i (up to 1000 cycles) — intensive iterative refinement "
            "while remaining more scalable than pairwise strategies."
        ),
        "arguments": ["--retree", "2", "--maxiterate", "1000"],
    },
    "nw-ns-2": {
        "description": "NW-NS-2 — progressive alignment without FFT approximation.",
        "arguments": ["--retree", "2", "--maxiterate", "0", "--nofft"],
    },
    "nw-ns-i-2": {
        "description": "NW-NS-i (2 cycles) — iterative refinement without FFT.",
        "arguments": ["--retree", "2", "--maxiterate", "2", "--nofft"],
    },
    "nw-ns-i-1000": {
        "description": (
            "NW-NS-i (up to 1000 cycles) — intensive iterative refinement "
            "without FFT."
        ),
        "arguments": ["--retree", "2", "--maxiterate", "1000", "--nofft"],
    },
    "l-ins-i": {
        "description": (
            "L-INS-i — high accuracy for locally alignable regions with "
            "variable flanks; generally for relatively small datasets."
        ),
        "arguments": ["--localpair", "--maxiterate", "1000"],
    },
    "g-ins-i": {
        "description": (
            "G-INS-i — high accuracy for globally alignable sequences of "
            "similar length; generally for relatively small datasets."
        ),
        "arguments": ["--globalpair", "--maxiterate", "1000"],
    },
    "e-ins-i": {
        "description": (
            "E-INS-i — high accuracy when conserved regions are separated by "
            "large difficult-to-align regions."
        ),
        "arguments": ["--ep", "0", "--genafpair", "--maxiterate", "1000"],
    },
    "parttree": {
        "description": (
            "NW-NS-PartTree-1 — scalable strategy intended for extremely "
            "large datasets."
        ),
        "arguments": [
            "--retree",
            "1",
            "--maxiterate",
            "0",
            "--nofft",
            "--parttree",
        ],
    },
}


# ---------------------------------------------------------------------------
# Final MAFFT alignment
# ---------------------------------------------------------------------------

def choose_mafft_strategy() -> str:
    choices = [
        (key, str(value["description"]))
        for key, value in MAFFT_STRATEGIES.items()
    ]
    return ask_required_choice("Final MAFFT alignment strategy", choices)


def warn_about_mafft_strategy(strategy: str, sequence_count: int) -> None:
    """Warn about unusually expensive or mismatched MAFFT choices."""

    if strategy in {"l-ins-i", "g-ins-i", "e-ins-i"} and sequence_count > 200:
        print(
            f"\nWarning: {strategy} is computationally expensive and is generally "
            f"used for relatively small datasets; this dataset has {sequence_count} sequences."
        )
        if not ask_yes_no("Continue with this MAFFT strategy?"):
            raise RuntimeError("MAFFT strategy rejected by the user.")

    if strategy == "parttree" and sequence_count < 10000:
        print(
            f"\nWarning: PartTree is intended for extremely large datasets, whereas "
            f"this dataset contains {sequence_count} sequences."
        )
        if not ask_yes_no("Continue with PartTree?"):
            raise RuntimeError("MAFFT PartTree strategy rejected by the user.")


def run_final_mafft(
    input_file: Path,
    output_file: Path,
    strategy: str,
    threads: int,
    mafft_executable: str,
) -> str:
    """
    Run the final multiple-sequence alignment and return the strategy reported
    by MAFFT when it can be detected from stderr.
    """

    mafft_arguments = list(MAFFT_STRATEGIES[strategy]["arguments"])
    command = [
        mafft_executable,
        *mafft_arguments,
        "--thread",
        str(threads),
        str(input_file),
    ]
    command_text = " ".join(command) + f" > {output_file}"
    log_line("COMMAND: " + command_text)

    if _run_management.VERBOSE:
        print("\nCommand:", command_text)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    stderr_lines: list[str] = []

    try:
        with output_file.open("w", encoding="utf-8") as handle:
            process = subprocess.Popen(
                command,
                stdout=handle,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
            )

            assert process.stderr is not None
            for line in process.stderr:
                stderr_lines.append(line)
                print(line, end="")

            return_code = process.wait()
    except OSError as error:
        raise RuntimeError(
            "Could not start MAFFT command: " + command_text
        ) from error

    if return_code != 0:
        raise RuntimeError(
            f"Command failed with exit code {return_code}: {command_text}"
        )

    stderr_text = "".join(stderr_lines)
    match = re.search(
        r"Strategy:\s*\n\s*([^\n]+)",
        stderr_text,
        flags=re.IGNORECASE,
    )

    if not match:
        return strategy

    reported = match.group(1).strip()
    # MAFFT may append a short explanation in parentheses.
    return reported.split("(", 1)[0].strip()




# ---------------------------------------------------------------------------
# Alignment QC and trimming
# ---------------------------------------------------------------------------

def read_alignment(file_path: Path):
    """Read and validate a FASTA multiple sequence alignment."""

    if not file_path.is_file():
        raise FileNotFoundError(f"Alignment file not found: {file_path}")

    alignment = AlignIO.read(file_path, "fasta")
    if len(alignment) == 0:
        raise ValueError("The alignment contains no sequences.")
    if alignment.get_alignment_length() == 0:
        raise ValueError("The alignment has zero columns.")
    return alignment


def calculate_alignment_statistics(
    alignment_file: Path,
) -> tuple[dict[str, float | int], pd.DataFrame]:
    """Calculate global, per-column and per-sequence gap/occupancy statistics."""

    alignment = read_alignment(alignment_file)
    n_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()

    column_gap_counts: list[int] = []
    for position in range(alignment_length):
        column_gap_counts.append(
            sum(
                str(record.seq[position]).upper() in GAP_CHARACTERS
                for record in alignment
            )
        )

    total_gaps = sum(column_gap_counts)
    total_characters = n_sequences * alignment_length
    mean_gap_content = 100.0 * total_gaps / total_characters
    gap_fractions = [count / n_sequences for count in column_gap_counts]

    columns_gt20 = sum(value > 0.20 for value in gap_fractions)
    columns_gt50 = sum(value > 0.50 for value in gap_fractions)
    fully_occupied = sum(value == 0 for value in gap_fractions)

    sequence_rows: list[dict[str, object]] = []
    sequence_occupancies: list[float] = []

    for record in alignment:
        sequence = str(record.seq).upper()
        gap_count = sum(char in GAP_CHARACTERS for char in sequence)
        occupancy = 1.0 - gap_count / alignment_length
        sequence_occupancies.append(occupancy)
        sequence_rows.append(
            {
                "sequence_id": record.id,
                "gap_count": gap_count,
                "gap_fraction": gap_count / alignment_length,
                "occupancy": occupancy,
            }
        )

    stats: dict[str, float | int] = {
        "number_of_sequences": n_sequences,
        "alignment_length": alignment_length,
        "total_gaps": total_gaps,
        "mean_gap_content_percent": mean_gap_content,
        "columns_gt20_gaps": columns_gt20,
        "columns_gt50_gaps": columns_gt50,
        "fully_occupied_columns": fully_occupied,
        "mean_sequence_occupancy_percent": 100.0
        * sum(sequence_occupancies)
        / len(sequence_occupancies),
        "lowest_sequence_occupancy_percent": 100.0 * min(sequence_occupancies),
        "sequences_lt90_occupancy": sum(value < 0.90 for value in sequence_occupancies),
        "sequences_lt80_occupancy": sum(value < 0.80 for value in sequence_occupancies),
    }

    return stats, pd.DataFrame(sequence_rows)


def trimming_recommendation(stats: dict[str, float | int]) -> tuple[str, str]:
    """Return a simple transparent heuristic assessment for trimming."""

    length = int(stats["alignment_length"])
    mean_gap = float(stats["mean_gap_content_percent"])
    frac_gt20 = int(stats["columns_gt20_gaps"]) / length
    frac_gt50 = int(stats["columns_gt50_gaps"]) / length
    lowest_occ = float(stats["lowest_sequence_occupancy_percent"])

    if mean_gap >= 15 or frac_gt50 >= 0.05 or lowest_occ < 70:
        return (
            "HIGHLY GAPPY",
            "Trimming should be considered, but first inspect whether the gaps "
            "represent partial/problematic sequences or genuine biological indels.",
        )

    if mean_gap >= 5 or frac_gt20 >= 0.05 or lowest_occ < 85:
        return (
            "MODERATELY GAPPY",
            "Trimming is optional. Review gappy columns and low-occupancy sequences "
            "before deciding which alignment to use for assay design.",
        )

    return (
        "GOOD",
        "Trimming is probably not necessary. A conservative cleanup can still be "
        "used if desired.",
    )


def print_alignment_statistics(
    stats: dict[str, float | int],
    *,
    title: str = "Alignment trimming checkpoint",
) -> None:
    """Display alignment statistics and a transparent recommendation."""

    status, recommendation = trimming_recommendation(stats)

    print(f"\n{title}")
    print("=" * len(title))
    print("\nAlignment statistics:")
    print("---------------------")
    print(f"Number of sequences:          {int(stats['number_of_sequences'])}")
    print(f"Alignment length:             {int(stats['alignment_length'])} positions")
    print(f"Mean gap content:             {float(stats['mean_gap_content_percent']):.2f} %")
    print(f"Columns with >20% gaps:       {int(stats['columns_gt20_gaps'])}")
    print(f"Columns with >50% gaps:       {int(stats['columns_gt50_gaps'])}")
    print(f"Fully occupied columns:       {int(stats['fully_occupied_columns'])}")
    print(
        "Mean sequence occupancy:     "
        f"{float(stats['mean_sequence_occupancy_percent']):.2f} %"
    )
    print(
        "Lowest sequence occupancy:   "
        f"{float(stats['lowest_sequence_occupancy_percent']):.2f} %"
    )
    print(
        "Sequences with <90% occupancy: "
        f"{int(stats['sequences_lt90_occupancy'])}"
    )
    print(
        "Sequences with <80% occupancy: "
        f"{int(stats['sequences_lt80_occupancy'])}"
    )
    print("\nHeuristic assessment:")
    print("---------------------")
    print(status)
    print(recommendation)
    print(
        "Note: this recommendation is a workflow heuristic, not a biological "
        "proof that columns should be removed."
    )


def save_alignment_statistics(
    stats: dict[str, float | int],
    per_sequence_df: pd.DataFrame,
    results_dir: Path,
    project_name: str,
    label: str,
) -> None:
    """Save alignment QC summaries."""

    results_dir.mkdir(parents=True, exist_ok=True)
    summary_file = results_dir / f"{project_name}_{label}_alignment_qc.txt"
    sequence_file = results_dir / f"{project_name}_{label}_sequence_occupancy.csv"

    status, recommendation = trimming_recommendation(stats)
    lines = [
        "ALIGNMENT QUALITY SUMMARY",
        "=" * 55,
        f"Dataset label: {label}",
        f"Number of sequences: {int(stats['number_of_sequences'])}",
        f"Alignment length: {int(stats['alignment_length'])}",
        f"Mean gap content: {float(stats['mean_gap_content_percent']):.4f} %",
        f"Columns with >20% gaps: {int(stats['columns_gt20_gaps'])}",
        f"Columns with >50% gaps: {int(stats['columns_gt50_gaps'])}",
        f"Fully occupied columns: {int(stats['fully_occupied_columns'])}",
        (
            "Mean sequence occupancy: "
            f"{float(stats['mean_sequence_occupancy_percent']):.4f} %"
        ),
        (
            "Lowest sequence occupancy: "
            f"{float(stats['lowest_sequence_occupancy_percent']):.4f} %"
        ),
        f"Sequences with <90% occupancy: {int(stats['sequences_lt90_occupancy'])}",
        f"Sequences with <80% occupancy: {int(stats['sequences_lt80_occupancy'])}",
        "",
        f"Heuristic assessment: {status}",
        recommendation,
        "",
        "This assessment is advisory; inspect biological indels, partial sequences,",
        "and terminal gaps before deciding whether to trim the alignment.",
    ]
    summary_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    per_sequence_df.to_csv(sequence_file, index=False)


def choose_trimming_strategy() -> str:
    """Choose a trimAl strategy or keep the original alignment."""

    return ask_required_choice(
        "trimAl strategy",
        [
            (
                "noallgaps",
                "Remove only all-gap columns — very conservative trimming "
                "[trimAl -noallgaps].",
            ),
            (
                "gappyout",
                "Gappyout — automatically remove excessively gappy columns "
                "[trimAl -gappyout].",
            ),
            (
                "automated1",
                "Automated1 — let trimAl heuristically select an automatic method "
                "[trimAl -automated1].",
            ),
            (
                "manual",
                "Manual gap threshold — define the minimum fraction of sequences "
                "without a gap in retained columns [trimAl -gt VALUE].",
            ),
            (
                "none",
                "Cancel trimming and keep the original MAFFT alignment.",
            ),
        ],
    )


def run_trimal(
    input_alignment: Path,
    output_alignment: Path,
    strategy: str,
    trimal_executable: str,
    manual_gap_threshold: float | None = None,
) -> None:
    """Run trimAl with one of the supported column-trimming modes."""

    command = [
        trimal_executable,
        "-in",
        str(input_alignment),
        "-out",
        str(output_alignment),
        "-fasta",
        "-keepheader",
        "-keepseqs",
    ]

    if strategy == "noallgaps":
        command.append("-noallgaps")
    elif strategy == "gappyout":
        command.append("-gappyout")
    elif strategy == "automated1":
        command.append("-automated1")
    elif strategy == "manual":
        if manual_gap_threshold is None:
            raise ValueError("A trimAl -gt value is required for manual trimming.")
        validate_fraction(
            manual_gap_threshold,
            "--trimal-gap-threshold",
            allow_zero=True,
        )
        command.extend(["-gt", str(manual_gap_threshold)])
    else:
        raise ValueError(f"Unsupported trimAl strategy: {strategy}")

    output_alignment.parent.mkdir(parents=True, exist_ok=True)
    run_command(command)
    read_alignment(output_alignment)


