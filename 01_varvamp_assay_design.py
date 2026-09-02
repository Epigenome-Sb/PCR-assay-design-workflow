#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""
01_varvamp_assay_design.py

Generic stepwise VarVAMP assay-design workflow for nucleotide sequence datasets.

Scientific workflow
-------------------
1. Validate the input nucleotide FASTA.
2. Optionally normalize sequence orientation with MAFFT.
3. For complete circular sequences, optionally normalize starts with MARS.
4. Handle redundancy with no filtering, SeqKit exact deduplication, or CD-HIT-EST.
5. Perform the final MAFFT multiple-sequence alignment.
6. Calculate alignment QC statistics and optionally trim with trimAl.
7. Analyse alignment conservation.
8. Design SINGLE, QPCR or TILED assays with VarVAMP.
9. Save VarVAMP outputs and visual summaries.
10. Write assay_design_manifest.json for reproducible handoff to
    02_in_silico_validation.py.

Interface
---------
Scientific decisions remain stepwise and explicit in interactive mode.

Normal terminal mode is concise.
Use --verbose to display complete external commands and detailed tables.
All commands are always recorded in results/<project>/workflow.log.

The workflow writes:
- workflow_summary.txt
- workflow.log
- assay_design_manifest.json
and copies the manifest into the selected varvamp_<mode>/ result directory.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from collections import Counter
from io import BytesIO
from math import log2
from pathlib import Path
from typing import Iterable

import matplotlib

# Allow execution on servers and in GitHub Codespaces without a display.
matplotlib.use("Agg")

import fitz
import matplotlib.pyplot as plt
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from PIL import Image


VALID_BASES = {"A", "C", "G", "T"}
GAP_CHARACTERS = {"-", "."}
ALLOWED_NUCLEOTIDES = set("ACGTURYSWKMBDHVN-.")


VERBOSE = False
WORKFLOW_LOG: Path | None = None


def log_line(message: str) -> None:
    """Append one line to the workflow log when logging is active."""
    if WORKFLOW_LOG is None:
        return
    WORKFLOW_LOG.parent.mkdir(parents=True, exist_ok=True)
    with WORKFLOW_LOG.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")


def compact_path(path: Path | str) -> str:
    """Prefer a path relative to the current working directory for display."""
    path_obj = Path(path)
    try:
        return str(path_obj.resolve().relative_to(Path.cwd().resolve()))
    except Exception:
        return str(path_obj)


def section(title: str) -> None:
    print(f"\n{title}")
    print("─" * len(title))


def progress(label: str, status: str = "OK") -> None:
    dots = "." * max(2, 52 - len(label))
    print(f"{label} {dots} {status}")


def executable_version(executable: str) -> str:
    """
    Best-effort capture of an external tool version.

    Version lookup failures never stop the scientific workflow.
    """
    candidates = (
        [executable, "--version"],
        [executable, "-version"],
        [executable, "version"],
    )

    for command in candidates:
        try:
            result = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=10,
            )
        except Exception:
            continue

        output = (result.stdout or result.stderr or "").strip()
        if output:
            return output.splitlines()[0][:300]

    return "unknown"


def write_project_manifest(
    path: Path,
    *,
    project_name: str,
    input_file: Path,
    workdir: Path,
    results_dir: Path,
    downstream_alignment: Path,
    initial_count: int,
    after_redundancy: int,
    args: argparse.Namespace,
    varvamp_results: Path | None,
) -> None:
    """
    Write a machine-readable handoff between assay design and validation.

    02_in_silico_validation.py can use this file to know which VarVAMP mode
    and result directory belong to this project.
    """
    selected_tools = {
        "mafft": executable_version(args.mafft_executable),
    }

    if args.topology == "circular":
        selected_tools["mars"] = executable_version(args.mars_executable)

    if args.redundancy == "seqkit":
        selected_tools["seqkit"] = executable_version(args.seqkit_executable)
    elif args.redundancy == "cdhit":
        selected_tools["cd-hit-est"] = executable_version(args.cdhit_executable)

    if args.trimming != "none":
        selected_tools["trimal"] = executable_version(args.trimal_executable)

    if not args.skip_varvamp:
        selected_tools["varvamp"] = executable_version(args.varvamp_executable)

    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project_name": project_name,
        "input_fasta": str(input_file.resolve()),
        "workdir": str(workdir.resolve()),
        "results_dir": str(results_dir.resolve()),
        "downstream_alignment": str(downstream_alignment.resolve()),
        "sequence_counts": {
            "initial_input": int(initial_count),
            "entering_final_mafft": int(after_redundancy),
        },
        "preprocessing": {
            "orientation": args.orientation,
            "topology": args.topology,
            "redundancy": args.redundancy,
            "cdhit_identity": args.identity,
            "cdhit_mode": args.cdhit_mode,
            "cdhit_strand": args.cdhit_strand,
            "mafft_strategy": args.mafft_strategy,
            "trimming": args.trimming,
            "trimal_gap_threshold": args.trimal_gap_threshold,
            "min_occupancy": args.min_occupancy,
            "min_major_frequency": args.min_major_frequency,
        },
        "varvamp": {
            "executed": not args.skip_varvamp,
            "mode": None if args.skip_varvamp else args.varvamp_mode,
            "result_dir": (
                None
                if varvamp_results is None
                else str(varvamp_results.resolve())
            ),
            "consensus_threshold": args.varvamp_threshold,
            "primer_ambiguity": args.primer_ambiguity,
            "probe_ambiguity": args.probe_ambiguity,
            "single_opt_length": args.single_opt_length,
            "single_max_length": args.single_max_length,
            "single_report_n": args.single_report_n,
            "tiled_opt_length": args.tiled_opt_length,
            "tiled_max_length": args.tiled_max_length,
            "tiled_overlap": args.tiled_overlap,
            "qpcr_test_n": args.qpcr_test_n,
            "qpcr_deltag": args.qpcr_deltag,
            "custom_config": (
                None
                if args.varvamp_config is None
                else str(args.varvamp_config.expanduser().resolve())
            ),
        },
        "tool_versions": selected_tools,
        "workflow_log": (
            None
            if WORKFLOW_LOG is None
            else str(WORKFLOW_LOG.resolve())
        ),
    }

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )


def write_workflow_summary(
    path: Path,
    *,
    project_name: str,
    input_file: Path,
    downstream_alignment: Path,
    results_dir: Path,
    initial_count: int,
    after_redundancy: int,
    args: argparse.Namespace,
    varvamp_results: Path | None,
) -> None:
    lines = [
        "VarVAMP ASSAY-DESIGN WORKFLOW SUMMARY",
        "=" * 55,
        f"Project: {project_name}",
        f"Input FASTA: {input_file}",
        f"Initial sequences: {initial_count}",
        f"Sequences entering final MAFFT: {after_redundancy}",
        f"Orientation: {args.orientation}",
        f"Topology: {args.topology}",
        f"Redundancy: {args.redundancy}",
        f"Final MAFFT strategy: {args.mafft_strategy}",
        f"Trimming: {args.trimming}",
        f"Downstream alignment: {downstream_alignment}",
        "",
        "VarVAMP:",
        f"Executed: {'yes' if not args.skip_varvamp else 'no'}",
        f"Mode: {args.varvamp_mode if not args.skip_varvamp else 'N/A'}",
        (
            f"Result directory: {varvamp_results}"
            if varvamp_results is not None
            else "Result directory: N/A"
        ),
        "",
        f"Results directory: {results_dir}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")




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
# Argument parsing and interactive helpers
# ---------------------------------------------------------------------------


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Nucleotide preprocessing, MAFFT alignment, alignment QC, "
            "conservation analysis and PCR/qPCR/tiled assay design with VarVAMP."
        )
    )

    parser.add_argument("--input", type=Path, default=None, help="Input FASTA file.")
    parser.add_argument("--project-name", type=str, default=None)
    parser.add_argument("--workdir", type=Path, default=None)
    parser.add_argument("--results", type=Path, default=None)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--verbose",
        action="store_true",
        help=(
            "Display complete external commands and detailed file paths. "
            "Commands are always saved to workflow.log."
        ),
    )

    parser.add_argument(
        "--orientation",
        choices=("keep", "adjust", "accurate"),
        default=None,
        help=(
            "Sequence-orientation strategy: keep, MAFFT --adjustdirection, or "
            "MAFFT --adjustdirectionaccurately."
        ),
    )
    parser.add_argument(
        "--topology",
        choices=("linear", "circular"),
        default=None,
        help="Sequence topology. MARS is run only for circular datasets.",
    )

    parser.add_argument(
        "--redundancy",
        choices=("none", "seqkit", "cdhit"),
        default=None,
        help=(
            "Redundancy strategy: keep all sequences, remove exact duplicates "
            "with SeqKit, or cluster with CD-HIT-EST."
        ),
    )
    # Backward-compatible flag from the previous workflow.
    parser.add_argument(
        "--skip-cdhit",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--identity", type=float, default=None)
    parser.add_argument(
        "--cdhit-mode", choices=("fast", "accurate"), default=None
    )
    parser.add_argument(
        "--cdhit-strand", choices=("both", "same"), default=None
    )

    parser.add_argument(
        "--mafft-strategy",
        choices=tuple(MAFFT_STRATEGIES),
        default=None,
        help="Final MAFFT alignment strategy.",
    )

    parser.add_argument(
        "--trimming",
        choices=("none", "noallgaps", "gappyout", "automated1", "manual"),
        default=None,
        help="Alignment trimming strategy after final MAFFT alignment.",
    )
    parser.add_argument(
        "--trimal-gap-threshold",
        type=float,
        default=None,
        help=(
            "trimAl -gt value for --trimming manual. This is the minimum "
            "fraction of sequences without a gap in a retained column."
        ),
    )

    parser.add_argument(
        "--min-occupancy",
        type=float,
        default=0.95,
        help="Minimum occupancy for conservation analysis (default: 0.95).",
    )
    parser.add_argument(
        "--min-major-frequency",
        type=float,
        default=0.95,
        help="Minimum major-base frequency for conservation (default: 0.95).",
    )

    parser.add_argument(
        "--varvamp-mode",
        choices=("single", "qpcr", "tiled"),
        default=None,
        help=(
            "VarVAMP mode: single (conventional PCR), qpcr, or tiled "
            "(full-genome amplicon sequencing)."
        ),
    )
    parser.add_argument("--varvamp-threshold", type=float, default=None)
    parser.add_argument("--primer-ambiguity", type=int, default=None)
    parser.add_argument("--probe-ambiguity", type=int, default=None)
    parser.add_argument("--skip-varvamp", action="store_true")

    # Optional mode-specific VarVAMP CLI parameters. When omitted, VarVAMP's
    # own documented defaults are left untouched.
    parser.add_argument("--single-opt-length", type=int, default=None)
    parser.add_argument("--single-max-length", type=int, default=None)
    parser.add_argument(
        "--single-report-n",
        type=str,
        default=None,
        help="VarVAMP single -n value; integer or 'inf'.",
    )
    parser.add_argument("--tiled-opt-length", type=int, default=None)
    parser.add_argument("--tiled-max-length", type=int, default=None)
    parser.add_argument("--tiled-overlap", type=int, default=None)
    parser.add_argument("--qpcr-test-n", type=int, default=None)
    parser.add_argument("--qpcr-deltag", type=float, default=None)
    parser.add_argument(
        "--varvamp-config",
        type=Path,
        default=None,
        help="Optional custom VarVAMP config file passed through VARVAMP_CONFIG.",
    )

    # Executable names/paths make the script easier to use in custom envs.
    parser.add_argument("--mafft-executable", default="mafft")
    parser.add_argument("--mars-executable", default="mars")
    parser.add_argument("--seqkit-executable", default="seqkit")
    parser.add_argument("--cdhit-executable", default="cd-hit-est")
    parser.add_argument("--trimal-executable", default="trimal")
    parser.add_argument("--varvamp-executable", default="varvamp")

    return parser.parse_args()


def ask_for_input_file() -> Path:
    """Interactively request the input FASTA file."""

    while True:
        try:
            raw_value = input("Name or path of the input FASTA file: ").strip()
        except EOFError as error:
            raise RuntimeError(
                "No interactive input is available. Use --input."
            ) from error

        raw_value = raw_value.strip("'\"")
        if not raw_value:
            print("Please enter a filename or file path.")
            continue

        candidate = Path(raw_value).expanduser()
        possible_paths = [candidate]
        if not candidate.is_absolute() and candidate.parent == Path("."):
            possible_paths.append(Path("data") / candidate)

        for possible_path in possible_paths:
            resolved = possible_path.resolve()
            if resolved.is_file():
                return resolved

        print(
            "File not found. Checked paths: "
            + ", ".join(str(path.resolve()) for path in possible_paths)
        )


def ask_required_choice(title: str, choices: list[tuple[str, str]]) -> str:
    """Display a numbered menu and require an explicit selection."""

    print(f"\n{title}")
    print("=" * len(title))
    for index, (_, description) in enumerate(choices, start=1):
        print(f"{index}. {description}")

    while True:
        try:
            raw_value = input("Select an option by number: ").strip()
        except EOFError as error:
            raise RuntimeError(
                "A required interactive choice is missing."
            ) from error

        if not raw_value:
            print("A selection is required.")
            continue

        try:
            selected_index = int(raw_value)
        except ValueError:
            print("Please enter one of the displayed numbers.")
            continue

        if 1 <= selected_index <= len(choices):
            return choices[selected_index - 1][0]

        print("The selected number is outside the available range.")


def ask_yes_no(question: str) -> bool:
    """Ask an explicit yes/no question without assuming a default."""

    while True:
        try:
            answer = input(f"{question} [y/n]: ").strip().lower()
        except EOFError as error:
            raise RuntimeError("A required yes/no answer is missing.") from error

        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        print("Please answer y or n.")


def ask_required_float(
    label: str,
    option_name: str,
    *,
    minimum: float,
    maximum: float,
) -> float:
    """Request a required floating-point value in a closed interval."""

    while True:
        try:
            raw_value = input(
                f"{label} ({option_name}) [{minimum}-{maximum}]: "
            ).strip()
        except EOFError as error:
            raise RuntimeError(f"Missing required parameter {option_name}.") from error

        if not raw_value:
            print("A value is required; pressing Enter alone is not accepted.")
            continue

        try:
            value = float(raw_value)
        except ValueError:
            print("Please enter a numeric value.")
            continue

        if not minimum <= value <= maximum:
            print(f"The value must be between {minimum} and {maximum}.")
            continue

        return value


def ask_required_non_negative_integer(label: str, option_name: str) -> int:
    """Request a required non-negative integer."""

    while True:
        try:
            raw_value = input(f"{label} ({option_name}): ").strip()
        except EOFError as error:
            raise RuntimeError(f"Missing required parameter {option_name}.") from error

        if not raw_value:
            print("A value is required; pressing Enter alone is not accepted.")
            continue

        try:
            value = int(raw_value)
        except ValueError:
            print("Please enter an integer such as 0, 1 or 2.")
            continue

        if value < 0:
            print("The value cannot be negative.")
            continue

        return value


def ask_positive_integer_with_default(label: str, default: int) -> int:
    """Request a positive integer while explicitly displaying a default."""

    while True:
        raw = input(f"{label} [default {default}]: ").strip()
        if not raw:
            return default
        try:
            value = int(raw)
        except ValueError:
            print("Please enter a positive integer.")
            continue
        if value <= 0:
            print("The value must be greater than zero.")
            continue
        return value


def ask_float_with_default(label: str, default: float) -> float:
    """Request a float while explicitly displaying a default."""

    while True:
        raw = input(f"{label} [default {default}]: ").strip()
        if not raw:
            return default
        try:
            return float(raw)
        except ValueError:
            print("Please enter a numeric value.")


# ---------------------------------------------------------------------------
# Validation and generic execution helpers
# ---------------------------------------------------------------------------


def sanitize_project_name(value: str) -> str:
    """Convert a project name into a safe file/directory identifier."""

    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    if not cleaned:
        raise ValueError("The project name is empty or invalid.")
    return cleaned


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


def resolve_executable(executable: str) -> str | None:
    """Return a usable executable path or None."""

    path = Path(executable).expanduser()
    if path.parent != Path(".") or path.is_absolute():
        return str(path.resolve()) if path.is_file() else None
    return shutil.which(executable)


def check_required_tools(tools: Iterable[tuple[str, str]]) -> None:
    """Check that selected external tools are available."""

    missing: list[str] = []
    for label, executable in tools:
        if resolve_executable(executable) is None:
            missing.append(f"{label} ({executable})")

    if missing:
        raise RuntimeError(
            "Required program(s) not found: " + ", ".join(missing) + "."
        )


def run_command(
    command: list[str],
    *,
    stdout_file: Path | None = None,
    env: dict[str, str] | None = None,
) -> None:
    """Run an external command, log it, and stop on failure."""
    command_text = " ".join(command)
    if stdout_file is not None:
        command_text += f" > {stdout_file}"

    log_line("COMMAND: " + command_text)

    if VERBOSE:
        print("\nCommand:", command_text)

    try:
        if stdout_file is None:
            subprocess.run(command, check=True, env=env)
        else:
            stdout_file.parent.mkdir(parents=True, exist_ok=True)
            with stdout_file.open("w", encoding="utf-8") as handle:
                subprocess.run(
                    command,
                    check=True,
                    stdout=handle,
                    env=env,
                )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"Command failed with exit code {error.returncode}: "
            + command_text
        ) from error


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


# ---------------------------------------------------------------------------
# Orientation and circular-genome preprocessing
# ---------------------------------------------------------------------------


def choose_orientation_strategy() -> str:
    return ask_required_choice(
        "Sequence orientation checkpoint",
        [
            (
                "keep",
                "Keep the original orientation — use when all sequences are "
                "known to be in the same direction.",
            ),
            (
                "adjust",
                "Check/correct orientation with MAFFT --adjustdirection "
                "(faster; suitable for most datasets).",
            ),
            (
                "accurate",
                "Check/correct orientation with MAFFT --adjustdirectionaccurately "
                "(slower DP-based method for more divergent sequences).",
            ),
        ],
    )


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

    if VERBOSE:
        print(f"Oriented FASTA: {oriented_file.resolve()}")
    print(f"Reverse-complemented sequences: {len(reversed_ids)}")
    return oriented_file, len(reversed_ids)


def choose_sequence_topology() -> str:
    return ask_required_choice(
        "Sequence topology checkpoint",
        [
            (
                "linear",
                "Linear sequences — do not perform circular start-position normalization.",
            ),
            (
                "circular",
                "Complete circular sequences — normalize cyclic starting positions with MARS.",
            ),
        ],
    )


def run_mars_rotation(
    input_file: Path,
    workdir: Path,
    project_name: str,
    mars_executable: str,
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
        ]
    )

    validate_nucleotide_fasta(rotated_file)
    if count_fasta_sequences(rotated_file) != count_fasta_sequences(input_file):
        raise RuntimeError(
            "MARS output contains a different number of sequences than its input."
        )

    if VERBOSE:
        print(f"MARS-rotated FASTA: {rotated_file.resolve()}")
    return rotated_file


# ---------------------------------------------------------------------------
# Redundancy reduction
# ---------------------------------------------------------------------------


def choose_redundancy_strategy() -> str:
    return ask_required_choice(
        "Sequence redundancy checkpoint",
        [
            ("none", "Keep all sequences — no redundancy reduction."),
            (
                "seqkit",
                "Remove exact duplicate sequences (100% identical) with SeqKit.",
            ),
            (
                "cdhit",
                "Cluster similar sequences with CD-HIT-EST using a user-defined "
                "identity threshold.",
            ),
        ],
    )


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
    if 0.75 <= identity < 0.80:
        return 4
    raise ValueError("CD-HIT-EST identity must be between 0.75 and 1.0.")


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
) -> None:
    """Run the final multiple sequence alignment."""

    mafft_arguments = list(MAFFT_STRATEGIES[strategy]["arguments"])
    run_command(
        [
            mafft_executable,
            *mafft_arguments,
            "--thread",
            str(threads),
            str(input_file),
        ],
        stdout_file=output_file,
    )


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

    run_command(command)
    read_alignment(output_alignment)


# ---------------------------------------------------------------------------
# Conservation analysis
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# VarVAMP design modes and result handling
# ---------------------------------------------------------------------------


def choose_varvamp_mode() -> str:
    return ask_required_choice(
        "Assay design checkpoint",
        [
            (
                "single",
                "Conventional PCR / individual amplicons — VarVAMP SINGLE mode.",
            ),
            (
                "qpcr",
                "qPCR primers with an optimized internal probe — VarVAMP QPCR mode.",
            ),
            (
                "tiled",
                "Overlapping tiled amplicons for full-genome sequencing — "
                "VarVAMP TILED mode.",
            ),
        ],
    )


def configure_interactive_varvamp_advanced(args: argparse.Namespace, mode: str) -> None:
    """Optionally collect documented mode-specific VarVAMP parameters."""

    if not ask_yes_no("Configure additional mode-specific VarVAMP parameters?"):
        return

    if mode == "single":
        args.single_opt_length = ask_positive_integer_with_default(
            "Optimal amplicon length (-ol)", 1000
        )
        args.single_max_length = ask_positive_integer_with_default(
            "Maximum amplicon length (-ml)", 1500
        )
        while True:
            raw = input("Number of top hits to report (-n) [default inf]: ").strip()
            if not raw:
                args.single_report_n = "inf"
                break
            if raw.lower() == "inf":
                args.single_report_n = "inf"
                break
            try:
                value = int(raw)
            except ValueError:
                print("Enter a positive integer or 'inf'.")
                continue
            if value <= 0:
                print("Enter a positive integer or 'inf'.")
                continue
            args.single_report_n = str(value)
            break

    elif mode == "tiled":
        args.tiled_opt_length = ask_positive_integer_with_default(
            "Optimal amplicon length (-ol)", 1000
        )
        args.tiled_max_length = ask_positive_integer_with_default(
            "Maximum amplicon length (-ml)", 1500
        )
        args.tiled_overlap = ask_positive_integer_with_default(
            "Minimum amplicon-insert overlap (-o)", 100
        )

    elif mode == "qpcr":
        args.qpcr_test_n = ask_positive_integer_with_default(
            "Number of top qPCR amplicons tested for secondary structures (-n)",
            50,
        )
        args.qpcr_deltag = ask_float_with_default(
            "Minimum deltaG cutoff (-d)", -3.0
        )


def build_varvamp_command(
    args: argparse.Namespace,
    mode: str,
    alignment_file: Path,
    output_dir: Path,
) -> list[str]:
    """Build a VarVAMP command using documented mode-specific arguments."""

    command = [
        args.varvamp_executable,
        mode,
        "-t",
        str(args.varvamp_threshold),
        "-a",
        str(args.primer_ambiguity),
        "-th",
        str(args.threads),
    ]

    if mode == "qpcr":
        command.extend(["-pa", str(args.probe_ambiguity)])
        if args.qpcr_test_n is not None:
            command.extend(["-n", str(args.qpcr_test_n)])
        if args.qpcr_deltag is not None:
            command.extend(["-d", str(args.qpcr_deltag)])

    elif mode == "single":
        if args.single_opt_length is not None:
            command.extend(["-ol", str(args.single_opt_length)])
        if args.single_max_length is not None:
            command.extend(["-ml", str(args.single_max_length)])
        if args.single_report_n is not None:
            command.extend(["-n", str(args.single_report_n)])

    elif mode == "tiled":
        if args.tiled_opt_length is not None:
            command.extend(["-ol", str(args.tiled_opt_length)])
        if args.tiled_max_length is not None:
            command.extend(["-ml", str(args.tiled_max_length)])
        if args.tiled_overlap is not None:
            command.extend(["-o", str(args.tiled_overlap)])

    command.extend([str(alignment_file), str(output_dir)])
    return command


def validate_varvamp_mode_parameters(args: argparse.Namespace, mode: str) -> None:
    """Validate common and mode-specific VarVAMP parameters."""

    if args.varvamp_threshold is None:
        raise RuntimeError("VarVAMP consensus threshold (-t) is required.")
    validate_fraction(args.varvamp_threshold, "--varvamp-threshold")

    if args.primer_ambiguity is None:
        raise RuntimeError("Primer ambiguity (-a) is required.")
    if args.primer_ambiguity < 0:
        raise ValueError("--primer-ambiguity cannot be negative.")

    if mode == "qpcr":
        if args.probe_ambiguity is None:
            raise RuntimeError("Probe ambiguity (-pa) is required for qPCR mode.")
        if args.probe_ambiguity < 0:
            raise ValueError("--probe-ambiguity cannot be negative.")

    if mode == "single":
        if (
            args.single_opt_length is not None
            and args.single_max_length is not None
            and args.single_max_length < args.single_opt_length
        ):
            raise ValueError(
                "--single-max-length must be >= --single-opt-length."
            )
        if args.single_report_n is not None:
            if args.single_report_n.lower() != "inf":
                try:
                    report_n = int(args.single_report_n)
                except ValueError as error:
                    raise ValueError(
                        "--single-report-n must be a positive integer or 'inf'."
                    ) from error
                if report_n <= 0:
                    raise ValueError(
                        "--single-report-n must be a positive integer or 'inf'."
                    )

    if mode == "tiled":
        if (
            args.tiled_opt_length is not None
            and args.tiled_max_length is not None
            and args.tiled_max_length < args.tiled_opt_length
        ):
            raise ValueError(
                "--tiled-max-length must be >= --tiled-opt-length."
            )
        if args.tiled_overlap is not None and args.tiled_overlap < 0:
            raise ValueError("--tiled-overlap cannot be negative.")

    if args.qpcr_test_n is not None and args.qpcr_test_n <= 0:
        raise ValueError("--qpcr-test-n must be greater than zero.")


def copy_varvamp_tree(varvamp_dir: Path, results_dir: Path, mode: str) -> Path:
    """Copy the complete VarVAMP result tree so no mode-specific output is lost."""

    destination = results_dir / f"varvamp_{mode}"
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(varvamp_dir, destination)
    return destination


def plot_bed_segments(
    bed_file: Path,
    output_file: Path,
    title: str,
    y_label: str,
    line_width: float,
) -> None:
    """Create a simple overview plot from a BED file."""

    if not bed_file.is_file():
        return

    try:
        data = pd.read_csv(bed_file, sep="\t", header=None, comment="#")
    except pd.errors.EmptyDataError:
        return

    if data.empty or data.shape[1] < 3:
        return

    plt.figure(figsize=(12, 4))
    for index, row in data.iterrows():
        plt.plot(
            [row.iloc[1], row.iloc[2]],
            [index, index],
            linewidth=line_width,
        )
    plt.xlabel("Alignment / consensus position")
    plt.ylabel(y_label)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def save_pdf_first_page_as_png(
    pdf_file: Path,
    output_file: Path,
    *,
    zoom: int = 4,
) -> None:
    """Convert the first page of a PDF to PNG."""

    if not pdf_file.is_file():
        return
    with fitz.open(pdf_file) as pdf:
        if len(pdf) == 0:
            return
        page = pdf[0]
        pixmap = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
    image = Image.open(BytesIO(pixmap.tobytes("png"))).copy()
    image.save(output_file)


def print_varvamp_tables(varvamp_results_dir: Path, mode: str) -> None:
    """
    Summarize principal VarVAMP tables.

    Normal mode shows row counts only; --verbose prints the complete tables.
    """
    if mode == "qpcr":
        tables = [
            ("qpcr_primers.tsv", "qPCR primers/probes", "\t"),
            ("qpcr_design.tsv", "qPCR designs", "\t"),
        ]
    else:
        tables = [
            ("primer.tsv", "Primers", "\t"),
            (
                "primer_to_amplicon_assignments.tabular",
                "Primer-to-amplicon assignments",
                "\t",
            ),
        ]

    section("VarVAMP output summary")

    for filename, label, separator in tables:
        path = varvamp_results_dir / filename
        if not path.is_file():
            continue

        try:
            table = pd.read_csv(path, sep=separator)
        except (pd.errors.ParserError, pd.errors.EmptyDataError):
            print(f"{label:<32}: present")
            continue

        print(f"{label:<32}: {len(table)} rows")

        if VERBOSE:
            print()
            print(table.to_string(index=False))
            print()

    print(f"Result directory               : {compact_path(varvamp_results_dir)}")


def create_varvamp_visual_summaries(varvamp_results_dir: Path) -> None:
    """Create optional PNG summaries from standard VarVAMP outputs."""

    plot_bed_segments(
        varvamp_results_dir / "amplicons.bed",
        varvamp_results_dir / "amplicons_overview.png",
        "Amplicon positions",
        "Amplicon",
        6,
    )
    plot_bed_segments(
        varvamp_results_dir / "primers.bed",
        varvamp_results_dir / "primers_overview.png",
        "Primer / probe positions",
        "Oligonucleotide",
        4,
    )
    save_pdf_first_page_as_png(
        varvamp_results_dir / "amplicon_plot.pdf",
        varvamp_results_dir / "amplicon_plot.png",
    )
    save_pdf_first_page_as_png(
        varvamp_results_dir / "per_base_mismatches.pdf",
        varvamp_results_dir / "per_base_mismatches.png",
    )


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------


def main() -> int:
    global VERBOSE, WORKFLOW_LOG

    args = parse_arguments()
    VERBOSE = bool(args.verbose)

    try:
        if args.threads <= 0:
            raise ValueError("--threads must be greater than zero.")

        validate_fraction(args.min_occupancy, "--min-occupancy")
        validate_fraction(args.min_major_frequency, "--min-major-frequency")

        if args.input is None:
            input_file = ask_for_input_file()
        else:
            input_file = args.input.expanduser().resolve()

        validate_nucleotide_fasta(input_file)

        project_name = sanitize_project_name(
            args.project_name if args.project_name else input_file.stem
        )
        workdir = (
            args.workdir.expanduser().resolve()
            if args.workdir is not None
            else (Path("work") / project_name).resolve()
        )
        results_dir = (
            args.results.expanduser().resolve()
            if args.results is not None
            else (Path("results") / project_name).resolve()
        )
        workdir.mkdir(parents=True, exist_ok=True)
        results_dir.mkdir(parents=True, exist_ok=True)

        WORKFLOW_LOG = results_dir / "workflow.log"
        WORKFLOW_LOG.write_text(
            "VarVAMP assay-design workflow log\n",
            encoding="utf-8",
        )

        print("╔══════════════════════════════════════════════╗")
        print("║       VarVAMP ASSAY DESIGN WORKFLOW         ║")
        print("╚══════════════════════════════════════════════╝")
        print(f"Project : {project_name}")
        print(f"Input   : {compact_path(input_file)}")
        print(f"Results : {compact_path(results_dir)}")
        if VERBOSE:
            print(f"Workdir : {workdir}")
            print(f"Log     : {WORKFLOW_LOG}")

        interactive = sys.stdin.isatty()

        # Interactive scientific choices are intentionally collected step by step.
        # Each tool is configured immediately before it is executed, instead of
        # asking for all workflow parameters at startup. Command-line values still
        # bypass the corresponding interactive question.

        # Manual trimAl parameters may be supplied in advance in non-interactive
        # mode, but interactive trimming remains deliberately deferred until after
        # alignment QC statistics have been displayed.
        if args.trimming == "manual" and args.trimal_gap_threshold is None and not interactive:
            raise RuntimeError(
                "--trimal-gap-threshold is required for --trimming manual."
            )

        # ------------------------------------------------------------------
        # 1. Orientation normalization
        # ------------------------------------------------------------------
        current_fasta = input_file
        initial_count = count_fasta_sequences(input_file)

        print("\n1. Sequence orientation checkpoint")
        print_fasta_statistics(
            calculate_fasta_statistics(current_fasta),
            title="Pre-orientation dataset statistics",
        )

        if args.orientation is None:
            if not interactive:
                raise RuntimeError("--orientation is required in non-interactive mode.")
            args.orientation = choose_orientation_strategy()

        if args.orientation not in {"keep", "adjust", "accurate"}:
            raise ValueError("Unsupported orientation strategy.")

        reversed_count = 0
        if args.orientation == "keep":
            print("Sequence orientation will be kept as provided.")
            print("Orientation changes:           0 (orientation check skipped)")
        else:
            check_required_tools([("MAFFT", args.mafft_executable)])
            print("Running sequence orientation normalization with MAFFT.")
            current_fasta, reversed_count = normalize_orientation_with_mafft(
                current_fasta,
                workdir,
                project_name,
                args.orientation,
                args.threads,
                args.mafft_executable,
            )
            print("\nOrientation result")
            print("==================")
            print(f"Sequences checked:             {initial_count}")
            print(
                f"Reverse-complemented:           {reversed_count} "
                f"({100.0 * reversed_count / initial_count if initial_count else 0.0:.2f} %)"
            )
            print(f"Orientation unchanged:         {initial_count - reversed_count}")

        print_fasta_statistics(
            calculate_fasta_statistics(current_fasta),
            title="Post-orientation dataset statistics",
        )

        # ------------------------------------------------------------------
        # 2. Circular start-position normalization
        # ------------------------------------------------------------------
        print("\n2. Sequence topology checkpoint")
        print_fasta_statistics(
            calculate_fasta_statistics(current_fasta),
            title="Pre-topology-normalization dataset statistics",
        )

        if args.topology is None:
            if not interactive:
                raise RuntimeError("--topology is required in non-interactive mode.")
            args.topology = choose_sequence_topology()

        if args.topology == "circular":
            check_required_tools([("MARS", args.mars_executable)])
            print("Running circular sequence rotation with MARS.")
            before_mars = current_fasta
            current_fasta = run_mars_rotation(
                current_fasta,
                workdir,
                project_name,
                args.mars_executable,
            )
            print_rotation_statistics(
                calculate_cyclic_rotation_statistics(before_mars, current_fasta)
            )
        elif args.topology == "linear":
            print("Circular start-position normalization skipped (linear sequences).")
            print("Cyclic start positions changed: 0 (MARS not applicable)")
        else:
            raise ValueError("--topology must be 'linear' or 'circular'.")

        print_fasta_statistics(
            calculate_fasta_statistics(current_fasta),
            title="Post-topology-normalization dataset statistics",
        )

        # ------------------------------------------------------------------
        # 3. Redundancy handling
        # ------------------------------------------------------------------
        print("\n3. Sequence redundancy checkpoint")
        print_fasta_statistics(
            calculate_fasta_statistics(current_fasta),
            title="Pre-redundancy dataset statistics",
        )
        if args.redundancy is None:
            if args.skip_cdhit:
                args.redundancy = "none"
            elif not interactive:
                raise RuntimeError("--redundancy is required in non-interactive mode.")
            else:
                args.redundancy = choose_redundancy_strategy()

        if args.redundancy == "cdhit":
            if args.identity is None:
                if not interactive:
                    raise RuntimeError("--identity is required for CD-HIT-EST.")
                args.identity = ask_required_float(
                    "CD-HIT-EST identity threshold",
                    "-c",
                    minimum=0.75,
                    maximum=1.0,
                )
            if not 0.75 <= args.identity <= 1.0:
                raise ValueError("--identity must be between 0.75 and 1.0.")

            if args.cdhit_mode is None:
                if not interactive:
                    raise RuntimeError("--cdhit-mode is required for CD-HIT-EST.")
                args.cdhit_mode = ask_required_choice(
                    "CD-HIT-EST clustering mode",
                    [
                        (
                            "fast",
                            "Fast mode (-g 0): assign to the first qualifying cluster.",
                        ),
                        (
                            "accurate",
                            "Accurate mode (-g 1): select the most similar qualifying representative.",
                        ),
                    ],
                )

            if args.cdhit_strand is None:
                if not interactive:
                    raise RuntimeError("--cdhit-strand is required for CD-HIT-EST.")
                args.cdhit_strand = ask_required_choice(
                    "CD-HIT-EST strand comparison",
                    [
                        ("both", "Compare both +/+ and +/- orientations (-r 1)."),
                        ("same", "Compare only the same orientation, +/+ (-r 0)."),
                    ],
                )

        redundancy_output = workdir / f"{project_name}_redundancy_filtered.fasta"
        before_redundancy = count_fasta_sequences(current_fasta)

        if args.redundancy == "none":
            shutil.copy2(current_fasta, redundancy_output)
            print("No redundancy reduction selected.")

        elif args.redundancy == "seqkit":
            check_required_tools([("SeqKit", args.seqkit_executable)])
            print("Running exact duplicate removal with SeqKit.")
            run_seqkit_deduplication(
                current_fasta,
                redundancy_output,
                args.seqkit_executable,
            )

        elif args.redundancy == "cdhit":
            check_required_tools([("CD-HIT-EST", args.cdhit_executable)])
            print("Running sequence clustering with CD-HIT-EST.")
            word_size = run_cdhit_clustering(
                current_fasta,
                redundancy_output,
                identity=args.identity,
                mode=args.cdhit_mode,
                strand=args.cdhit_strand,
                threads=args.threads,
                cdhit_executable=args.cdhit_executable,
            )
            print(f"Automatically selected CD-HIT-EST word size: {word_size}")

        else:
            raise RuntimeError(f"Unknown redundancy strategy: {args.redundancy}")

        after_redundancy = count_fasta_sequences(redundancy_output)
        removed_redundancy = before_redundancy - after_redundancy
        reduction = (
            100.0 * removed_redundancy / before_redundancy
            if before_redundancy
            else 0.0
        )
        retained_percentage = (
            100.0 * after_redundancy / before_redundancy
            if before_redundancy
            else 0.0
        )
        print("\nRedundancy reduction statistics")
        print("===============================")
        print(f"Strategy:                       {args.redundancy}")
        print(f"Sequences before:               {before_redundancy}")
        print(f"Sequences removed:              {removed_redundancy} ({reduction:.2f} %)")
        print(f"Sequences retained:             {after_redundancy} ({retained_percentage:.2f} %)")
        if args.redundancy == "cdhit":
            print(f"CD-HIT-EST cluster representatives: {after_redundancy}")

        print_fasta_statistics(
            calculate_fasta_statistics(redundancy_output),
            title="Post-redundancy dataset statistics",
        )

        # ------------------------------------------------------------------
        # 4. Final MAFFT multiple sequence alignment
        # ------------------------------------------------------------------
        print("\n4. Final MAFFT alignment checkpoint")
        if args.mafft_strategy is None:
            if not interactive:
                raise RuntimeError(
                    "--mafft-strategy is required in non-interactive mode."
                )
            args.mafft_strategy = choose_mafft_strategy()

        if args.mafft_strategy not in MAFFT_STRATEGIES:
            raise ValueError(f"Unsupported MAFFT strategy: {args.mafft_strategy}")

        if interactive:
            warn_about_mafft_strategy(args.mafft_strategy, after_redundancy)

        check_required_tools([("MAFFT", args.mafft_executable)])
        final_alignment = workdir / f"{project_name}_alignment.fasta"
        print(f"Running final MAFFT alignment: {args.mafft_strategy}")
        run_final_mafft(
            redundancy_output,
            final_alignment,
            args.mafft_strategy,
            args.threads,
            args.mafft_executable,
        )
        read_alignment(final_alignment)
        shutil.copy2(
            final_alignment,
            results_dir / f"{project_name}_alignment.fasta",
        )

        # ------------------------------------------------------------------
        # 5. Alignment QC and trimming checkpoint
        # ------------------------------------------------------------------
        print("\n5. Alignment QC and trimming checkpoint")
        pretrim_stats, pretrim_sequences = calculate_alignment_statistics(
            final_alignment
        )
        print_alignment_statistics(pretrim_stats)
        save_alignment_statistics(
            pretrim_stats,
            pretrim_sequences,
            results_dir,
            project_name,
            "pretrim",
        )

        if args.trimming is None:
            # Only possible in interactive mode; non-interactive mode was checked.
            args.trimming = choose_trimming_strategy()

        if args.trimming == "manual" and args.trimal_gap_threshold is None:
            print("\nManual trimAl gap threshold")
            print("============================")
            print(
                "The -gt value is the minimum fraction of sequences without a gap "
                "in a retained column. Examples: 0.90 = >=90% occupancy; "
                "1.00 = no gaps allowed; 0.01 is extremely permissive."
            )
            args.trimal_gap_threshold = ask_required_float(
                "trimAl gap threshold",
                "-gt",
                minimum=0.0,
                maximum=1.0,
            )
            if args.trimal_gap_threshold < 0.50:
                print(
                    "Warning: this is a permissive threshold; many highly gapped "
                    "columns may be retained."
                )
                if not ask_yes_no(
                    f"Continue with -gt {args.trimal_gap_threshold}?"
                ):
                    raise RuntimeError("Manual trimAl threshold rejected by user.")

        downstream_alignment = final_alignment

        if args.trimming != "none":
            check_required_tools([("trimAl", args.trimal_executable)])
            trimmed_alignment = workdir / f"{project_name}_alignment_trimmed.fasta"
            print(f"\nRunning trimAl strategy: {args.trimming}")
            run_trimal(
                final_alignment,
                trimmed_alignment,
                args.trimming,
                args.trimal_executable,
                args.trimal_gap_threshold,
            )
            downstream_alignment = trimmed_alignment

            shutil.copy2(
                trimmed_alignment,
                results_dir / f"{project_name}_alignment_trimmed.fasta",
            )

            posttrim_stats, posttrim_sequences = calculate_alignment_statistics(
                trimmed_alignment
            )
            print_alignment_statistics(
                posttrim_stats,
                title="Post-trimming alignment statistics",
            )
            save_alignment_statistics(
                posttrim_stats,
                posttrim_sequences,
                results_dir,
                project_name,
                "posttrim",
            )

            original_length = int(pretrim_stats["alignment_length"])
            trimmed_length = int(posttrim_stats["alignment_length"])
            print(
                f"Columns removed by trimming: {original_length - trimmed_length} "
                f"({100.0 * (original_length - trimmed_length) / original_length:.2f} %)"
            )
        else:
            print("Keeping the complete final MAFFT alignment.")

        # ------------------------------------------------------------------
        # 6. Conservation analysis
        # ------------------------------------------------------------------
        print("\n6. Conservation analysis")
        perform_conservation_analysis(
            downstream_alignment,
            results_dir,
            project_name,
            args.min_occupancy,
            args.min_major_frequency,
        )

        # ------------------------------------------------------------------
        # 7. VarVAMP assay design
        # ------------------------------------------------------------------
        varvamp_results: Path | None = None

        if not args.skip_varvamp:
            print("\n7. VarVAMP assay design checkpoint")

            if args.varvamp_mode is None:
                if not interactive:
                    raise RuntimeError(
                        "--varvamp-mode is required in non-interactive mode."
                    )
                args.varvamp_mode = choose_varvamp_mode()

            if args.varvamp_threshold is None:
                if not interactive:
                    raise RuntimeError("--varvamp-threshold is required.")
                args.varvamp_threshold = ask_required_float(
                    "VarVAMP consensus threshold",
                    "-t",
                    minimum=0.01,
                    maximum=1.0,
                )

            if args.primer_ambiguity is None:
                if not interactive:
                    raise RuntimeError("--primer-ambiguity is required.")
                args.primer_ambiguity = ask_required_non_negative_integer(
                    "Maximum ambiguous bases in each primer",
                    "-a",
                )

            if args.varvamp_mode == "qpcr" and args.probe_ambiguity is None:
                if not interactive:
                    raise RuntimeError(
                        "--probe-ambiguity is required for qPCR mode."
                    )
                args.probe_ambiguity = ask_required_non_negative_integer(
                    "Maximum ambiguous bases in the probe",
                    "-pa",
                )

            if interactive:
                configure_interactive_varvamp_advanced(args, args.varvamp_mode)

            validate_varvamp_mode_parameters(args, args.varvamp_mode)
            check_required_tools([("VarVAMP", args.varvamp_executable)])
            print(f"Running assay design with VarVAMP ({args.varvamp_mode}).")

            varvamp_workdir = workdir / f"{project_name}_varvamp_{args.varvamp_mode}"
            if varvamp_workdir.exists():
                shutil.rmtree(varvamp_workdir)

            command = build_varvamp_command(
                args,
                args.varvamp_mode,
                downstream_alignment,
                varvamp_workdir,
            )

            command_env = os.environ.copy()
            if args.varvamp_config is not None:
                config_path = args.varvamp_config.expanduser().resolve()
                if not config_path.is_file():
                    raise FileNotFoundError(
                        f"VarVAMP config file not found: {config_path}"
                    )
                command_env["VARVAMP_CONFIG"] = str(config_path)

            run_command(command, env=command_env)
            varvamp_results = copy_varvamp_tree(
                varvamp_workdir,
                results_dir,
                args.varvamp_mode,
            )
            print_varvamp_tables(varvamp_results, args.varvamp_mode)
            create_varvamp_visual_summaries(varvamp_results)
        else:
            print("\n7. VarVAMP assay design skipped")

        manifest_file = results_dir / "assay_design_manifest.json"
        summary_file = results_dir / "workflow_summary.txt"

        write_project_manifest(
            manifest_file,
            project_name=project_name,
            input_file=input_file,
            workdir=workdir,
            results_dir=results_dir,
            downstream_alignment=downstream_alignment,
            initial_count=initial_count,
            after_redundancy=after_redundancy,
            args=args,
            varvamp_results=varvamp_results,
        )

        write_workflow_summary(
            summary_file,
            project_name=project_name,
            input_file=input_file,
            downstream_alignment=downstream_alignment,
            results_dir=results_dir,
            initial_count=initial_count,
            after_redundancy=after_redundancy,
            args=args,
            varvamp_results=varvamp_results,
        )

        # Put a copy of the design manifest inside the selected VarVAMP
        # result directory so validation can identify the exact parent project.
        if varvamp_results is not None:
            shutil.copy2(
                manifest_file,
                varvamp_results / "assay_design_manifest.json",
            )

        section("Workflow completed")
        print(f"Initial sequences        : {initial_count}")
        print(f"Final MAFFT sequences    : {after_redundancy}")
        print(f"Downstream alignment     : {compact_path(downstream_alignment)}")
        if varvamp_results is not None:
            print(f"VarVAMP mode             : {args.varvamp_mode.upper()}")
            print(f"VarVAMP results          : {compact_path(varvamp_results)}")
        else:
            print("VarVAMP                  : skipped")
        print(f"Project manifest         : {compact_path(manifest_file)}")
        print(f"Workflow summary         : {compact_path(summary_file)}")
        print(f"Workflow log             : {compact_path(WORKFLOW_LOG)}")
        return 0

    except (
        FileNotFoundError,
        ValueError,
        RuntimeError,
        pd.errors.ParserError,
    ) as error:
        log_line("ERROR: " + str(error))
        print(f"\nERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())