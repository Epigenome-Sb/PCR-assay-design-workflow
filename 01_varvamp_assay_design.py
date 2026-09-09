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
All commands are recorded in the run-specific workflow log.

Default storage architecture
----------------------------
Input design FASTA:
    data/design/

Temporary/reconstructible files:
    work/<project>/design/<run_id>/

Permanent scientific results:
    results/<project>/design/<run_id>/

Each run is isolated by a unique run_id so previous analyses are not
overwritten. The run contains preprocessing, alignment, conservation,
VarVAMP, configuration, manifest, summary and log outputs.
"""

from __future__ import annotations

import argparse
import ast
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


def available_cpu_threads() -> int:
    """Return the logical CPU threads available to this workflow process."""

    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


def automatic_thread_count(cpu_threads: int) -> int:
    """Keep two logical CPU threads free when possible."""

    if cpu_threads <= 2:
        return 1
    return max(1, cpu_threads - 2)


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


def portable_manifest_path(path: Path | str) -> str:
    """
    Serialize a path for machine-readable manifests.

    Paths located inside the repository/current working directory are stored
    relative to that directory so the same manifest can be reused after moving
    the project, including between native execution and Docker (/app).

    Truly external paths remain absolute because there is no safe repository-
    relative representation for them.
    """
    resolved = Path(path).expanduser().resolve()
    repository_root = Path.cwd().resolve()

    try:
        return resolved.relative_to(repository_root).as_posix()
    except ValueError:
        return str(resolved)


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
    run_id: str,
    input_file: Path,
    workdir: Path,
    results_dir: Path,
    downstream_alignment: Path,
    initial_count: int,
    after_redundancy: int,
    args: argparse.Namespace,
    varvamp_results: Path | None,
    varvamp_attempts: list[dict[str, object]] | None = None,
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
        "schema_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "workflow_stage": "design",
        "path_base": "repository_root",
        "path_policy": (
            "repository-relative when inside the project; "
            "absolute only for external paths"
        ),
        "project_name": project_name,
        "run_id": run_id,
        "input_role": "design_database",
        "input_fasta": portable_manifest_path(input_file),
        "input_sha256": sha256_file(input_file),
        "workdir": portable_manifest_path(workdir),
        "results_dir": portable_manifest_path(results_dir),
        "result_directories": {
            "preprocessing": portable_manifest_path(results_dir / "preprocessing"),
            "alignment": portable_manifest_path(results_dir / "alignment"),
            "conservation": portable_manifest_path(results_dir / "conservation"),
            "config": portable_manifest_path(results_dir / "config"),
            "varvamp": (
                None
                if varvamp_results is None
                else portable_manifest_path(varvamp_results)
            ),
        },
        "downstream_alignment": portable_manifest_path(downstream_alignment),
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
            "threads": int(args.threads),
        },
        "varvamp": {
            "executed": bool(varvamp_attempts)
            if varvamp_attempts is not None
            else not args.skip_varvamp,
            "mode": None if args.skip_varvamp else args.varvamp_mode,
            "result_dir": (
                None
                if varvamp_results is None
                else portable_manifest_path(varvamp_results)
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
                else portable_manifest_path(args.varvamp_config)
            ),
            "attempts": varvamp_attempts or [],
        },
        "tool_versions": selected_tools,
        "workflow_log": (
            None
            if WORKFLOW_LOG is None
            else portable_manifest_path(WORKFLOW_LOG)
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
    run_id: str,
    input_file: Path,
    downstream_alignment: Path,
    workdir: Path,
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
        f"Run ID: {run_id}",
        "Workflow stage: design",
        f"Design FASTA: {compact_path(input_file)}",
        f"Initial sequences: {initial_count}",
        f"Sequences entering final MAFFT: {after_redundancy}",
        f"Orientation: {args.orientation}",
        f"Topology: {args.topology}",
        f"Redundancy: {args.redundancy}",
        f"Final MAFFT strategy: {args.mafft_strategy}",
        f"Threads: {args.threads}",
        f"Trimming: {args.trimming}",
        f"Downstream alignment: {compact_path(downstream_alignment)}",
        "",
        "VarVAMP:",
        f"Executed: {'yes' if not args.skip_varvamp else 'no'}",
        f"Mode: {args.varvamp_mode if not args.skip_varvamp else 'N/A'}",
        (
            f"Result directory: {compact_path(varvamp_results)}"
            if varvamp_results is not None
            else "Result directory: N/A"
        ),
        "",
        f"Work directory: {compact_path(workdir)}",
        f"Results directory: {compact_path(results_dir)}",
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
    parser.add_argument(
        "--run-id",
        type=str,
        default=None,
        help=(
            "Optional run identifier. By default a timestamp such as "
            "20260906_164500 is generated automatically."
        ),
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=None,
        help=(
            "Optional override for this run work directory. Default: "
            "work/<project>/design/<run_id>."
        ),
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=None,
        help=(
            "Optional override for this run results directory. Default: "
            "results/<project>/design/<run_id>."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help=(
            "Number of CPU threads used by supported tools. If omitted, the "
            "workflow detects the logical CPU threads available to the process "
            "and keeps two threads free when possible."
        ),
    )
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
    parser.add_argument("--qpcr-deltag", type=int, default=None)
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
    """Interactively select the Workflow 01 design FASTA."""

    design_dir = Path("data") / "design"
    candidates: list[Path] = []
    if design_dir.is_dir():
        for pattern in ("*.fasta", "*.fa", "*.fna", "*.fas"):
            candidates.extend(design_dir.glob(pattern))
        candidates = sorted({path.resolve() for path in candidates if path.is_file()})

    print("\nDesign database")
    print("===============")
    print("Workflow 01 uses the FASTA stored in data/design/ for assay design.")

    if candidates:
        for index, path in enumerate(candidates, start=1):
            print(f"{index}. {compact_path(path)}")
        print(f"{len(candidates) + 1}. Enter another FASTA path")

        while True:
            try:
                raw = input("Select the design FASTA by number: ").strip()
            except EOFError as error:
                raise RuntimeError(
                    "No interactive input is available. Use --input."
                ) from error

            try:
                selected = int(raw)
            except ValueError:
                print("Please enter one of the displayed numbers.")
                continue

            if 1 <= selected <= len(candidates):
                return candidates[selected - 1]
            if selected == len(candidates) + 1:
                break
            print("Selection outside the available range.")

    while True:
        try:
            raw_value = input("Name or path of the design FASTA file: ").strip()
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
            possible_paths.append(design_dir / candidate)
            # Compatibility fallback for older repository layouts.
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


def ask_integer_with_default(
    label: str,
    default: int,
    *,
    minimum: int | None = None,
    maximum: int | None = None,
) -> int:
    """Request an integer while explicitly displaying and validating a default."""

    while True:
        raw = input(f"{label} [default {default}]: ").strip()
        if not raw:
            value = default
        else:
            try:
                value = int(raw)
            except ValueError:
                print("Please enter an integer (for example -3, 0, 50).")
                continue

        if minimum is not None and value < minimum:
            print(f"The value must be >= {minimum}.")
            continue
        if maximum is not None and value > maximum:
            print(f"The value must be <= {maximum}.")
            continue
        return value


def ask_non_negative_integer_with_default(label: str, default: int) -> int:
    """Request a non-negative integer with a default."""
    return ask_integer_with_default(label, default, minimum=0)


# ---------------------------------------------------------------------------
# Validation and generic execution helpers
# ---------------------------------------------------------------------------


def sanitize_project_name(value: str) -> str:
    """Convert a project name into a safe file/directory identifier."""

    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    if not cleaned:
        raise ValueError("The project name is empty or invalid.")
    return cleaned


def sanitize_run_id(value: str) -> str:
    """Validate/sanitize a run identifier used as a directory name."""

    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    if not cleaned:
        raise ValueError("The run ID is empty or invalid.")
    return cleaned


def generate_run_id() -> str:
    """Generate a human-readable local timestamp for one workflow execution."""

    return datetime.now().strftime("%Y%m%d_%H%M%S")


def resolve_run_directories(
    args: argparse.Namespace,
    project_name: str,
) -> tuple[str, Path, Path]:
    """
    Resolve isolated work/results directories for one assay-design run.

    The default layout is:
      work/<project>/design/<run_id>/
      results/<project>/design/<run_id>/

    When the automatically generated timestamp already exists, a numeric
    suffix is added so an earlier run is never overwritten. Explicit
    --workdir/--results values remain exact per-run overrides.
    """

    requested_run_id = (
        sanitize_run_id(args.run_id)
        if args.run_id is not None
        else generate_run_id()
    )

    if args.workdir is not None or args.results is not None:
        workdir = (
            args.workdir.expanduser().resolve()
            if args.workdir is not None
            else (Path("work") / project_name / "design" / requested_run_id).resolve()
        )
        results_dir = (
            args.results.expanduser().resolve()
            if args.results is not None
            else (Path("results") / project_name / "design" / requested_run_id).resolve()
        )
        return requested_run_id, workdir, results_dir

    candidate = requested_run_id
    suffix = 1
    while True:
        workdir = (Path("work") / project_name / "design" / candidate).resolve()
        results_dir = (Path("results") / project_name / "design" / candidate).resolve()
        if not workdir.exists() and not results_dir.exists():
            return candidate, workdir, results_dir
        suffix += 1
        candidate = f"{requested_run_id}_{suffix:02d}"


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file for provenance tracking."""

    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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

    if VERBOSE:
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


def ask_single_report_n(default: str = "inf") -> str:
    """Request VarVAMP SINGLE -n as a positive integer or 'inf'."""
    while True:
        raw = input(
            f"Number of top hits to report (-n) [default {default}]: "
        ).strip()
        if not raw:
            return default
        if raw.lower() == "inf":
            return "inf"
        try:
            value = int(raw)
        except ValueError:
            print("Enter a positive integer or 'inf'.")
            continue
        if value <= 0:
            print("Enter a positive integer or 'inf'.")
            continue
        return str(value)


def configure_interactive_varvamp_advanced(
    args: argparse.Namespace,
    mode: str,
    *,
    force: bool = False,
) -> None:
    """
    Collect documented mode-specific VarVAMP CLI parameters.

    Values are validated before VarVAMP is launched so type errors such as
    passing '-3.0' to an integer-only option cannot reach VarVAMP.
    """

    if not force and not ask_yes_no(
        "Configure additional mode-specific VarVAMP parameters?"
    ):
        return

    if mode == "single":
        while True:
            opt_length = ask_positive_integer_with_default(
                "Optimal amplicon length (-ol)",
                args.single_opt_length or 1000,
            )
            max_length = ask_positive_integer_with_default(
                "Maximum amplicon length (-ml)",
                args.single_max_length or 1500,
            )
            if max_length < opt_length:
                print(
                    "Maximum amplicon length (-ml) must be >= "
                    "optimal amplicon length (-ol)."
                )
                continue
            args.single_opt_length = opt_length
            args.single_max_length = max_length
            break

        args.single_report_n = ask_single_report_n(
            args.single_report_n or "inf"
        )

    elif mode == "tiled":
        while True:
            opt_length = ask_positive_integer_with_default(
                "Optimal amplicon length (-ol)",
                args.tiled_opt_length or 1000,
            )
            max_length = ask_positive_integer_with_default(
                "Maximum amplicon length (-ml)",
                args.tiled_max_length or 1500,
            )
            if max_length < opt_length:
                print(
                    "Maximum amplicon length (-ml) must be >= "
                    "optimal amplicon length (-ol)."
                )
                continue
            args.tiled_opt_length = opt_length
            args.tiled_max_length = max_length
            break

        args.tiled_overlap = ask_non_negative_integer_with_default(
            "Minimum amplicon-insert overlap (-o)",
            args.tiled_overlap if args.tiled_overlap is not None else 100,
        )

    elif mode == "qpcr":
        args.qpcr_test_n = ask_positive_integer_with_default(
            "Number of top qPCR amplicons tested for secondary structures (-n)",
            args.qpcr_test_n or 50,
        )
        args.qpcr_deltag = ask_integer_with_default(
            "Minimum deltaG cutoff (-d)",
            args.qpcr_deltag if args.qpcr_deltag is not None else -3,
        )


VARVAMP_ADVANCED_CONFIG_DEFAULTS: dict[str, dict[str, object]] = {
    "single": {
        "PRIMER_SIZES": (18, 25, 21),
        "PRIMER_TMP": (56, 63, 60),
        "PRIMER_GC_RANGE": (30, 75, 50),
    },
    "tiled": {
        "PRIMER_SIZES": (18, 25, 21),
        "PRIMER_TMP": (56, 63, 60),
        "PRIMER_GC_RANGE": (30, 75, 50),
    },
    "qpcr": {
        "PRIMER_SIZES": (18, 25, 21),
        "PRIMER_TMP": (56, 63, 60),
        "PRIMER_GC_RANGE": (30, 75, 50),
        "QPROBE_SIZES": (18, 30, 25),
        "QPROBE_TMP": (60, 72, 67),
        "QPROBE_GC_RANGE": (35, 85, 60),
        "QPRIMER_DIFF": 3,
        "QPROBE_TEMP_DIFF": (2, 12),
        "QPROBE_DISTANCE": (4, 25),
        "QAMPLICON_LENGTH": (70, 200),
        "QAMPLICON_GC": (40, 60),
    },
}


VARVAMP_ADVANCED_CONFIG_ORDER: dict[str, tuple[str, ...]] = {
    "single": (
        "PRIMER_SIZES",
        "PRIMER_TMP",
        "PRIMER_GC_RANGE",
    ),
    "tiled": (
        "PRIMER_SIZES",
        "PRIMER_TMP",
        "PRIMER_GC_RANGE",
    ),
    "qpcr": (
        "PRIMER_SIZES",
        "PRIMER_TMP",
        "PRIMER_GC_RANGE",
        "QPROBE_SIZES",
        "QPROBE_TMP",
        "QPROBE_GC_RANGE",
        "QPRIMER_DIFF",
        "QPROBE_TEMP_DIFF",
        "QPROBE_DISTANCE",
        "QAMPLICON_LENGTH",
        "QAMPLICON_GC",
    ),
}


def render_custom_varvamp_config(mode: str, values: dict[str, object]) -> str:
    """Render a mode-specific custom VarVAMP Python configuration."""
    lines = [
        f"# Custom VarVAMP {mode.upper()} configuration.",
        "# Created interactively by 01_varvamp_assay_design.py.",
        "# This file is used only when explicitly selected for this workflow run.",
        "",
    ]

    for name in VARVAMP_ADVANCED_CONFIG_ORDER[mode]:
        lines.append(f"{name} = {values[name]!r}")
        if name in {
            "PRIMER_GC_RANGE",
            "QPROBE_GC_RANGE",
            "QPRIMER_DIFF",
        }:
            lines.append("")
        elif name == "QPROBE_DISTANCE":
            lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def _literal_assignments_from_python(path: Path) -> dict[str, object]:
    """Safely read simple literal assignments from a Python config file."""
    source = path.read_text(encoding="utf-8")
    try:
        tree = ast.parse(source, filename=str(path))
    except SyntaxError as error:
        raise ValueError(
            f"Python syntax error in {path.name}: {error.msg} "
            f"(line {error.lineno})."
        ) from error

    values: dict[str, object] = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            continue
        try:
            values[target.id] = ast.literal_eval(node.value)
        except (ValueError, TypeError):
            # Unknown/non-literal VarVAMP settings are left for VarVAMP itself.
            continue
    return values


def _validate_min_max_opt(
    name: str,
    value: object,
    *,
    positive: bool = False,
    percentage: bool = False,
) -> None:
    """Validate a (min, max, opt) tuple used by VarVAMP."""
    if not isinstance(value, tuple) or len(value) != 3:
        raise ValueError(f"{name} must be a 3-item tuple: (min, max, opt).")
    if not all(isinstance(item, (int, float)) for item in value):
        raise ValueError(f"{name} values must be numeric.")

    minimum, maximum, optimum = value
    if minimum > optimum or optimum > maximum:
        raise ValueError(
            f"{name} must satisfy min <= opt <= max; found {value}."
        )
    if positive and minimum <= 0:
        raise ValueError(f"{name} values must be greater than zero.")
    if percentage and (minimum < 0 or maximum > 100):
        raise ValueError(f"{name} values must remain between 0 and 100.")


def _validate_min_max(
    name: str,
    value: object,
    *,
    non_negative: bool = False,
    positive: bool = False,
    percentage: bool = False,
) -> None:
    """Validate a (min, max) tuple used by VarVAMP."""
    if not isinstance(value, tuple) or len(value) != 2:
        raise ValueError(f"{name} must be a 2-item tuple: (min, max).")
    if not all(isinstance(item, (int, float)) for item in value):
        raise ValueError(f"{name} values must be numeric.")

    minimum, maximum = value
    if minimum > maximum:
        raise ValueError(f"{name} must satisfy min <= max; found {value}.")
    if non_negative and minimum < 0:
        raise ValueError(f"{name} values cannot be negative.")
    if positive and minimum <= 0:
        raise ValueError(f"{name} values must be greater than zero.")
    if percentage and (minimum < 0 or maximum > 100):
        raise ValueError(f"{name} values must remain between 0 and 100.")


def validate_custom_varvamp_config(path: Path, mode: str) -> None:
    """
    Validate Python syntax and the known core settings written by this workflow.

    Unknown VarVAMP settings are not rejected, which lets an experienced user
    add other valid VarVAMP configuration variables manually.
    """
    if not path.is_file():
        raise FileNotFoundError(f"VarVAMP config file not found: {path}")

    values = _literal_assignments_from_python(path)

    for name in ("PRIMER_SIZES", "PRIMER_TMP", "PRIMER_GC_RANGE"):
        if name not in values:
            continue
        if name == "PRIMER_SIZES":
            _validate_min_max_opt(name, values[name], positive=True)
        elif name == "PRIMER_GC_RANGE":
            _validate_min_max_opt(name, values[name], percentage=True)
        else:
            _validate_min_max_opt(name, values[name])

    if mode == "qpcr":
        if "QPROBE_SIZES" in values:
            _validate_min_max_opt(
                "QPROBE_SIZES",
                values["QPROBE_SIZES"],
                positive=True,
            )
        if "QPROBE_TMP" in values:
            _validate_min_max_opt("QPROBE_TMP", values["QPROBE_TMP"])
        if "QPROBE_GC_RANGE" in values:
            _validate_min_max_opt(
                "QPROBE_GC_RANGE",
                values["QPROBE_GC_RANGE"],
                percentage=True,
            )
        if "QPRIMER_DIFF" in values:
            value = values["QPRIMER_DIFF"]
            if not isinstance(value, (int, float)) or value < 0:
                raise ValueError("QPRIMER_DIFF must be a non-negative number.")
        if "QPROBE_TEMP_DIFF" in values:
            _validate_min_max(
                "QPROBE_TEMP_DIFF",
                values["QPROBE_TEMP_DIFF"],
                non_negative=True,
            )
        if "QPROBE_DISTANCE" in values:
            _validate_min_max(
                "QPROBE_DISTANCE",
                values["QPROBE_DISTANCE"],
                non_negative=True,
            )
        if "QAMPLICON_LENGTH" in values:
            _validate_min_max(
                "QAMPLICON_LENGTH",
                values["QAMPLICON_LENGTH"],
                positive=True,
            )
        if "QAMPLICON_GC" in values:
            _validate_min_max(
                "QAMPLICON_GC",
                values["QAMPLICON_GC"],
                percentage=True,
            )


def stage_varvamp_config_for_run(
    source: Path,
    config_dir: Path,
    mode: str,
    attempt_number: int,
) -> Path:
    """Copy an externally selected config into the permanent run results."""

    source = source.expanduser().resolve()
    validate_custom_varvamp_config(source, mode)
    config_dir.mkdir(parents=True, exist_ok=True)

    try:
        source.relative_to(config_dir.resolve())
        return source
    except ValueError:
        pass

    destination = (
        config_dir
        / f"varvamp_{mode}_selected_attempt_{attempt_number:02d}.py"
    )
    if destination.exists():
        raise RuntimeError(
            "Refusing to overwrite a staged VarVAMP configuration: "
            f"{destination}"
        )
    shutil.copy2(source, destination)
    validate_custom_varvamp_config(destination, mode)
    return destination.resolve()


def next_custom_varvamp_config_path(
    config_dir: Path,
    mode: str,
) -> Path:
    """Return a new attempt-specific config filename without overwriting older files."""
    attempt = 1
    while True:
        candidate = (
            config_dir
            / f"varvamp_{mode}_config_attempt_{attempt}.py"
        )
        if not candidate.exists():
            return candidate
        attempt += 1


def ask_integer_keep_default(
    label: str,
    default: int,
    *,
    minimum: int | None = None,
    maximum: int | None = None,
) -> int:
    """Ask for an integer, keeping the displayed value when Enter is pressed."""
    while True:
        raw = input(f"  {label} [{default}]: ").strip()
        if not raw:
            return default
        try:
            value = int(raw)
        except ValueError:
            print("  Please enter an integer, or press Enter to keep the current value.")
            continue
        if minimum is not None and value < minimum:
            print(f"  The value must be >= {minimum}.")
            continue
        if maximum is not None and value > maximum:
            print(f"  The value must be <= {maximum}.")
            continue
        return value


def edit_min_max_opt_value(
    name: str,
    current: tuple[int, int, int],
    *,
    positive: bool = False,
    percentage: bool = False,
) -> tuple[int, int, int]:
    """Edit a (minimum, maximum, optimum) setting directly in the terminal."""
    while True:
        minimum_limit = 1 if positive else (0 if percentage else None)
        maximum_limit = 100 if percentage else None

        print(f"\n{name} = {current}")
        minimum = ask_integer_keep_default(
            "Minimum",
            current[0],
            minimum=minimum_limit,
            maximum=maximum_limit,
        )
        maximum = ask_integer_keep_default(
            "Maximum",
            current[1],
            minimum=minimum_limit,
            maximum=maximum_limit,
        )
        optimum = ask_integer_keep_default(
            "Optimal",
            current[2],
            minimum=minimum_limit,
            maximum=maximum_limit,
        )
        candidate = (minimum, maximum, optimum)

        try:
            _validate_min_max_opt(
                name,
                candidate,
                positive=positive,
                percentage=percentage,
            )
        except ValueError as error:
            print(f"  Invalid value: {error}")
            print("  Please enter this parameter again.")
            continue
        return candidate


def edit_min_max_value(
    name: str,
    current: tuple[int, int],
    *,
    non_negative: bool = False,
    positive: bool = False,
    percentage: bool = False,
) -> tuple[int, int]:
    """Edit a (minimum, maximum) setting directly in the terminal."""
    while True:
        minimum_limit = None
        if positive:
            minimum_limit = 1
        elif non_negative or percentage:
            minimum_limit = 0
        maximum_limit = 100 if percentage else None

        print(f"\n{name} = {current}")
        minimum = ask_integer_keep_default(
            "Minimum",
            current[0],
            minimum=minimum_limit,
            maximum=maximum_limit,
        )
        maximum = ask_integer_keep_default(
            "Maximum",
            current[1],
            minimum=minimum_limit,
            maximum=maximum_limit,
        )
        candidate = (minimum, maximum)

        try:
            _validate_min_max(
                name,
                candidate,
                non_negative=non_negative,
                positive=positive,
                percentage=percentage,
            )
        except ValueError as error:
            print(f"  Invalid value: {error}")
            print("  Please enter this parameter again.")
            continue
        return candidate


def edit_varvamp_config_in_terminal(
    mode: str,
    initial_values: dict[str, object] | None = None,
) -> dict[str, object]:
    """Display default/current settings and let the user modify them in-terminal."""
    values = dict(VARVAMP_ADVANCED_CONFIG_DEFAULTS[mode])
    if initial_values:
        for name in VARVAMP_ADVANCED_CONFIG_ORDER[mode]:
            if name in initial_values:
                values[name] = initial_values[name]

    section(f"Edit advanced {mode.upper()} VarVAMP configuration")
    print("The current values are shown in brackets.")
    print("Press Enter without typing anything to keep a value unchanged.")
    print("The configuration is saved only after you confirm it at the end.")

    values["PRIMER_SIZES"] = edit_min_max_opt_value(
        "PRIMER_SIZES",
        values["PRIMER_SIZES"],
        positive=True,
    )
    values["PRIMER_TMP"] = edit_min_max_opt_value(
        "PRIMER_TMP",
        values["PRIMER_TMP"],
    )
    values["PRIMER_GC_RANGE"] = edit_min_max_opt_value(
        "PRIMER_GC_RANGE",
        values["PRIMER_GC_RANGE"],
        percentage=True,
    )

    if mode == "qpcr":
        values["QPROBE_SIZES"] = edit_min_max_opt_value(
            "QPROBE_SIZES",
            values["QPROBE_SIZES"],
            positive=True,
        )
        values["QPROBE_TMP"] = edit_min_max_opt_value(
            "QPROBE_TMP",
            values["QPROBE_TMP"],
        )
        values["QPROBE_GC_RANGE"] = edit_min_max_opt_value(
            "QPROBE_GC_RANGE",
            values["QPROBE_GC_RANGE"],
            percentage=True,
        )

        print(f"\nQPRIMER_DIFF = {values['QPRIMER_DIFF']}")
        values["QPRIMER_DIFF"] = ask_integer_keep_default(
            "Value",
            int(values["QPRIMER_DIFF"]),
            minimum=0,
        )
        values["QPROBE_TEMP_DIFF"] = edit_min_max_value(
            "QPROBE_TEMP_DIFF",
            values["QPROBE_TEMP_DIFF"],
            non_negative=True,
        )
        values["QPROBE_DISTANCE"] = edit_min_max_value(
            "QPROBE_DISTANCE",
            values["QPROBE_DISTANCE"],
            non_negative=True,
        )
        values["QAMPLICON_LENGTH"] = edit_min_max_value(
            "QAMPLICON_LENGTH",
            values["QAMPLICON_LENGTH"],
            positive=True,
        )
        values["QAMPLICON_GC"] = edit_min_max_value(
            "QAMPLICON_GC",
            values["QAMPLICON_GC"],
            percentage=True,
        )

    return values


def configure_custom_varvamp_file(
    args: argparse.Namespace,
    mode: str,
    config_dir: Path,
) -> bool:
    """
    Explicitly create/edit or select a custom VarVAMP Python configuration.

    New custom configurations are edited directly in the terminal. No
    previously created custom configuration is auto-detected or auto-loaded.
    """
    choice = ask_required_choice(
        f"Advanced {mode.upper()} VarVAMP configuration",
        [
            (
                "create",
                "Create/edit a new mode-specific configuration directly in the terminal.",
            ),
            (
                "existing",
                "Use an existing Python configuration file by entering its path.",
            ),
            (
                "return",
                "Return without changing the configuration.",
            ),
        ],
    )

    if choice == "return":
        return False

    if choice == "create":
        values = dict(VARVAMP_ADVANCED_CONFIG_DEFAULTS[mode])

        while True:
            section(f"Default {mode.upper()} VarVAMP configuration")
            print(render_custom_varvamp_config(mode, values).rstrip())
            print()
            print(
                "You can now modify these values directly in the terminal. "
                "Press Enter at any prompt to keep the displayed value."
            )

            values = edit_varvamp_config_in_terminal(mode, values)

            section("Configuration preview")
            preview = render_custom_varvamp_config(mode, values)
            print(preview.rstrip())

            if ask_yes_no("Save this configuration and use it for this VarVAMP attempt?"):
                config_dir.mkdir(parents=True, exist_ok=True)
                path = next_custom_varvamp_config_path(config_dir, mode)
                path.write_text(preview, encoding="utf-8")
                validate_custom_varvamp_config(path, mode)
                args.varvamp_config = path

                section("Custom VarVAMP configuration saved")
                print(f"File: {compact_path(path)}")
                print(
                    "This file will NOT be loaded automatically in a future workflow run."
                )
                return True

            if not ask_yes_no("Modify the configuration again?"):
                return False

    while True:
        raw = input("Path to the existing VarVAMP config.py: ").strip().strip("'\"")
        if not raw:
            print("A file path is required.")
            continue

        path = Path(raw).expanduser()
        if not path.is_absolute():
            path = (Path.cwd() / path).resolve()
        else:
            path = path.resolve()

        try:
            validate_custom_varvamp_config(path, mode)
        except (FileNotFoundError, ValueError) as error:
            print(f"Configuration error: {error}")
            if ask_yes_no("Try another configuration file?"):
                continue
            return False

        args.varvamp_config = path
        return True

def build_varvamp_command(
    args: argparse.Namespace,
    mode: str,
    alignment_file: Path,
    output_dir: Path,
) -> list[str]:
    """Build a VarVAMP command using validated mode-specific arguments."""

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
            command.extend(["-d", str(int(args.qpcr_deltag))])

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
    """Validate common and mode-specific VarVAMP parameters before execution."""

    if args.varvamp_threshold is None:
        raise RuntimeError("VarVAMP consensus threshold (-t) is required.")
    validate_fraction(args.varvamp_threshold, "--varvamp-threshold")

    if args.primer_ambiguity is None:
        raise RuntimeError("Primer ambiguity (-a) is required.")
    if not isinstance(args.primer_ambiguity, int) or args.primer_ambiguity < 0:
        raise ValueError("--primer-ambiguity must be a non-negative integer.")

    if mode == "qpcr":
        if args.probe_ambiguity is None:
            raise RuntimeError("Probe ambiguity (-pa) is required for qPCR mode.")
        if not isinstance(args.probe_ambiguity, int) or args.probe_ambiguity < 0:
            raise ValueError("--probe-ambiguity must be a non-negative integer.")

        if args.qpcr_test_n is not None:
            if not isinstance(args.qpcr_test_n, int) or args.qpcr_test_n <= 0:
                raise ValueError("--qpcr-test-n must be a positive integer.")

        if args.qpcr_deltag is not None and not isinstance(args.qpcr_deltag, int):
            raise ValueError(
                "--qpcr-deltag must be an integer, for example -3 (not -3.0)."
            )

    if mode == "single":
        for name, value in (
            ("--single-opt-length", args.single_opt_length),
            ("--single-max-length", args.single_max_length),
        ):
            if value is not None and (
                not isinstance(value, int) or value <= 0
            ):
                raise ValueError(f"{name} must be a positive integer.")

        if (
            args.single_opt_length is not None
            and args.single_max_length is not None
            and args.single_max_length < args.single_opt_length
        ):
            raise ValueError(
                "--single-max-length must be >= --single-opt-length."
            )

        if args.single_report_n is not None:
            if str(args.single_report_n).lower() != "inf":
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
        for name, value in (
            ("--tiled-opt-length", args.tiled_opt_length),
            ("--tiled-max-length", args.tiled_max_length),
        ):
            if value is not None and (
                not isinstance(value, int) or value <= 0
            ):
                raise ValueError(f"{name} must be a positive integer.")

        if (
            args.tiled_opt_length is not None
            and args.tiled_max_length is not None
            and args.tiled_max_length < args.tiled_opt_length
        ):
            raise ValueError(
                "--tiled-max-length must be >= --tiled-opt-length."
            )

        if args.tiled_overlap is not None and (
            not isinstance(args.tiled_overlap, int)
            or args.tiled_overlap < 0
        ):
            raise ValueError("--tiled-overlap must be a non-negative integer.")

    if args.varvamp_config is not None:
        validate_custom_varvamp_config(
            args.varvamp_config.expanduser().resolve(),
            mode,
        )


def run_varvamp_command(
    command: list[str],
    *,
    env: dict[str, str],
) -> tuple[int, str]:
    """
    Run VarVAMP while streaming its combined stdout/stderr to the terminal.

    Unlike run_command(), a non-zero return code is returned to the caller so
    interactive mode can adjust parameters and retry VarVAMP without restarting
    orientation, MARS, clustering, MAFFT, trimming or conservation analysis.
    """
    command_text = " ".join(command)
    log_line("COMMAND: " + command_text)

    if VERBOSE:
        print("\nCommand:", command_text)

    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=env,
        )
    except OSError as error:
        raise RuntimeError(
            "Could not start VarVAMP command: " + command_text
        ) from error

    output_lines: list[str] = []
    assert process.stdout is not None

    for line in process.stdout:
        output_lines.append(line)
        print(line, end="")
        log_line("VARVAMP: " + line.rstrip())

    return_code = process.wait()
    return return_code, "".join(output_lines)


def parse_varvamp_diagnostics(output: str) -> dict[str, object]:
    """Extract useful counts and the reported error from VarVAMP terminal output."""
    diagnostics: dict[str, object] = {}

    patterns = {
        "forward_primers": r"(\d+)\s+fw\s+and\s+\d+\s+rv potential primers",
        "reverse_primers": r"\d+\s+fw\s+and\s+(\d+)\s+rv potential primers",
        "potential_probes": r"(\d+)\s+potential qPCR probes",
        "unique_qpcr_amplicons": r"(\d+)\s+unique amplicons with internal probe",
        "delta_g_schemes": (
            r"(\d+)\s+non-overlapping qPCR schemes that passed deltaG cutoff"
        ),
    }

    for key, pattern in patterns.items():
        matches = re.findall(pattern, output, flags=re.IGNORECASE)
        if matches:
            diagnostics[key] = int(matches[-1])

    error_lines = []
    for line in output.splitlines():
        stripped = line.strip()
        if stripped.upper().startswith("ERROR:"):
            error_lines.append(stripped[6:].strip())

    if error_lines:
        diagnostics["reason"] = error_lines[-1]
    else:
        nonempty = [line.strip() for line in output.splitlines() if line.strip()]
        diagnostics["reason"] = (
            nonempty[-1] if nonempty else "VarVAMP exited with a non-zero status."
        )

    return diagnostics


def print_varvamp_failure(
    args: argparse.Namespace,
    mode: str,
    diagnostics: dict[str, object],
) -> None:
    """Display a mode-aware failure summary before offering a retry."""
    title = f"VarVAMP did not produce a {mode.upper()} scheme"
    print("\n" + title)
    print("=" * len(title))
    print("\nReason reported by VarVAMP:")
    print(str(diagnostics.get("reason", "Unknown VarVAMP error.")))

    print("\nCurrent parameters:")
    print(f"  Mode                 : {mode.upper()}")
    print(f"  Consensus threshold  : {args.varvamp_threshold}")
    print(f"  Primer ambiguity     : {args.primer_ambiguity}")
    if mode == "qpcr":
        print(f"  Probe ambiguity      : {args.probe_ambiguity}")
        if args.qpcr_test_n is not None:
            print(f"  qPCR test count (-n) : {args.qpcr_test_n}")
        if args.qpcr_deltag is not None:
            print(f"  deltaG cutoff (-d)   : {args.qpcr_deltag}")

    if mode == "qpcr" and any(
        key in diagnostics
        for key in (
            "forward_primers",
            "reverse_primers",
            "potential_probes",
            "unique_qpcr_amplicons",
            "delta_g_schemes",
        )
    ):
        print("\nVarVAMP nevertheless found:")
        if "forward_primers" in diagnostics:
            print(f"  Forward primers      : {diagnostics['forward_primers']}")
        if "reverse_primers" in diagnostics:
            print(f"  Reverse primers      : {diagnostics['reverse_primers']}")
        if "potential_probes" in diagnostics:
            print(f"  Potential probes     : {diagnostics['potential_probes']}")
        if "unique_qpcr_amplicons" in diagnostics:
            print(
                "  Probe-containing amps : "
                f"{diagnostics['unique_qpcr_amplicons']}"
            )
        if "delta_g_schemes" in diagnostics:
            print(f"  deltaG-passing schemes: {diagnostics['delta_g_schemes']}")

    print("\nPossible adjustment:")
    print("  - lower or otherwise adjust consensus threshold (-t)")
    print("  - increase primer ambiguity (-a)")
    if mode == "qpcr":
        print("  - increase probe ambiguity (-pa)")
    print("  - change mode-specific parameters or use a custom config")


def ensure_mode_required_interactive_parameters(
    args: argparse.Namespace,
    mode: str,
) -> None:
    """Collect required parameters after switching VarVAMP mode."""
    if mode == "qpcr" and args.probe_ambiguity is None:
        if not sys.stdin.isatty():
            raise RuntimeError(
                "--probe-ambiguity is required for qPCR mode "
                "in non-interactive execution."
            )
        args.probe_ambiguity = ask_required_non_negative_integer(
            "Maximum ambiguous bases in the probe",
            "-pa",
        )


def adjust_varvamp_after_failure(
    args: argparse.Namespace,
    mode: str,
    config_dir: Path,
) -> str:
    """
    Adjust parameters after a failed VarVAMP attempt.

    Returns:
        "retry"       retry VarVAMP
        "mode_changed" retry after changing mode
        "stop"        stop VarVAMP design
    """
    if mode == "qpcr":
        choices = [
            ("threshold", "Change consensus threshold (-t)."),
            ("primer", "Change primer ambiguity (-a)."),
            ("probe", "Change probe ambiguity (-pa)."),
            ("several", "Change several VarVAMP parameters."),
            ("advanced", "Configure advanced qPCR parameters (Python config)."),
            ("mode", "Change VarVAMP assay mode."),
            ("stop", "Stop the workflow."),
        ]
    else:
        choices = [
            ("threshold", "Change consensus threshold (-t)."),
            ("primer", "Change primer ambiguity (-a)."),
            (
                "mode_cli",
                f"Change additional {mode.upper()} command-line parameters.",
            ),
            ("several", "Change several VarVAMP parameters."),
            (
                "advanced",
                f"Configure advanced {mode.upper()} parameters (Python config).",
            ),
            ("mode", "Change VarVAMP assay mode."),
            ("stop", "Stop the workflow."),
        ]

    choice = ask_required_choice("What would you like to do?", choices)

    if choice == "threshold":
        args.varvamp_threshold = ask_required_float(
            "VarVAMP consensus threshold",
            "-t",
            minimum=0.01,
            maximum=1.0,
        )
        return "retry"

    if choice == "primer":
        args.primer_ambiguity = ask_required_non_negative_integer(
            "Maximum ambiguous bases in each primer",
            "-a",
        )
        return "retry"

    if choice == "probe":
        args.probe_ambiguity = ask_required_non_negative_integer(
            "Maximum ambiguous bases in the probe",
            "-pa",
        )
        return "retry"

    if choice == "mode_cli":
        configure_interactive_varvamp_advanced(
            args,
            mode,
            force=True,
        )
        return "retry"

    if choice == "several":
        args.varvamp_threshold = ask_required_float(
            "VarVAMP consensus threshold",
            "-t",
            minimum=0.01,
            maximum=1.0,
        )
        args.primer_ambiguity = ask_required_non_negative_integer(
            "Maximum ambiguous bases in each primer",
            "-a",
        )
        if mode == "qpcr":
            args.probe_ambiguity = ask_required_non_negative_integer(
                "Maximum ambiguous bases in the probe",
                "-pa",
            )
        if ask_yes_no(
            "Also configure additional mode-specific command-line parameters?"
        ):
            configure_interactive_varvamp_advanced(
                args,
                mode,
                force=True,
            )
        return "retry"

    if choice == "advanced":
        changed = configure_custom_varvamp_file(
            args,
            mode,
            config_dir,
        )
        if changed:
            return "retry"
        return adjust_varvamp_after_failure(
            args,
            mode,
            config_dir,
        )

    if choice == "mode":
        args.varvamp_mode = choose_varvamp_mode()
        # A mode change returns to VarVAMP's normal/default Python config unless
        # the user explicitly selects a custom config again.
        args.varvamp_config = None
        ensure_mode_required_interactive_parameters(
            args,
            args.varvamp_mode,
        )
        if ask_yes_no(
            "Configure additional mode-specific VarVAMP parameters for the new mode?"
        ):
            configure_interactive_varvamp_advanced(
                args,
                args.varvamp_mode,
                force=True,
            )
        return "mode_changed"

    return "stop"


def copy_varvamp_tree(varvamp_dir: Path, results_dir: Path, mode: str) -> Path:
    """Copy the complete successful VarVAMP tree into this isolated run."""

    destination = results_dir / "varvamp"
    if destination.exists():
        raise RuntimeError(
            "The VarVAMP result directory already exists inside this run: "
            f"{destination}. Refusing to overwrite an existing result."
        )
    shutil.copytree(varvamp_dir, destination)
    (destination / "assay_mode.txt").write_text(
        mode.upper() + "\n",
        encoding="utf-8",
    )
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
    Display VarVAMP text results directly in the terminal after a successful run.

    PDF and image files are intentionally excluded. Text-based scientific output
    such as TSV, TABULAR, BED, CSV, TXT and FASTA is displayed.
    """
    text_extensions = {
        ".tsv",
        ".tabular",
        ".bed",
        ".txt",
        ".csv",
        ".fasta",
        ".fa",
        ".fna",
    }

    preferred_by_mode = {
        "qpcr": [
            "qpcr_primers.tsv",
            "qpcr_design.tsv",
            "oligos.fasta",
            "primers.bed",
            "amplicons.bed",
        ],
        "single": [
            "primer.tsv",
            "primer_to_amplicon_assignments.tabular",
            "primers.bed",
            "amplicons.bed",
        ],
        "tiled": [
            "primer.tsv",
            "primer_to_amplicon_assignments.tabular",
            "primers.bed",
            "amplicons.bed",
        ],
    }

    files = [
        path
        for path in varvamp_results_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in text_extensions
    ]

    preferred_names = preferred_by_mode.get(mode, [])
    preferred_rank = {
        name: index
        for index, name in enumerate(preferred_names)
    }

    files.sort(
        key=lambda path: (
            preferred_rank.get(path.name, len(preferred_rank) + 1),
            str(path.relative_to(varvamp_results_dir)),
        )
    )

    section(f"VarVAMP {mode.upper()} text results")

    if not files:
        print("No terminal-readable VarVAMP result files were found.")
        print(f"Result directory: {compact_path(varvamp_results_dir)}")
        return

    for path in files:
        relative = path.relative_to(varvamp_results_dir)
        title = str(relative)
        print(f"\n{title}")
        print("-" * len(title))

        try:
            content = path.read_text(
                encoding="utf-8",
                errors="replace",
            ).rstrip()
        except OSError as error:
            print(f"[Could not read file: {error}]")
            continue

        if content:
            print(content)
        else:
            print("[empty file]")

    print(
        f"\nResult directory: {compact_path(varvamp_results_dir)}"
    )


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
        cpu_threads = available_cpu_threads()

        if args.threads is None:
            args.threads = automatic_thread_count(cpu_threads)
            thread_selection = "automatic"
        else:
            thread_selection = "manual"

        if args.threads <= 0:
            raise ValueError("--threads must be greater than zero.")

        if args.threads > cpu_threads:
            print(
                f"Warning: --threads={args.threads} exceeds the "
                f"{cpu_threads} logical CPU threads currently available."
            )

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
        run_id, workdir, results_dir = resolve_run_directories(
            args,
            project_name,
        )

        # Run-specific working areas contain reconstructible/intermediate files.
        preprocessing_workdir = workdir / "preprocessing"
        alignment_workdir = workdir / "alignment"
        varvamp_work_root = workdir / "varvamp"

        # Run-specific result areas contain permanent scientific outputs.
        preprocessing_results_dir = results_dir / "preprocessing"
        alignment_results_dir = results_dir / "alignment"
        conservation_results_dir = results_dir / "conservation"
        config_results_dir = results_dir / "config"

        for directory in (
            preprocessing_workdir,
            alignment_workdir,
            varvamp_work_root,
            preprocessing_results_dir,
            alignment_results_dir,
            conservation_results_dir,
            config_results_dir,
        ):
            directory.mkdir(parents=True, exist_ok=True)

        WORKFLOW_LOG = results_dir / "workflow.log"
        WORKFLOW_LOG.write_text(
            "VarVAMP assay-design workflow log\n"
            f"Project: {project_name}\n"
            f"Run ID: {run_id}\n"
            f"Design FASTA: {compact_path(input_file)}\n"
            f"CPU threads available: {cpu_threads}\n"
            f"Workflow threads: {args.threads} ({thread_selection})\n",
            encoding="utf-8",
        )

        print("╔══════════════════════════════════════════════╗")
        print("║       VarVAMP ASSAY DESIGN WORKFLOW         ║")
        print("╚══════════════════════════════════════════════╝")
        print(f"Project : {project_name}")
        print(f"Run ID  : {run_id}")
        print(f"Input   : {compact_path(input_file)}")
        print(f"Work    : {compact_path(workdir)}")
        print(f"Results : {compact_path(results_dir)}")
        print(f"CPU     : {cpu_threads} logical threads available")
        print(
            f"Threads : {args.threads} used by supported tools "
            f"({thread_selection})"
        )
        if VERBOSE:
            print(f"Log     : {WORKFLOW_LOG}")

        interactive = sys.stdin.isatty()

        if not interactive and args.trimming is None:
            raise RuntimeError(
                "--trimming is required in non-interactive mode."
            )

        if (
            args.trimming == "manual"
            and args.trimal_gap_threshold is None
            and not interactive
        ):
            raise RuntimeError(
                "--trimal-gap-threshold is required for --trimming manual."
            )

        # ------------------------------------------------------------------
        # Input overview
        # ------------------------------------------------------------------
        current_fasta = input_file
        initial_count = count_fasta_sequences(input_file)

        print_fasta_statistics(
            calculate_fasta_statistics(input_file),
            title="Input dataset statistics",
        )

        # ------------------------------------------------------------------
        # 1. Orientation normalization
        # ------------------------------------------------------------------
        print("\n1. Sequence orientation checkpoint")

        if args.orientation is None:
            if not interactive:
                raise RuntimeError(
                    "--orientation is required in non-interactive mode."
                )
            args.orientation = choose_orientation_strategy()

        if args.orientation not in {"keep", "adjust", "accurate"}:
            raise ValueError("Unsupported orientation strategy.")

        reversed_count = 0

        if args.orientation == "keep":
            section("Orientation result")
            print("Method                    : keep original orientation")
            print("MAFFT orientation check   : skipped")
            print(f"Sequences retained        : {initial_count}")
            print("Sequences changed         : 0")
        else:
            check_required_tools([("MAFFT", args.mafft_executable)])
            print("Running sequence orientation normalization with MAFFT.")
            current_fasta, reversed_count = normalize_orientation_with_mafft(
                current_fasta,
                preprocessing_workdir,
                project_name,
                args.orientation,
                args.threads,
                args.mafft_executable,
            )

            section("Orientation result")
            print(
                "Method                    : "
                + (
                    "MAFFT --adjustdirection"
                    if args.orientation == "adjust"
                    else "MAFFT --adjustdirectionaccurately"
                )
            )
            print(f"Sequences checked         : {initial_count}")
            print(
                f"Reverse-complemented      : {reversed_count} "
                f"({100.0 * reversed_count / initial_count if initial_count else 0.0:.2f} %)"
            )
            print(
                f"Orientation unchanged     : {initial_count - reversed_count}"
            )

        if args.orientation != "keep":
            oriented_result = (
                preprocessing_results_dir / f"{project_name}_oriented.fasta"
            )
            shutil.copy2(current_fasta, oriented_result)
            orientation_report = (
                preprocessing_workdir / f"{project_name}_orientation_report.txt"
            )
            if orientation_report.is_file():
                shutil.copy2(
                    orientation_report,
                    preprocessing_results_dir / orientation_report.name,
                )

        # ------------------------------------------------------------------
        # 2. Circular start-position normalization
        # ------------------------------------------------------------------
        print("\n2. Sequence topology checkpoint")

        if args.topology is None:
            if not interactive:
                raise RuntimeError(
                    "--topology is required in non-interactive mode."
                )
            args.topology = choose_sequence_topology()

        if args.topology == "circular":
            check_required_tools([("MARS", args.mars_executable)])
            print("Running circular sequence rotation with MARS.")
            before_mars = current_fasta
            current_fasta = run_mars_rotation(
                current_fasta,
                preprocessing_workdir,
                project_name,
                args.mars_executable,
                args.threads,
            )
            print_rotation_statistics(
                calculate_cyclic_rotation_statistics(
                    before_mars,
                    current_fasta,
                )
            )
        elif args.topology == "linear":
            section("Topology result")
            print("Topology                  : linear")
            print("MARS                      : not applicable")
            print("Cyclic starts changed     : 0")
        else:
            raise ValueError("--topology must be 'linear' or 'circular'.")

        if args.topology == "circular":
            shutil.copy2(
                current_fasta,
                preprocessing_results_dir / f"{project_name}_rotated.fasta",
            )

        # ------------------------------------------------------------------
        # 3. Redundancy handling
        # ------------------------------------------------------------------
        print("\n3. Sequence redundancy checkpoint")

        if args.redundancy is None:
            if args.skip_cdhit:
                args.redundancy = "none"
            elif not interactive:
                raise RuntimeError(
                    "--redundancy is required in non-interactive mode."
                )
            else:
                args.redundancy = choose_redundancy_strategy()

        if args.redundancy == "cdhit":
            if args.identity is None:
                if not interactive:
                    raise RuntimeError(
                        "--identity is required for CD-HIT-EST."
                    )
                args.identity = ask_required_float(
                    "CD-HIT-EST identity threshold",
                    "-c",
                    minimum=0.80,
                    maximum=1.0,
                )

            if not 0.80 <= args.identity <= 1.0:
                raise ValueError(
                    "--identity must be between 0.80 and 1.0."
                )

            if args.cdhit_mode is None:
                if not interactive:
                    raise RuntimeError(
                        "--cdhit-mode is required for CD-HIT-EST."
                    )
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
                    raise RuntimeError(
                        "--cdhit-strand is required for CD-HIT-EST."
                    )
                args.cdhit_strand = ask_required_choice(
                    "CD-HIT-EST strand comparison",
                    [
                        (
                            "both",
                            "Compare both +/+ and +/- orientations (-r 1).",
                        ),
                        (
                            "same",
                            "Compare only the same orientation, +/+ (-r 0).",
                        ),
                    ],
                )

        redundancy_output = (
            preprocessing_workdir / f"{project_name}_redundancy_filtered.fasta"
        )
        before_redundancy = count_fasta_sequences(current_fasta)
        word_size: int | None = None

        if args.redundancy == "none":
            shutil.copy2(current_fasta, redundancy_output)

        elif args.redundancy == "seqkit":
            check_required_tools([("SeqKit", args.seqkit_executable)])
            print("Running exact duplicate removal with SeqKit.")
            run_seqkit_deduplication(
                current_fasta,
                redundancy_output,
                args.seqkit_executable,
            )

        elif args.redundancy == "cdhit":
            check_required_tools(
                [("CD-HIT-EST", args.cdhit_executable)]
            )
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

        else:
            raise RuntimeError(
                f"Unknown redundancy strategy: {args.redundancy}"
            )

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

        section("Redundancy reduction result")
        print(f"Strategy                  : {args.redundancy}")
        if args.redundancy == "seqkit":
            print("Matching criterion         : exact sequence identity")
        elif args.redundancy == "cdhit":
            print(f"Identity threshold (-c)    : {args.identity}")
            print(f"Word size (-n)             : {word_size}")
            print(
                "Clustering mode (-g)       : "
                f"{args.cdhit_mode} "
                f"({'1' if args.cdhit_mode == 'accurate' else '0'})"
            )
            print(
                "Strand comparison (-r)     : "
                f"{args.cdhit_strand} "
                f"({'1' if args.cdhit_strand == 'both' else '0'})"
            )

        print(f"Sequences before          : {before_redundancy}")
        print(
            f"Sequences removed         : {removed_redundancy} "
            f"({reduction:.2f} %)"
        )
        print(
            f"Sequences retained        : {after_redundancy} "
            f"({retained_percentage:.2f} %)"
        )
        if args.redundancy == "cdhit":
            print(
                f"Cluster representatives   : {after_redundancy}"
            )

        # Keep the biologically meaningful preprocessing output for provenance.
        redundancy_result = (
            preprocessing_results_dir
            / f"{project_name}_redundancy_filtered.fasta"
        )
        shutil.copy2(redundancy_output, redundancy_result)
        if args.redundancy == "cdhit":
            cluster_file = Path(str(redundancy_output) + ".clstr")
            if cluster_file.is_file():
                shutil.copy2(
                    cluster_file,
                    preprocessing_results_dir / cluster_file.name,
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
            raise ValueError(
                f"Unsupported MAFFT strategy: {args.mafft_strategy}"
            )

        if interactive:
            warn_about_mafft_strategy(
                args.mafft_strategy,
                after_redundancy,
            )

        check_required_tools([("MAFFT", args.mafft_executable)])
        final_alignment = alignment_workdir / f"{project_name}_alignment.fasta"

        print(
            f"Running final MAFFT alignment: {args.mafft_strategy} "
            f"({args.threads} threads)",
            flush=True,
        )
        mafft_reported_strategy = run_final_mafft(
            redundancy_output,
            final_alignment,
            args.mafft_strategy,
            args.threads,
            args.mafft_executable,
        )
        read_alignment(final_alignment)

        final_alignment_result = (
            alignment_results_dir / f"{project_name}_alignment.fasta"
        )
        shutil.copy2(final_alignment, final_alignment_result)
        downstream_alignment_result = final_alignment_result

        pretrim_stats, pretrim_sequences = calculate_alignment_statistics(
            final_alignment
        )

        section("MAFFT alignment result")
        print(f"Sequences aligned         : {after_redundancy}")
        print(f"Requested strategy        : {args.mafft_strategy}")
        print(f"MAFFT reported strategy   : {mafft_reported_strategy}")
        print(
            f"Alignment length          : "
            f"{int(pretrim_stats['alignment_length'])} positions"
        )
        print(
            f"Mean gap content          : "
            f"{float(pretrim_stats['mean_gap_content_percent']):.2f} %"
        )
        print(
            f"Fully occupied columns    : "
            f"{int(pretrim_stats['fully_occupied_columns'])}"
        )

        # ------------------------------------------------------------------
        # 5. Alignment QC and trimming checkpoint
        # ------------------------------------------------------------------
        print("\n5. Alignment QC and trimming checkpoint")

        print_alignment_statistics(pretrim_stats)
        save_alignment_statistics(
            pretrim_stats,
            pretrim_sequences,
            alignment_results_dir,
            project_name,
            "pretrim",
        )

        if args.trimming is None:
            args.trimming = choose_trimming_strategy()

        if (
            args.trimming == "manual"
            and args.trimal_gap_threshold is None
        ):
            print("\nManual trimAl gap threshold")
            print("============================")
            print(
                "The -gt value is the minimum fraction of sequences without "
                "a gap in a retained column. Examples: 0.90 = >=90% "
                "occupancy; 1.00 = no gaps allowed."
            )
            args.trimal_gap_threshold = ask_required_float(
                "trimAl gap threshold",
                "-gt",
                minimum=0.0,
                maximum=1.0,
            )

            if args.trimal_gap_threshold < 0.50:
                print(
                    "Warning: this is a permissive threshold; many highly "
                    "gapped columns may be retained."
                )
                if not ask_yes_no(
                    f"Continue with -gt {args.trimal_gap_threshold}?"
                ):
                    raise RuntimeError(
                        "Manual trimAl threshold rejected by user."
                    )

        downstream_alignment = final_alignment

        if args.trimming != "none":
            check_required_tools(
                [("trimAl", args.trimal_executable)]
            )
            trimmed_alignment = (
                alignment_workdir / f"{project_name}_alignment_trimmed.fasta"
            )

            print(f"\nRunning trimAl strategy: {args.trimming}")
            run_trimal(
                final_alignment,
                trimmed_alignment,
                args.trimming,
                args.trimal_executable,
                args.trimal_gap_threshold,
            )
            downstream_alignment = trimmed_alignment

            trimmed_alignment_result = (
                alignment_results_dir
                / f"{project_name}_alignment_trimmed.fasta"
            )
            shutil.copy2(trimmed_alignment, trimmed_alignment_result)
            downstream_alignment_result = trimmed_alignment_result

            posttrim_stats, posttrim_sequences = (
                calculate_alignment_statistics(trimmed_alignment)
            )

            save_alignment_statistics(
                posttrim_stats,
                posttrim_sequences,
                alignment_results_dir,
                project_name,
                "posttrim",
            )

            original_length = int(pretrim_stats["alignment_length"])
            trimmed_length = int(posttrim_stats["alignment_length"])
            removed_columns = original_length - trimmed_length

            section("trimAl effect")
            print(f"Strategy                  : {args.trimming}")
            if args.trimming == "manual":
                print(
                    f"Gap threshold (-gt)       : "
                    f"{args.trimal_gap_threshold}"
                )
            print(f"Columns before            : {original_length}")
            print(f"Columns after             : {trimmed_length}")
            print(
                f"Columns removed           : {removed_columns} "
                f"({100.0 * removed_columns / original_length:.2f} %)"
            )
            print(
                f"Mean gap content before   : "
                f"{float(pretrim_stats['mean_gap_content_percent']):.2f} %"
            )
            print(
                f"Mean gap content after    : "
                f"{float(posttrim_stats['mean_gap_content_percent']):.2f} %"
            )

            print_alignment_statistics(
                posttrim_stats,
                title="Post-trimming alignment statistics",
            )
        else:
            section("trimAl result")
            print("Trimming                  : skipped")
            print("Alignment used downstream : complete MAFFT alignment")

        # ------------------------------------------------------------------
        # 6. Conservation analysis
        # ------------------------------------------------------------------
        print("\n6. Conservation analysis")
        perform_conservation_analysis(
            downstream_alignment,
            conservation_results_dir,
            project_name,
            args.min_occupancy,
            args.min_major_frequency,
        )

        # ------------------------------------------------------------------
        # 7. VarVAMP assay design with retry-on-failure
        # ------------------------------------------------------------------
        varvamp_results: Path | None = None
        varvamp_attempts: list[dict[str, object]] = []

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
                    raise RuntimeError(
                        "--varvamp-threshold is required."
                    )
                args.varvamp_threshold = ask_required_float(
                    "VarVAMP consensus threshold",
                    "-t",
                    minimum=0.01,
                    maximum=1.0,
                )

            if args.primer_ambiguity is None:
                if not interactive:
                    raise RuntimeError(
                        "--primer-ambiguity is required."
                    )
                args.primer_ambiguity = (
                    ask_required_non_negative_integer(
                        "Maximum ambiguous bases in each primer",
                        "-a",
                    )
                )

            ensure_mode_required_interactive_parameters(
                args,
                args.varvamp_mode,
            )

            if interactive:
                configure_interactive_varvamp_advanced(
                    args,
                    args.varvamp_mode,
                )

            check_required_tools(
                [("VarVAMP", args.varvamp_executable)]
            )

            attempt_number = 0

            while True:
                attempt_number += 1
                mode = args.varvamp_mode

                ensure_mode_required_interactive_parameters(
                    args,
                    mode,
                )
                validate_varvamp_mode_parameters(args, mode)

                varvamp_workdir = (
                    varvamp_work_root
                    / f"attempt_{attempt_number:02d}_{mode}"
                )
                if varvamp_workdir.exists():
                    raise RuntimeError(
                        "Unexpected existing VarVAMP attempt directory: "
                        f"{varvamp_workdir}. Refusing to overwrite it."
                    )

                command = build_varvamp_command(
                    args,
                    mode,
                    downstream_alignment,
                    varvamp_workdir,
                )

                command_env = os.environ.copy()

                # Every new script execution starts with VarVAMP's normal
                # configuration unless a custom config is explicitly selected.
                command_env.pop("VARVAMP_CONFIG", None)

                config_path: Path | None = None
                if args.varvamp_config is not None:
                    config_path = stage_varvamp_config_for_run(
                        args.varvamp_config,
                        config_results_dir,
                        mode,
                        attempt_number,
                    )
                    # Keep the run-local copy as the canonical provenance path.
                    args.varvamp_config = config_path
                    command_env["VARVAMP_CONFIG"] = str(config_path)

                section(f"VarVAMP attempt {attempt_number}")
                print(f"Mode                      : {mode.upper()}")
                print(
                    f"Consensus threshold (-t)  : "
                    f"{args.varvamp_threshold}"
                )
                print(
                    f"Primer ambiguity (-a)     : "
                    f"{args.primer_ambiguity}"
                )
                if mode == "qpcr":
                    print(
                        f"Probe ambiguity (-pa)     : "
                        f"{args.probe_ambiguity}"
                    )
                print(
                    "Python configuration      : "
                    + (
                        compact_path(config_path)
                        if config_path is not None
                        else "VarVAMP default"
                    )
                )

                return_code, varvamp_output = run_varvamp_command(
                    command,
                    env=command_env,
                )

                diagnostics = parse_varvamp_diagnostics(
                    varvamp_output
                )

                attempt_record: dict[str, object] = {
                    "attempt": attempt_number,
                    "mode": mode,
                    "consensus_threshold": args.varvamp_threshold,
                    "primer_ambiguity": args.primer_ambiguity,
                    "probe_ambiguity": (
                        args.probe_ambiguity
                        if mode == "qpcr"
                        else None
                    ),
                    "single_opt_length": args.single_opt_length,
                    "single_max_length": args.single_max_length,
                    "single_report_n": args.single_report_n,
                    "tiled_opt_length": args.tiled_opt_length,
                    "tiled_max_length": args.tiled_max_length,
                    "tiled_overlap": args.tiled_overlap,
                    "qpcr_test_n": args.qpcr_test_n,
                    "qpcr_deltag": args.qpcr_deltag,
                    "custom_config": (
                        portable_manifest_path(config_path)
                        if config_path is not None
                        else None
                    ),
                    "return_code": return_code,
                    "status": (
                        "success"
                        if return_code == 0
                        else "failed"
                    ),
                }

                if return_code != 0:
                    attempt_record["reason"] = diagnostics.get(
                        "reason"
                    )

                varvamp_attempts.append(attempt_record)
                log_line(
                    "VARVAMP_ATTEMPT: "
                    + json.dumps(
                        attempt_record,
                        sort_keys=True,
                    )
                )

                if return_code == 0:
                    varvamp_results = copy_varvamp_tree(
                        varvamp_workdir,
                        results_dir,
                        mode,
                    )
                    print_varvamp_tables(
                        varvamp_results,
                        mode,
                    )
                    create_varvamp_visual_summaries(
                        varvamp_results
                    )
                    break

                if not interactive:
                    raise RuntimeError(
                        "VarVAMP failed in non-interactive mode: "
                        + str(
                            diagnostics.get(
                                "reason",
                                "unknown error",
                            )
                        )
                    )

                print_varvamp_failure(
                    args,
                    mode,
                    diagnostics,
                )

                action = adjust_varvamp_after_failure(
                    args,
                    mode,
                    config_results_dir,
                )

                if action == "stop":
                    print(
                        "VarVAMP design stopped by the user. "
                        "Previous workflow steps remain saved."
                    )
                    args.skip_varvamp = True
                    break

                # "retry" and "mode_changed" both return here and rerun only
                # the VarVAMP stage. All previous workflow products are reused.
                print("\nRetrying VarVAMP only...")
                print("Orientation              : already completed")
                print("Topology normalization   : already completed")
                print("Redundancy handling      : already completed")
                print("Final MAFFT              : already completed")
                print("trimAl / alignment QC    : already completed")
                print("Conservation analysis    : already completed")

        else:
            print("\n7. VarVAMP assay design skipped")

        manifest_file = results_dir / "assay_design_manifest.json"
        summary_file = results_dir / "workflow_summary.txt"

        write_project_manifest(
            manifest_file,
            project_name=project_name,
            run_id=run_id,
            input_file=input_file,
            workdir=workdir,
            results_dir=results_dir,
            downstream_alignment=downstream_alignment_result,
            initial_count=initial_count,
            after_redundancy=after_redundancy,
            args=args,
            varvamp_results=varvamp_results,
            varvamp_attempts=varvamp_attempts,
        )

        write_workflow_summary(
            summary_file,
            project_name=project_name,
            run_id=run_id,
            input_file=input_file,
            downstream_alignment=downstream_alignment_result,
            workdir=workdir,
            results_dir=results_dir,
            initial_count=initial_count,
            after_redundancy=after_redundancy,
            args=args,
            varvamp_results=varvamp_results,
        )

        if varvamp_results is not None:
            shutil.copy2(
                manifest_file,
                varvamp_results
                / "assay_design_manifest.json",
            )

        section("Workflow completed")
        print(f"Run ID                   : {run_id}")
        print(f"Initial sequences        : {initial_count}")
        print(
            f"Final MAFFT sequences    : {after_redundancy}"
        )
        print(
            f"Downstream alignment     : "
            f"{compact_path(downstream_alignment_result)}"
        )

        if varvamp_results is not None:
            print(
                f"VarVAMP mode             : "
                f"{args.varvamp_mode.upper()}"
            )
            print(
                f"VarVAMP results          : "
                f"{compact_path(varvamp_results)}"
            )
        elif args.skip_varvamp:
            print("VarVAMP                  : no final design")
        else:
            print("VarVAMP                  : skipped")

        print(
            f"Project manifest         : "
            f"{compact_path(manifest_file)}"
        )
        print(
            f"Workflow summary         : "
            f"{compact_path(summary_file)}"
        )
        print(
            f"Workflow log             : "
            f"{compact_path(WORKFLOW_LOG)}"
        )
        print(f"Permanent run results    : {compact_path(results_dir)}")
        print(f"Reconstructible work     : {compact_path(workdir)}")
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