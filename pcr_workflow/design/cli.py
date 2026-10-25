# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import argparse

from ..core.paths import compact_path
from .constants import MAFFT_STRATEGIES

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
