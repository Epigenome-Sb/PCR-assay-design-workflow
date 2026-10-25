# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from math import ceil
from pathlib import Path

from ..core.paths import compact_path
from .database import (
    apply_sequence_length_filter,
    create_mfeprimer_index,
    find_candidate_target_fastas,
    mfeprimer_index_files,
    mfeprimer_index_status,
)
from .fasta import fasta_statistics
from .runtime import progress

def ask_yes_no(question: str) -> bool:
    while True:
        value = input(f"{question} [y/n]: ").strip().lower()
        if value in {"y", "yes"}:
            return True
        if value in {"n", "no"}:
            return False
        print("Please answer y or n.")
def ask_choice(title: str, choices: list[str]) -> int:
    print(f"\n{title}")
    print("=" * len(title))
    for idx, item in enumerate(choices, start=1):
        print(f"{idx}. {item}")

    while True:
        raw = input("Select an option by number: ").strip()
        try:
            selected = int(raw)
        except ValueError:
            print("Please enter one of the displayed numbers.")
            continue
        if 1 <= selected <= len(choices):
            return selected
        print("Selection outside the available range.")
def ask_int(label: str, minimum: int, maximum: int) -> int:
    while True:
        raw = input(f"{label} [{minimum}-{maximum}]: ").strip()
        try:
            value = int(raw)
        except ValueError:
            print("Please enter an integer.")
            continue
        if minimum <= value <= maximum:
            return value
        print(f"Value must be between {minimum} and {maximum}.")
def ask_float(label: str, minimum: float, maximum: float) -> float:
    while True:
        raw = input(f"{label} [{minimum}-{maximum}]: ").strip()
        try:
            value = float(raw)
        except ValueError:
            print("Please enter a numeric value.")
            continue
        if minimum <= value <= maximum:
            return value
        print(f"Value must be between {minimum} and {maximum}.")
def choose_sequence_length_filter(
    database: Path,
    database_workdir: Path,
) -> tuple[Path, str, dict[str, float | int]]:
    """
    Decide which sequence-length set enters MFEprimer preparation.

    The "near-complete" option is an operational length criterion relative
    to the longest sequence in the selected database; it does not infer
    biological completeness from annotation.
    """
    before_stats = fasta_statistics(database)
    max_length = int(before_stats["max_length"])

    print("\n3. Sequence-length filtering")
    print("============================")
    print(
        f"Current database: {before_stats['sequences']} sequences; "
        f"length range {before_stats['min_length']}-{before_stats['max_length']} nt."
    )

    mode = ask_choice(
        "How should the validation database be prepared?",
        [
            "Keep all sequences.",
            "Keep complete / nearly complete sequences using a relative-length threshold.",
            "Apply a custom minimum sequence length.",
        ],
    )

    threshold: int | None = None
    strategy: str

    if mode == 1:
        strategy = "keep all sequences"
        retained_database = database
    elif mode == 2:
        print(
            "\nNear-complete filtering is defined here by sequence length "
            "relative to the longest sequence in this database."
        )
        relative_mode = ask_choice(
            "Near-complete length threshold",
            [
                "Keep sequences >=90% of the maximum sequence length.",
                "Keep sequences >=95% of the maximum sequence length.",
                "Enter a custom percentage of the maximum sequence length.",
            ],
        )
        if relative_mode == 1:
            percentage = 90.0
        elif relative_mode == 2:
            percentage = 95.0
        else:
            percentage = ask_float(
                "Minimum relative length (% of maximum)",
                50.0,
                100.0,
            )
        threshold = ceil(max_length * percentage / 100.0)
        strategy = (
            f"near-complete length >= {percentage:g}% of maximum "
            f"({threshold} nt)"
        )
    else:
        threshold = ask_int("Minimum sequence length (nt)", 1, max_length)
        strategy = f"minimum sequence length >= {threshold} nt"

    return apply_sequence_length_filter(
        database, database_workdir, threshold, strategy, before_stats
    )
def choose_target_database(
    requested_database: Path | None = None,
) -> Path:
    """Select the independent validation FASTA, normally from data/validation/."""
    print("\n2. Validation database")
    print("======================")

    if requested_database is not None:
        candidate = requested_database.expanduser().resolve()
        if not candidate.is_file():
            raise FileNotFoundError(f"Validation FASTA not found: {candidate}")
        return candidate

    candidates = find_candidate_target_fastas()
    if candidates:
        for idx, path in enumerate(candidates, start=1):
            print(f"{idx}. {compact_path(path)}")
        print(f"{len(candidates) + 1}. Enter another FASTA path")

        while True:
            raw = input("Select the validation database by number: ").strip()
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
        raw = input("Name or path of the validation FASTA: ").strip().strip("'\"")
        if not raw:
            print("A FASTA path is required.")
            continue
        candidate = Path(raw).expanduser()
        paths = [candidate]
        if not candidate.is_absolute():
            paths.extend(
                [
                    Path.cwd() / candidate,
                    Path.cwd() / "data" / "validation" / candidate,
                ]
            )
        for path in paths:
            if path.is_file():
                return path.resolve()
        print("FASTA file not found. Preferred location: data/validation/.")
def choose_schemes(schemes: dict[str, dict[str, str]]) -> list[str]:
    """Select one, several, or all VarVAMP schemes."""
    scheme_names = sorted(schemes)

    print("\nAssays to validate")
    print("------------------")
    for idx, scheme in enumerate(scheme_names, start=1):
        probe_text = (
            f" | PROBE {len(schemes[scheme]['PROBE'])} nt"
            if "PROBE" in schemes[scheme]
            else ""
        )
        print(
            f"{idx}. {scheme}: LEFT {len(schemes[scheme]['LEFT'])} nt | "
            f"RIGHT {len(schemes[scheme]['RIGHT'])} nt{probe_text}"
        )
    print(f"{len(scheme_names) + 1}. ALL")

    while True:
        raw = input("Select one/multiple numbers (e.g. 1,3) or ALL: ").strip()
        try:
            numbers = [int(item.strip()) for item in raw.split(",") if item.strip()]
        except ValueError:
            print("Use numbers separated by commas, e.g. 1,2.")
            continue
        if not numbers:
            print("At least one selection is required.")
            continue

        all_option = len(scheme_names) + 1
        if all_option in numbers:
            if len(numbers) > 1:
                print("Choose ALL alone, or specific assay numbers.")
                continue
            return scheme_names

        if any(number < 1 or number > len(scheme_names) for number in numbers):
            print("One or more selected numbers are outside the available range.")
            continue

        selected: list[str] = []
        for number in numbers:
            scheme = scheme_names[number - 1]
            if scheme not in selected:
                selected.append(scheme)
        return selected
def choose_mfeprimer_parameters() -> dict[str, float | int] | None:
    print("\n6. Validation parameters")
    print("========================")
    mode = ask_choice(
        "MFEprimer search filters",
        [
            "Use defaults (0-2000 bp; Tm >=30 °C).",
            "Set custom amplicon size and Tm thresholds.",
        ],
    )
    if mode == 1:
        return None

    min_size = ask_int("Minimum predicted amplicon size (-s)", 0, 100000)
    max_size = ask_int("Maximum predicted amplicon size (-S)", min_size, 100000)
    tm_cutoff = ask_float("Minimum Tm cutoff (-t)", 0.0, 100.0)
    return {
        "min_size": min_size,
        "max_size": max_size,
        "tm_cutoff": tm_cutoff,
    }
def choose_coverage_amplicon_filter() -> tuple[int | None, int | None]:
    """
    Choose the size range used to count an amplicon as a target-coverage hit.

    This filter is applied AFTER MFEprimer has generated its report, so
    non-target/self-priming products remain visible in the QC output.
    """
    mode = ask_choice(
        "Target-amplicon extraction filter",
        [
            "qPCR range: keep target amplicons from 70 to 250 bp.",
            "Enter a custom target-amplicon size range.",
            "Do not apply an additional size filter.",
        ],
    )

    if mode == 1:
        return 70, 250

    if mode == 2:
        min_size = ask_int(
            "Minimum target amplicon size",
            0,
            100000,
        )
        max_size = ask_int(
            "Maximum target amplicon size",
            min_size,
            100000,
        )
        return min_size, max_size

    return None, None
def index_target_database(database: Path) -> dict[str, object]:
    """
    Ensure that the formatted validation database has an MFEprimer index and
    explicitly report whether it was created or reused.
    """
    print("\n5. MFEprimer database indexing")
    print("==============================")

    indexed, detail = mfeprimer_index_status(database)
    existing_files = mfeprimer_index_files(database)

    if indexed:
        print("Existing MFEprimer index detected.")
        print("Status   : EXISTING INDEX REUSED")
        print(f"Database : {compact_path(database)}")
        print(f"Layout   : {detail}")
        print("\nIndex files detected")
        print("--------------------")
        for path in existing_files:
            print(f"✓ {compact_path(path)}")
        print("\nNo re-indexing was necessary.")
        return {
            "status": "REUSED",
            "detail": detail,
            "files": existing_files,
        }

    print("No complete MFEprimer index was detected for this formatted FASTA.")
    if not ask_yes_no("Index this database with MFEprimer now?"):
        raise RuntimeError(
            "Coverage validation requires an indexed MFEprimer database."
        )

    index_mode = ask_choice(
        "MFEprimer indexing parameters",
        [
            "Use MFEprimer default indexing parameters.",
            "Set a custom k-mer seed length (-k).",
        ],
    )

    k_value = None
    if index_mode == 2:
        k_value = ask_int("MFEprimer index k-mer length (-k)", 9, 15)

    progress("Indexing database with MFEprimer", "RUN")
    result = create_mfeprimer_index(database, k_value)
    progress("Indexing database with MFEprimer", "OK")
    detail = result["detail"]
    generated_files = result["files"]

    print("\nIndexing result")
    print("---------------")
    print("Status   : INDEX CREATED SUCCESSFULLY")
    print(f"Database : {compact_path(database)}")
    print(f"Layout   : {detail}")

    print("\nGenerated index files")
    print("---------------------")
    for path in generated_files:
        print(f"✓ {compact_path(path)}")

    print("\nMFEprimer database is ready for validation.")

    return result
