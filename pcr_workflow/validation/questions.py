# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""Interactive questions and validation parameters for Workflow 02."""

from __future__ import annotations


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
