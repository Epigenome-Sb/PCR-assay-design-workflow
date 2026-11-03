
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""Interactive questions and numeric input validation for Workflow 01."""

from __future__ import annotations


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
