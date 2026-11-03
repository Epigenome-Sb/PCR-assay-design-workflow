# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""VarVAMP mode selection, parameter validation, and custom configuration."""

from __future__ import annotations

import argparse
import ast
import shutil
from pathlib import Path

from .questions import (
    ask_required_choice,
    ask_yes_no,
    ask_positive_integer_with_default,
    ask_integer_with_default,
    ask_non_negative_integer_with_default,
    ask_integer_keep_default,
    validate_fraction,
)
from .run_management import compact_path, section


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
        # VarVAMP uses infinity as the internal default, but its CLI -n
        # accepts only integers. Omit -n to request the unlimited default.
        if (
            args.single_report_n is not None
            and str(args.single_report_n).lower() != "inf"
        ):
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
