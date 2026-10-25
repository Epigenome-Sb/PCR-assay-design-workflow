# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import argparse
import sys

from ..core.paths import compact_path
from ..core.runtime import section
from .cli import ask_integer_with_default, ask_non_negative_integer_with_default, ask_positive_integer_with_default, ask_required_choice, ask_required_float, ask_required_non_negative_integer, ask_yes_no
from .constants import MAFFT_STRATEGIES, VARVAMP_ADVANCED_CONFIG_DEFAULTS, VARVAMP_ADVANCED_CONFIG_ORDER
from .varvamp_config import _validate_min_max, _validate_min_max_opt, next_custom_varvamp_config_path, render_custom_varvamp_config, validate_custom_varvamp_config

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
def progress(label: str, status: str = "OK") -> None:
    dots = "." * max(2, 52 - len(label))
    print(f"{label} {dots} {status}")
