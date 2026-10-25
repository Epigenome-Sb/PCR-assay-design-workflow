# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import ast
import shutil

from .constants import VARVAMP_ADVANCED_CONFIG_ORDER

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
