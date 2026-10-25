# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import json

from ..core.paths import compact_path, resolve_manifest_path, sanitize_identifier
from .ui import ask_yes_no

def load_design_manifest(manifest_path: Path) -> dict[str, object]:
    """Load and validate the Workflow 01 handoff manifest."""
    manifest_path = manifest_path.expanduser().resolve()
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Design manifest not found: {manifest_path}")

    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(
            f"Invalid JSON in design manifest {manifest_path}: {error}"
        ) from error

    if manifest.get("workflow_stage") not in {None, "design"}:
        raise ValueError(
            f"Manifest is not an assay-design manifest: {manifest_path}"
        )

    project_name = manifest.get("project_name")
    if not isinstance(project_name, str) or not project_name.strip():
        raise ValueError("Design manifest does not contain a valid project_name.")

    varvamp = manifest.get("varvamp")
    if not isinstance(varvamp, dict):
        raise ValueError("Design manifest does not contain a VarVAMP section.")

    mode = str(varvamp.get("mode") or "").lower()
    if mode not in {"single", "qpcr", "tiled"}:
        raise ValueError(
            "Design manifest does not identify a completed SINGLE/QPCR/TILED run."
        )

    result_value = varvamp.get("result_dir")
    if not result_value:
        result_directories = manifest.get("result_directories")
        if isinstance(result_directories, dict):
            result_value = result_directories.get("varvamp")
    if not result_value:
        raise ValueError("Design manifest does not contain a VarVAMP result directory.")

    result_dir = resolve_manifest_path(str(result_value), manifest_path)
    if not result_dir.is_dir():
        raise FileNotFoundError(
            "VarVAMP result directory referenced by the design manifest does not "
            f"exist: {result_dir}"
        )

    manifest["_manifest_path"] = str(manifest_path)
    manifest["_resolved_varvamp_result_dir"] = str(result_dir)
    manifest["_resolved_mode"] = mode
    return manifest
def find_design_manifests() -> list[Path]:
    """Find run-level Workflow 01 manifests without selecting nested copies."""
    root = Path.cwd() / "results"
    if not root.is_dir():
        return []
    candidates = [
        path.resolve()
        for path in root.glob("*/design/*/assay_design_manifest.json")
        if path.is_file()
    ]
    return sorted(
        set(candidates),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
def choose_design_run(
    requested_manifest: Path | None,
) -> tuple[str, Path, Path, dict[str, object], str, str]:
    """
    Select a Workflow 01 run through assay_design_manifest.json.

    Returns mode, VarVAMP result directory, manifest path, manifest data,
    project name and design run ID.
    """
    if requested_manifest is not None:
        paths = [requested_manifest.expanduser().resolve()]
    else:
        paths = find_design_manifests()
        if not paths:
            raise FileNotFoundError(
                "No assay_design_manifest.json was found under "
                "results/<project>/design/<run_id>/. Run Workflow 01 first or "
                "provide --design-manifest."
            )

    valid: list[tuple[Path, dict[str, object]]] = []
    errors: list[str] = []
    for path in paths:
        try:
            valid.append((path, load_design_manifest(path)))
        except (FileNotFoundError, ValueError) as error:
            errors.append(f"{path}: {error}")

    if not valid:
        raise RuntimeError(
            "No usable Workflow 01 design manifest was found.\n" + "\n".join(errors)
        )

    print("\n1. Assay-design run")
    print("===================")

    if requested_manifest is None:
        for idx, (path, manifest) in enumerate(valid, start=1):
            mode = str(manifest["_resolved_mode"]).upper()
            project = str(manifest["project_name"])
            design_run = str(manifest.get("run_id") or path.parent.name)
            design_input = Path(str(manifest.get("input_fasta", "unknown"))).name
            print(
                f"{idx}. {project} | {mode} | design run {design_run} | "
                f"input {design_input}"
            )

        while True:
            raw = input("Select the assay-design run by number: ").strip()
            try:
                selected = int(raw)
            except ValueError:
                print("Please enter one of the displayed numbers.")
                continue
            if 1 <= selected <= len(valid):
                manifest_path, manifest = valid[selected - 1]
                break
            print("Selection outside the available range.")
    else:
        manifest_path, manifest = valid[0]

    mode = str(manifest["_resolved_mode"])
    result_dir = Path(str(manifest["_resolved_varvamp_result_dir"]))
    project_name = sanitize_identifier(str(manifest["project_name"]), "Project name")
    design_run_id = sanitize_identifier(
        str(manifest.get("run_id") or manifest_path.parent.name),
        "Design run ID",
    )

    print(f"Project     : {project_name}")
    print(f"Design run  : {design_run_id}")
    print(f"Mode        : {mode.upper()}")
    print(f"VarVAMP dir : {compact_path(result_dir)}")
    print(f"Manifest    : {compact_path(manifest_path)}")

    if requested_manifest is None and not ask_yes_no("Use this design run?"):
        raise RuntimeError("Assay-design run selection cancelled.")

    return mode, result_dir, manifest_path, manifest, project_name, design_run_id
def find_varvamp_result_directories() -> list[tuple[str, Path]]:
    """Compatibility helper backed by Workflow 01 run manifests."""
    found: list[tuple[str, Path]] = []
    for manifest_path in find_design_manifests():
        try:
            manifest = load_design_manifest(manifest_path)
        except (FileNotFoundError, ValueError):
            continue
        found.append(
            (
                str(manifest["_resolved_mode"]),
                Path(str(manifest["_resolved_varvamp_result_dir"])),
            )
        )
    return found
def choose_varvamp_result() -> tuple[str, Path]:
    """Compatibility wrapper; new code should use choose_design_run()."""
    mode, result_dir, _, _, _, _ = choose_design_run(None)
    return mode, result_dir
