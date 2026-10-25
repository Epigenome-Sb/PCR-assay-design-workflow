# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse
import re

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


def sanitize_identifier(value: str, label: str) -> str:
    """Convert a project/run label into a safe directory identifier."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    if not cleaned:
        raise ValueError(f"{label} is empty or invalid.")
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
def portable_project_path(path: Path | str) -> str:
    """
    Serialize paths for manifests and provenance tables.

    Files located inside the current project are stored relative to the
    repository root. Truly external files remain absolute.
    """
    resolved = Path(path).expanduser().resolve()
    repository_root = Path.cwd().resolve()

    try:
        return resolved.relative_to(repository_root).as_posix()
    except ValueError:
        return str(resolved)
def optional_portable_project_path(value: object) -> str:
    """Serialize an optional path while preserving empty/not-applicable values."""
    if value is None or str(value) == "":
        return ""
    return portable_project_path(Path(str(value)))
def _legacy_project_relative_candidate(path: Path) -> Path | None:
    """
    Recover a project-relative suffix from an old absolute path.

    This is mainly for manifests created before portable paths were introduced,
    for example:
        /app/results/<project>/design/<run_id>/varvamp
        /home/user/repo/results/<project>/design/<run_id>/varvamp
    """
    parts = path.parts
    for marker in ("results", "work", "data"):
        if marker in parts:
            index = parts.index(marker)
            return Path(*parts[index:])
    return None
def resolve_manifest_path(value: str | Path, manifest_path: Path) -> Path:
    """
    Resolve a path stored in a Workflow 01 manifest.

    Resolution order:
    1. Existing absolute path.
    2. Legacy absolute project path remapped under the current repository root.
    3. Repository-relative path (schema v2 portable path policy).
    4. Legacy manifest-directory-relative path.

    This makes a design run portable between native execution and Docker while
    retaining compatibility with older manifests that stored /app/... paths.
    """
    repository_root = Path.cwd().resolve()
    path = Path(value).expanduser()

    if path.is_absolute():
        resolved = path.resolve()
        if resolved.exists():
            return resolved

        legacy_relative = _legacy_project_relative_candidate(resolved)
        if legacy_relative is not None:
            remapped = (repository_root / legacy_relative).resolve()
            if remapped.exists():
                return remapped

        return resolved

    repository_candidate = (repository_root / path).resolve()
    if repository_candidate.exists():
        return repository_candidate

    manifest_candidate = (manifest_path.parent / path).resolve()
    if manifest_candidate.exists():
        return manifest_candidate

    return repository_candidate
