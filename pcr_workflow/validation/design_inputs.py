# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""Read Workflow 01 outputs and prepare VarVAMP assays for validation."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from . import run_management as _run_management
from .fasta import IUPAC_DNA, read_fasta, write_fasta
from .questions import ask_yes_no
from .run_management import compact_path, progress, sanitize_identifier, section


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


def _pick_column(fieldnames: list[str], tokens: tuple[str, ...]) -> str | None:
    lowered = {name: name.lower() for name in fieldnames}
    for token in tokens:
        for name, lower in lowered.items():
            if lower == token:
                return name
    for token in tokens:
        for name, lower in lowered.items():
            if token in lower:
                return name
    return None


def read_varvamp_primer_table(primer_tsv: Path) -> dict[str, str]:
    if not primer_tsv.is_file():
        raise FileNotFoundError(f"VarVAMP primer table not found: {primer_tsv}")
    with primer_tsv.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise RuntimeError(f"No header was found in {primer_tsv.name}.")
        id_col = _pick_column(fieldnames, ("primer", "primer_name", "name", "id"))
        seq_col = _pick_column(
            fieldnames,
            ("sequence", "ambiguous_sequence", "primer_sequence", "seq"),
        )
        if seq_col is None:
            raise RuntimeError(
                f"Could not identify the primer-sequence column in {primer_tsv.name}. "
                f"Columns found: {', '.join(fieldnames)}"
            )
        rows = list(reader)
    primers: dict[str, str] = {}
    for row_index, row in enumerate(rows, 1):
        seq = str(row.get(seq_col, "")).strip().upper()
        if not seq or not set(seq) <= IUPAC_DNA:
            continue
        primer_id = str(row.get(id_col, "")).strip() if id_col else ""
        if not primer_id:
            primer_id = f"primer_{row_index}"
        primers[primer_id] = seq
    if not primers:
        raise RuntimeError(f"No usable primer sequences were parsed from {primer_tsv}.")
    return primers


def read_varvamp_assignments(
    assignment_file: Path,
    primer_ids: set[str],
) -> list[tuple[str, str, str]]:
    if not assignment_file.is_file():
        raise FileNotFoundError(
            "VarVAMP primer-to-amplicon assignment file not found: "
            f"{assignment_file}"
        )
    pairs: list[tuple[str, str, str]] = []
    seen: set[tuple[str, str]] = set()
    with assignment_file.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            matched = [
                pid for pid in primer_ids
                if re.search(
                    rf"(?<![A-Za-z0-9_.-]){re.escape(pid)}(?![A-Za-z0-9_.-])",
                    line,
                )
            ]
            matched = sorted(set(matched), key=lambda pid: line.find(pid))
            if len(matched) < 2:
                continue
            left_id, right_id = matched[:2]
            if (left_id, right_id) in seen:
                continue
            seen.add((left_id, right_id))
            first = line.split("\t")[0].strip()
            if first and first not in {left_id, right_id} and first.lower() not in {"amplicon", "primer"}:
                scheme = re.sub(r"[^A-Za-z0-9._-]+", "_", first).strip("_")
            else:
                scheme = f"amplicon_{len(pairs)+1}"
            pairs.append((scheme or f"amplicon_{len(pairs)+1}", left_id, right_id))
    if not pairs:
        raise RuntimeError(
            "No primer pairs could be recovered from "
            f"{assignment_file.name}."
        )
    return pairs



def _resolve_varvamp_file(result_dir: Path, *filenames: str) -> Path:
    """Find a native VarVAMP output, accepting current and legacy filenames."""
    found = [result_dir / name for name in filenames if (result_dir / name).is_file()]
    if not found:
        raise FileNotFoundError(
            f"VarVAMP output not found in {result_dir}. "
            f"Expected one of: {', '.join(filenames)}."
        )
    if len(found) > 1:
        first_content = found[0].read_bytes()
        if any(path.read_bytes() != first_content for path in found[1:]):
            raise RuntimeError(
                "Conflicting VarVAMP output files: "
                + ", ".join(str(path) for path in found)
            )
    return found[0]


def prepare_single_or_tiled_files(
    result_dir: Path,
    assay_inputs_dir: Path,
    mode: str,
) -> dict[str, dict[str, str]]:
    primer_file = _resolve_varvamp_file(result_dir, "primers.tsv", "primer.tsv")
    assignment_file = _resolve_varvamp_file(
        result_dir,
        "primer_to_amplicon_assignment.tabular",
        "primer_to_amplicon_assignments.tabular",
    )
    primers = read_varvamp_primer_table(primer_file)
    pairs = read_varvamp_assignments(assignment_file, set(primers))
    primer_dir = assay_inputs_dir / "primer_pairs"
    scheme_root = assay_inputs_dir / "schemes"
    primer_dir.mkdir(parents=True, exist_ok=True)
    scheme_root.mkdir(parents=True, exist_ok=True)
    schemes: dict[str, dict[str, str]] = {}
    for index, (raw_scheme, left_id, right_id) in enumerate(pairs, 1):
        scheme = re.sub(r"[^A-Za-z0-9._-]+", "_", raw_scheme).strip("_")
        if not scheme or scheme in schemes:
            scheme = f"amplicon_{index}"
        left_seq, right_seq = primers[left_id], primers[right_id]
        schemes[scheme] = {"LEFT": left_seq, "RIGHT": right_seq}
        records = [(f"{scheme}_LEFT", left_seq), (f"{scheme}_RIGHT", right_seq)]
        write_fasta(primer_dir / f"{scheme}_primers.fasta", records)
        write_fasta(scheme_root / scheme / "primers.fasta", records)
    section("VarVAMP assay preparation")
    print(f"Mode         : {mode.upper()}")
    print(f"Primer pairs : {len(schemes)}")
    progress("Preparing primer-pair validation files")
    return schemes


def find_qpcr_oligos(result_dir: Path) -> Path:
    candidate = result_dir / "oligos.fasta"
    if candidate.is_file():
        return candidate.resolve()
    matches = sorted(result_dir.glob("*oligo*.fasta"))
    if matches:
        return matches[0].resolve()
    raise FileNotFoundError(f"No oligos.fasta was found inside {result_dir}.")


def find_oligos_fasta() -> Path:
    cwd = Path.cwd()

    preferred_names = (
        "oligos.fasta",
        "oligos.fa",
        "qpcr_oligos.fasta",
        "varvamp_oligos.fasta",
    )

    for name in preferred_names:
        candidate = cwd / name
        if candidate.is_file():
            return candidate.resolve()

    candidates: list[Path] = []

    for root in (cwd / "results", cwd / "work"):
        if not root.is_dir():
            continue
        for pattern in ("*oligo*.fasta", "*oligo*.fa"):
            candidates.extend(root.rglob(pattern))

    candidates = sorted({p.resolve() for p in candidates if p.is_file()})

    if not candidates:
        raise FileNotFoundError(
            "No oligos FASTA file was found automatically.\n"
            "Use a Workflow 01 assay-design manifest that points to the VarVAMP result directory."
        )

    if len(candidates) == 1:
        return candidates[0]

    print("\nSeveral oligos FASTA files were found:")
    for idx, path in enumerate(candidates, start=1):
        try:
            display = path.relative_to(cwd)
        except ValueError:
            display = path
        print(f"{idx}. {display}")

    while True:
        raw = input("Select the oligos FASTA file by number: ").strip()
        try:
            idx = int(raw)
        except ValueError:
            print("Please enter one of the displayed numbers.")
            continue

        if 1 <= idx <= len(candidates):
            return candidates[idx - 1]

        print("Selection outside the available range.")


def parse_varvamp_header(header: str) -> tuple[str, str]:
    match = re.fullmatch(
        r"(.+)_(LEFT|PROBE|RIGHT)",
        header,
        flags=re.IGNORECASE,
    )
    if not match:
        raise ValueError(
            f"Unrecognized VarVAMP header: {header!r}. "
            "Expected _LEFT, _PROBE or _RIGHT."
        )

    return match.group(1), match.group(2).upper()


def group_schemes(
    records: list[tuple[str, str]],
) -> dict[str, dict[str, str]]:
    schemes: dict[str, dict[str, str]] = {}

    for header, sequence in records:
        scheme, kind = parse_varvamp_header(header)
        schemes.setdefault(scheme, {})

        if kind in schemes[scheme]:
            raise ValueError(f"Duplicate {kind} for scheme {scheme}")

        schemes[scheme][kind] = sequence

    return schemes


def prepare_varvamp_files(
    oligos_file: Path,
    assay_inputs_dir: Path,
) -> dict[str, dict[str, str]]:
    records = read_fasta(oligos_file)
    schemes = group_schemes(records)

    primer_dir = assay_inputs_dir / "primer_pairs"
    probe_dir = assay_inputs_dir / "probes"
    scheme_root = assay_inputs_dir / "schemes"

    primer_dir.mkdir(parents=True, exist_ok=True)
    probe_dir.mkdir(parents=True, exist_ok=True)
    scheme_root.mkdir(parents=True, exist_ok=True)

    all_primers: list[tuple[str, str]] = []
    all_probes: list[tuple[str, str]] = []
    manifest_rows: list[list[object]] = []

    section("VarVAMP assay preparation")
    print("Mode         : QPCR")
    print(f"Schemes      : {len(schemes)}")
    if _run_management.VERBOSE:
        print(f"Oligos FASTA : {oligos_file}")

    for scheme in sorted(schemes):
        oligos = schemes[scheme]
        missing = [kind for kind in ("LEFT", "PROBE", "RIGHT") if kind not in oligos]
        if missing:
            raise ValueError(
                f"Scheme {scheme} is incomplete. Missing: {', '.join(missing)}"
            )

        pair_records = [
            (f"{scheme}_LEFT", oligos["LEFT"]),
            (f"{scheme}_RIGHT", oligos["RIGHT"]),
        ]
        probe_records = [(f"{scheme}_PROBE", oligos["PROBE"])]

        write_fasta(primer_dir / f"{scheme}_primers.fasta", pair_records)
        write_fasta(probe_dir / f"{scheme}_probe.fasta", probe_records)
        all_primers.extend(pair_records)
        all_probes.extend(probe_records)

        one_dir = scheme_root / scheme
        write_fasta(one_dir / "primers.fasta", pair_records)
        write_fasta(one_dir / "probe.fasta", probe_records)
        with (one_dir / "scheme.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["scheme", "oligo_type", "sequence", "length"])
            for kind in ("LEFT", "PROBE", "RIGHT"):
                seq = oligos[kind]
                writer.writerow([scheme, kind, seq, len(seq)])

        manifest_rows.append(
            [
                scheme,
                oligos["LEFT"], oligos["PROBE"], oligos["RIGHT"],
                len(oligos["LEFT"]), len(oligos["PROBE"]), len(oligos["RIGHT"]),
            ]
        )
        print(
            f"  {scheme}: LEFT {len(oligos['LEFT'])} nt | "
            f"PROBE {len(oligos['PROBE'])} nt | RIGHT {len(oligos['RIGHT'])} nt"
        )

    write_fasta(primer_dir / "all_primers.fasta", all_primers)
    write_fasta(probe_dir / "all_probes.fasta", all_probes)
    with (assay_inputs_dir / "manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "scheme", "left_sequence", "probe_sequence", "right_sequence",
                "left_length", "probe_length", "right_length",
            ]
        )
        writer.writerows(manifest_rows)

    return schemes
