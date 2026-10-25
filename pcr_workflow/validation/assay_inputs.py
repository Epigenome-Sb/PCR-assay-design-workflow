# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import csv
import re

from ..core import runtime
from ..core.runtime import section
from .constants import IUPAC_DNA
from .fasta import read_fasta, write_fasta
from .runtime import progress

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
                "Could not identify the primer-sequence column in primer.tsv. "
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
            "primer_to_amplicon_assignments.tabular."
        )
    return pairs
def prepare_single_or_tiled_files(
    result_dir: Path,
    assay_inputs_dir: Path,
    mode: str,
) -> dict[str, dict[str, str]]:
    primers = read_varvamp_primer_table(result_dir / "primer.tsv")
    pairs = read_varvamp_assignments(
        result_dir / "primer_to_amplicon_assignments.tabular",
        set(primers),
    )
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
    if runtime.VERBOSE:
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
