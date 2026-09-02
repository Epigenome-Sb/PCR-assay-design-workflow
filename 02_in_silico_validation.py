#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""
02_in_silico_validation.py

Mode-aware in silico validation for VarVAMP assay designs.

Branches
--------
SINGLE: validate assigned primer pairs with MFEprimer.
QPCR: validate LEFT/RIGHT with MFEprimer, extract valid amplicons, expand any
existing IUPAC ambiguity in the VarVAMP probe into concrete A/C/G/T variants,
and search those variants with blastn-short. No new ambiguity is proposed.
TILED: validate assigned primer pairs with MFEprimer; no probe BLAST.

Normal terminal mode is concise. Use --verbose to show commands and detailed
paths. All external commands are logged to validation_inputs/validation.log.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

IUPAC_DNA = set("ACGTRYSWKMBDHVN")

VERBOSE = False
VALIDATION_LOG: Path | None = None

def parse_cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Mode-aware in silico validation of VarVAMP SINGLE, QPCR or TILED "
            "assay designs with MFEprimer and, for QPCR, BLAST+ probe analysis."
        )
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show external commands and detailed paths in the terminal.",
    )
    return parser.parse_args()

def log_line(message: str) -> None:
    if VALIDATION_LOG is None:
        return
    VALIDATION_LOG.parent.mkdir(parents=True, exist_ok=True)
    with VALIDATION_LOG.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")

def compact_path(path: Path | str) -> str:
    path_obj = Path(path)
    try:
        return str(path_obj.resolve().relative_to(Path.cwd().resolve()))
    except Exception:
        return str(path_obj)

def section(title: str) -> None:
    print(f"\n{title}")
    print("─" * len(title))

def progress(label: str, status: str = "OK") -> None:
    dots = "." * max(2, 48 - len(label))
    print(f"{label} {dots} {status}")


# ---------------------------------------------------------------------------
# Generic terminal helpers
# ---------------------------------------------------------------------------

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


def run_command(command: list[str]) -> None:
    command_text = " ".join(command)
    log_line("COMMAND: " + command_text)
    if VERBOSE:
        print("\nCommand:", command_text)
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"Command failed with exit code {error.returncode}: " + command_text
        ) from error


def run_command_to_file(command: list[str], output_file: Path) -> None:
    command_text = " ".join(command) + " > " + str(output_file)
    log_line("COMMAND: " + command_text)
    if VERBOSE:
        print("\nCommand:", command_text)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with output_file.open("w", encoding="utf-8") as handle:
            subprocess.run(command, check=True, stdout=handle)
    except subprocess.CalledProcessError as error:
        output_file.unlink(missing_ok=True)
        raise RuntimeError(
            f"Command failed with exit code {error.returncode}: " + command_text
        ) from error


# ---------------------------------------------------------------------------
# FASTA handling
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> list[tuple[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    records: list[tuple[str, str]] = []
    current: str | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal current, seq_parts
        if current is None:
            return

        sequence = "".join(seq_parts).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"Empty sequence for FASTA record: {current}")

        invalid = sorted(set(sequence) - IUPAC_DNA)
        if invalid:
            raise ValueError(
                f"Invalid nucleotide character(s) in {current}: "
                + ", ".join(invalid)
            )

        records.append((current, sequence))

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">"):
                flush()
                current = line[1:].split()[0]
                seq_parts = []
            else:
                if current is None:
                    raise ValueError(
                        "Sequence encountered before the first FASTA header."
                    )
                seq_parts.append(line)

    flush()

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def fasta_statistics(path: Path) -> dict[str, float | int]:
    records = read_fasta(path)
    lengths = [len(seq) for _, seq in records]
    ambiguous = sum(
        sum(base not in {"A", "C", "G", "T"} for base in seq)
        for _, seq in records
    )
    return {
        "sequences": len(records),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": sum(lengths) / len(lengths),
        "ambiguous": ambiguous,
    }


def print_fasta_statistics(path: Path, title: str) -> None:
    stats = fasta_statistics(path)
    print(f"\n{title}")
    print("=" * len(title))
    print(f"FASTA:                 {path}")
    print(f"Sequences:             {stats['sequences']}")
    print(
        f"Sequence length range: {stats['min_length']}-"
        f"{stats['max_length']} nt"
    )
    print(f"Mean sequence length:  {stats['mean_length']:.2f} nt")
    print(f"Ambiguous bases:       {stats['ambiguous']}")



# ---------------------------------------------------------------------------
# VarVAMP result-mode detection
# ---------------------------------------------------------------------------

def find_varvamp_result_directories() -> list[tuple[str, Path]]:
    found: list[tuple[str, Path]] = []
    root = Path.cwd() / "results"
    if not root.is_dir():
        return found
    for mode in ("single", "qpcr", "tiled"):
        for directory in root.rglob(f"varvamp_{mode}"):
            if directory.is_dir():
                found.append((mode, directory.resolve()))
    unique = []
    seen = set()
    for mode, directory in sorted(found, key=lambda x: (x[0], str(x[1]))):
        key = (mode, str(directory))
        if key not in seen:
            seen.add(key)
            unique.append((mode, directory))
    return unique


def choose_varvamp_result() -> tuple[str, Path]:
    candidates = find_varvamp_result_directories()
    if not candidates:
        raise FileNotFoundError(
            "No VarVAMP result directory was found under results/. Expected "
            "varvamp_single, varvamp_qpcr or varvamp_tiled."
        )
    if len(candidates) == 1:
        mode, directory = candidates[0]
        section("VarVAMP design detected")
        print(f"Mode       : {mode.upper()}")
        print(f"Result dir : {compact_path(directory)}")
        if ask_yes_no("Use this VarVAMP design for validation?"):
            return mode, directory
        raise RuntimeError("VarVAMP design selection cancelled.")
    section("VarVAMP designs detected")
    for idx, (mode, directory) in enumerate(candidates, 1):
        print(f"{idx}. {mode.upper():6s}  {compact_path(directory)}")
    while True:
        raw = input("Select the VarVAMP design by number: ").strip()
        try:
            idx = int(raw)
        except ValueError:
            print("Please enter one of the displayed numbers.")
            continue
        if 1 <= idx <= len(candidates):
            return candidates[idx - 1]
        print("Selection outside the available range.")


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
    validation_root: Path,
    mode: str,
) -> dict[str, dict[str, str]]:
    primers = read_varvamp_primer_table(result_dir / "primer.tsv")
    pairs = read_varvamp_assignments(
        result_dir / "primer_to_amplicon_assignments.tabular",
        set(primers),
    )
    primer_dir = validation_root / "primer_pairs"
    scheme_root = validation_root / "schemes"
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

# ---------------------------------------------------------------------------
# VarVAMP oligo preparation
# ---------------------------------------------------------------------------

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
            "Place oligos.fasta in the current directory or inside results/ or work/."
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
    root: Path,
) -> dict[str, dict[str, str]]:
    records = read_fasta(oligos_file)
    schemes = group_schemes(records)

    primer_dir = root / "primer_pairs"
    probe_dir = root / "probes"
    scheme_root = root / "schemes"

    primer_dir.mkdir(parents=True, exist_ok=True)
    probe_dir.mkdir(parents=True, exist_ok=True)
    scheme_root.mkdir(parents=True, exist_ok=True)

    all_primers: list[tuple[str, str]] = []
    all_probes: list[tuple[str, str]] = []
    manifest_rows: list[list[object]] = []

    section("VarVAMP assay preparation")
    print("Mode         : QPCR")
    print(f"Schemes      : {len(schemes)}")
    if VERBOSE:
        print(f"Oligos FASTA : {oligos_file}")

    for scheme in sorted(schemes):
        oligos = schemes[scheme]
        missing = [
            kind
            for kind in ("LEFT", "PROBE", "RIGHT")
            if kind not in oligos
        ]

        if missing:
            raise ValueError(
                f"Scheme {scheme} is incomplete. Missing: {', '.join(missing)}"
            )

        pair_records = [
            (f"{scheme}_LEFT", oligos["LEFT"]),
            (f"{scheme}_RIGHT", oligos["RIGHT"]),
        ]
        probe_records = [
            (f"{scheme}_PROBE", oligos["PROBE"]),
        ]

        write_fasta(
            primer_dir / f"{scheme}_primers.fasta",
            pair_records,
        )
        write_fasta(
            probe_dir / f"{scheme}_probe.fasta",
            probe_records,
        )

        all_primers.extend(pair_records)
        all_probes.extend(probe_records)

        one_dir = scheme_root / scheme
        write_fasta(one_dir / "primers.fasta", pair_records)
        write_fasta(one_dir / "probe.fasta", probe_records)

        with (one_dir / "scheme.tsv").open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                ["scheme", "oligo_type", "sequence", "length"]
            )
            for kind in ("LEFT", "PROBE", "RIGHT"):
                seq = oligos[kind]
                writer.writerow([scheme, kind, seq, len(seq)])

        manifest_rows.append(
            [
                scheme,
                oligos["LEFT"],
                oligos["PROBE"],
                oligos["RIGHT"],
                len(oligos["LEFT"]),
                len(oligos["PROBE"]),
                len(oligos["RIGHT"]),
            ]
        )

        print(
            f"  {scheme}: LEFT {len(oligos['LEFT'])} nt | "
            f"PROBE {len(oligos['PROBE'])} nt | "
            f"RIGHT {len(oligos['RIGHT'])} nt"
        )

    write_fasta(primer_dir / "all_primers.fasta", all_primers)
    write_fasta(probe_dir / "all_probes.fasta", all_probes)

    with (root / "manifest.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "scheme",
                "left_sequence",
                "probe_sequence",
                "right_sequence",
                "left_length",
                "probe_length",
                "right_length",
            ]
        )
        writer.writerows(manifest_rows)

    return schemes


# ---------------------------------------------------------------------------
# Target database selection and MFEprimer indexing
# ---------------------------------------------------------------------------

def is_inside(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def find_candidate_target_fastas(validation_root: Path) -> list[Path]:
    cwd = Path.cwd()

    excluded_roots = [
        validation_root.resolve(),
        (cwd / "work").resolve(),
        (cwd / "results").resolve(),
    ]

    candidates: list[Path] = []

    # Current directory + common data directories.
    search_roots = [cwd]
    for name in ("data", "database", "databases", "validation_data"):
        candidate = cwd / name
        if candidate.is_dir():
            search_roots.append(candidate)

    patterns = ("*.fasta", "*.fa", "*.fna")

    for root in search_roots:
        if not root.exists():
            continue

        if root == cwd:
            iterator = []
            for pattern in patterns:
                iterator.extend(root.glob(pattern))
        else:
            iterator = []
            for pattern in patterns:
                iterator.extend(root.rglob(pattern))

        for path in iterator:
            if not path.is_file():
                continue

            resolved = path.resolve()

            if any(is_inside(resolved, excluded) for excluded in excluded_roots):
                continue

            lower_name = path.name.lower()
            if any(
                token in lower_name
                for token in ("oligo", "primer", "probe", "consensus")
            ):
                continue

            candidates.append(resolved)

    return sorted(set(candidates))


def choose_target_database(validation_root: Path) -> Path:
    candidates = find_candidate_target_fastas(validation_root)

    print("\n2. Target coverage database")
    print("==========================")

    if candidates:
        print(
            "FASTA databases detected locally "
            "(generated primer/probe files are excluded):"
        )
        for idx, path in enumerate(candidates, start=1):
            try:
                display = path.relative_to(Path.cwd())
            except ValueError:
                display = path
            print(f"{idx}. {display}")

        print(f"{len(candidates) + 1}. Enter another FASTA path")

        while True:
            raw = input("Select the target database by number: ").strip()
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
        raw = input("Name or path of the target validation FASTA: ").strip()
        raw = raw.strip("'\"")
        candidate = Path(raw).expanduser()

        paths = [candidate]
        if not candidate.is_absolute():
            paths.extend(
                [
                    Path.cwd() / candidate,
                    Path.cwd() / "data" / candidate,
                    Path.cwd() / "database" / candidate,
                    Path.cwd() / "databases" / candidate,
                ]
            )

        for path in paths:
            if path.is_file():
                return path.resolve()

        print("FASTA file not found.")


def mfeprimer_index_status(database: Path) -> tuple[bool, str]:
    """
    Support current MFEprimer 4.x index and legacy 3.x index layouts.

    Current 4.5.x: <db>.primerqc.bin
    Legacy 3.x:    <db>.primerqc + <db>.primerqc.fai
    """
    current_binary = Path(str(database) + ".primerqc.bin")
    legacy_index = Path(str(database) + ".primerqc")
    legacy_fai = Path(str(database) + ".primerqc.fai")

    if current_binary.is_file():
        return True, "current binary index (.primerqc.bin)"

    if legacy_index.is_file() and legacy_fai.is_file():
        return True, "legacy MFEprimer index (.primerqc + .primerqc.fai)"

    return False, "not indexed"


def ensure_validation_tools_available(assay_mode: str) -> None:
    """Check tools required for the selected validation branch."""
    tools = ["seqkit", "mfeprimer"]
    if assay_mode == "qpcr":
        tools.append("blastn")
    missing = [tool for tool in tools if shutil.which(tool) is None]

    if missing:
        raise RuntimeError(
            "Required validation tool(s) not found in PATH: "
            + ", ".join(missing)
            + ". Install/activate them before running coverage validation."
        )


def format_target_database_for_mfeprimer(
    database: Path,
    validation_root: Path,
) -> Path:
    """
    Reformat the target FASTA so that every sequence is written on one line.

    This preserves FASTA headers and nucleotide/IUPAC characters and only
    changes line wrapping, using:

        seqkit seq -w 0 input.fasta > input_fixed.fasta
    """
    database_dir = validation_root / "databases"
    database_dir.mkdir(parents=True, exist_ok=True)

    safe_stem = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        database.stem,
    )
    fixed_database = database_dir / f"{safe_stem}_fixed.fasta"

    print("\n3. MFEprimer FASTA formatting")
    print("=============================")
    print(f"Original database:  {database}")
    print(f"Formatted database: {fixed_database}")
    print(
        "Formatting rule: one complete sequence per FASTA line "
        "(SeqKit -w 0)."
    )
    print(
        "Headers and nucleotide/IUPAC characters are preserved."
    )

    run_command_to_file(
        [
            "seqkit",
            "seq",
            "-w",
            "0",
            str(database),
        ],
        fixed_database,
    )

    original_stats = fasta_statistics(database)
    fixed_stats = fasta_statistics(fixed_database)

    print("\nFASTA formatting validation")
    print("===========================")
    print(
        f"Sequences before: {original_stats['sequences']} | "
        f"after: {fixed_stats['sequences']}"
    )
    print(
        f"Mean length before: {original_stats['mean_length']:.2f} nt | "
        f"after: {fixed_stats['mean_length']:.2f} nt"
    )
    print(
        f"Ambiguous bases before: {original_stats['ambiguous']} | "
        f"after: {fixed_stats['ambiguous']}"
    )

    if (
        original_stats["sequences"] != fixed_stats["sequences"]
        or original_stats["min_length"] != fixed_stats["min_length"]
        or original_stats["max_length"] != fixed_stats["max_length"]
        or original_stats["ambiguous"] != fixed_stats["ambiguous"]
    ):
        raise RuntimeError(
            "SeqKit formatting changed FASTA content unexpectedly. "
            "The formatted database will not be used."
        )

    print("FASTA content preserved; only sequence line wrapping was changed.")

    return fixed_database


def index_target_database(database: Path) -> None:
    indexed, detail = mfeprimer_index_status(database)

    print("\nMFEprimer database index")
    print("========================")
    print(f"Database: {database}")
    print(f"Status:   {detail}")

    if indexed:
        print("The target database is already indexed.")
        return

    if not ask_yes_no("Index this target database with MFEprimer now?"):
        raise RuntimeError(
            "Coverage validation requires an indexed MFEprimer database."
        )

    index_mode = ask_choice(
        "MFEprimer indexing strategy",
        [
            "Use MFEprimer default indexing parameters.",
            (
                "Use a custom k-mer seed length (-k). "
                "Mainly useful for very large/repeat-rich genomes."
            ),
        ],
    )

    command = ["mfeprimer", "index", "-i", str(database)]

    if index_mode == 2:
        k_value = ask_int(
            "MFEprimer index k-mer length (-k)",
            9,
            15,
        )
        command.extend(["-k", str(k_value)])

    run_command(command)

    indexed, detail = mfeprimer_index_status(database)
    print(f"\nIndex status after command: {detail}")

    if not indexed:
        print(
            "Warning: index files were not recognized by this script. "
            "MFEprimer may use a different index layout in your installed version."
        )


# ---------------------------------------------------------------------------
# Pair selection and MFEprimer execution
# ---------------------------------------------------------------------------

def choose_schemes(schemes: dict[str, dict[str, str]]) -> list[str]:
    scheme_names = sorted(schemes)

    print("\n5. Primer-pair selection")
    print("========================")

    for idx, scheme in enumerate(scheme_names, start=1):
        probe_text = (
            f" / PROBE {len(schemes[scheme]['PROBE'])} nt"
            if "PROBE" in schemes[scheme]
            else ""
        )
        print(
            f"{idx}. {scheme} "
            f"(LEFT {len(schemes[scheme]['LEFT'])} nt / "
            f"RIGHT {len(schemes[scheme]['RIGHT'])} nt{probe_text})"
        )

    print(f"{len(scheme_names) + 1}. Test ALL primer pairs")

    while True:
        raw = input(
            "Select one/multiple numbers (e.g. 1,3) or ALL option: "
        ).strip()

        try:
            numbers = [
                int(item.strip())
                for item in raw.split(",")
                if item.strip()
            ]
        except ValueError:
            print("Use numbers separated by commas, e.g. 1,2.")
            continue

        if not numbers:
            print("At least one selection is required.")
            continue

        all_option = len(scheme_names) + 1

        if all_option in numbers:
            if len(numbers) > 1:
                print("Choose ALL alone, or specific pair numbers.")
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
    mode = ask_choice(
        "MFEprimer amplicon-filter parameters",
        [
            (
                "Use MFEprimer defaults "
                "(min product 0 bp, max product 2000 bp, Tm cutoff 30 °C)."
            ),
            "Set custom amplicon size and Tm thresholds.",
        ],
    )

    if mode == 1:
        return None

    min_size = ask_int("Minimum predicted amplicon size (-s)", 0, 100000)
    max_size = ask_int(
        "Maximum predicted amplicon size (-S)",
        min_size,
        100000,
    )
    tm_cutoff = ask_float(
        "Minimum amplicon/primer Tm cutoff (-t)",
        0.0,
        100.0,
    )

    return {
        "min_size": min_size,
        "max_size": max_size,
        "tm_cutoff": tm_cutoff,
    }



# ---------------------------------------------------------------------------
# MFEprimer 3.x report parsing and target-amplicon extraction
# ---------------------------------------------------------------------------

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


def normalize_primer_identifier(value: str) -> str:
    """
    Normalize an identifier for robust comparisons.

    Examples:
        varVAMP_0      -> varvamp0
        VARVAMP-0      -> varvamp0
        varVAMP_0_LEFT -> varvamp0left
    """
    return re.sub(
        r"[^A-Za-z0-9]+",
        "",
        value,
    ).lower()


def parse_primer_identity(
    primer_name: str,
) -> tuple[str | None, str | None]:
    """
    Extract the VarVAMP scheme and LEFT/RIGHT role from a MFEprimer
    primer identifier.

    This intentionally tolerates harmless decorations added around the
    FASTA identifier, for example:
        varVAMP_0_LEFT
        varVAMP_0_LEFT(+)
        varVAMP_0-LEFT
        varVAMP_0_LEFT_1

    LEFT/RIGHT must still appear as a real token, not as part of another
    ordinary word.
    """
    raw = primer_name.strip().lstrip(">")

    role_match = re.search(
        r"(^|[^A-Za-z0-9])(LEFT|RIGHT)(?=$|[^A-Za-z0-9])",
        raw,
        flags=re.IGNORECASE,
    )

    if not role_match:
        return None, None

    role = role_match.group(2).upper()

    # Everything before the LEFT/RIGHT token is treated as the scheme
    # component. Remove punctuation/separators surrounding the token.
    scheme = raw[:role_match.start(2)]
    scheme = re.sub(
        r"[^A-Za-z0-9]+$",
        "",
        scheme,
    ).strip()

    if not scheme:
        return None, role

    return scheme, role


def same_scheme(
    observed: str | None,
    expected: str,
) -> bool:
    """Compare scheme identifiers while ignoring separators and case."""
    if observed is None:
        return False

    return (
        normalize_primer_identifier(observed)
        == normalize_primer_identifier(expected)
    )


def classify_mfeprimer_product(
    primer_1: str,
    primer_2: str,
    expected_scheme: str,
) -> tuple[str, str, str]:
    """
    Classify a predicted MFEprimer product.

    TARGET:
        one LEFT and one RIGHT primer from the selected VarVAMP scheme.

    SELF_PRIMING:
        LEFT+LEFT or RIGHT+RIGHT from the selected scheme.

    CROSS_SCHEME:
        LEFT/RIGHT roles are recognized, but the primers do not belong to
        the same scheme.

    OTHER:
        at least one raw primer identifier cannot be interpreted safely.
    """
    scheme_1, role_1 = parse_primer_identity(primer_1)
    scheme_2, role_2 = parse_primer_identity(primer_2)

    if role_1 is None or role_2 is None:
        return (
            "OTHER",
            role_1 or "UNKNOWN",
            role_2 or "UNKNOWN",
        )

    # If roles were recognized but either scheme prefix was not recoverable,
    # do not guess.
    if scheme_1 is None or scheme_2 is None:
        return "OTHER", role_1, role_2

    if normalize_primer_identifier(
        scheme_1
    ) != normalize_primer_identifier(
        scheme_2
    ):
        return "CROSS_SCHEME", role_1, role_2

    if not same_scheme(
        scheme_1,
        expected_scheme,
    ):
        return "OTHER", role_1, role_2

    roles = {role_1, role_2}

    if roles == {"LEFT", "RIGHT"}:
        return "TARGET", role_1, role_2

    if (
        role_1 == role_2
        and role_1 in {"LEFT", "RIGHT"}
    ):
        return "SELF_PRIMING", role_1, role_2

    return "OTHER", role_1, role_2



def parse_mfeprimer3_report(report_file: Path) -> dict[str, object]:
    """
    Parse a legacy MFEprimer 3.x specificity text report.

    The parser deliberately avoids fixed-width column positions.

    It uses two independent structures:
    1. the AmpID/HitID/Size summary table;
    2. FASTA-like blocks beginning with:
           >Amp_N PRIMER1 + PRIMER2 ==> HitID

    The sequence block is considered authoritative for extraction, while
    the summary table is used for consistency checks.
    """
    if not report_file.is_file():
        raise FileNotFoundError(
            f"MFEprimer report not found: {report_file}"
        )

    lines = report_file.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    expected_count: int | None = None

    for line in lines:
        match = re.search(
            r"Descriptions\s+of\s+\[\s*(\d+)\s*\]\s+potential\s+amplicons",
            line,
            flags=re.IGNORECASE,
        )
        if match:
            expected_count = int(match.group(1))
            break

    # Parse the summary table without relying on fixed-width spacing.
    summary_by_ampid: dict[int, dict[str, object]] = {}
    in_table = False

    for line in lines:
        stripped = line.strip()

        if re.match(
            r"^AmpID\s+HitID\s+Size\b",
            stripped,
            flags=re.IGNORECASE,
        ):
            in_table = True
            continue

        if in_table and re.match(
            r"^Amplicon\s+details\b",
            stripped,
            flags=re.IGNORECASE,
        ):
            break

        if not in_table:
            continue

        row = re.match(
            r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+",
            line,
        )
        if not row:
            continue

        amp_id = int(row.group(1))
        hit_id = row.group(2)
        size = int(row.group(3))

        if amp_id in summary_by_ampid:
            raise RuntimeError(
                f"Duplicate AmpID {amp_id} in MFEprimer summary table."
            )

        summary_by_ampid[amp_id] = {
            "hit_id": hit_id,
            "size": size,
        }

    # Parse FASTA-like amplicon blocks.
    header_re = re.compile(
        r"^>Amp_(\d+)\s+(\S+)\s+\+\s+(\S+)\s+==>\s+(\S+)"
    )
    sequence_re = re.compile(
        r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$"
    )

    parsed_blocks: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    sequence_parts: list[str] = []

    def flush_current() -> None:
        nonlocal current, sequence_parts

        if current is None:
            return

        sequence = "".join(sequence_parts).upper()

        if not sequence:
            raise RuntimeError(
                f"No amplicon sequence found for Amp_{current['amp_id']}."
            )

        current["sequence"] = sequence
        current["sequence_length"] = len(sequence)
        parsed_blocks.append(current)

        current = None
        sequence_parts = []

    for raw_line in lines:
        stripped = raw_line.strip()

        header_match = header_re.match(stripped)

        if header_match:
            flush_current()

            current = {
                "amp_id": int(header_match.group(1)),
                "primer_1": header_match.group(2),
                "primer_2": header_match.group(3),
                "hit_id": header_match.group(4),
            }
            sequence_parts = []
            continue

        if current is None:
            continue

        if stripped and sequence_re.fullmatch(stripped):
            sequence_parts.append(stripped)
            continue

        # A blank line or a non-sequence line terminates the FASTA block.
        if sequence_parts:
            flush_current()

    flush_current()

    # If no amplicons were reported, allow a clean empty result.
    if expected_count == 0 and not parsed_blocks:
        return {
            "expected_count": 0,
            "summary_by_ampid": summary_by_ampid,
            "amplicons": [],
        }

    if not parsed_blocks:
        raise RuntimeError(
            "No '>Amp_N ... ==> HitID' amplicon sequence blocks were found "
            "in the MFEprimer report. This parser expects the legacy "
            "MFEprimer 3.x text report."
        )

    amp_ids = [int(item["amp_id"]) for item in parsed_blocks]

    if len(amp_ids) != len(set(amp_ids)):
        raise RuntimeError(
            "Duplicate AmpID values were found in the MFEprimer sequence blocks."
        )

    if (
        expected_count is not None
        and len(parsed_blocks) != expected_count
    ):
        raise RuntimeError(
            "MFEprimer report consistency check failed: "
            f"report announces {expected_count} potential amplicons, "
            f"but {len(parsed_blocks)} sequence blocks were parsed."
        )

    # Cross-check HitID and Size between the table and sequence blocks.
    for item in parsed_blocks:
        amp_id = int(item["amp_id"])
        table_item = summary_by_ampid.get(amp_id)

        if table_item is None:
            raise RuntimeError(
                f"Amp_{amp_id} is present in amplicon details but missing "
                "from the AmpID/HitID/Size table."
            )

        if str(table_item["hit_id"]) != str(item["hit_id"]):
            raise RuntimeError(
                f"HitID mismatch for Amp_{amp_id}: "
                f"table={table_item['hit_id']} vs "
                f"details={item['hit_id']}."
            )

        reported_size = int(table_item["size"])
        sequence_length = int(item["sequence_length"])

        if reported_size != sequence_length:
            raise RuntimeError(
                f"Size mismatch for Amp_{amp_id}: "
                f"table reports {reported_size} bp but extracted sequence "
                f"contains {sequence_length} nt."
            )

        item["size"] = reported_size

    return {
        "expected_count": expected_count,
        "summary_by_ampid": summary_by_ampid,
        "amplicons": parsed_blocks,
    }


def analyse_and_extract_mfeprimer_amplicons(
    report_file: Path,
    scheme: str,
    scheme_dir: Path,
    total_target_sequences: int,
    min_target_size: int | None,
    max_target_size: int | None,
) -> dict[str, object]:
    """
    Classify MFEprimer products, extract valid target amplicons and calculate
    unique primer-pair target coverage.
    """
    parsed = parse_mfeprimer3_report(report_file)
    amplicons = list(parsed["amplicons"])

    all_rows: list[dict[str, object]] = []
    valid_rows: list[dict[str, object]] = []
    non_target_rows: list[dict[str, object]] = []

    for item in amplicons:
        classification, role_1, role_2 = classify_mfeprimer_product(
            str(item["primer_1"]),
            str(item["primer_2"]),
            scheme,
        )

        size = int(item["size"])

        if min_target_size is None or max_target_size is None:
            size_pass = True
        else:
            size_pass = min_target_size <= size <= max_target_size

        coverage_eligible = (
            classification == "TARGET"
            and size_pass
        )

        row = {
            "AmpID": int(item["amp_id"]),
            "HitID": str(item["hit_id"]),
            "Primer1": str(item["primer_1"]),
            "Primer2": str(item["primer_2"]),
            "Role1": role_1,
            "Role2": role_2,
            "ProductClass": classification,
            "Size": size,
            "SizePass": "YES" if size_pass else "NO",
            "CoverageEligible": "YES" if coverage_eligible else "NO",
            "Sequence": str(item["sequence"]).upper(),
        }

        all_rows.append(row)

        if coverage_eligible:
            valid_rows.append(row)

        if classification != "TARGET":
            non_target_rows.append(row)

    # Stable first-occurrence order, without double-counting a HitID.
    unique_positive_hitids = list(
        dict.fromkeys(
            str(row["HitID"])
            for row in valid_rows
        )
    )

    valid_amplicons = len(valid_rows)
    unique_positive_targets = len(unique_positive_hitids)
    duplicate_target_products = (
        valid_amplicons - unique_positive_targets
    )

    if total_target_sequences <= 0:
        raise ValueError(
            "Total target sequence count must be greater than zero."
        )

    coverage_pct = (
        100.0
        * unique_positive_targets
        / total_target_sequences
    )

    # Output files.
    all_tsv = scheme_dir / "amplicons_all.tsv"
    valid_tsv = scheme_dir / "amplicons_valid.tsv"
    valid_fasta = scheme_dir / "amplicons_valid.fasta"
    hitids_file = scheme_dir / "positive_hitids.txt"
    non_target_tsv = scheme_dir / "non_target_products.tsv"
    summary_file = scheme_dir / "coverage_summary.txt"

    tsv_columns = [
        "AmpID",
        "HitID",
        "Primer1",
        "Primer2",
        "Role1",
        "Role2",
        "ProductClass",
        "Size",
        "SizePass",
        "CoverageEligible",
    ]

    def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
        with path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=tsv_columns,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)

    write_rows(all_tsv, all_rows)
    write_rows(valid_tsv, valid_rows)
    write_rows(non_target_tsv, non_target_rows)

    with valid_fasta.open("w", encoding="utf-8") as handle:
        for row in valid_rows:
            pair = f"{row['Role1']}+{row['Role2']}"
            handle.write(
                f">Amp_{row['AmpID']}|HitID={row['HitID']}|"
                f"Scheme={scheme}|Pair={pair}|Size={row['Size']}\n"
            )
            sequence = str(row["Sequence"])
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    hitids_file.write_text(
        "".join(f"{hit_id}\n" for hit_id in unique_positive_hitids),
        encoding="utf-8",
    )

    class_counts: dict[str, int] = {}
    for row in all_rows:
        key = str(row["ProductClass"])
        class_counts[key] = class_counts.get(key, 0) + 1

    unrecognized_examples: list[str] = []
    seen_unrecognized: set[tuple[str, str, str]] = set()

    for row in all_rows:
        product_class = str(row["ProductClass"])

        if product_class not in {"OTHER", "CROSS_SCHEME"}:
            continue

        key = (
            str(row["Primer1"]),
            str(row["Primer2"]),
            product_class,
        )

        if key in seen_unrecognized:
            continue

        seen_unrecognized.add(key)
        unrecognized_examples.append(
            f"{product_class}: "
            f"{row['Primer1']} + {row['Primer2']}"
        )

        if len(unrecognized_examples) >= 10:
            break

    if min_target_size is None or max_target_size is None:
        size_filter_text = "none"
    else:
        size_filter_text = f"{min_target_size}-{max_target_size} bp"

    summary_lines = [
        f"MFEprimer target-coverage summary — {scheme}",
        "=" * (38 + len(scheme)),
        "",
        f"MFEprimer report: {report_file}",
        f"Total validation sequences: {total_target_sequences}",
        f"Potential amplicons parsed: {len(all_rows)}",
        f"TARGET products (LEFT+RIGHT or RIGHT+LEFT): "
        f"{class_counts.get('TARGET', 0)}",
        f"SELF_PRIMING products (LEFT+LEFT or RIGHT+RIGHT): "
        f"{class_counts.get('SELF_PRIMING', 0)}",
        f"CROSS_SCHEME products: {class_counts.get('CROSS_SCHEME', 0)}",
        f"OTHER products: {class_counts.get('OTHER', 0)}",
        f"Target amplicon size filter: {size_filter_text}",
        f"Valid target amplicons after filter: {valid_amplicons}",
        f"Unique positive HitIDs: {unique_positive_targets}",
        f"Additional valid amplicons on already-positive HitIDs: "
        f"{duplicate_target_products}",
        (
            "Primer-pair target coverage: "
            f"{unique_positive_targets} / {total_target_sequences} "
            f"= {coverage_pct:.2f} %"
        ),
        "",
        "Valid target definition:",
        "- same VarVAMP scheme;",
        "- exactly one LEFT and one RIGHT primer;",
        "- LEFT+RIGHT and RIGHT+LEFT are both accepted;",
        "- optional target-amplicon size filter passed.",
    ]

    if unrecognized_examples:
        summary_lines.extend(
            [
                "",
                "Examples of products not recognized as TARGET/SELF_PRIMING:",
                *[
                    f"- {example}"
                    for example in unrecognized_examples
                ],
            ]
        )

    summary_file.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    section(f"MFEprimer results — {scheme}")
    print(f"Validation sequences       {total_target_sequences}")
    print(f"Potential products         {len(all_rows)}")
    print(f"Valid LEFT/RIGHT products  {valid_amplicons}")
    print(f"Self-priming products      {class_counts.get('SELF_PRIMING', 0)}")
    print(f"Unique PCR-positive HitIDs {unique_positive_targets}")
    print(f"PCR coverage               {coverage_pct:.2f} %")
    if class_counts.get("CROSS_SCHEME", 0):
        print(f"Cross-scheme products      {class_counts.get('CROSS_SCHEME', 0)}")
    if class_counts.get("OTHER", 0):
        print(f"Unclassified products      {class_counts.get('OTHER', 0)}")
    if unrecognized_examples and VERBOSE:
        print("\nUnrecognized primer-name examples:")
        for example in unrecognized_examples[:5]:
            print(f"- {example}")
    if VERBOSE:
        print(f"Valid amplicons TSV        {valid_tsv}")
        print(f"Valid amplicons FASTA      {valid_fasta}")
        print(f"Positive HitIDs            {hitids_file}")
        print(f"Non-target products        {non_target_tsv}")
        print(f"Coverage summary           {summary_file}")

    return {
        "scheme": scheme,
        "potential_amplicons": len(all_rows),
        "target_products": class_counts.get("TARGET", 0),
        "self_priming_products": class_counts.get("SELF_PRIMING", 0),
        "valid_target_amplicons": valid_amplicons,
        "unique_positive_hitids": unique_positive_targets,
        "coverage_pct": coverage_pct,
        "valid_tsv": valid_tsv,
        "valid_fasta": valid_fasta,
        "positive_hitids": hitids_file,
        "summary": summary_file,
        "valid_rows": valid_rows,
    }



# ---------------------------------------------------------------------------
# BLAST+ analysis of the original VarVAMP probe
# ---------------------------------------------------------------------------

PROBE_IUPAC_TO_BASES = {
    "A": ("A",),
    "C": ("C",),
    "G": ("G",),
    "T": ("T",),
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T"),
}


def expand_iupac_probe(
    sequence: str,
    max_variants: int = 256,
) -> list[str]:
    """
    Expand an already-degenerate VarVAMP probe into concrete A/C/G/T queries.

    This is NOT probe optimization. It only allows BLAST to evaluate an
    existing IUPAC probe because short-query BLAST works best with concrete
    nucleotide queries.
    """
    sequence = sequence.upper()

    invalid = sorted(set(sequence) - set(PROBE_IUPAC_TO_BASES))
    if invalid:
        raise ValueError(
            "Unsupported probe character(s): " + ", ".join(invalid)
        )

    degeneracy = 1
    for symbol in sequence:
        degeneracy *= len(PROBE_IUPAC_TO_BASES[symbol])

    if degeneracy > max_variants:
        raise RuntimeError(
            f"Probe degeneracy is {degeneracy}, which would require more "
            f"than {max_variants} concrete BLAST queries. "
            "This BLAST-only version stops instead of silently truncating "
            "the probe variants."
        )

    variants = [""]

    for symbol in sequence:
        variants = [
            prefix + base
            for prefix in variants
            for base in PROBE_IUPAC_TO_BASES[symbol]
        ]

    return variants


def write_probe_blast_queries(
    scheme: str,
    probe: str,
    output_file: Path,
) -> list[str]:
    """Write one or more concrete probe variants for blastn-short."""
    variants = expand_iupac_probe(probe)

    with output_file.open("w", encoding="utf-8") as handle:
        for idx, variant in enumerate(variants, start=1):
            handle.write(
                f">{scheme}_PROBE_variant_{idx}\n"
                f"{variant}\n"
            )

    return variants


def parse_amplicon_subject_id(subject_id: str) -> int | None:
    """
    Extract AmpID from headers written by amplicons_valid.fasta.

    Example:
        Amp_2|HitID=OR613708.1|Scheme=varVAMP_0|...
        -> 2
    """
    match = re.match(r"^Amp_(\d+)(?:\||$)", subject_id)

    if not match:
        return None

    return int(match.group(1))


def run_probe_blast(
    scheme: str,
    probe: str,
    valid_amplicons_fasta: Path,
    valid_rows: list[dict[str, object]],
    scheme_dir: Path,
) -> tuple[Path, Path, int]:
    """
    Run blastn-short directly against amplicons_valid.fasta using -subject.

    No second BLAST database is required. The raw tabular BLAST output is
    preserved for later manual/scientific analysis.
    """
    query_file = scheme_dir / "probe_blast_queries.fasta"
    variants = write_probe_blast_queries(
        scheme,
        probe,
        query_file,
    )

    raw_file = scheme_dir / "probe_blast_raw.tsv"

    # Keep enough target sequences even when a future dataset contains
    # substantially more than BLAST's default number of reported targets.
    max_targets = max(1000, len(valid_rows) * 2)

    outfmt = (
        "6 "
        "qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qseq sseq"
    )

    command = [
        "blastn",
        "-task",
        "blastn-short",
        "-query",
        str(query_file),
        "-subject",
        str(valid_amplicons_fasta),
        "-strand",
        "both",
        "-dust",
        "no",
        "-evalue",
        "1000",
        "-max_target_seqs",
        str(max_targets),
        "-max_hsps",
        "20",
        "-outfmt",
        outfmt,
        "-out",
        str(raw_file),
    ]

    ambiguous_positions = [
        f"{idx}:{base}"
        for idx, base in enumerate(probe.upper(), start=1)
        if base not in {"A", "C", "G", "T"}
    ]

    section(f"Probe BLAST — {scheme}")
    print(f"Original probe    : {probe}")
    print(f"Length            : {len(probe)} nt")
    print(
        "IUPAC ambiguity   : "
        + (", ".join(ambiguous_positions) if ambiguous_positions else "none")
    )
    print(f"Concrete variants : {len(variants)}")
    print(f"Target amplicons  : {len(valid_rows)}")
    if VERBOSE:
        print("BLAST task        : blastn-short")
        print("Strands           : both")
        print("DUST              : disabled")
        print("E-value           : 1000")

    run_command(command)

    return raw_file, query_file, len(variants)


def read_probe_blast_raw(
    raw_file: Path,
) -> list[dict[str, object]]:
    """Parse the custom BLAST outfmt 6 produced by run_probe_blast()."""
    columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qseq",
        "sseq",
    ]

    hits: list[dict[str, object]] = []

    if not raw_file.is_file():
        raise FileNotFoundError(
            f"BLAST output file not found: {raw_file}"
        )

    with raw_file.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) != len(columns):
                raise RuntimeError(
                    f"Unexpected BLAST field count on line {line_number}: "
                    f"expected {len(columns)}, found {len(fields)}."
                )

            row = dict(zip(columns, fields))

            amp_id = parse_amplicon_subject_id(
                str(row["sseqid"])
            )

            row.update(
                {
                    "AmpID": amp_id,
                    "pident": float(row["pident"]),
                    "length": int(row["length"]),
                    "mismatch": int(row["mismatch"]),
                    "gapopen": int(row["gapopen"]),
                    "qstart": int(row["qstart"]),
                    "qend": int(row["qend"]),
                    "sstart": int(row["sstart"]),
                    "send": int(row["send"]),
                    "evalue": float(row["evalue"]),
                    "bitscore": float(row["bitscore"]),
                }
            )

            hits.append(row)

    return hits


def classify_full_length_probe_hit(
    hit: dict[str, object],
    probe_length: int,
) -> tuple[bool, str]:
    """
    Define the transparent BLAST classes used by this workflow.

    Full-length ungapped means:
    - alignment starts at query base 1;
    - alignment ends at the complete probe length;
    - aligned length equals probe length;
    - zero gap openings.

    No biological PASS/FAIL label is imposed.
    """
    full_length = (
        int(hit["qstart"]) == 1
        and int(hit["qend"]) == probe_length
        and int(hit["length"]) == probe_length
        and int(hit["gapopen"]) == 0
    )

    if not full_length:
        return False, "PARTIAL_OR_GAPPED"

    mismatches = int(hit["mismatch"])

    if mismatches == 0:
        return True, "0_MISMATCH"

    if mismatches == 1:
        return True, "1_MISMATCH"

    if mismatches == 2:
        return True, "2_MISMATCHES"

    return True, "GE3_MISMATCHES"



def build_blast_match_line(
    query_aligned: str,
    subject_aligned: str,
) -> str:
    """
    Create a BLAST-like visual match line.

    |  exact nucleotide identity
    .  aligned non-identical bases
       gap
    """
    symbols: list[str] = []

    for query_base, subject_base in zip(
        query_aligned.upper(),
        subject_aligned.upper(),
    ):
        if query_base == "-" or subject_base == "-":
            symbols.append(" ")
        elif query_base == subject_base:
            symbols.append("|")
        else:
            symbols.append(".")

    return "".join(symbols)


def extract_full_length_mismatches(
    row: dict[str, object],
) -> list[dict[str, object]]:
    """
    Extract mismatch positions from a full-length ungapped BLAST alignment.

    ProbePosition is 1-based relative to the concrete BLAST query variant.
    SubjectCoordinate respects the BLAST subject orientation.
    """
    if str(row["FullLengthUngapped"]) != "YES":
        return []

    query = str(row["QueryAligned"]).upper()
    subject = str(row["SubjectAligned"]).upper()

    if len(query) != len(subject):
        return []

    subject_start = int(row["SubjectStart"])
    subject_end = int(row["SubjectEnd"])
    subject_step = 1 if subject_end >= subject_start else -1
    subject_coord = subject_start

    mismatch_rows: list[dict[str, object]] = []

    for probe_position, (query_base, subject_base) in enumerate(
        zip(query, subject),
        start=1,
    ):
        if query_base == "-" or subject_base == "-":
            if subject_base != "-":
                subject_coord += subject_step
            continue

        current_subject_coord = subject_coord

        if query_base != subject_base:
            mismatch_rows.append(
                {
                    "AmpID": int(row["AmpID"]),
                    "HitID": str(row["HitID"]),
                    "BlastCategory": str(row["BlastCategory"]),
                    "QueryVariant": str(row["QueryVariant"]),
                    "ProbePosition": probe_position,
                    "ProbeBase": query_base,
                    "TargetBase": subject_base,
                    "SubjectCoordinate": current_subject_coord,
                    "SubjectStrand": (
                        "+"
                        if subject_step == 1
                        else "-"
                    ),
                }
            )

        subject_coord += subject_step

    return mismatch_rows


def format_probe_blast_alignment_block(
    row: dict[str, object],
) -> str:
    """Create one human-readable BLAST alignment block."""
    category = str(row["BlastCategory"])
    hit_id = str(row["HitID"])
    amp_id = int(row["AmpID"])
    amp_size = int(row["AmpliconSize"])

    lines = [
        "=" * 78,
        f"HitID: {hit_id} | AmpID: {amp_id} | Amplicon size: {amp_size} bp",
        f"Category: {category}",
    ]

    if category == "NO_HIT":
        lines.extend(
            [
                "No BLAST HSP was reported for this amplicon.",
                "",
            ]
        )
        return "\n".join(lines)

    qseq = str(row["QueryAligned"])
    sseq = str(row["SubjectAligned"])
    match_line = build_blast_match_line(qseq, sseq)

    sstart = int(row["SubjectStart"])
    send = int(row["SubjectEnd"])
    strand = "+" if send >= sstart else "-"

    lines.extend(
        [
            (
                f"Identity: {row['PercentIdentity']}% | "
                f"Alignment: {row['AlignmentLength']} nt | "
                f"Mismatches: {row['Mismatches']} | "
                f"Gap openings: {row['GapOpenings']}"
            ),
            (
                f"Query variant: {row['QueryVariant']} | "
                f"Query coordinates: {row['QueryStart']}..{row['QueryEnd']}"
            ),
            (
                f"Subject coordinates: {sstart}..{send} | "
                f"Subject strand: {strand}"
            ),
            f"E-value: {row['Evalue']} | Bit score: {row['BitScore']}",
            "",
            f"Probe  {str(row['QueryStart']).rjust(4)}  {qseq}  {row['QueryEnd']}",
            f"             {match_line}",
            f"Target {str(sstart).rjust(4)}  {sseq}  {send}",
        ]
    )

    mismatch_rows = extract_full_length_mismatches(row)

    if str(row["FullLengthUngapped"]) == "YES":
        if mismatch_rows:
            mismatch_text = ", ".join(
                (
                    f"pos{item['ProbePosition']} "
                    f"{item['ProbeBase']}>{item['TargetBase']} "
                    f"(target:{item['SubjectCoordinate']})"
                )
                for item in mismatch_rows
            )
            lines.extend(
                [
                    "",
                    f"Mismatch positions: {mismatch_text}",
                ]
            )
        else:
            lines.extend(
                [
                    "",
                    "Mismatch positions: none",
                ]
            )
    else:
        lines.extend(
            [
                "",
                (
                    "Note: this is not a full-length ungapped alignment; "
                    "mismatch-position interpretation is not used."
                ),
            ]
        )

    lines.append("")
    return "\n".join(lines)


def write_probe_blast_alignment_reports(
    amplicon_summary: list[dict[str, object]],
    hitid_summary: list[dict[str, object]],
    scheme_dir: Path,
) -> tuple[Path, Path, Path]:
    """
    Write human-readable alignments and a machine-readable mismatch table.
    """
    category_order = {
        "0_MISMATCH": 0,
        "1_MISMATCH": 1,
        "2_MISMATCHES": 2,
        "GE3_MISMATCHES": 3,
        "PARTIAL_ONLY": 4,
        "NO_HIT": 5,
    }

    def sort_rows(
        rows: list[dict[str, object]],
    ) -> list[dict[str, object]]:
        return sorted(
            rows,
            key=lambda row: (
                category_order.get(
                    str(row["BlastCategory"]),
                    99,
                ),
                str(row["HitID"]),
                int(row["AmpID"]),
            ),
        )

    amplicon_report = (
        scheme_dir / "probe_blast_alignments_by_amplicon.txt"
    )
    hitid_report = (
        scheme_dir / "probe_blast_alignments_by_hitid.txt"
    )
    mismatch_file = (
        scheme_dir / "probe_blast_mismatch_positions.tsv"
    )

    with amplicon_report.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "BLAST probe alignments — one block per valid MFEprimer amplicon\n"
        )
        handle.write(
            "Symbols: | exact identity, . mismatch, blank = gap\n\n"
        )

        for row in sort_rows(amplicon_summary):
            handle.write(
                format_probe_blast_alignment_block(row)
            )

    with hitid_report.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "BLAST probe alignments — best amplicon per unique HitID\n"
        )
        handle.write(
            "Symbols: | exact identity, . mismatch, blank = gap\n\n"
        )

        for row in sort_rows(hitid_summary):
            handle.write(
                format_probe_blast_alignment_block(row)
            )

    mismatch_rows: list[dict[str, object]] = []

    for row in hitid_summary:
        mismatch_rows.extend(
            extract_full_length_mismatches(row)
        )

    mismatch_columns = [
        "AmpID",
        "HitID",
        "BlastCategory",
        "QueryVariant",
        "ProbePosition",
        "ProbeBase",
        "TargetBase",
        "SubjectCoordinate",
        "SubjectStrand",
    ]

    with mismatch_file.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=mismatch_columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(mismatch_rows)

    return (
        amplicon_report,
        hitid_report,
        mismatch_file,
    )




def summarize_probe_variants(
    query_file: Path,
    raw_hits: list[dict[str, object]],
    valid_rows: list[dict[str, object]],
    probe_length: int,
    scheme_dir: Path,
    category_rank: dict[str, int],
) -> dict[str, object]:
    """
    Analyse every concrete IUPAC-expanded probe variant independently.

    Statistics are calculated at unique HitID level:
      - each concrete probe variant is evaluated against every PCR-positive HitID;
      - if one HitID has multiple valid MFEprimer amplicons, the best result for
        that concrete variant is retained once;
      - variants are NOT merged in this table.

    A separate best-of-all-variants summary is still produced elsewhere and is
    used for the final qPCR coverage.
    """
    query_records = read_fasta(query_file)

    variants = [
        {
            "query_id": header,
            "sequence": sequence,
        }
        for header, sequence in query_records
    ]

    valid_by_ampid = {
        int(row["AmpID"]): row
        for row in valid_rows
    }

    # Stable first-occurrence order of unique biological PCR-positive targets.
    pcr_positive_hitids = list(
        dict.fromkeys(
            str(row["HitID"])
            for row in valid_rows
        )
    )

    # Hits have already been classified as full-length / partial by the caller.
    hits_by_variant_amp: dict[
        tuple[str, int],
        list[dict[str, object]],
    ] = {}

    for hit in raw_hits:
        amp_id = hit.get("AmpID")
        if amp_id is None:
            continue

        amp_id = int(amp_id)
        if amp_id not in valid_by_ampid:
            continue

        query_id = str(hit["qseqid"])
        hits_by_variant_amp.setdefault(
            (query_id, amp_id),
            [],
        ).append(hit)

    def choose_best_amp_result(
        query_id: str,
        amp_id: int,
    ) -> dict[str, object]:
        amp_hits = hits_by_variant_amp.get(
            (query_id, amp_id),
            [],
        )

        full_hits = [
            hit
            for hit in amp_hits
            if hit.get("FullLengthUngapped") == "YES"
        ]

        if full_hits:
            full_hits.sort(
                key=lambda hit: (
                    int(hit["mismatch"]),
                    -float(hit["pident"]),
                    -float(hit["bitscore"]),
                    float(hit["evalue"]),
                )
            )
            best = full_hits[0]
            return {
                "category": str(best["BlastCategory"]),
                "mismatches": int(best["mismatch"]),
                "alignment_length": int(best["length"]),
                "bitscore": float(best["bitscore"]),
                "amp_id": amp_id,
            }

        if amp_hits:
            amp_hits.sort(
                key=lambda hit: (
                    -int(hit["length"]),
                    int(hit["mismatch"]),
                    -float(hit["bitscore"]),
                    float(hit["evalue"]),
                )
            )
            best = amp_hits[0]
            return {
                "category": "PARTIAL_ONLY",
                "mismatches": int(best["mismatch"]),
                "alignment_length": int(best["length"]),
                "bitscore": float(best["bitscore"]),
                "amp_id": amp_id,
            }

        return {
            "category": "NO_HIT",
            "mismatches": 999999,
            "alignment_length": 0,
            "bitscore": 0.0,
            "amp_id": amp_id,
        }

    def result_rank(result: dict[str, object]) -> tuple:
        return (
            category_rank[str(result["category"])],
            int(result["mismatches"]),
            -int(result["alignment_length"]),
            -float(result["bitscore"]),
            int(result["amp_id"]),
        )

    # Map each HitID to all its valid amplicons.
    ampids_by_hitid: dict[str, list[int]] = {}
    for row in valid_rows:
        ampids_by_hitid.setdefault(
            str(row["HitID"]),
            [],
        ).append(int(row["AmpID"]))

    matrix: dict[str, dict[str, str]] = {
        hit_id: {}
        for hit_id in pcr_positive_hitids
    }

    stats_rows: list[dict[str, object]] = []

    for variant in variants:
        query_id = str(variant["query_id"])
        sequence = str(variant["sequence"])

        counts = {
            category: 0
            for category in category_rank
        }

        for hit_id in pcr_positive_hitids:
            amp_results = [
                choose_best_amp_result(
                    query_id,
                    amp_id,
                )
                for amp_id in ampids_by_hitid[hit_id]
            ]

            best = min(
                amp_results,
                key=result_rank,
            )
            category = str(best["category"])

            counts[category] += 1
            matrix[hit_id][query_id] = category

        full_length = (
            counts["0_MISMATCH"]
            + counts["1_MISMATCH"]
            + counts["2_MISMATCHES"]
            + counts["GE3_MISMATCHES"]
        )

        denominator = len(pcr_positive_hitids)

        stats_rows.append(
            {
                "QueryVariant": query_id,
                "ConcreteSequence": sequence,
                "PCRPositiveHitIDs": denominator,
                "Exact0Mismatch": counts["0_MISMATCH"],
                "OneMismatch": counts["1_MISMATCH"],
                "TwoMismatches": counts["2_MISMATCHES"],
                "GE3Mismatches": counts["GE3_MISMATCHES"],
                "PartialOnly": counts["PARTIAL_ONLY"],
                "NoHit": counts["NO_HIT"],
                "FullLengthSite": full_length,
                "ExactPercentPCRPositive": (
                    100.0 * counts["0_MISMATCH"] / denominator
                    if denominator
                    else 0.0
                ),
                "LE1MismatchPercentPCRPositive": (
                    100.0
                    * (
                        counts["0_MISMATCH"]
                        + counts["1_MISMATCH"]
                    )
                    / denominator
                    if denominator
                    else 0.0
                ),
                "LE2MismatchPercentPCRPositive": (
                    100.0
                    * (
                        counts["0_MISMATCH"]
                        + counts["1_MISMATCH"]
                        + counts["2_MISMATCHES"]
                    )
                    / denominator
                    if denominator
                    else 0.0
                ),
                "FullLengthPercentPCRPositive": (
                    100.0 * full_length / denominator
                    if denominator
                    else 0.0
                ),
            }
        )

    stats_file = scheme_dir / "probe_variant_statistics.tsv"
    matrix_file = scheme_dir / "probe_variant_hitid_matrix.tsv"

    stats_columns = [
        "QueryVariant",
        "ConcreteSequence",
        "PCRPositiveHitIDs",
        "Exact0Mismatch",
        "OneMismatch",
        "TwoMismatches",
        "GE3Mismatches",
        "PartialOnly",
        "NoHit",
        "FullLengthSite",
        "ExactPercentPCRPositive",
        "LE1MismatchPercentPCRPositive",
        "LE2MismatchPercentPCRPositive",
        "FullLengthPercentPCRPositive",
    ]

    with stats_file.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=stats_columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(stats_rows)

    variant_ids = [
        str(variant["query_id"])
        for variant in variants
    ]

    with matrix_file.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )
        writer.writerow(
            ["HitID", *variant_ids, "BestAcrossAllVariants"]
        )

        for hit_id in pcr_positive_hitids:
            categories = [
                matrix[hit_id][variant_id]
                for variant_id in variant_ids
            ]

            best_category = min(
                categories,
                key=lambda category: category_rank[category],
            )

            writer.writerow(
                [
                    hit_id,
                    *categories,
                    best_category,
                ]
            )

    # Display this table when expansion generated more than one concrete probe.
    if len(stats_rows) > 1:
        print("\nConcrete probe variant statistics — unique PCR-positive HitIDs")
        print(
            f"{'Variant':<25} {'Sequence':<{probe_length + 2}} "
            f"{'Exact':>6} {'1MM':>6} {'2MM':>6} "
            f"{'>=3':>6} {'Partial':>8} {'No hit':>7} {'Full':>6}"
        )
        print(
            "-" * (
                25 + probe_length + 2
                + 6 + 6 + 6 + 6 + 8 + 7 + 6 + 8
            )
        )

        for row in stats_rows:
            print(
                f"{str(row['QueryVariant']):<25} "
                f"{str(row['ConcreteSequence']):<{probe_length + 2}} "
                f"{int(row['Exact0Mismatch']):6d} "
                f"{int(row['OneMismatch']):6d} "
                f"{int(row['TwoMismatches']):6d} "
                f"{int(row['GE3Mismatches']):6d} "
                f"{int(row['PartialOnly']):8d} "
                f"{int(row['NoHit']):7d} "
                f"{int(row['FullLengthSite']):6d}"
            )

        print(
            "\nNote: each row is one concrete version of the original "
            "ambiguous probe; targets are unique HitIDs."
        )

    return {
        "stats_rows": stats_rows,
        "stats_file": stats_file,
        "matrix_file": matrix_file,
    }



def summarize_probe_blast(
    scheme: str,
    probe: str,
    raw_file: Path,
    query_file: Path,
    concrete_query_count: int,
    valid_rows: list[dict[str, object]],
    scheme_dir: Path,
    total_target_sequences: int,
) -> dict[str, object]:
    """
    Produce BLAST-only summaries.

    There is intentionally:
    - no IUPAC optimization;
    - no ambiguity proposal;
    - no PASS/FAIL probe decision.
    """
    probe_length = len(probe)
    raw_hits = read_probe_blast_raw(raw_file)

    valid_by_ampid = {
        int(row["AmpID"]): row
        for row in valid_rows
    }

    hits_by_ampid: dict[int, list[dict[str, object]]] = {}

    for hit in raw_hits:
        amp_id = hit["AmpID"]

        if amp_id is None:
            continue

        amp_id = int(amp_id)

        if amp_id not in valid_by_ampid:
            continue

        full_length, category = classify_full_length_probe_hit(
            hit,
            probe_length,
        )

        hit["FullLengthUngapped"] = (
            "YES" if full_length else "NO"
        )
        hit["BlastCategory"] = category

        hits_by_ampid.setdefault(amp_id, []).append(hit)

    category_rank = {
        "0_MISMATCH": 0,
        "1_MISMATCH": 1,
        "2_MISMATCHES": 2,
        "GE3_MISMATCHES": 3,
        "PARTIAL_ONLY": 4,
        "NO_HIT": 5,
    }

    variant_analysis = summarize_probe_variants(
        query_file,
        raw_hits,
        valid_rows,
        probe_length,
        scheme_dir,
        category_rank,
    )

    amplicon_summary: list[dict[str, object]] = []

    for amp_id in sorted(valid_by_ampid):
        source = valid_by_ampid[amp_id]
        amp_hits = hits_by_ampid.get(amp_id, [])

        full_hits = [
            hit
            for hit in amp_hits
            if hit.get("FullLengthUngapped") == "YES"
        ]

        if full_hits:
            full_hits.sort(
                key=lambda hit: (
                    int(hit["mismatch"]),
                    -float(hit["pident"]),
                    -float(hit["bitscore"]),
                    float(hit["evalue"]),
                )
            )
            best = full_hits[0]
            category = str(best["BlastCategory"])
        elif amp_hits:
            amp_hits.sort(
                key=lambda hit: (
                    -int(hit["length"]),
                    int(hit["mismatch"]),
                    -float(hit["bitscore"]),
                    float(hit["evalue"]),
                )
            )
            best = amp_hits[0]
            category = "PARTIAL_ONLY"
        else:
            best = None
            category = "NO_HIT"

        if best is None:
            row = {
                "AmpID": amp_id,
                "HitID": str(source["HitID"]),
                "AmpliconSize": int(source["Size"]),
                "BlastCategory": category,
                "FullLengthUngapped": "NO",
                "QueryVariant": "",
                "PercentIdentity": "",
                "AlignmentLength": "",
                "Mismatches": "",
                "GapOpenings": "",
                "QueryStart": "",
                "QueryEnd": "",
                "SubjectStart": "",
                "SubjectEnd": "",
                "Evalue": "",
                "BitScore": "",
                "QueryAligned": "",
                "SubjectAligned": "",
                "RawHSPCount": 0,
            }
        else:
            row = {
                "AmpID": amp_id,
                "HitID": str(source["HitID"]),
                "AmpliconSize": int(source["Size"]),
                "BlastCategory": category,
                "FullLengthUngapped": str(
                    best.get("FullLengthUngapped", "NO")
                ),
                "QueryVariant": str(best["qseqid"]),
                "PercentIdentity": float(best["pident"]),
                "AlignmentLength": int(best["length"]),
                "Mismatches": int(best["mismatch"]),
                "GapOpenings": int(best["gapopen"]),
                "QueryStart": int(best["qstart"]),
                "QueryEnd": int(best["qend"]),
                "SubjectStart": int(best["sstart"]),
                "SubjectEnd": int(best["send"]),
                "Evalue": float(best["evalue"]),
                "BitScore": float(best["bitscore"]),
                "QueryAligned": str(best["qseq"]),
                "SubjectAligned": str(best["sseq"]),
                "RawHSPCount": len(amp_hits),
            }

        amplicon_summary.append(row)

    # One best amplicon result per unique target HitID.
    hitid_best: dict[str, dict[str, object]] = {}

    def result_rank(row: dict[str, object]) -> tuple:
        category = str(row["BlastCategory"])

        mismatch_value = (
            int(row["Mismatches"])
            if row["Mismatches"] != ""
            else 999999
        )

        alignment_length = (
            int(row["AlignmentLength"])
            if row["AlignmentLength"] != ""
            else 0
        )

        bitscore = (
            float(row["BitScore"])
            if row["BitScore"] != ""
            else 0.0
        )

        return (
            category_rank[category],
            mismatch_value,
            -alignment_length,
            -bitscore,
            int(row["AmpID"]),
        )

    for row in amplicon_summary:
        hit_id = str(row["HitID"])

        if (
            hit_id not in hitid_best
            or result_rank(row) < result_rank(hitid_best[hit_id])
        ):
            hitid_best[hit_id] = row

    hitid_summary = list(hitid_best.values())
    hitid_summary.sort(
        key=lambda row: str(row["HitID"])
    )

    (
        alignment_amplicon_file,
        alignment_hitid_file,
        mismatch_positions_file,
    ) = write_probe_blast_alignment_reports(
        amplicon_summary,
        hitid_summary,
        scheme_dir,
    )

    counts = {
        category: 0
        for category in category_rank
    }

    for row in hitid_summary:
        counts[str(row["BlastCategory"])] += 1

    amplicon_counts = {
        category: 0
        for category in category_rank
    }

    for row in amplicon_summary:
        amplicon_counts[str(row["BlastCategory"])] += 1

    full_length_amplicons = (
        amplicon_counts["0_MISMATCH"]
        + amplicon_counts["1_MISMATCH"]
        + amplicon_counts["2_MISMATCHES"]
        + amplicon_counts["GE3_MISMATCHES"]
    )

    exact_probe_amplicons = amplicon_counts["0_MISMATCH"]

    full_length_hitids = (
        counts["0_MISMATCH"]
        + counts["1_MISMATCH"]
        + counts["2_MISMATCHES"]
        + counts["GE3_MISMATCHES"]
    )

    exact_probe_hitids = counts["0_MISMATCH"]
    pcr_positive_hitids = len(hitid_summary)
    total_valid_amplicons = len(amplicon_summary)

    # -------------------- write processed outputs --------------------
    columns = [
        "AmpID",
        "HitID",
        "AmpliconSize",
        "BlastCategory",
        "FullLengthUngapped",
        "QueryVariant",
        "PercentIdentity",
        "AlignmentLength",
        "Mismatches",
        "GapOpenings",
        "QueryStart",
        "QueryEnd",
        "SubjectStart",
        "SubjectEnd",
        "Evalue",
        "BitScore",
        "QueryAligned",
        "SubjectAligned",
        "RawHSPCount",
    ]

    amplicon_file = scheme_dir / "probe_blast_amplicon_summary.tsv"
    hitid_file = scheme_dir / "probe_blast_hitid_summary.tsv"
    summary_file = scheme_dir / "probe_blast_summary.txt"

    def write_table(
        path: Path,
        rows: list[dict[str, object]],
    ) -> None:
        with path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=columns,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)

    write_table(amplicon_file, amplicon_summary)
    write_table(hitid_file, hitid_summary)

    summary_lines = [
        f"BLAST probe summary — {scheme}",
        "=" * (22 + len(scheme)),
        "",
        f"Original VarVAMP probe: {probe}",
        f"Probe length: {probe_length} nt",
        f"Concrete BLAST query variants: {concrete_query_count}",
        f"Valid MFEprimer target amplicons: {len(valid_rows)}",
        f"Unique PCR-positive HitIDs: {pcr_positive_hitids}",
        f"Raw BLAST HSPs: {len(raw_hits)}",
        "",
        "AMPlicon-level probe statistics:",
        f"Valid target amplicons analysed: {total_valid_amplicons}",
        (
            "Exact probe sequence (0 mismatch): "
            f"{exact_probe_amplicons} / {total_valid_amplicons} "
            f"= {(100.0 * exact_probe_amplicons / total_valid_amplicons if total_valid_amplicons else 0.0):.2f} %"
        ),
        (
            "Full-length probe alignments (any mismatch count): "
            f"{full_length_amplicons} / {total_valid_amplicons} "
            f"= {(100.0 * full_length_amplicons / total_valid_amplicons if total_valid_amplicons else 0.0):.2f} %"
        ),
        f"0 mismatch: {amplicon_counts['0_MISMATCH']}",
        f"1 mismatch: {amplicon_counts['1_MISMATCH']}",
        f"2 mismatches: {amplicon_counts['2_MISMATCHES']}",
        f">=3 mismatches: {amplicon_counts['GE3_MISMATCHES']}",
        f"Partial/gapped hit only: {amplicon_counts['PARTIAL_ONLY']}",
        f"No BLAST hit: {amplicon_counts['NO_HIT']}",
        "",
        "HITID-level probe statistics (unique PCR-positive targets):",
        f"Unique PCR-positive HitIDs analysed: {pcr_positive_hitids}",
        (
            "Exact probe sequence (0 mismatch): "
            f"{exact_probe_hitids} / {pcr_positive_hitids} "
            f"= {(100.0 * exact_probe_hitids / pcr_positive_hitids if pcr_positive_hitids else 0.0):.2f} %"
        ),
        (
            "Full-length probe alignments (any mismatch count): "
            f"{full_length_hitids} / {pcr_positive_hitids} "
            f"= {(100.0 * full_length_hitids / pcr_positive_hitids if pcr_positive_hitids else 0.0):.2f} %"
        ),
        f"0 mismatch: {counts['0_MISMATCH']}",
        f"1 mismatch: {counts['1_MISMATCH']}",
        f"2 mismatches: {counts['2_MISMATCHES']}",
        f">=3 mismatches: {counts['GE3_MISMATCHES']}",
        f"Partial/gapped hit only: {counts['PARTIAL_ONLY']}",
        f"No BLAST hit: {counts['NO_HIT']}",
        "",
        (
            "Full-length ungapped BLAST matches among PCR-positive HitIDs: "
            f"{full_length_hitids} / {pcr_positive_hitids} "
            f"= "
            f"{(100.0 * full_length_hitids / pcr_positive_hitids if pcr_positive_hitids else 0.0):.2f} %"
        ),
        (
            "Full-length ungapped probe matches relative to the complete "
            f"validation database: {full_length_hitids} / "
            f"{total_target_sequences} = "
            f"{(100.0 * full_length_hitids / total_target_sequences if total_target_sequences else 0.0):.2f} %"
        ),
        "",
        "Full-length ungapped definition:",
        f"- query starts at 1 and ends at {probe_length};",
        f"- alignment length = {probe_length} nt;",
        "- zero gap openings.",
        "",
        "This version reports BLAST results only.",
        "No probe optimization or IUPAC substitution is proposed.",
        f"Raw BLAST output: {raw_file}",
        f"BLAST query FASTA: {query_file}",
        f"Visual alignments by amplicon: {alignment_amplicon_file}",
        f"Visual alignments by HitID: {alignment_hitid_file}",
        f"Mismatch positions TSV: {mismatch_positions_file}",
        f"Per-variant statistics TSV: {variant_analysis['stats_file']}",
        f"HitID x probe-variant matrix TSV: {variant_analysis['matrix_file']}",
    ]

    summary_file.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    section(f"Probe results — {scheme}")

    def _pct(value: int, total: int) -> float:
        return 100.0 * value / total if total else 0.0

    print("Biological coverage — unique PCR-positive targets")
    print(f"  PCR-positive HitIDs   {pcr_positive_hitids}")
    print(
        f"  Exact probe           {counts['0_MISMATCH']:>4}  "
        f"{_pct(counts['0_MISMATCH'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  1 mismatch            {counts['1_MISMATCH']:>4}  "
        f"{_pct(counts['1_MISMATCH'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  2 mismatches          {counts['2_MISMATCHES']:>4}  "
        f"{_pct(counts['2_MISMATCHES'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  >=3 mismatches        {counts['GE3_MISMATCHES']:>4}  "
        f"{_pct(counts['GE3_MISMATCHES'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  Partial/gapped        {counts['PARTIAL_ONLY']:>4}  "
        f"{_pct(counts['PARTIAL_ONLY'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  No hit                {counts['NO_HIT']:>4}  "
        f"{_pct(counts['NO_HIT'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  Full-length site      {full_length_hitids:>4}  "
        f"{_pct(full_length_hitids, pcr_positive_hitids):6.2f} %"
    )

    print("\nTechnical amplicon statistics")
    print(f"  Valid amplicons       {total_valid_amplicons}")
    print(f"  Exact probe           {amplicon_counts['0_MISMATCH']}")
    print(f"  1 mismatch            {amplicon_counts['1_MISMATCH']}")
    print(f"  2 mismatches          {amplicon_counts['2_MISMATCHES']}")
    print(f"  >=3 mismatches        {amplicon_counts['GE3_MISMATCHES']}")
    print(f"  Partial/gapped        {amplicon_counts['PARTIAL_ONLY']}")
    print(f"  No hit                {amplicon_counts['NO_HIT']}")

    if VERBOSE:
        print(f"\nRaw BLAST               {raw_file}")
        print(f"Amplicon summary         {amplicon_file}")
        print(f"HitID summary            {hitid_file}")
        print(f"Alignment by amplicon    {alignment_amplicon_file}")
        print(f"Alignment by HitID       {alignment_hitid_file}")
        print(f"Mismatch positions       {mismatch_positions_file}")
        print(f"Variant statistics       {variant_analysis['stats_file']}")
        print(f"HitID x variant matrix   {variant_analysis['matrix_file']}")
        print(f"BLAST summary            {summary_file}")

    return {
        "raw_hsps": len(raw_hits),
        "amplicon_total": total_valid_amplicons,
        "amplicon_exact_0": exact_probe_amplicons,
        "amplicon_full_length": full_length_amplicons,
        "amplicon_1_mismatch": amplicon_counts["1_MISMATCH"],
        "amplicon_2_mismatch": amplicon_counts["2_MISMATCHES"],
        "amplicon_ge3_mismatch": amplicon_counts["GE3_MISMATCHES"],
        "exact_0": counts["0_MISMATCH"],
        "mismatch_1": counts["1_MISMATCH"],
        "mismatch_2": counts["2_MISMATCHES"],
        "mismatch_ge3": counts["GE3_MISMATCHES"],
        "partial_only": counts["PARTIAL_ONLY"],
        "no_hit": counts["NO_HIT"],
        "full_length_hitids": full_length_hitids,
        "pcr_positive_hitids": pcr_positive_hitids,
        "amplicon_summary": amplicon_file,
        "hitid_summary": hitid_file,
        "blast_summary": summary_file,
        "alignment_amplicon_file": alignment_amplicon_file,
        "alignment_hitid_file": alignment_hitid_file,
        "mismatch_positions_file": mismatch_positions_file,
        "variant_statistics_file": variant_analysis["stats_file"],
        "variant_hitid_matrix_file": variant_analysis["matrix_file"],
    }


def analyse_probe_with_blast(
    scheme: str,
    probe: str,
    valid_amplicons_fasta: Path,
    valid_rows: list[dict[str, object]],
    scheme_dir: Path,
    total_target_sequences: int,
) -> dict[str, object]:
    """Run and summarize BLAST analysis for the original VarVAMP probe."""
    raw_file, query_file, variant_count = run_probe_blast(
        scheme,
        probe,
        valid_amplicons_fasta,
        valid_rows,
        scheme_dir,
    )

    return summarize_probe_blast(
        scheme,
        probe,
        raw_file,
        query_file,
        variant_count,
        valid_rows,
        scheme_dir,
        total_target_sequences,
    )


def run_mfeprimer_pairs(
    schemes: dict[str, dict[str, str]],
    selected_schemes: list[str],
    validation_root: Path,
    database: Path,
    parameters: dict[str, float | int] | None,
    total_target_sequences: int,
    min_target_size: int | None,
    max_target_size: int | None,
    assay_mode: str,
) -> None:
    db_name = re.sub(r"[^A-Za-z0-9._-]+", "_", database.stem)
    coverage_root = validation_root / "mfeprimer_coverage" / db_name
    coverage_root.mkdir(parents=True, exist_ok=True)

    run_manifest: list[list[str]] = []
    coverage_manifest: list[list[object]] = []

    print("\n6. MFEprimer target-coverage testing")
    print("===================================")
    print(f"Target database: {database}")
    print(f"Primer pairs selected: {len(selected_schemes)}")

    for scheme in selected_schemes:
        pair_file = (
            validation_root
            / "primer_pairs"
            / f"{scheme}_primers.fasta"
        )

        scheme_dir = coverage_root / scheme
        scheme_dir.mkdir(parents=True, exist_ok=True)

        # Since MFEprimer >=4.2 uses TSV output by default and -o is mandatory,
        # use a clean prefix. Current versions create <prefix>.spec.tsv.
        output_prefix = scheme_dir / scheme

        command = [
            "mfeprimer",
            "spec",
            "-i",
            str(pair_file),
            "-d",
            str(database),
            "-o",
            str(output_prefix),
        ]

        if parameters is not None:
            command.extend(
                [
                    "-s",
                    str(parameters["min_size"]),
                    "-S",
                    str(parameters["max_size"]),
                    "-t",
                    str(parameters["tm_cutoff"]),
                ]
            )

        print(f"\nTesting {scheme}")
        print("-" * (8 + len(scheme)))

        run_command(command)

        expected_tsv = Path(str(output_prefix) + ".spec.tsv")
        legacy_output = output_prefix

        if legacy_output.is_file():
            # MFEprimer 3.x legacy text report. This is the report format
            # from which this workflow extracts amplicon sequences.
            result_path = legacy_output
        elif expected_tsv.is_file():
            result_path = expected_tsv
        else:
            raise RuntimeError(
                "MFEprimer completed but its output file could not be found "
                f"for scheme {scheme}."
            )

        run_manifest.append(
            [
                scheme,
                str(pair_file),
                str(database),
                str(result_path),
            ]
        )

        if result_path == expected_tsv:
            raise RuntimeError(
                "A .spec.tsv result was produced, but amplicon-sequence "
                "extraction in this version of the workflow is implemented "
                "for the legacy MFEprimer 3.x text report. "
                "Your previous tests produced the supported text format."
            )

        analysis = analyse_and_extract_mfeprimer_amplicons(
            result_path,
            scheme,
            scheme_dir,
            total_target_sequences,
            min_target_size,
            max_target_size,
        )

        if assay_mode != "qpcr":
            blast_analysis = {
                "amplicon_total": 0,
                "amplicon_exact_0": 0,
                "amplicon_full_length": 0,
                "amplicon_1_mismatch": 0,
                "amplicon_2_mismatch": 0,
                "amplicon_ge3_mismatch": 0,
                "exact_0": 0,
                "mismatch_1": 0,
                "mismatch_2": 0,
                "mismatch_ge3": 0,
                "partial_only": 0,
                "no_hit": 0,
                "full_length_hitids": 0,
                "hitid_summary": "",
                "blast_summary": "",
                "alignment_hitid_file": "",
                "mismatch_positions_file": "",
                "variant_statistics_file": "",
                "variant_hitid_matrix_file": "",
            }
        elif analysis["valid_target_amplicons"] == 0:
            print(
                f"\nBLAST probe analysis skipped — {scheme}"
            )
            print(
                "Reason: MFEprimer produced no amplicons classified "
                "as valid TARGET products."
            )
            print(
                "Inspect non_target_products.tsv and coverage_summary.txt "
                "for the raw primer identifiers before interpreting this "
                "as true 0 % biological coverage."
            )

            blast_analysis = {
                "amplicon_total": 0,
                "amplicon_exact_0": 0,
                "amplicon_full_length": 0,
                "amplicon_1_mismatch": 0,
                "amplicon_2_mismatch": 0,
                "amplicon_ge3_mismatch": 0,
                "exact_0": 0,
                "mismatch_1": 0,
                "mismatch_2": 0,
                "mismatch_ge3": 0,
                "partial_only": 0,
                "no_hit": 0,
                "full_length_hitids": 0,
                "hitid_summary": "",
                "blast_summary": "",
                "alignment_hitid_file": "",
                "mismatch_positions_file": "",
                "variant_statistics_file": "",
                "variant_hitid_matrix_file": "",
            }
        else:
            blast_analysis = analyse_probe_with_blast(
                scheme,
                schemes[scheme]["PROBE"],
                analysis["valid_fasta"],
                analysis["valid_rows"],
                scheme_dir,
                total_target_sequences,
            )

        coverage_manifest.append(
            [
                assay_mode,
                scheme,
                analysis["potential_amplicons"],
                analysis["target_products"],
                analysis["self_priming_products"],
                analysis["valid_target_amplicons"],
                analysis["unique_positive_hitids"],
                f"{analysis['coverage_pct']:.4f}",
                blast_analysis["amplicon_total"],
                blast_analysis["amplicon_exact_0"],
                blast_analysis["amplicon_full_length"],
                blast_analysis["amplicon_1_mismatch"],
                blast_analysis["amplicon_2_mismatch"],
                blast_analysis["amplicon_ge3_mismatch"],
                blast_analysis["exact_0"],
                blast_analysis["mismatch_1"],
                blast_analysis["mismatch_2"],
                blast_analysis["mismatch_ge3"],
                blast_analysis["partial_only"],
                blast_analysis["no_hit"],
                blast_analysis["full_length_hitids"],
                str(analysis["valid_fasta"]),
                str(analysis["positive_hitids"]),
                str(analysis["summary"]),
                str(blast_analysis["hitid_summary"]),
                str(blast_analysis["blast_summary"]),
                str(blast_analysis["alignment_hitid_file"]),
                str(blast_analysis["mismatch_positions_file"]),
                (
                    f"{(100.0 * blast_analysis['exact_0'] / total_target_sequences if total_target_sequences else 0.0):.4f}"
                    if assay_mode == "qpcr"
                    else ""
                ),
                (
                    f"{(100.0 * (blast_analysis['exact_0'] + blast_analysis['mismatch_1']) / total_target_sequences if total_target_sequences else 0.0):.4f}"
                    if assay_mode == "qpcr"
                    else ""
                ),
                (
                    f"{(100.0 * (blast_analysis['exact_0'] + blast_analysis['mismatch_1'] + blast_analysis['mismatch_2']) / total_target_sequences if total_target_sequences else 0.0):.4f}"
                    if assay_mode == "qpcr"
                    else ""
                ),
                (
                    f"{(100.0 * blast_analysis['full_length_hitids'] / total_target_sequences if total_target_sequences else 0.0):.4f}"
                    if assay_mode == "qpcr"
                    else ""
                ),
                str(blast_analysis["variant_statistics_file"]),
                str(blast_analysis["variant_hitid_matrix_file"]),
            ]
        )

    section("Final assay comparison")
    if assay_mode == "qpcr":
        print(
            "All percentages below use the COMPLETE validation database "
            f"(N={total_target_sequences}) as denominator."
        )
        print(
            f"{'Scheme':<16} {'PCR':>8} {'PCR+probe':>11} "
            f"{'PCR+probe':>11} {'PCR+probe':>11} {'PCR+full':>10}"
        )
        print(
            f"{'':<16} {'':>8} {'exact':>11} "
            f"{'<=1 MM':>11} {'<=2 MM':>11} {'probe site':>10}"
        )
        print("-" * 73)

        for row in coverage_manifest:
            scheme_name = str(row[1])
            pcr_pct = float(row[7])

            exact = int(row[14])
            mm1 = int(row[15])
            mm2 = int(row[16])
            full_length = int(row[20])

            exact_global = (
                100.0 * exact / total_target_sequences
                if total_target_sequences
                else 0.0
            )
            le1_global = (
                100.0 * (exact + mm1) / total_target_sequences
                if total_target_sequences
                else 0.0
            )
            le2_global = (
                100.0 * (exact + mm1 + mm2) / total_target_sequences
                if total_target_sequences
                else 0.0
            )
            full_global = (
                100.0 * full_length / total_target_sequences
                if total_target_sequences
                else 0.0
            )

            print(
                f"{scheme_name:<16} "
                f"{pcr_pct:7.2f}% "
                f"{exact_global:10.2f}% "
                f"{le1_global:10.2f}% "
                f"{le2_global:10.2f}% "
                f"{full_global:9.2f}%"
            )

        print(
            "\nDefinitions: PCR = LEFT+RIGHT; all 'PCR+probe' columns require "
            "a PCR-positive target plus the indicated probe criterion."
        )
    else:
        print(
            f"{'Scheme':<24} {'PCR coverage':>14} "
            f"{'Unique targets':>16}"
        )
        print("-" * 58)
        for row in coverage_manifest:
            print(
                f"{str(row[1]):<24} {float(row[7]):13.2f}% "
                f"{int(row[6]):16d}"
            )

    with (coverage_root / "mfeprimer_runs.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "scheme",
                "primer_pair_fasta",
                "target_database",
                "mfeprimer_result",
            ]
        )
        writer.writerows(run_manifest)

    coverage_manifest_file = coverage_root / "coverage_results.tsv"

    with coverage_manifest_file.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "mode",
                "scheme",
                "potential_amplicons",
                "target_products",
                "self_priming_products",
                "valid_target_amplicons",
                "unique_positive_hitids",
                "primer_pair_coverage_percent",
                "probe_blast_total_valid_amplicons",
                "probe_blast_exact_0_mismatch_amplicons",
                "probe_blast_full_length_amplicons",
                "probe_blast_1_mismatch_amplicons",
                "probe_blast_2_mismatch_amplicons",
                "probe_blast_ge3_mismatch_amplicons",
                "probe_blast_0_mismatch_hitids",
                "probe_blast_1_mismatch_hitids",
                "probe_blast_2_mismatch_hitids",
                "probe_blast_ge3_mismatch_hitids",
                "probe_blast_partial_only_hitids",
                "probe_blast_no_hit_hitids",
                "probe_blast_full_length_ungapped_hitids",
                "valid_amplicons_fasta",
                "positive_hitids_file",
                "coverage_summary",
                "probe_blast_hitid_summary",
                "probe_blast_summary",
                "probe_blast_visual_alignment_by_hitid",
                "probe_blast_mismatch_positions",
                "qpcr_exact_global_percent",
                "qpcr_le1_mismatch_global_percent",
                "qpcr_le2_mismatch_global_percent",
                "qpcr_full_probe_site_global_percent",
                "probe_variant_statistics",
                "probe_variant_hitid_matrix",
            ]
        )
        writer.writerows(coverage_manifest)

    section("Validation completed")
    print(f"Mode    : {assay_mode.upper()}")
    print(f"Results : {compact_path(coverage_root)}")
    print(f"Summary : {compact_path(coverage_manifest_file)}")
    if VALIDATION_LOG is not None:
        print(f"Log     : {compact_path(VALIDATION_LOG)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    global VERBOSE, VALIDATION_LOG

    cli = parse_cli()
    VERBOSE = bool(cli.verbose)

    try:
        validation_root = (Path.cwd() / "validation_inputs").resolve()
        validation_root.mkdir(parents=True, exist_ok=True)
        VALIDATION_LOG = validation_root / "validation.log"
        VALIDATION_LOG.write_text(
            "VarVAMP in silico validation log\n",
            encoding="utf-8",
        )

        print("╔══════════════════════════════════════════════╗")
        print("║       VarVAMP IN SILICO VALIDATION          ║")
        print("╚══════════════════════════════════════════════╝")

        assay_mode, result_dir = choose_varvamp_result()

        if assay_mode == "qpcr":
            oligos_file = find_qpcr_oligos(result_dir)
            schemes = prepare_varvamp_files(
                oligos_file,
                validation_root,
            )
        elif assay_mode in {"single", "tiled"}:
            schemes = prepare_single_or_tiled_files(
                result_dir,
                validation_root,
                assay_mode,
            )
        else:
            raise RuntimeError(f"Unsupported VarVAMP mode: {assay_mode}")

        progress("Preparing assay files")

        if not ask_yes_no(
            "Continue to target-coverage validation with MFEprimer?"
        ):
            print("Validation cancelled.")
            return 0

        ensure_validation_tools_available(assay_mode)

        original_database = choose_target_database(validation_root)
        stats = fasta_statistics(original_database)

        section("Validation database")
        print(f"FASTA       : {original_database.name}")
        print(f"Sequences   : {stats['sequences']}")
        print(
            f"Length      : {stats['min_length']}-"
            f"{stats['max_length']} nt"
        )
        print(f"Mean length : {stats['mean_length']:.2f} nt")
        print(f"Ambiguous   : {stats['ambiguous']}")

        if not ask_yes_no(
            "Use this FASTA as the target coverage database?"
        ):
            print("Target coverage validation cancelled.")
            return 0

        database = format_target_database_for_mfeprimer(
            original_database,
            validation_root,
        )
        index_target_database(database)

        selected_schemes = choose_schemes(schemes)
        parameters = choose_mfeprimer_parameters()
        min_target_size, max_target_size = choose_coverage_amplicon_filter()
        total_target_sequences = int(stats["sequences"])

        section("Validation plan")
        print(f"Mode              : {assay_mode.upper()}")
        print(f"Database          : {original_database.name}")
        print(f"Sequences         : {total_target_sequences}")
        print("Selected assays   : " + ", ".join(selected_schemes))
        if parameters is None:
            print("MFEprimer filters : defaults")
        else:
            print(
                "MFEprimer filters : "
                f"{parameters['min_size']}-{parameters['max_size']} bp; "
                f"Tm >= {parameters['tm_cutoff']} °C"
            )
        if min_target_size is None or max_target_size is None:
            print("Coverage size     : no additional filter")
        else:
            print(
                "Coverage size     : "
                f"{min_target_size}-{max_target_size} bp"
            )
        if assay_mode == "qpcr":
            print(
                "Probe analysis    : blastn-short on all concrete "
                "IUPAC probe variants"
            )
        else:
            print("Probe analysis    : not applicable")

        if not ask_yes_no("Run validation now?"):
            print("Validation cancelled.")
            return 0

        progress("Starting MFEprimer validation", "RUN")

        run_mfeprimer_pairs(
            schemes,
            selected_schemes,
            validation_root,
            database,
            parameters,
            total_target_sequences,
            min_target_size,
            max_target_size,
            assay_mode,
        )
        return 0

    except (
        FileNotFoundError,
        ValueError,
        RuntimeError,
        KeyboardInterrupt,
    ) as error:
        log_line("ERROR: " + str(error))
        print(f"\nERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())