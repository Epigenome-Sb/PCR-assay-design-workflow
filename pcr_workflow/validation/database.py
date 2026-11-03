# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""Validation FASTA selection, length filtering and MFEprimer database preparation."""

from __future__ import annotations

import re
import shutil
from math import ceil
from pathlib import Path

from .fasta import (
    fasta_statistics,
    read_fasta_with_full_headers,
    write_fasta_with_full_headers,
)
from .questions import ask_choice, ask_float, ask_int, ask_yes_no
from .run_management import (
    compact_path,
    progress,
    run_command,
    run_command_to_file,
)


def choose_sequence_length_filter(
    database: Path,
    database_workdir: Path,
) -> tuple[Path, str, dict[str, float | int]]:
    """
    Decide which sequence-length set enters MFEprimer preparation.

    The "near-complete" option is an operational length criterion relative
    to the longest sequence in the selected database; it does not infer
    biological completeness from annotation.
    """
    before_stats = fasta_statistics(database)
    max_length = int(before_stats["max_length"])

    print("\n3. Sequence-length filtering")
    print("============================")
    print(
        f"Current database: {before_stats['sequences']} sequences; "
        f"length range {before_stats['min_length']}-{before_stats['max_length']} nt."
    )

    mode = ask_choice(
        "How should the validation database be prepared?",
        [
            "Keep all sequences.",
            "Keep complete / nearly complete sequences using a relative-length threshold.",
            "Apply a custom minimum sequence length.",
        ],
    )

    threshold: int | None = None
    strategy: str

    if mode == 1:
        strategy = "keep all sequences"
        retained_database = database
    elif mode == 2:
        print(
            "\nNear-complete filtering is defined here by sequence length "
            "relative to the longest sequence in this database."
        )
        relative_mode = ask_choice(
            "Near-complete length threshold",
            [
                "Keep sequences >=90% of the maximum sequence length.",
                "Keep sequences >=95% of the maximum sequence length.",
                "Enter a custom percentage of the maximum sequence length.",
            ],
        )
        if relative_mode == 1:
            percentage = 90.0
        elif relative_mode == 2:
            percentage = 95.0
        else:
            percentage = ask_float(
                "Minimum relative length (% of maximum)",
                50.0,
                100.0,
            )
        threshold = ceil(max_length * percentage / 100.0)
        strategy = (
            f"near-complete length >= {percentage:g}% of maximum "
            f"({threshold} nt)"
        )
    else:
        threshold = ask_int("Minimum sequence length (nt)", 1, max_length)
        strategy = f"minimum sequence length >= {threshold} nt"

    if threshold is not None:
        source_records = read_fasta_with_full_headers(database)
        retained_records = [
            (header, sequence)
            for header, sequence in source_records
            if len(sequence) >= threshold
        ]
        if not retained_records:
            raise RuntimeError(
                "The selected length filter removed every sequence. "
                "Choose a less restrictive threshold."
            )
        length_filter_dir = database_workdir / "length_filter"
        length_filter_dir.mkdir(parents=True, exist_ok=True)
        safe_stem = re.sub(r"[^A-Za-z0-9._-]+", "_", database.stem)
        retained_database = length_filter_dir / f"{safe_stem}_min{threshold}.fasta"
        write_fasta_with_full_headers(retained_database, retained_records)

    after_stats = fasta_statistics(retained_database)
    before_n = int(before_stats["sequences"])
    after_n = int(after_stats["sequences"])
    removed = before_n - after_n
    retained_pct = 100.0 * after_n / before_n if before_n else 0.0
    removed_pct = 100.0 * removed / before_n if before_n else 0.0

    print("\nSequence filtering results")
    print("--------------------------")
    print(f"Strategy           : {strategy}")
    print(f"Sequences before   : {before_n}")
    print(f"Sequences retained : {after_n} ({retained_pct:.2f} %)")
    print(f"Sequences removed  : {removed} ({removed_pct:.2f} %)")
    print(f"Retained length    : {after_stats['min_length']}-{after_stats['max_length']} nt")
    print(f"Mean length        : {after_stats['mean_length']:.2f} nt")
    print(f"Median length      : {after_stats['median_length']:.2f} nt")
    print(f"Ambiguous bases    : {after_stats['ambiguous']}")
    if threshold is not None:
        print(f"Filtered FASTA     : {compact_path(retained_database)}")
    else:
        print("Filtered FASTA     : not created (all sequences retained)")

    return retained_database.resolve(), strategy, after_stats


def is_inside(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def find_candidate_target_fastas(validation_root: Path | None = None) -> list[Path]:
    """Return user validation FASTA files from data/validation only."""
    root = (Path.cwd() / "data" / "validation").resolve()
    if not root.is_dir():
        return []
    candidates: list[Path] = []
    for pattern in ("*.fasta", "*.fa", "*.fna"):
        candidates.extend(path.resolve() for path in root.rglob(pattern) if path.is_file())
    return sorted(set(candidates))


def choose_target_database(
    requested_database: Path | None = None,
) -> Path:
    """Select the independent validation FASTA, normally from data/validation/."""
    print("\n2. Validation database")
    print("======================")

    if requested_database is not None:
        candidate = requested_database.expanduser().resolve()
        if not candidate.is_file():
            raise FileNotFoundError(f"Validation FASTA not found: {candidate}")
        return candidate

    candidates = find_candidate_target_fastas()
    if candidates:
        for idx, path in enumerate(candidates, start=1):
            print(f"{idx}. {compact_path(path)}")
        print(f"{len(candidates) + 1}. Enter another FASTA path")

        while True:
            raw = input("Select the validation database by number: ").strip()
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
        raw = input("Name or path of the validation FASTA: ").strip().strip("'\"")
        if not raw:
            print("A FASTA path is required.")
            continue
        candidate = Path(raw).expanduser()
        paths = [candidate]
        if not candidate.is_absolute():
            paths.extend(
                [
                    Path.cwd() / candidate,
                    Path.cwd() / "data" / "validation" / candidate,
                ]
            )
        for path in paths:
            if path.is_file():
                return path.resolve()
        print("FASTA file not found. Preferred location: data/validation/.")


def mfeprimer_index_files(database: Path) -> list[Path]:
    """Return recognized MFEprimer index files that currently exist."""
    candidates = [
        Path(str(database) + ".primerqc.bin"),
        Path(str(database) + ".primerqc"),
        Path(str(database) + ".primerqc.fai"),
    ]
    return [path for path in candidates if path.is_file()]


def mfeprimer_index_status(database: Path) -> tuple[bool, str]:
    """
    Recognize common MFEprimer index layouts.

    This function only recognizes index files. The downstream report parser
    still intentionally requires the legacy MFEprimer 3.x full text report.
    """
    binary_index = Path(str(database) + ".primerqc.bin")
    legacy_index = Path(str(database) + ".primerqc")
    legacy_fai = Path(str(database) + ".primerqc.fai")

    if binary_index.is_file():
        return True, "binary index (.primerqc.bin)"
    if legacy_index.is_file() and legacy_fai.is_file():
        return True, "legacy index (.primerqc + .primerqc.fai)"
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
    database_workdir: Path,
    threads: int,
) -> Path:
    """
    Write a one-sequence-per-line FASTA and verify that complete headers and
    nucleotide content are preserved. Stale adjacent indexes are removed when
    the formatted FASTA changes.
    """
    mfeprimer_dir = database_workdir / "mfeprimer"
    mfeprimer_dir.mkdir(parents=True, exist_ok=True)
    safe_stem = re.sub(r"[^A-Za-z0-9._-]+", "_", database.stem)
    fixed_database = mfeprimer_dir / f"{safe_stem}_fixed.fasta"
    temporary = mfeprimer_dir / f".{safe_stem}_fixed.tmp.fasta"

    print("\n4. MFEprimer database formatting")
    print("================================")
    progress("Formatting database with SeqKit", "RUN")
    seqkit_threads = min(max(1, int(threads)), 4)
    run_command_to_file(
        [
            "seqkit",
            "seq",
            "-j",
            str(seqkit_threads),
            "-w",
            "0",
            str(database),
        ],
        temporary,
    )

    source_records = read_fasta_with_full_headers(database)
    formatted_records = read_fasta_with_full_headers(temporary)
    headers_preserved = [h for h, _ in source_records] == [h for h, _ in formatted_records]
    sequences_preserved = [s for _, s in source_records] == [s for _, s in formatted_records]
    if not headers_preserved or not sequences_preserved:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            "SeqKit formatting changed FASTA headers or sequence content unexpectedly. "
            "The formatted database will not be used."
        )

    reused_existing = fixed_database.is_file() and fixed_database.read_bytes() == temporary.read_bytes()
    if reused_existing:
        temporary.unlink(missing_ok=True)
    else:
        # Prevent stale indexes from being associated with changed FASTA content.
        for index_file in mfeprimer_index_files(fixed_database):
            index_file.unlink(missing_ok=True)
        temporary.replace(fixed_database)

    source_stats = fasta_statistics(database)
    fixed_stats = fasta_statistics(fixed_database)
    progress("Formatting database with SeqKit", "OK")

    print("\nFormatting validation")
    print("---------------------")
    print(f"Input FASTA        : {compact_path(database)}")
    print(f"Formatted FASTA    : {compact_path(fixed_database)}")
    print(f"Sequences          : {source_stats['sequences']} -> {fixed_stats['sequences']}")
    print("Headers preserved  : YES")
    print("Sequences preserved: YES")
    print("IUPAC preserved    : YES")
    print("Modification       : FASTA line wrapping only")
    print(
        "Status             : "
        + (
            "EXISTING FORMATTED FASTA REUSED"
            if reused_existing
            else "FORMATTED FASTA CREATED/UPDATED"
        )
    )
    return fixed_database.resolve()


def index_target_database(database: Path) -> dict[str, object]:
    """
    Ensure that the formatted validation database has an MFEprimer index and
    explicitly report whether it was created or reused.
    """
    print("\n5. MFEprimer database indexing")
    print("==============================")

    indexed, detail = mfeprimer_index_status(database)
    existing_files = mfeprimer_index_files(database)

    if indexed:
        print("Existing MFEprimer index detected.")
        print("Status   : EXISTING INDEX REUSED")
        print(f"Database : {compact_path(database)}")
        print(f"Layout   : {detail}")
        print("\nIndex files detected")
        print("--------------------")
        for path in existing_files:
            print(f"✓ {compact_path(path)}")
        print("\nNo re-indexing was necessary.")
        return {
            "status": "REUSED",
            "detail": detail,
            "files": existing_files,
        }

    print("No complete MFEprimer index was detected for this formatted FASTA.")
    if not ask_yes_no("Index this database with MFEprimer now?"):
        raise RuntimeError(
            "Coverage validation requires an indexed MFEprimer database."
        )

    index_mode = ask_choice(
        "MFEprimer indexing parameters",
        [
            "Use MFEprimer default indexing parameters.",
            "Set a custom k-mer seed length (-k).",
        ],
    )

    command = ["mfeprimer", "index", "-i", str(database)]
    if index_mode == 2:
        k_value = ask_int("MFEprimer index k-mer length (-k)", 9, 15)
        command.extend(["-k", str(k_value)])

    progress("Indexing database with MFEprimer", "RUN")
    run_command(command)
    progress("Indexing database with MFEprimer", "OK")

    indexed, detail = mfeprimer_index_status(database)
    generated_files = mfeprimer_index_files(database)

    if not indexed:
        raise RuntimeError(
            "MFEprimer finished indexing, but the expected index layout "
            "could not be recognized. Check the installed MFEprimer version."
        )

    print("\nIndexing result")
    print("---------------")
    print("Status   : INDEX CREATED SUCCESSFULLY")
    print(f"Database : {compact_path(database)}")
    print(f"Layout   : {detail}")

    print("\nGenerated index files")
    print("---------------------")
    for path in generated_files:
        print(f"✓ {compact_path(path)}")

    print("\nMFEprimer database is ready for validation.")

    return {
        "status": "CREATED",
        "detail": detail,
        "files": generated_files,
    }
