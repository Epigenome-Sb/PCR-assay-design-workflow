# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import re
import shutil

from ..core.paths import compact_path
from .commands import run_command, run_command_to_file
from .fasta import (
    fasta_statistics,
    read_fasta_with_full_headers,
    write_fasta_with_full_headers,
)
from .runtime import progress


def apply_sequence_length_filter(
    database: Path,
    database_workdir: Path,
    threshold: int | None,
    strategy: str,
    before_stats: dict[str, float | int],
) -> tuple[Path, str, dict[str, float | int]]:
    """Write the selected length-filtered FASTA and report retained counts."""
    retained_database = database
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
        database_workdir.mkdir(parents=True, exist_ok=True)
        safe_stem = re.sub(r"[^A-Za-z0-9._-]+", "_", database.stem)
        retained_database = database_workdir / f"{safe_stem}_min{threshold}.fasta"
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
    database_workdir.mkdir(parents=True, exist_ok=True)
    safe_stem = re.sub(r"[^A-Za-z0-9._-]+", "_", database.stem)
    fixed_database = database_workdir / f"{safe_stem}_fixed.fasta"
    temporary = database_workdir / f".{safe_stem}_fixed.tmp.fasta"

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


def create_mfeprimer_index(
    database: Path,
    k_value: int | None = None,
) -> dict[str, object]:
    """Create and verify an MFEprimer index for a formatted FASTA."""
    command = ["mfeprimer", "index", "-i", str(database)]
    if k_value is not None:
        command.extend(["-k", str(k_value)])
    run_command(command)

    indexed, detail = mfeprimer_index_status(database)
    generated_files = mfeprimer_index_files(database)
    if not indexed:
        raise RuntimeError(
            "MFEprimer finished indexing, but the expected index layout "
            "could not be recognized. Check the installed MFEprimer version."
        )
    return {
        "status": "CREATED",
        "detail": detail,
        "files": generated_files,
    }
