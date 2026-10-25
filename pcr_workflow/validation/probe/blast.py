# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import re

from ...core import runtime
from ...core.runtime import section
from ..commands import run_command
from .iupac import expand_iupac_probe

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
    threads: int,
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
        "-num_threads",
        str(threads),
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
    if runtime.VERBOSE:
        print("BLAST task        : blastn-short")
        print("Strands           : both")
        print("DUST              : disabled")
        print("E-value           : 1000")
        print(f"Threads           : {threads}")

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
