# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import argparse

def parse_cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Mode-aware in silico validation of VarVAMP SINGLE, QPCR or TILED "
            "assay designs with MFEprimer and, for QPCR, BLAST+ probe analysis."
        )
    )
    parser.add_argument(
        "--design-manifest",
        type=Path,
        default=None,
        help=(
            "Workflow 01 assay_design_manifest.json. When omitted, completed "
            "design runs under results/<project>/design/<run_id>/ are listed."
        ),
    )
    parser.add_argument(
        "--validation-db",
        type=Path,
        default=None,
        help=(
            "Validation FASTA. When omitted, FASTA files under data/validation/ "
            "are listed."
        ),
    )
    parser.add_argument(
        "--run-id",
        type=str,
        default=None,
        help="Optional validation run identifier; default is a local timestamp.",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=None,
        help=(
            "Optional exact work directory override. Default: "
            "work/<project>/validation/<run_id>/"
        ),
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=None,
        help=(
            "Optional exact results directory override. Default: "
            "results/<project>/validation/<run_id>/"
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help=(
            "CPU threads for tools that expose a thread option. If omitted, "
            "the workflow detects the logical CPU threads available to the "
            "process and keeps two threads free when possible."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show external commands and detailed paths in the terminal.",
    )
    return parser.parse_args()
