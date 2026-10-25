# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

import argparse
import sys

def main() -> int:
    parser = argparse.ArgumentParser(
        prog="python -m pcr_workflow",
        description="PCR assay design and in silico validation workflows.",
    )
    parser.add_argument("workflow", choices=("design", "validate"))
    forwarded = sys.argv[1:]
    if not forwarded:
        parser.print_help()
        return 2
    if forwarded[0] in {"-h", "--help"}:
        parser.parse_args(forwarded)
        return 0

    args = parser.parse_args(forwarded[:1])
    sys.argv = [sys.argv[0], *forwarded[1:]]

    if args.workflow == "design":
        from .design.pipeline import main as run
    else:
        from .validation.pipeline import main as run
    return run()
