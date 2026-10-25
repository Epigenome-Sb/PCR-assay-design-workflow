#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""Compatibility launcher for the modular PCR assay design workflow."""

from pcr_workflow.design.pipeline import main

if __name__ == "__main__":
    raise SystemExit(main())
