# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

def progress(label: str, status: str = "OK") -> None:
    dots = "." * max(2, 48 - len(label))
    print(f"{label} {dots} {status}")
