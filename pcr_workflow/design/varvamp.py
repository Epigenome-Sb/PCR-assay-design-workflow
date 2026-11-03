# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""Run VarVAMP, inspect failures, and display mode-specific outputs."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

from . import run_management as _run_management
from .run_management import compact_path, log_line, section
from .questions import (
    ask_required_choice,
    ask_required_float,
    ask_required_non_negative_integer,
    ask_yes_no,
)
from .varvamp_config import (
    choose_varvamp_mode,
    configure_custom_varvamp_file,
    configure_interactive_varvamp_advanced,
)


def run_varvamp_command(
    command: list[str],
    *,
    env: dict[str, str],
) -> tuple[int, str]:
    """
    Run VarVAMP while streaming its combined stdout/stderr to the terminal.

    Unlike run_command(), a non-zero return code is returned to the caller so
    interactive mode can adjust parameters and retry VarVAMP without restarting
    orientation, MARS, clustering, MAFFT, trimming or conservation analysis.
    """
    command_text = " ".join(command)
    log_line("COMMAND: " + command_text)

    if _run_management.VERBOSE:
        print("\nCommand:", command_text)

    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=env,
        )
    except OSError as error:
        raise RuntimeError(
            "Could not start VarVAMP command: " + command_text
        ) from error

    output_lines: list[str] = []
    assert process.stdout is not None

    for line in process.stdout:
        output_lines.append(line)
        print(line, end="")
        log_line("VARVAMP: " + line.rstrip())

    return_code = process.wait()
    return return_code, "".join(output_lines)


def parse_varvamp_diagnostics(output: str) -> dict[str, object]:
    """Extract useful counts and the reported error from VarVAMP terminal output."""
    diagnostics: dict[str, object] = {}

    patterns = {
        "forward_primers": r"(\d+)\s+fw\s+and\s+\d+\s+rv potential primers",
        "reverse_primers": r"\d+\s+fw\s+and\s+(\d+)\s+rv potential primers",
        "potential_probes": r"(\d+)\s+potential qPCR probes",
        "unique_qpcr_amplicons": r"(\d+)\s+unique amplicons with internal probe",
        "delta_g_schemes": (
            r"(\d+)\s+non-overlapping qPCR schemes that passed deltaG cutoff"
        ),
    }

    for key, pattern in patterns.items():
        matches = re.findall(pattern, output, flags=re.IGNORECASE)
        if matches:
            diagnostics[key] = int(matches[-1])

    error_lines = []
    for line in output.splitlines():
        stripped = line.strip()
        if stripped.upper().startswith("ERROR:"):
            error_lines.append(stripped[6:].strip())

    if error_lines:
        diagnostics["reason"] = error_lines[-1]
    else:
        nonempty = [line.strip() for line in output.splitlines() if line.strip()]
        diagnostics["reason"] = (
            nonempty[-1] if nonempty else "VarVAMP exited with a non-zero status."
        )

    return diagnostics


def print_varvamp_failure(
    args: argparse.Namespace,
    mode: str,
    diagnostics: dict[str, object],
) -> None:
    """Display a mode-aware failure summary before offering a retry."""
    title = f"VarVAMP did not produce a {mode.upper()} scheme"
    print("\n" + title)
    print("=" * len(title))
    print("\nReason reported by VarVAMP:")
    print(str(diagnostics.get("reason", "Unknown VarVAMP error.")))

    print("\nCurrent parameters:")
    print(f"  Mode                 : {mode.upper()}")
    print(f"  Consensus threshold  : {args.varvamp_threshold}")
    print(f"  Primer ambiguity     : {args.primer_ambiguity}")
    if mode == "qpcr":
        print(f"  Probe ambiguity      : {args.probe_ambiguity}")
        if args.qpcr_test_n is not None:
            print(f"  qPCR test count (-n) : {args.qpcr_test_n}")
        if args.qpcr_deltag is not None:
            print(f"  deltaG cutoff (-d)   : {args.qpcr_deltag}")

    if mode == "qpcr" and any(
        key in diagnostics
        for key in (
            "forward_primers",
            "reverse_primers",
            "potential_probes",
            "unique_qpcr_amplicons",
            "delta_g_schemes",
        )
    ):
        print("\nVarVAMP nevertheless found:")
        if "forward_primers" in diagnostics:
            print(f"  Forward primers      : {diagnostics['forward_primers']}")
        if "reverse_primers" in diagnostics:
            print(f"  Reverse primers      : {diagnostics['reverse_primers']}")
        if "potential_probes" in diagnostics:
            print(f"  Potential probes     : {diagnostics['potential_probes']}")
        if "unique_qpcr_amplicons" in diagnostics:
            print(
                "  Probe-containing amps : "
                f"{diagnostics['unique_qpcr_amplicons']}"
            )
        if "delta_g_schemes" in diagnostics:
            print(f"  deltaG-passing schemes: {diagnostics['delta_g_schemes']}")

    print("\nPossible adjustment:")
    print("  - lower or otherwise adjust consensus threshold (-t)")
    print("  - increase primer ambiguity (-a)")
    if mode == "qpcr":
        print("  - increase probe ambiguity (-pa)")
    print("  - change mode-specific parameters or use a custom config")


def ensure_mode_required_interactive_parameters(
    args: argparse.Namespace,
    mode: str,
) -> None:
    """Collect required parameters after switching VarVAMP mode."""
    if mode == "qpcr" and args.probe_ambiguity is None:
        if not sys.stdin.isatty():
            raise RuntimeError(
                "--probe-ambiguity is required for qPCR mode "
                "in non-interactive execution."
            )
        args.probe_ambiguity = ask_required_non_negative_integer(
            "Maximum ambiguous bases in the probe",
            "-pa",
        )


def adjust_varvamp_after_failure(
    args: argparse.Namespace,
    mode: str,
    config_dir: Path,
) -> str:
    """
    Adjust parameters after a failed VarVAMP attempt.

    Returns:
        "retry"       retry VarVAMP
        "mode_changed" retry after changing mode
        "stop"        stop VarVAMP design
    """
    if mode == "qpcr":
        choices = [
            ("threshold", "Change consensus threshold (-t)."),
            ("primer", "Change primer ambiguity (-a)."),
            ("probe", "Change probe ambiguity (-pa)."),
            ("several", "Change several VarVAMP parameters."),
            ("advanced", "Configure advanced qPCR parameters (Python config)."),
            ("mode", "Change VarVAMP assay mode."),
            ("stop", "Stop the workflow."),
        ]
    else:
        choices = [
            ("threshold", "Change consensus threshold (-t)."),
            ("primer", "Change primer ambiguity (-a)."),
            (
                "mode_cli",
                f"Change additional {mode.upper()} command-line parameters.",
            ),
            ("several", "Change several VarVAMP parameters."),
            (
                "advanced",
                f"Configure advanced {mode.upper()} parameters (Python config).",
            ),
            ("mode", "Change VarVAMP assay mode."),
            ("stop", "Stop the workflow."),
        ]

    choice = ask_required_choice("What would you like to do?", choices)

    if choice == "threshold":
        args.varvamp_threshold = ask_required_float(
            "VarVAMP consensus threshold",
            "-t",
            minimum=0.01,
            maximum=1.0,
        )
        return "retry"

    if choice == "primer":
        args.primer_ambiguity = ask_required_non_negative_integer(
            "Maximum ambiguous bases in each primer",
            "-a",
        )
        return "retry"

    if choice == "probe":
        args.probe_ambiguity = ask_required_non_negative_integer(
            "Maximum ambiguous bases in the probe",
            "-pa",
        )
        return "retry"

    if choice == "mode_cli":
        configure_interactive_varvamp_advanced(
            args,
            mode,
            force=True,
        )
        return "retry"

    if choice == "several":
        args.varvamp_threshold = ask_required_float(
            "VarVAMP consensus threshold",
            "-t",
            minimum=0.01,
            maximum=1.0,
        )
        args.primer_ambiguity = ask_required_non_negative_integer(
            "Maximum ambiguous bases in each primer",
            "-a",
        )
        if mode == "qpcr":
            args.probe_ambiguity = ask_required_non_negative_integer(
                "Maximum ambiguous bases in the probe",
                "-pa",
            )
        if ask_yes_no(
            "Also configure additional mode-specific command-line parameters?"
        ):
            configure_interactive_varvamp_advanced(
                args,
                mode,
                force=True,
            )
        return "retry"

    if choice == "advanced":
        changed = configure_custom_varvamp_file(
            args,
            mode,
            config_dir,
        )
        if changed:
            return "retry"
        return adjust_varvamp_after_failure(
            args,
            mode,
            config_dir,
        )

    if choice == "mode":
        args.varvamp_mode = choose_varvamp_mode()
        # A mode change returns to VarVAMP's normal/default Python config unless
        # the user explicitly selects a custom config again.
        args.varvamp_config = None
        ensure_mode_required_interactive_parameters(
            args,
            args.varvamp_mode,
        )
        if ask_yes_no(
            "Configure additional mode-specific VarVAMP parameters for the new mode?"
        ):
            configure_interactive_varvamp_advanced(
                args,
                args.varvamp_mode,
                force=True,
            )
        return "mode_changed"

    return "stop"


def copy_varvamp_tree(varvamp_dir: Path, results_dir: Path, mode: str) -> Path:
    """Copy the complete successful VarVAMP tree into this isolated run."""

    destination = results_dir / "varvamp"
    if destination.exists():
        raise RuntimeError(
            "The VarVAMP result directory already exists inside this run: "
            f"{destination}. Refusing to overwrite an existing result."
        )
    shutil.copytree(varvamp_dir, destination)
    (destination / "assay_mode.txt").write_text(
        mode.upper() + "\n",
        encoding="utf-8",
    )
    return destination


def _varvamp_result_file(directory: Path, *filenames: str) -> Path | None:
    """Find the first available name variant of a VarVAMP output."""
    for filename in filenames:
        matches = sorted(path for path in directory.rglob(filename) if path.is_file())
        if matches:
            return matches[0]
    return None


def _varvamp_tsv(path: Path | None) -> list[dict[str, str]]:
    """Keep native VarVAMP values; normalize only column names for display."""
    if path is None:
        return []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return [
            {key.strip().lower(): (value or "").strip()
             for key, value in row.items() if key is not None}
            for row in csv.DictReader(stream, delimiter="\t")
        ]


def _varvamp_tiled_coverage(directory: Path) -> str:
    """Read reported coverage; never infer it from consensus/BED coordinates."""
    pattern = re.compile(
        r"(\d+(?:\.\d+)?)\s*%\s*total coverage with\s*\d+\s*amplicons?",
        flags=re.IGNORECASE,
    )
    for path in (directory / "varvamp_log.txt", directory.parent / "workflow.log"):
        if path.is_file():
            matches = pattern.findall(path.read_text(encoding="utf-8", errors="replace"))
            if matches:
                return f"{matches[-1]} % (VarVAMP reported)"
    return "not available in saved logs"


def print_varvamp_tables(varvamp_results_dir: Path, mode: str) -> None:
    """Display a consistent, non-destructive report for SINGLE, QPCR and TILED."""
    mode = mode.lower()
    directory = Path(varvamp_results_dir)
    section(f"VarVAMP {mode.upper()} design results")

    if mode not in {"single", "qpcr", "tiled"}:
        print(f"Unknown VarVAMP mode: {mode}")
        return
    if not directory.is_dir():
        print(f"Result directory not found: {compact_path(directory)}")
        return

    if mode == "qpcr":
        primers_path = _varvamp_result_file(directory, "qpcr_primers.tsv")
        design_path = _varvamp_result_file(directory, "qpcr_design.tsv")
        fasta_paths = [_varvamp_result_file(directory, "oligos.fasta")]
        assignments_path = None
    else:
        primers_path = _varvamp_result_file(directory, "primers.tsv", "primer.tsv")
        design_path = None
        assignments_path = _varvamp_result_file(
            directory,
            "primer_to_amplicon_assignment.tabular",
            "primer_to_amplicon_assignments.tabular",
        )
        fasta_paths = (
            sorted(directory.rglob("primers_pool_*.fasta")) if mode == "tiled"
            else [_varvamp_result_file(directory, "primers.fasta")]
        )
    fasta_paths = [path for path in fasta_paths if path and path.is_file()]

    # ONE shared structure for the three modes, using actual VarVAMP fields.
    assays: dict[str, dict[str, object]] = {}
    for row in _varvamp_tsv(design_path):
        name = row.get("qpcr_scheme", "")
        if name:
            assays[name] = {
                "length": row.get("length", ""), "start": row.get("start", ""),
                "stop": row.get("stop", ""), "pool": "",
                "penalty": row.get("penalty", ""),
                "deltag": row.get("deltag", ""),
                "offtargets": row.get("off_target_amplicons", ""),
                "oligos": [],
            }

    for row in _varvamp_tsv(primers_path):
        primer_name = row.get("primer_name", "")
        oligo_type = row.get("oligo_type", "") or primer_name.rsplit("_", 1)[-1]
        name = (
            row.get("qpcr_scheme") or row.get("amplicon_name")
            or row.get("amlicon_name")  # VarVAMP 1.3.2 spelling.
            or (primer_name.rsplit("_", 1)[0] if "_" in primer_name else "")
        )
        if not name:
            continue
        assay = assays.setdefault(name, {
            "length": "", "start": "", "stop": "", "pool": "",
            "penalty": "", "deltag": "", "offtargets": "", "oligos": [],
        })
        if mode != "qpcr":
            assay["length"] = assay["length"] or row.get("amplicon_length", "")
            assay["pool"] = assay["pool"] or row.get("pool", "")
            assay["offtargets"] = assay["offtargets"] or row.get("off_target_amplicons", "")
        assay["oligos"].append({
            "type": oligo_type.upper(), "seq": row.get("seq", ""),
            "size": row.get("size", ""), "start": row.get("start", ""),
            "stop": row.get("stop", ""), "tm": row.get("mean_temp", ""),
            "gc": row.get("mean_gc", ""), "penalty": row.get("penalty", ""),
        })

    if not assays:
        print("No assay data found. Check expected TSV files in this output directory.")
    else:
        if mode != "qpcr":
            for assay in assays.values():
                starts = [int(o["start"]) for o in assay["oligos"] if o["start"].isdigit()]
                stops = [int(o["stop"]) for o in assay["oligos"] if o["stop"].isdigit()]
                if starts and stops:
                    assay["start"] = str(min(starts))
                    assay["stop"] = str(max(stops))

        def display(value: object) -> str:
            return str(value) if value not in ("", None) else "—"

        oligos = [oligo for assay in assays.values() for oligo in assay["oligos"]]
        primer_count = sum(o["type"] in {"LEFT", "RIGHT"} for o in oligos)
        probe_count = sum(o["type"] == "PROBE" for o in oligos)
        pools = sorted({str(assay["pool"]) for assay in assays.values() if assay["pool"]})

        print("\nDesign summary")
        print("--------------")
        print(f"  Assays / amplicons  : {len(assays)}")
        print(f"  Primer oligos       : {primer_count}")
        print(f"  Probe oligos        : {probe_count}")
        print(f"  Pools               : {', '.join(pools) if pools else 'N/A'}")
        if mode == "tiled":
            print(f"  Amplicon coverage   : {_varvamp_tiled_coverage(directory)}")

        for index, (name, assay) in enumerate(assays.items(), 1):
            print(f"\nAssay {index} — {name}")
            print("-" * (len(name) + len(str(index)) + 9))
            print(f"  Amplicon length     : {display(assay['length'])} bp")
            print(f"  Amplicon positions  : {display(assay['start'])}–{display(assay['stop'])}")
            print(f"  Pool                : {display(assay['pool'])}")
            print(f"  Scheme penalty      : {display(assay['penalty'])}")
            print(f"  Scheme deltaG       : {display(assay['deltag'])}")
            print(f"  Off-target amplicons: {display(assay['offtargets'])}")
            print("  Oligonucleotides:")
            types = {"LEFT": 0, "PROBE": 1, "RIGHT": 2}
            for oligo in sorted(assay["oligos"], key=lambda o: types.get(o["type"], 3)):
                print(
                    f"    {oligo['type']:<5} {display(oligo['start'])}–{display(oligo['stop'])}"
                    f" | {display(oligo['size'])} nt | Tm {display(oligo['tm'])} °C"
                    f" | GC {display(oligo['gc'])} % | penalty {display(oligo['penalty'])}"
                )
                print(f"      5'-{display(oligo['seq'])}-3'")

    print("\nPrincipal output files")
    print("----------------------")
    for path in (design_path, primers_path, assignments_path, *fasta_paths):
        if path:
            print(f"  {path.relative_to(directory)}")
    print(f"\nComplete VarVAMP result directory: {compact_path(directory)}")
