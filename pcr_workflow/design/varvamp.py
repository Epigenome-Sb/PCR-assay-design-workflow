# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import argparse
import re
import subprocess

from ..core import runtime
from ..core.runtime import log_line
from .fasta import validate_fraction
from .varvamp_config import validate_custom_varvamp_config

def build_varvamp_command(
    args: argparse.Namespace,
    mode: str,
    alignment_file: Path,
    output_dir: Path,
) -> list[str]:
    """Build a VarVAMP command using validated mode-specific arguments."""

    command = [
        args.varvamp_executable,
        mode,
        "-t",
        str(args.varvamp_threshold),
        "-a",
        str(args.primer_ambiguity),
        "-th",
        str(args.threads),
    ]

    if mode == "qpcr":
        command.extend(["-pa", str(args.probe_ambiguity)])
        if args.qpcr_test_n is not None:
            command.extend(["-n", str(args.qpcr_test_n)])
        if args.qpcr_deltag is not None:
            command.extend(["-d", str(int(args.qpcr_deltag))])

    elif mode == "single":
        if args.single_opt_length is not None:
            command.extend(["-ol", str(args.single_opt_length)])
        if args.single_max_length is not None:
            command.extend(["-ml", str(args.single_max_length)])
        if args.single_report_n is not None:
            command.extend(["-n", str(args.single_report_n)])

    elif mode == "tiled":
        if args.tiled_opt_length is not None:
            command.extend(["-ol", str(args.tiled_opt_length)])
        if args.tiled_max_length is not None:
            command.extend(["-ml", str(args.tiled_max_length)])
        if args.tiled_overlap is not None:
            command.extend(["-o", str(args.tiled_overlap)])

    command.extend([str(alignment_file), str(output_dir)])
    return command
def validate_varvamp_mode_parameters(args: argparse.Namespace, mode: str) -> None:
    """Validate common and mode-specific VarVAMP parameters before execution."""

    if args.varvamp_threshold is None:
        raise RuntimeError("VarVAMP consensus threshold (-t) is required.")
    validate_fraction(args.varvamp_threshold, "--varvamp-threshold")

    if args.primer_ambiguity is None:
        raise RuntimeError("Primer ambiguity (-a) is required.")
    if not isinstance(args.primer_ambiguity, int) or args.primer_ambiguity < 0:
        raise ValueError("--primer-ambiguity must be a non-negative integer.")

    if mode == "qpcr":
        if args.probe_ambiguity is None:
            raise RuntimeError("Probe ambiguity (-pa) is required for qPCR mode.")
        if not isinstance(args.probe_ambiguity, int) or args.probe_ambiguity < 0:
            raise ValueError("--probe-ambiguity must be a non-negative integer.")

        if args.qpcr_test_n is not None:
            if not isinstance(args.qpcr_test_n, int) or args.qpcr_test_n <= 0:
                raise ValueError("--qpcr-test-n must be a positive integer.")

        if args.qpcr_deltag is not None and not isinstance(args.qpcr_deltag, int):
            raise ValueError(
                "--qpcr-deltag must be an integer, for example -3 (not -3.0)."
            )

    if mode == "single":
        for name, value in (
            ("--single-opt-length", args.single_opt_length),
            ("--single-max-length", args.single_max_length),
        ):
            if value is not None and (
                not isinstance(value, int) or value <= 0
            ):
                raise ValueError(f"{name} must be a positive integer.")

        if (
            args.single_opt_length is not None
            and args.single_max_length is not None
            and args.single_max_length < args.single_opt_length
        ):
            raise ValueError(
                "--single-max-length must be >= --single-opt-length."
            )

        if args.single_report_n is not None:
            if str(args.single_report_n).lower() != "inf":
                try:
                    report_n = int(args.single_report_n)
                except ValueError as error:
                    raise ValueError(
                        "--single-report-n must be a positive integer or 'inf'."
                    ) from error
                if report_n <= 0:
                    raise ValueError(
                        "--single-report-n must be a positive integer or 'inf'."
                    )

    if mode == "tiled":
        for name, value in (
            ("--tiled-opt-length", args.tiled_opt_length),
            ("--tiled-max-length", args.tiled_max_length),
        ):
            if value is not None and (
                not isinstance(value, int) or value <= 0
            ):
                raise ValueError(f"{name} must be a positive integer.")

        if (
            args.tiled_opt_length is not None
            and args.tiled_max_length is not None
            and args.tiled_max_length < args.tiled_opt_length
        ):
            raise ValueError(
                "--tiled-max-length must be >= --tiled-opt-length."
            )

        if args.tiled_overlap is not None and (
            not isinstance(args.tiled_overlap, int)
            or args.tiled_overlap < 0
        ):
            raise ValueError("--tiled-overlap must be a non-negative integer.")

    if args.varvamp_config is not None:
        validate_custom_varvamp_config(
            args.varvamp_config.expanduser().resolve(),
            mode,
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

    if runtime.VERBOSE:
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
