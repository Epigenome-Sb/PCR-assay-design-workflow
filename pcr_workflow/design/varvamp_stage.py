# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import json
import os

from ..core.paths import compact_path, portable_manifest_path
from ..core.runtime import log_line, section
from .cli import ask_required_float, ask_required_non_negative_integer
from .commands import check_required_tools
from .reporting import copy_varvamp_tree, create_varvamp_visual_summaries, print_varvamp_tables
from .ui import adjust_varvamp_after_failure, choose_varvamp_mode, configure_interactive_varvamp_advanced, ensure_mode_required_interactive_parameters
from .varvamp import build_varvamp_command, parse_varvamp_diagnostics, print_varvamp_failure, run_varvamp_command, validate_varvamp_mode_parameters
from .varvamp_config import stage_varvamp_config_for_run

def run_varvamp_stage(
    args,
    interactive: bool,
    downstream_alignment: Path,
    varvamp_work_root: Path,
    config_results_dir: Path,
    results_dir: Path,
) -> tuple[Path | None, list[dict[str, object]]]:
    # ------------------------------------------------------------------
    # 7. VarVAMP assay design with retry-on-failure
    # ------------------------------------------------------------------
    varvamp_results: Path | None = None
    varvamp_attempts: list[dict[str, object]] = []

    if not args.skip_varvamp:
        print("\n7. VarVAMP assay design checkpoint")

        if args.varvamp_mode is None:
            if not interactive:
                raise RuntimeError(
                    "--varvamp-mode is required in non-interactive mode."
                )
            args.varvamp_mode = choose_varvamp_mode()

        if args.varvamp_threshold is None:
            if not interactive:
                raise RuntimeError(
                    "--varvamp-threshold is required."
                )
            args.varvamp_threshold = ask_required_float(
                "VarVAMP consensus threshold",
                "-t",
                minimum=0.01,
                maximum=1.0,
            )

        if args.primer_ambiguity is None:
            if not interactive:
                raise RuntimeError(
                    "--primer-ambiguity is required."
                )
            args.primer_ambiguity = (
                ask_required_non_negative_integer(
                    "Maximum ambiguous bases in each primer",
                    "-a",
                )
            )

        ensure_mode_required_interactive_parameters(
            args,
            args.varvamp_mode,
        )

        if interactive:
            configure_interactive_varvamp_advanced(
                args,
                args.varvamp_mode,
            )

        check_required_tools(
            [("VarVAMP", args.varvamp_executable)]
        )

        attempt_number = 0

        while True:
            attempt_number += 1
            mode = args.varvamp_mode

            ensure_mode_required_interactive_parameters(
                args,
                mode,
            )
            validate_varvamp_mode_parameters(args, mode)

            varvamp_workdir = (
                varvamp_work_root
                / f"attempt_{attempt_number:02d}_{mode}"
            )
            if varvamp_workdir.exists():
                raise RuntimeError(
                    "Unexpected existing VarVAMP attempt directory: "
                    f"{varvamp_workdir}. Refusing to overwrite it."
                )

            command = build_varvamp_command(
                args,
                mode,
                downstream_alignment,
                varvamp_workdir,
            )

            command_env = os.environ.copy()

            # Every new script execution starts with VarVAMP's normal
            # configuration unless a custom config is explicitly selected.
            command_env.pop("VARVAMP_CONFIG", None)

            config_path: Path | None = None
            if args.varvamp_config is not None:
                config_path = stage_varvamp_config_for_run(
                    args.varvamp_config,
                    config_results_dir,
                    mode,
                    attempt_number,
                )
                # Keep the run-local copy as the canonical provenance path.
                args.varvamp_config = config_path
                command_env["VARVAMP_CONFIG"] = str(config_path)

            section(f"VarVAMP attempt {attempt_number}")
            print(f"Mode                      : {mode.upper()}")
            print(
                f"Consensus threshold (-t)  : "
                f"{args.varvamp_threshold}"
            )
            print(
                f"Primer ambiguity (-a)     : "
                f"{args.primer_ambiguity}"
            )
            if mode == "qpcr":
                print(
                    f"Probe ambiguity (-pa)     : "
                    f"{args.probe_ambiguity}"
                )
            print(
                "Python configuration      : "
                + (
                    compact_path(config_path)
                    if config_path is not None
                    else "VarVAMP default"
                )
            )

            return_code, varvamp_output = run_varvamp_command(
                command,
                env=command_env,
            )

            diagnostics = parse_varvamp_diagnostics(
                varvamp_output
            )

            attempt_record: dict[str, object] = {
                "attempt": attempt_number,
                "mode": mode,
                "consensus_threshold": args.varvamp_threshold,
                "primer_ambiguity": args.primer_ambiguity,
                "probe_ambiguity": (
                    args.probe_ambiguity
                    if mode == "qpcr"
                    else None
                ),
                "single_opt_length": args.single_opt_length,
                "single_max_length": args.single_max_length,
                "single_report_n": args.single_report_n,
                "tiled_opt_length": args.tiled_opt_length,
                "tiled_max_length": args.tiled_max_length,
                "tiled_overlap": args.tiled_overlap,
                "qpcr_test_n": args.qpcr_test_n,
                "qpcr_deltag": args.qpcr_deltag,
                "custom_config": (
                    portable_manifest_path(config_path)
                    if config_path is not None
                    else None
                ),
                "return_code": return_code,
                "status": (
                    "success"
                    if return_code == 0
                    else "failed"
                ),
            }

            if return_code != 0:
                attempt_record["reason"] = diagnostics.get(
                    "reason"
                )

            varvamp_attempts.append(attempt_record)
            log_line(
                "VARVAMP_ATTEMPT: "
                + json.dumps(
                    attempt_record,
                    sort_keys=True,
                )
            )

            if return_code == 0:
                varvamp_results = copy_varvamp_tree(
                    varvamp_workdir,
                    results_dir,
                    mode,
                )
                print_varvamp_tables(
                    varvamp_results,
                    mode,
                )
                create_varvamp_visual_summaries(
                    varvamp_results
                )
                break

            if not interactive:
                raise RuntimeError(
                    "VarVAMP failed in non-interactive mode: "
                    + str(
                        diagnostics.get(
                            "reason",
                            "unknown error",
                        )
                    )
                )

            print_varvamp_failure(
                args,
                mode,
                diagnostics,
            )

            action = adjust_varvamp_after_failure(
                args,
                mode,
                config_results_dir,
            )

            if action == "stop":
                print(
                    "VarVAMP design stopped by the user. "
                    "Previous workflow steps remain saved."
                )
                args.skip_varvamp = True
                break

            # "retry" and "mode_changed" both return here and rerun only
            # the VarVAMP stage. All previous workflow products are reused.
            print("\nRetrying VarVAMP only...")
            print("Orientation              : already completed")
            print("Topology normalization   : already completed")
            print("Redundancy handling      : already completed")
            print("Final MAFFT              : already completed")
            print("trimAl / alignment QC    : already completed")
            print("Conservation analysis    : already completed")

    else:
        print("\n7. VarVAMP assay design skipped")
    return varvamp_results, varvamp_attempts
