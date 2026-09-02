# Changelog

All notable changes to this project will be documented in this file.

## [1.2.0] - 2026-09-02

### Added

- Added `02_in_silico_validation.py` as a dedicated in silico validation workflow.
- Added automatic detection and validation of VarVAMP `SINGLE`, `QPCR`, and `TILED` assay modes.
- Added MFEprimer-based primer-pair target-coverage analysis.
- Added classification of valid target products, self-priming products, cross-scheme products, and unclassified products.
- Added extraction of valid MFEprimer target amplicons.
- Added unique-HitID primer coverage statistics to distinguish biological target coverage from technical amplicon counts.
- Added local BLAST+ `blastn-short` probe validation for QPCR assays.
- Added automatic expansion of IUPAC-degenerate qPCR probes into all concrete A/C/G/T-compatible variants.
- Added independent match/mismatch statistics for each concrete probe variant.
- Added `probe_variant_statistics.tsv`.
- Added `probe_variant_hitid_matrix.tsv`.
- Added best-across-all-probe-variants classification for each PCR-positive HitID.
- Added probe mismatch categories:
  - exact match;
  - 1 mismatch;
  - 2 mismatches;
  - 3 or more mismatches;
  - partial/gapped alignment;
  - no hit.
- Added full-length ungapped probe-site statistics.
- Added global qPCR coverage statistics using the complete validation dataset as denominator:
  - PCR;
  - PCR + exact probe;
  - PCR + probe with <=1 mismatch;
  - PCR + probe with <=2 mismatches;
  - PCR + full-length probe site.
- Added visual probe/target alignment reports.
- Added mismatch-position tables for probe analysis.
- Added `validation.log`.
- Added `--verbose` mode for detailed validation commands and file paths.
- Added `assay_design_manifest.json` to Workflow 01 for reproducible handoff to Workflow 02.
- Added `workflow_summary.txt`.
- Added `workflow.log`.
- Added software-version capture in the assay-design manifest.
- Added sequence-orientation normalization options:
  - keep original orientation;
  - MAFFT `--adjustdirection`;
  - MAFFT `--adjustdirectionaccurately`.
- Added sequence-topology selection:
  - linear;
  - circular.
- Added optional MARS cyclic start-position normalization for complete circular sequences.
- Added SeqKit exact sequence deduplication.
- Added explicit redundancy strategy selection:
  - none;
  - SeqKit;
  - CD-HIT-EST.
- Added alignment quality-control statistics before trimming.
- Added optional trimAl trimming with:
  - `-noallgaps`;
  - `-gappyout`;
  - `-automated1`;
  - manual `-gt` threshold.
- Added post-trimming alignment quality-control statistics.
- Added support for all three VarVAMP assay modes:
  - `single`;
  - `qpcr`;
  - `tiled`.
- Added mode-specific VarVAMP parameter handling.
- Added compact VarVAMP output summaries in normal terminal mode.
- Added `--verbose` support to Workflow 01.

### Changed

- Renamed the primary assay-design workflow to `01_varvamp_assay_design.py`.
- Reorganized the project into two explicit workflow stages:
  1. assay design;
  2. in silico validation.
- Expanded the project scope from qPCR-only design to VarVAMP `SINGLE`, `QPCR`, and `TILED` assay design.
- Improved terminal output to emphasize scientific results while hiding long external commands by default.
- External commands are now recorded in log files even when not displayed in the terminal.
- VarVAMP result tables are summarized instead of being fully printed in normal mode.
- Final qPCR comparison now clearly distinguishes primer-only PCR coverage from primer-plus-probe coverage.
- qPCR probe percentages in the final comparison now use the complete validation database as denominator.
- Probe-level statistics now distinguish:
  - per-concrete-variant performance;
  - best result across all concrete variants.
- Updated output organization to include validation intermediates, summaries, and machine-readable reports.
- Updated `README.md` to document both workflows.
- Updated `environment.yml` to include dependencies required by both workflows.
- Updated `CITATION.cff` for version 1.2.0 and the expanded project scope.
- Updated `CONTRIBUTING.md` to include tests and reporting requirements for both workflows.
- Updated `LICENSE` formatting while retaining GNU GPL v3 or later.

### Improved

- Reproducibility through machine-readable workflow metadata.
- Traceability of scientific parameter selections.
- Separation between assay design and assay validation.
- Distinction between predicted amplicon counts and unique biological target coverage.
- Handling of degenerate qPCR probes.
- Interpretation of probe mismatch statistics.
- Transparency of qPCR coverage denominators.
- Portability through relative project paths and explicit executable handling.
- Error handling for missing tools, invalid inputs, unsupported parameters, and empty validation results.
- Scientific reporting of MFEprimer and BLAST results.
- Readability of command-line output.

## [1.1.0] - 2026-07-25

### Added

- Interactive CD-HIT-EST parameter checkpoints.
- Optional CD-HIT-EST dereplication.
- Explicit CD-HIT-EST identity threshold selection.
- CD-HIT-EST fast (`-g 0`) and accurate (`-g 1`) clustering modes.
- CD-HIT-EST strand comparison selection using `-r`.
- Automatic CD-HIT-EST word-size selection according to the chosen identity threshold.
- Interactive MAFFT alignment strategy selection.
- Support for the following MAFFT strategies:
  - Auto;
  - FFT-NS-1;
  - FFT-NS-2;
  - FFT-NS-i with 2 refinement cycles;
  - FFT-NS-i with up to 1000 refinement cycles;
  - NW-NS-2;
  - NW-NS-i with 2 refinement cycles;
  - NW-NS-i with up to 1000 refinement cycles;
  - L-INS-i;
  - G-INS-i;
  - E-INS-i;
  - NW-NS-PartTree-1.
- Warnings for computationally intensive MAFFT strategies when the dataset size is outside their typical use range.
- Interactive VarVAMP parameter checkpoints.
- Explicit VarVAMP consensus threshold (`-t`) selection.
- Explicit primer ambiguity (`-a`) selection.
- Explicit probe ambiguity (`-pa`) selection.
- Final scientific parameter summary before workflow execution.
- Explicit user confirmation before starting the analysis.
- Support for fully non-interactive execution through command-line arguments.

### Changed

- Removed fixed interactive defaults for the main CD-HIT-EST and VarVAMP scientific parameters.
- MAFFT is no longer restricted to the `--auto` strategy.
- CD-HIT-EST execution now exposes clustering mode and strand comparison options.
- The workflow now requires explicit scientific parameter selection during interactive execution.
- Documentation updated to reflect the new parameter checkpoints and MAFFT strategies.
- Repository structure updated to remove the notebook component.
- Conda environment simplified by removing JupyterLab and IPython kernel dependencies.

### Improved

- Reproducibility by requiring users to explicitly review the scientific parameters used for each analysis.
- Input validation for interactive numerical parameters.
- Error handling for invalid parameter selections.
- Transparency of CD-HIT-EST, MAFFT, and VarVAMP execution settings.

## [1.0.0] - 2026-07-22

### Added

- Generic nucleotide FASTA input selection.
- CD-HIT-EST sequence dereplication.
- MAFFT multiple sequence alignment.
- Per-position conservation analysis.
- Conservative consensus generation.
- VarVAMP qPCR primer and probe design.
- Result collection and graphical exports.
- English command-line interface and documentation.