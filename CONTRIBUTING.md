# Contributing

Contributions are welcome through GitHub issues and pull requests.

This repository contains two main workflows:

- `01_varvamp_assay_design.py` — nucleotide preprocessing, alignment, conservation analysis, and VarVAMP assay design;
- `02_in_silico_validation.py` — in silico primer-pair and qPCR probe validation.

## Reporting a problem

Please include:

- operating system;
- Python version;
- workflow/script concerned (`01` or `02`);
- complete command used;
- scientific parameters selected;
- full error message;
- relevant software versions;
- a minimal non-confidential nucleotide FASTA example when possible.

Depending on the workflow, useful software-version information may include:

- MAFFT;
- MARS;
- SeqKit;
- CD-HIT-EST;
- trimAl;
- VarVAMP;
- MFEprimer;
- BLAST+.

For Workflow 01 interactive runs, please also report relevant selections such as:

- sequence-orientation strategy;
- sequence topology;
- redundancy strategy;
- CD-HIT-EST identity threshold (`-c`), if used;
- CD-HIT-EST clustering mode (`-g`), if used;
- CD-HIT-EST strand comparison (`-r`), if used;
- MAFFT alignment strategy;
- trimAl strategy and threshold, if used;
- VarVAMP mode (`single`, `qpcr`, or `tiled`);
- VarVAMP consensus threshold (`-t`);
- primer ambiguity (`-a`);
- probe ambiguity (`-pa`) for qPCR;
- mode-specific VarVAMP parameters.

For Workflow 02, please also report:

- detected VarVAMP mode;
- validation FASTA used;
- MFEprimer settings;
- target-amplicon extraction filter;
- selected assay(s);
- whether the qPCR probe contained IUPAC ambiguity;
- number of concrete probe variants generated;
- BLAST-related error messages, when applicable.

When possible, include:

```text
results/<project>/workflow.log
```

for Workflow 01, or:

```text
validation_inputs/validation.log
```

for Workflow 02.

Do not upload confidential, identifiable, restricted, sensitive, or unpublished
sequence data to a public GitHub issue.

## Pull requests

Before submitting a pull request:

1. create a focused branch;
2. keep changes limited to one clear purpose;
3. preserve command-line compatibility when possible;
4. update the README, CHANGELOG, and other documentation for user-visible changes;
5. verify Python syntax;
6. test the affected workflow on a small nucleotide FASTA dataset;
7. test interactive input validation when interactive behavior is modified;
8. test non-interactive execution when command-line behavior is modified;
9. verify that generated files are written to the expected directories;
10. avoid introducing dataset-specific assumptions.

### Workflow 01 tests

When modifying `01_varvamp_assay_design.py`, test the relevant affected steps.

Depending on the change, this may include:

- orientation kept unchanged;
- MAFFT orientation correction;
- linear topology;
- circular topology with MARS;
- no redundancy reduction;
- SeqKit exact deduplication;
- CD-HIT-EST clustering;
- at least one MAFFT alignment strategy;
- alignment QC;
- trimAl enabled and disabled;
- conservation analysis;
- VarVAMP `single` mode;
- VarVAMP `qpcr` mode;
- VarVAMP `tiled` mode;
- `--skip-varvamp`;
- generation of `assay_design_manifest.json`;
- generation of `workflow_summary.txt`;
- generation of `workflow.log`.

When MAFFT-related code is modified, testing more than one MAFFT strategy is encouraged.

When CD-HIT-EST-related code is modified, test at least one clustering mode and the relevant strand option.

When VarVAMP-related code is modified, test the affected assay mode and parameter handling.

### Workflow 02 tests

When modifying `02_in_silico_validation.py`, test the relevant affected branches.

Depending on the change, this may include:

- VarVAMP mode detection;
- target FASTA selection;
- SeqKit FASTA reformatting;
- MFEprimer database indexing;
- MFEprimer target-coverage analysis;
- self-priming classification;
- unique-HitID coverage calculation;
- qPCR probe BLAST with `blastn-short`;
- a non-degenerate probe;
- an IUPAC-degenerate probe;
- concrete-probe expansion;
- per-probe-variant statistics;
- HitID × probe-variant matrix generation;
- best-across-all-variants classification;
- global PCR + probe coverage calculation;
- `--verbose` output;
- validation log generation.

For degenerate probes, verify that each concrete version is analysed separately
and that final biological coverage is deduplicated by unique target HitID.

## Scientific parameter changes

Changes affecting scientific parameters must be documented clearly.

Contributors should explain:

- which parameter was added or modified;
- which external software option it corresponds to;
- accepted values;
- whether it affects interactive execution, command-line execution, or both;
- whether existing analyses may produce different results after the change.

Scientific parameters should not be silently changed without corresponding
documentation.

Defaults inherited from third-party software should not be redefined without
a clear scientific and technical justification.

## Code quality

Contributions should:

- remain compatible with nucleotide FASTA input;
- provide clear error messages for invalid input;
- avoid hard-coded user- or dataset-specific paths;
- avoid unnecessary dependencies;
- preserve reproducibility whenever possible;
- preserve unique sequence identifiers when biologically appropriate;
- distinguish technical amplicon counts from unique biological target counts;
- avoid presenting computational mismatch categories as experimental diagnostic performance.

Before committing Python changes, syntax can be checked with:

```bash
python -m py_compile 01_varvamp_assay_design.py
python -m py_compile 02_in_silico_validation.py
```

Command-line options can be checked with:

```bash
python 01_varvamp_assay_design.py --help
python 02_in_silico_validation.py --help
```

## Reproducibility

New features should preserve or improve the ability to reproduce an analysis.

When relevant, contributors should ensure that:

- commands are recorded in the appropriate log file;
- selected scientific parameters are saved;
- software versions are recorded when possible;
- output files have stable and informative names;
- Workflow 01 continues to generate `assay_design_manifest.json`;
- Workflow 02 reports the validation denominator explicitly;
- qPCR probe results distinguish per-variant statistics from best-of-all-variants coverage.

## Data and privacy

Example datasets used for testing should be public, synthetic, or otherwise
authorized for redistribution.

Do not commit:

- confidential sequence data;
- personally identifiable data;
- restricted datasets;
- unpublished data without authorization;
- credentials or API keys;
- institution-specific private paths;
- large generated intermediate files.

Generated `work/`, validation intermediates, and large analysis outputs should
only be committed when they are intentionally provided as small documented
examples.

## Third-party software

Contributions must not redistribute third-party software unless its license
explicitly permits redistribution in the intended form.

This project interfaces with tools including:

- MAFFT;
- MARS;
- SeqKit;
- CD-HIT/CD-HIT-EST;
- trimAl;
- VarVAMP;
- MFEprimer;
- NCBI BLAST+.

Changes involving these tools should respect their respective licenses and
document any new installation requirements.

## License

By contributing, you agree that your contribution may be distributed under
the repository license:

**GNU General Public License v3.0 or later (GPL-3.0-or-later).**