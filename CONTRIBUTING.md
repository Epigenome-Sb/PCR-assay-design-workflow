# Contributing

Contributions are welcome through GitHub issues and pull requests.

This repository contains two main workflows:

- `pcr_workflow/design/` — nucleotide preprocessing, alignment, conservation analysis, and VarVAMP assay design;
- `pcr_workflow/validation/` — in silico primer-pair validation and, for QPCR designs, qPCR probe validation.

Run the workflows with `python3 -m pcr_workflow design` or
`python3 -m pcr_workflow validate`. The numbered root scripts remain as
compatibility launchers. Shared runtime, path, and command helpers live in
`pcr_workflow/core/`.

The project supports VarVAMP `SINGLE`, `QPCR`, and `TILED` assay modes.

---

## Reporting a problem

Please include:

- operating system;
- execution method: native Conda or Docker;
- Python version;
- workflow/script concerned (`01` or `02`);
- complete command used;
- scientific parameters selected;
- full error message;
- relevant software versions;
- a minimal non-confidential nucleotide FASTA example when possible.

For Docker-related problems, please also include:

```bash
docker --version
docker compose version
docker compose config
