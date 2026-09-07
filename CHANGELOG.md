# Changelog

All notable changes to this project will be documented in this file.

## [1.3.0] - 2026-09-07

### Added

- Added a reproducible Docker environment for the complete workflow.
- Added `Dockerfile`.
- Added `compose.yaml` with dedicated `design` and `validation` services.
- Added `.dockerignore`.
- Added automatic installation of MARS during Docker image construction.
- Added automatic installation of MFEprimer 3.1.0 during Docker image construction.
- Added Docker execution for Workflow 01:

  ```bash
  docker compose run --rm design