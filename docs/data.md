# Data Conventions

Generated CSV files are written under `data/` by default, controlled by the `saveFolder` value in `_params.json`.

Recommended conventions:

- Treat `data/` as generated output unless a file is deliberately promoted to a small reference dataset.
- Store large production outputs outside git, or archive them with a DOI-backed data repository.
- Record the parameter file used for any published run alongside the resulting data.
- Use descriptive subdirectories such as `data/raw/`, `data/processed/`, or domain-specific labels when running batches.
