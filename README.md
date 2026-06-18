# One-Dimensional Mixed-Field Ising Model

This repository contains C++ and Python code for analysis of a mixed-field Ising model. We consider the Hamiltonian

$$H = J \sum_{i=1}^N \sigma^z_i \sigma^z_{i+1} + h_x \sum_{i=1}^N \sigma^x_i + h_z \sum_{i=1}^N \sigma^z_i,$$

with periodic boundary conditions. In the absence of the longitudinal field ($h_z=0$), this model reduces to the transverse-field Ising chain, which is integrable. A longitudinal field generically breaks integrability, and the resulting mixed-field Ising model exhibits chaotic behavior.

## Repository Layout

- `src/`: C++ simulation drivers and shared Ising routines.
- `include/`: C++ headers and Eigen type aliases.
- `python/`: Python analysis and exact-diagonalization scripts.
- `configs/`: Example parameter files for reproducible runs.
- `data/`: Generated CSV outputs. Large/generated data should stay out of git unless a small reference dataset is intentionally curated.
- `docs/`: Notes about data and workflow conventions.
- `build/`: Generated object files and executables from `make`.

## Dependencies

C++:

- A C++17 compiler with OpenMP support.
- Eigen headers at `../Libraries/eigen` relative to this repo.
- nlohmann/json headers at `../Libraries/json/include` relative to this repo.

Python:

```bash
pip install -r requirements.txt
```

## Build

```bash
make
```

Executables are written to `build/bin/`, with object files in `build/obj/`.

Clean generated build files with:

```bash
make clean
```

## Run

The C++ drivers read parameters from `_params.json` in the repository root. Start from `configs/params.example.json` when creating a new run configuration.

Examples:

```bash
./build/bin/cyclicPerturbation.exe 10 0.1
./build/bin/cyclicPerturbationFast.exe 10 0.1
./build/bin/spectralFunctionKrylov.exe 100
```

Python analysis scripts live in `python/` and assume generated CSV files are available under `data/`.

## Citation

If you use this software, cite it using the metadata in `CITATION.cff`.
