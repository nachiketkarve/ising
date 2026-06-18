# One-Dimensional Mixed-Field Ising Model

This repository contains C++ and Python code for analysis of a mixed-field Ising model. We consider the Hamiltonian

$$H = J \sum_{i=1}^N \sigma^z_i \sigma^z_{i+1} + h_x \sum_{i=1}^N \sigma^x_i + h_z \sum_{i=1}^N \sigma^z_i,$$

with periodic boundary conditions. In the absence of the longitudinal field ($h_z=0$), this model reduces to the transverse-field Ising chain, which is integrable. A longitudinal field generically breaks integrability, and the resulting mixed-field Ising model exhibits chaotic behavior.

## Speed-Fisher Information

The chaotic/non-chaotic nature of the system is captured by the so called "speed-Fisher information". This quantity can be estimated by measuring thermodynamic drag coefficient under a cyclic perturbation of the system. Concretely, in a Gibbs state, the speed-Fisher information is given by

$$\mathscr{I}_{\bar{v}}(\mu) = \frac{2}{k_B T \bar{v}^2}\langle\Delta E\rangle,$$

where $\bar{v}$ is the average speed of the protocol and $\langle\Delta E\rangle$ is the average heating induced by the perturbation. Specifically, here we consider the perturbation $\lambda(t) V$, where

$$V = \sum_{i=1}^N \sigma^x_i$$

and

$$\lambda(t) = \begin{cases}
\frac{\pi\bar{v}}{2\mu}\sin^2 \mu t, & 0 \leq t \leq \frac{2\pi}{\mu} \\ 
0, & \text{otherwise.}
\end{cases}$$

These measurements are performed by the `cyclicPerturbation` applications.

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
