# ns_vortex_degradation_3d.py

## 3D Decaying Homogeneous Isotropic Turbulence with Scale-Adaptive Vorticity Graceful Degradation

---

## Table of Contents

1. [Overview](#overview)
2. [Physics & Governing Equations](#physics--governing-equations)
3. [Numerical Methods](#numerical-methods)
4. [Installation & Dependencies](#installation--dependencies)
5. [Usage](#usage)
6. [Configuration & Presets](#configuration--presets)
7. [Output Files](#output-files)
8. [Key Mechanisms & Constants](#key-mechanisms--constants)
9. [Code Architecture](#code-architecture)
10. [Current Limitations](#current-limitations)
11. [Future Enhancements](#future-enhancements)

---

## Overview

`ns_vortex_degradation_3d.py` is a fully self-contained Python solver for **3D decaying homogeneous isotropic turbulence (HIT)** with a novel **scale-adaptive vorticity graceful degradation** model. It is the complete 3D elevation of an original 1D solver, preserving the core philosophy while extending all physics, diagnostics, and output capabilities to three dimensions.

### Core Philosophy

The central innovation is that the sink term **−Λ(ℓ[ω])·ω** is not a fixed regulariser — it is the **primary predictive output**: the spatially resolved degradation rate field `D(x,y,z,t) = Λ·|ω|`. At every grid point, the local vortex geometry (via the instantaneous dominant scale `ℓ`) selects which physical mechanism governs energy dissipation. This yields a physics-informed, interpretable diagnostic for regularity risk in turbulent flows.

---

## Physics & Governing Equations

### Vorticity Transport (3D)

```
Dω/Dt = g(|ω|, ℓ) · (ω·∇)u  +  ν∇²ω  −  Λ(ℓ[ω]) · ω
```

| Term | Description |
|---|---|
| `g(|ω|, ℓ) · (ω·∇)u` | Scale-limited vortex stretching |
| `ν∇²ω` | Viscous diffusion (exact spectral) |
| `−Λ(ℓ[ω]) · ω` | Scale-adaptive degradation (primary output) |

### Local Quantities Computed On-the-Fly

| Quantity | Formula / Method |
|---|---|
| `ℓ(x,y,z,t)` | 3D multi-level Gaussian scale decomposition |
| `η` | Kolmogorov scale: `(ν³/ε)^(1/4)`, where `ε = ν|ω|²` |
| `Re_local` | `|u|·ℓ/ν` |
| Helicity density | `h = u·ω / (|u||ω|)` (normalised, real) |
| Vortex stretching | Full `(ω·∇)u` via central differences |
| `D(x,y,z,t)` | `Λ · |ω|` — primary degradation rate field |

### Degradation Mechanisms

Six physical mechanisms compete at each grid point, selected by local scale `ℓ` and flow conditions:

| Mechanism | Active When | Lambda Formula |
|---|---|---|
| Viscous floor | `ℓ < η` | `2ν / ℓ²` |
| Vortex reconnection | `η ≤ ℓ < 100η` | `α_r · (ν/ℓ²) · √Re_local` |
| Inertial cascade | Default (inertial range) | `α_c · (ε/ℓ²)^(1/3)` |
| Helicity brake | Applied globally (modifier) | Reduces Λ when `h > 0.7` |
| Coriolis quasi-2D | `Ro < 1`, `ℓ > 0.1` | `α_Co / (ℓ·Ro)` |
| Stratification shedding | `Fr < 1`, `ℓ > 0.05` | `α_St / (ℓ·Fr²)` |

---

## Numerical Methods

| Aspect | Method |
|---|---|
| Velocity reconstruction | Exact Biot-Savart projection (FFT, spectral) |
| Diffusion | Exact spectral: `−k²ω̂` |
| Nonlinear terms | Real-space, fully vectorised (NumPy) |
| Time integration | 4th-order Runge-Kutta (RK4) |
| Domain | Periodic box `[0, 2π]³` |
| Default grid | 64³ (~minutes on a laptop) |
| Scale extraction | 3D multi-level Gaussian decomposition (wavelet-inspired, 6 levels) |
| Hölder exponent | Local log-log regression over multi-scale increments (subsampled) |

---

## Installation & Dependencies

The solver is **100% self-contained** and requires only standard scientific Python libraries:

```bash
pip install numpy scipy matplotlib
```

| Library | Version | Purpose |
|---|---|---|
| `numpy` | ≥ 1.20 | Array operations, FFT, linear algebra |
| `scipy` | ≥ 1.7 | `gaussian_filter` for scale extraction |
| `matplotlib` | ≥ 3.4 | Visualisation and plot output |

No GPU, MPI, or compiled extensions are required.

---

## Usage

```bash
# Default: Low-Re decaying turbulence (Re_λ ≈ 50), 64³ grid
python ns_vortex_degradation_3d.py

# Higher Re (Re_λ ≈ 120)
python ns_vortex_degradation_3d.py --preset high

# Run the minimal test suite (16³ grid, 50 steps)
python ns_vortex_degradation_3d.py --test

# Custom parameters
python ns_vortex_degradation_3d.py --N 128 --nu 0.005 --n_steps 5000 --dt 0.003

# Custom output directory
python ns_vortex_degradation_3d.py --outdir my_results/
```

### Command-Line Arguments

| Argument | Default | Description |
|---|---|---|
| `--preset` | `low` | Simulation preset: `low` (Re_λ≈50) or `high` (Re_λ≈120) |
| `--test` | flag | Run minimal test suite and exit |
| `--N` | `64` | Grid resolution per dimension (N³ total) |
| `--nu` | preset | Kinematic viscosity |
| `--n_steps` | preset | Number of time steps |
| `--dt` | preset | Time step size |
| `--outdir` | `output_3d` | Output directory for plots and data |

---

## Configuration & Presets

Two built-in presets are provided and can be overridden via CLI flags or the Python API:

### Low Reynolds (`cfg_low_reynolds`)
```
Re_λ ≈ 50 | ν = 0.01 | dt = 0.005 | n_steps = 3000 | ω_crit = 15.0
```

### High Reynolds (`cfg_high_reynolds`)
```
Re_λ ≈ 120 | ν = 0.0025 | dt = 0.002 | n_steps = 6000 | ω_crit = 40.0
```

### Shared Base Config

| Parameter | Default | Description |
|---|---|---|
| `N` | 64 | Grid points per dimension |
| `L` | 2π | Domain size |
| `n_save` | 150 | Save diagnostics every N steps |
| `Ro` | None | Rossby number (enable Coriolis) |
| `Fr` | None | Froude number (enable stratification) |
| `use_helicity` | True | Enable helicity brake modifier |
| `degrade_sharp` | 2.0 | Sharpness of stretching limiter |

---

## Output Files

All outputs are written to `./output_3d/` (or custom `--outdir`):

| File | Description |
|---|---|
| `time_evolution.png` | Global scalar diagnostics vs time: enstrophy, peak vorticity, BKM integral, risk fraction, mechanism fractions, peak degradation rate |
| `final_slices.png` | Mid-plane (z = N/2) slice maps of `|ω|`, `D`, `ℓ`, mechanism, risk, and Hölder exponent at final time |
| `degradation_map.png` | Space-time evolution of D, risk, etc. on mid-plane |
| `mechanism_fractions.png` | Volume fraction of each mechanism vs time |
| `enstrophy_evolution.png` | `Ω(t)`, `𝒫(t)`, and BKM integral on semi-log axes |
| `final_fields.npy` | Full 3D NumPy arrays: `omega`, `D_mid`, `mechanism` (for post-processing) |
| `final_summary.txt` | Summary statistics |

### BKM Criterion
The Beale-Kato-Majda (BKM) integral `∫₀ᵀ ‖ω‖_∞ dt` is tracked throughout the simulation. If it exceeds the threshold `BKM_THRESHOLD = 50.0`, this signals a potential regularity breakdown.

---

## Key Mechanisms & Constants

```python
ALPHA_RECONNECT   = 0.18   # Reconnection rate coefficient
ALPHA_CASCADE     = 0.40   # Inertial cascade coefficient
ALPHA_CORIOLIS    = 0.05   # Coriolis mechanism coefficient
ALPHA_STRAT       = 0.12   # Stratification shedding coefficient
HELICITY_EXPONENT = 0.85   # Helicity brake exponent
BKM_THRESHOLD     = 50.0   # Regularity risk threshold
```

---

## Code Architecture

```
ns_vortex_degradation_3d.py
│
├── Global constants & mechanism labels
├── Simulation presets (cfg_low_reynolds, cfg_high_reynolds)
│
├── Diagnostics
│   ├── compute_enstrophy()          Ω = ½∫|ω|² dV
│   ├── compute_palinstrophy()       𝒫 = ½∫|∇ω|² dV
│   └── compute_max_vorticity()      ‖ω‖_∞
│
├── Scale Analysis
│   ├── extract_dominant_scale_3d()  Multi-level Gaussian decomposition
│   └── compute_holder_exponent_3d() Local Hölder exponent via log-log fit
│
├── Spectral Helpers
│   ├── get_wavenumbers()            3D FFT wavenumber arrays
│   └── velocity_from_vorticity()    Exact Biot-Savart in Fourier space
│
├── Degradation Physics
│   ├── _lambda_viscous()
│   ├── _lambda_reconnection()
│   ├── _lambda_cascade()
│   ├── _lambda_coriolis()
│   ├── _lambda_stratification()
│   ├── stretching_limiter()         Scale-dependent g(|ω|, ℓ)
│   └── compute_degradation_field_3d() Full mechanism selection + risk map
│
├── Nonlinear Terms
│   └── compute_vorticity_stretching()  Full (ω·∇)u via central differences
│
├── Time Integration
│   └── run_simulation()             RK4 loop + diagnostics collection
│
├── Post-processing
│   └── save_all_plots()             All output figures
│
└── CLI
    ├── parse_args()
    ├── run_tests()
    └── main()
```

---

## Current Limitations

### Numerical

- **No de-aliasing:** The pseudo-spectral approach applies FFTs without a 2/3 dealiasing rule (e.g., truncating the upper 1/3 of wavenumbers). At higher Reynolds numbers or coarser grids, aliasing errors from nonlinear term computation will accumulate and can corrupt the solution.

- **Fixed time step (no CFL control):** The time step `dt` is set manually. There is no adaptive stepping or CFL stability check. Inappropriate `dt` choices can cause numerical instability, especially at higher Re or finer grids, and the solver will not warn or abort automatically.

- **Real-space gradient computation:** Vortex stretching `(ω·∇)u` is computed via `numpy.gradient` (second-order central differences) rather than spectrally. This introduces `O(dx²)` truncation error in the nonlinear terms, inconsistent with the spectral accuracy of the diffusion and pressure terms.

- **Hölder exponent is a 64³ bottleneck:** The `compute_holder_exponent_3d()` function uses a triple nested Python `for` loop over grid points (even when subsampled every 2 points), making it prohibitively expensive for large grids. At 64³ it is called only at final time, but it will not scale to 128³ or beyond.

- **Mid-plane only for diagnostics:** The risk map, degradation rate, and mechanism fields saved in history are computed only on the mid-plane slice (`z = N//2`), not the full 3D volume. The global risk fraction is therefore a 2D proxy, not a true volumetric statistic.

### Physics

- **No rotation or stratification dynamics beyond Λ:** Rossby (`Ro`) and Froude (`Fr`) number parameters only affect the degradation rate `Λ`. There is no actual centrifugal force or buoyancy body force added to the momentum equation — the background rotation or density stratification is not physically coupled into the flow evolution.

- **Helicity mechanism labelling is incomplete:** The `MECH_HELICAL` mechanism label is assigned in the mechanism array, but helicity is never selected as the dominant mechanism in the `np.select` call — it only overrides the cascade label post-hoc. Helicity therefore has no dedicated `_lambda_helical()` computation; its effect enters only as a scalar multiplier on whichever mechanism was selected.

- **Constant kinematic viscosity:** The viscosity `ν` is spatially and temporally uniform. There is no support for temperature-dependent viscosity, turbulent eddy viscosity (LES sub-grid model), or hyperviscosity.

- **Fixed random seed:** The initial condition is always generated with `numpy.random.seed(42)`, making every run with the same parameters deterministic. There is no facility to run ensembles or vary initial conditions without modifying the source code.

- **Empirical degradation constants are uncalibrated:** The alpha coefficients (`ALPHA_RECONNECT`, `ALPHA_CASCADE`, etc.) are hard-coded constants that have not been validated against DNS benchmarks or experimental data. Results are qualitatively illustrative, not quantitatively predictive.

### Software & Usability

- **No checkpointing or restart capability:** If a long simulation is interrupted, it cannot be resumed. The full run must be restarted from scratch.

- **No 3D volume visualisation:** Outputs are limited to 2D mid-plane slices. There is no support for isosurface rendering, volume rendering, or interactive 3D exploration of the vorticity or degradation fields.

- **Minimal test coverage:** The test suite (`run_tests()`) runs a single 16³, 50-step simulation and checks only that `bkm_final > 0` and `len(history['times']) > 0`. There are no unit tests for individual physics routines, no convergence tests, and no regression baselines.

- **Memory scales as O(N³):** All fields (`omega`, `u`, `ell`, `D_field`, etc.) are kept as dense NumPy arrays. At 128³ this requires ~2 GB RAM; at 256³ ~16 GB. There is no memory-efficient streaming or out-of-core processing.

- **History grows unboundedly:** Mid-plane slices of `omega`, `D`, `ell`, `mechanism`, and `risk` are appended to Python lists at every save step. For long, high-resolution runs this history can itself consume significant RAM.

---

## Future Enhancements

### Numerical Accuracy

- **Implement 2/3-rule de-aliasing** by zeroing the upper third of wavenumber modes after each FFT operation. This is the standard fix for pseudo-spectral turbulence solvers and prevents nonlinear aliasing errors from polluting inertial range statistics.

- **Spectral vortex stretching:** Replace `numpy.gradient`-based computation of `(ω·∇)u` with a fully spectral calculation (multiply in Fourier space by `iK`, then transform). This removes the `O(dx²)` finite-difference error and makes the nonlinear terms consistent with the rest of the spectral scheme.

- **Adaptive time-stepping with CFL control:** Add a CFL condition `dt ≤ C · dx / max|u|` evaluated each step, with automatic `dt` reduction and a maximum `dt_max` guard. This prevents blow-up at high Re without requiring manual tuning.

- **Implicit or IMEX time integration:** Move the stiff diffusion term `ν∇²ω` to an implicit or exponential integrator (e.g., integrating factor RK4), allowing much larger time steps at high Re without stability penalty.

### Physics Extensions

- **Physical rotation and stratification:** Add a Coriolis body force `−2Ω × u` and/or buoyancy term `−(ρ'/ρ₀)g ê_z` to the momentum equation when `Ro` or `Fr` are set. This would make the rotation/stratification effects physically consistent rather than only modifying `Λ`.

- **Dedicated helicity mechanism:** Implement a `_lambda_helical()` function with its own scaling law (e.g., based on relative helicity suppression of the energy cascade) and integrate it properly into the mechanism selection logic with its own condition.

- **Turbulent sub-grid model (LES):** Add a Smagorinsky or dynamic Smagorinsky eddy-viscosity model for the sub-grid stress tensor. This would allow physically meaningful simulations at coarser resolutions by modelling unresolved small-scale dissipation.

- **Ensemble runs and statistics:** Expose the random seed as a CLI parameter (`--seed`). Add a batch mode that runs N realisations and computes ensemble-averaged diagnostics (mean, variance of enstrophy decay, BKM statistics, etc.).

### Performance & Scalability

- **GPU acceleration with CuPy:** Replace `numpy` and `scipy` with `cupy` and `cupyx.scipy` for GPU-accelerated FFTs and Gaussian filtering. This could yield 10–100× speedup for 128³–256³ grids.

- **Parallel CPU computation with mpi4py:** Decompose the domain in slab geometry along the z-axis and distribute FFTs across MPI ranks using `mpi4py` and `mpi4py-fft`. This enables weak scaling to high-resolution grids on HPC clusters.

- **Vectorise the Hölder exponent computation:** Rewrite `compute_holder_exponent_3d()` using vectorised NumPy operations over all grid points simultaneously, avoiding the triple Python `for` loop. Alternatively, implement this as a rolling-window convolution.

- **Chunked/streaming history saving:** Write mid-plane snapshots to disk (HDF5 or NumPy `.npy`) at each save step rather than accumulating them in RAM. This makes long runs feasible and supports post-hoc analysis of large datasets.

### Diagnostics & Visualisation

- **3D volumetric outputs:** Export full 3D fields to VTK or HDF5 format for rendering in ParaView or VisIt. Add isosurface plots of `|ω|` and `D` using `matplotlib`'s 3D axes or `pyvista`.

- **Full 3D risk statistics:** Compute risk fraction, mechanism volume fractions, and mean degradation rate as true 3D volume averages rather than mid-plane proxies. This requires computing the full 3D degradation field each save step (currently only done for the mid-plane).

- **Energy spectrum `E(k)`:** Compute the spherically averaged 3D kinetic energy spectrum `E(k)` at each save step and plot it alongside the Kolmogorov `k^{-5/3}` reference line. This is a standard validation diagnostic for turbulence solvers.

- **Inertial range diagnostics:** Track the compensated spectrum `k^{5/3} E(k)`, the integral length scale `L_int = (3π/4) ∫ E(k)/k dk / ∫ E(k) dk`, and the Taylor microscale `λ = √(10νk / ε)` over time.

### Code Quality & Testing

- **Expanded test suite:** Add unit tests for each physics function (`_lambda_viscous`, `_lambda_reconnection`, `compute_vorticity_stretching`, etc.) with known analytical inputs and expected outputs. Add a convergence test verifying that the enstrophy decay rate approaches the expected `dΩ/dt = −2ν𝒫` identity.

- **Checkpointing and restarts:** Save the full simulation state (omega, BKM integral, history, config, step number) to a checkpoint `.npz` file at regular intervals and on graceful exit. Add a `--resume <checkpoint>` CLI flag.

- **Configuration via YAML/JSON:** Accept simulation parameters from a configuration file rather than only CLI flags. This makes it easier to reproduce, share, and version-control specific runs.

- **Type annotations and docstrings:** Add NumPy-style docstrings with `Parameters`, `Returns`, and `Notes` sections to all functions. Add Python type annotations for improved IDE support and static analysis.
