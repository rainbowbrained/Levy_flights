# Revisiting the Levy Foraging Hypothesis in 2D with Reproducible Simulation

Master's thesis by **Anita Toleutaeva**, Skoltech, 2026.

Supervisor: **Vladimir Palyulin** (Associate Professor, PhD).

Department: Advanced Computational Science.

## Abstract

Random search arises whenever an agent must locate sparse targets with little or no sensory guidance. The Levy foraging hypothesis proposes that scale-free motion with power-law relocation lengths can be near-optimal in two-dimensional search. This work develops a transparent, reproducible 2D simulation framework for step-by-step random search (Levy walk at constant unit speed) on an infinite plane with Poisson-distributed targets. Within a unified framework, we compare four post-detection strategies -- no-interrupt, teleport, and Viswanathan restart in both destructive and non-destructive modes -- sweeping the Levy exponent mu in [1, 3]. Extensive Monte Carlo simulations across sparsity levels lambda/r_v in [10, 2000] yield an empirical regime map of the optimal exponent mu*.

**Note on terminology:** The simulation implements a *Levy walk* (constant-speed motion), not a Levy flight (instantaneous jumps). Following standard usage in the foraging literature, the term "Levy flight" appears throughout the thesis as shorthand. See Section 2.1 of the thesis for a detailed discussion.

## Repository structure

```
simulation/              C99 source code for Monte Carlo simulations
  levy_viswanathan.c       Viswanathan destructive & non-destructive strategies
  levy_no_interrupt.c      No-interrupt strategy
  levy_teleport.c          Teleport strategy
  levy_batch_infinite_v2.c Combined batch runner (all strategies)
  levy.c                   Standalone single-strategy runner

scripts/                 Python plotting and bash experiment-runner scripts
  plot_full_study.py       Generates all thesis figures (14 plots)
  run_mega_experiments.sh  Phase diagram sweep (8 sparsity x 4 strategies)
  run_mega_sensitivity.sh  Truncation parameter sweeps
  run_mega_convergence.sh  Convergence studies (flights and repetitions)
  run_all_mega.sh          Master script running all experiments

data/                    Raw experiment results (CSV)
  phase_diagram_mega/      Phase diagram sweep: 32 files
  sensitivity_mega/        Truncation sensitivity: 13 files
  convergence_mega/        Convergence study: 28 files
  mfpt_large/              Mean first-passage time: 8 files

figures/                 Generated figures (PNG, PDF, EPS)
  full_study_mega/         Final thesis figures

thesis/                  LaTeX thesis source (lualatex + bibtex)
  main.tex                 Master document
  skthesis.cls             Skoltech thesis class
  Bibliography.bib         References
  chapters/                Individual chapter .tex files
  images/                  Static images and generated figure copies
  fonts/                   Times New Roman fonts for Skoltech template

slides/                  Beamer presentation source (pdflatex + bibtex)
  sample.tex               Presentation source
  reference.bib            Slide references
  figs/                    Slide figures
```

## Software environment

All results were produced in the following environment:

| Component   | Version                                   |
|-------------|-------------------------------------------|
| Compiler    | gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0 |
| C standard  | C99 (`-std=c99 -O2`)                      |
| OS          | Ubuntu 22.04 LTS, Linux 6.9.3, x86_64     |
| Python      | 3.10.12                                    |
| Matplotlib  | 3.10.1                                     |
| NumPy       | 1.26.4                                     |
| Pandas      | 2.2.3                                      |

## Random seeds and determinism

The simulation uses a deterministic seeding scheme:

- **Default master seed:** 42 (overridable via `--seed <int>`)
- **Per-run seed:** derived as `master_seed XOR (mu_index << 32) XOR run_index`
- **World seed:** derived from the per-run seed via a bijective hash (`mix64`)
- **Chunk seeds:** derived deterministically from the world seed and chunk coordinates, so the Poisson target field is independent of traversal order
- **PRNG:** Mersenne Twister (mt19937), 64-bit

Given the same master seed and parameters, the simulation produces bit-identical output across runs.

## Reproducing the results from scratch

### 1. Compile the simulation

```bash
cd simulation
gcc -std=c99 -O2 -o levy_viswanathan levy_viswanathan.c -lm
gcc -std=c99 -O2 -o levy_no_interrupt levy_no_interrupt.c -lm
gcc -std=c99 -O2 -o levy_teleport levy_teleport.c -lm
```

### 2. Run the experiments

```bash
cd ../scripts

# Phase diagrams: 8 sparsity levels x 4 strategies (~ 30 CPU-hours)
bash run_mega_experiments.sh

# Truncation sensitivity sweeps
bash run_mega_sensitivity.sh

# Convergence studies (flights and repetitions)
bash run_mega_convergence.sh
```

Or run everything at once:
```bash
bash run_all_mega.sh
```

### 3. Generate all figures

```bash
python3 plot_full_study.py --mega
```

This produces all 14 thesis figures in `figures/full_study_mega/`.

### 4. Build the thesis

```bash
cd ../thesis
lualatex main && bibtex main && lualatex main && lualatex main
```

### 5. Build the slides

```bash
cd ../slides
pdflatex sample && bibtex sample && pdflatex sample && pdflatex sample
```

## Simulation parameters

| Parameter           | Symbol             | Default value       |
|---------------------|--------------------|---------------------|
| Vision radius       | r_v                | 0.5                 |
| Min flight length   | l_min              | r_v = 0.5           |
| Max flight length   | l_max              | lambda              |
| Scan segment length | dx                 | 50                  |
| Chunk side          | chunk_size         | 2000                |
| Teleport distance   | l_c                | 5 r_v = 2.5         |
| Levy exponent range | mu                 | 1.0 -- 3.0          |
| Levy exponent step  | delta_mu           | 0.1                 |
| Master seed         | seed               | 42                  |

Density-dependent parameters (flights per run, number of runs) vary with lambda/r_v; see Table 5.1 in the thesis for exact values.

## Key results

- **Sparse targets** (lambda/r_v >= 500): mu* ~ 2 for all four strategies, consistent with the Levy foraging hypothesis
- **Dense targets** (lambda/r_v <= 20): mu* -> 1 (ballistic) for three of four strategies
- **Non-destructive trapping:** revisit ratios ~ 10^3 at high density create spurious mu* shifts above 2
- **MFPT vs throughput:** MFPT-optimal mu ~ 1.5 < throughput-optimal mu* ~ 1.7-2.0
- **Truncation robustness:** standard choice l_min = r_v, l_max = lambda is robust for l_max/lambda in [0.5, 10]

## CSV data format

Each CSV file contains one row per mu value with columns:

```
rho, rv, mu, lmin, lmax, dx, chunk_size, remove_prob, flights, seed, reps,
unique_captures_global, visits_global, total_distance_global,
normalized_rate_unique_global, normalized_rate_visits_global,
mfpt_mean, mfpt_std, first_capture_time_mean, first_capture_time_std,
normalized_rate_unique_mean, normalized_rate_unique_std,
normalized_rate_visits_mean, normalized_rate_visits_std
```

## License

This repository accompanies a master's thesis at Skoltech. The code is provided for reproducibility and academic use.
