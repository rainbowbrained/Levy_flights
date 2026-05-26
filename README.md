# Revisiting the Lévy Foraging Hypothesis in 2D with Reproducible Simulation

Master's thesis by Anita Toleutaeva, Skoltech, 2026.

Supervisor: Vladimir Palyulin.

## Repository structure

```
simulation/          C99 source code for Monte Carlo simulations
scripts/             Python plotting and bash experiment-runner scripts
data/                Experiment results (CSV)
  convergence_mega/  Convergence study data
  phase_diagram_mega/ Phase diagram sweep data
  sensitivity_mega/  Truncation sensitivity data
  mfpt_large/        Mean first-passage time data
figures/             Generated figures (PNG, PDF, EPS)
thesis/              LaTeX source for the thesis (lualatex + bibtex)
slides/              Beamer presentation source (pdflatex + bibtex)
```

## Building the simulation

```bash
cd simulation
gcc -O2 -o levy_batch_infinite_v2 levy_batch_infinite_v2.c -lm
```

## Building the thesis

```bash
cd thesis
lualatex main && bibtex main && lualatex main && lualatex main
```

## Building the slides

```bash
cd slides
pdflatex sample && bibtex sample && pdflatex sample && pdflatex sample
```
