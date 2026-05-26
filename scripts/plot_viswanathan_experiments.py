#!/usr/bin/env python3
"""
Visualize quick Viswanathan experiments:
  1) Destructive vs non-destructive (remove_prob=1 vs 0)
  2) Density sweep: dense (λ/rv=20), baseline (λ/rv=100), sparse (λ/rv=500)
  3) Fine μ-step around the optimum
  4) Summary: all experiments on one canvas (unique-capture, global pooled)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

outdir = Path('figures/viswanathan_experiments')
outdir.mkdir(parents=True, exist_ok=True)


def load(csv):
    return pd.read_csv(csv)


# ── Load all datasets ──────────────────────────────────────────────
d_destr  = load('results_visw_destructive.csv')     # remove_prob=1, ρ=0.02
d_base   = load('results_viswanathan_100reps.csv')   # remove_prob=0, ρ=0.02, 100 reps
d_sparse = load('results_visw_sparse500.csv')        # ρ=0.002, λ/rv=500
d_dense  = load('results_visw_dense.csv')            # ρ=0.1,   λ/rv=20
d_fine   = load('results_visw_fine_mu.csv')           # ρ=0.02,  μ step=0.05

# Also load the two other strategies for comparison
d_noint  = load('results_no_interrupt_large.csv')
d_telep  = load('results_teleport_large.csv')


def get_lambda(df):
    return df['lambda'].iloc[0] if 'lambda' in df.columns else 1.0 / (2 * df['rho'].iloc[0] * df['rv'].iloc[0])


# ══════════════════════════════════════════════════════════════════
# Plot 1: Destructive vs non-destructive
# ══════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# Left: unique captures (global pooled)
ax = axes[0]
for df, label, color, ls in [
    (d_destr,  'Destructive (remove=1)',     '#d62728', '-'),
    (d_base,   'Non-destructive (remove=0)', '#2ca02c', '-'),
    (d_noint,  'No interrupt (remove=0)',     '#1f77b4', '--'),
]:
    mu = df['mu'].values
    y  = df['normalized_rate_unique_global'].values
    v  = np.isfinite(y)
    ax.plot(mu[v], y[v], 'o-', color=color, label=label,
            markersize=4, linewidth=2, linestyle=ls, alpha=0.85)
ax.set_xlabel('μ', fontsize=12)
ax.set_ylabel('λη (unique, pooled)', fontsize=12)
ax.set_title('Unique captures', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Right: total visits (global pooled)
ax = axes[1]
for df, label, color, ls in [
    (d_destr,  'Destructive',     '#d62728', '-'),
    (d_base,   'Non-destructive', '#2ca02c', '-'),
    (d_noint,  'No interrupt',    '#1f77b4', '--'),
]:
    mu = df['mu'].values
    y  = df['normalized_rate_visits_global'].values
    v  = np.isfinite(y)
    ax.plot(mu[v], y[v], 'o-', color=color, label=label,
            markersize=4, linewidth=2, linestyle=ls, alpha=0.85)
ax.set_xlabel('μ', fontsize=12)
ax.set_ylabel('λη (visits, pooled)', fontsize=12)
ax.set_title('Total visits (including revisits)', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.suptitle('Viswanathan: Destructive vs Non-destructive Search', fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
p = outdir / 'destructive_vs_nondestructive.png'
plt.savefig(p, dpi=150, bbox_inches='tight'); plt.close()
print(f'  {p}')


# ══════════════════════════════════════════════════════════════════
# Plot 2: Density sweep (unique captures, global)
# ══════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

ax = axes[0]
for df, label, color in [
    (d_dense,  f'Dense  ρ=0.1  (λ/rv={get_lambda(d_dense)/d_dense["rv"].iloc[0]:.0f})', '#d62728'),
    (d_base,   f'Base   ρ=0.02 (λ/rv={get_lambda(d_base)/d_base["rv"].iloc[0]:.0f})',   '#2ca02c'),
    (d_sparse, f'Sparse ρ=0.002 (λ/rv={get_lambda(d_sparse)/d_sparse["rv"].iloc[0]:.0f})', '#1f77b4'),
]:
    mu = df['mu'].values
    y  = df['normalized_rate_unique_global'].values
    v  = np.isfinite(y)
    ax.plot(mu[v], y[v], 'o-', color=color, label=label,
            markersize=4, linewidth=2, alpha=0.85)
ax.set_xlabel('μ', fontsize=12)
ax.set_ylabel('λη (unique, pooled)', fontsize=12)
ax.set_title('Unique captures vs density', fontsize=13, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

ax = axes[1]
for df, label, color in [
    (d_dense,  'Dense λ/rv=20',   '#d62728'),
    (d_base,   'Base λ/rv=100',   '#2ca02c'),
    (d_sparse, 'Sparse λ/rv=500', '#1f77b4'),
]:
    mu = df['mu'].values
    y  = df['normalized_rate_visits_global'].values
    v  = np.isfinite(y)
    ax.plot(mu[v], y[v], 'o-', color=color, label=label,
            markersize=4, linewidth=2, alpha=0.85)
ax.set_xlabel('μ', fontsize=12)
ax.set_ylabel('λη (visits, pooled)', fontsize=12)
ax.set_title('Total visits vs density', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.suptitle('Viswanathan Strategy: Effect of Target Density', fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
p = outdir / 'density_sweep.png'
plt.savefig(p, dpi=150, bbox_inches='tight'); plt.close()
print(f'  {p}')


# ══════════════════════════════════════════════════════════════════
# Plot 3: Fine μ-step (unique global) + baseline coarse + destructive
# ══════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 6))

# Baseline (non-destructive, coarse μ)
mu = d_base['mu'].values
y  = d_base['normalized_rate_unique_global'].values
v  = np.isfinite(y)
ax.plot(mu[v], y[v], 's-', color='#2ca02c', label='Non-destr. (Δμ=0.1, 100 reps)',
        markersize=5, linewidth=2, alpha=0.7)

# Fine step (non-destructive)
mu = d_fine['mu'].values
y  = d_fine['normalized_rate_unique_global'].values
v  = np.isfinite(y)
ax.plot(mu[v], y[v], 'o-', color='#9467bd', label='Non-destr. (Δμ=0.05, 80 reps)',
        markersize=5, linewidth=2.5, alpha=0.9)

# Destructive
mu = d_destr['mu'].values
y  = d_destr['normalized_rate_unique_global'].values
v  = np.isfinite(y)
ax.plot(mu[v], y[v], '^-', color='#d62728', label='Destructive (Δμ=0.1, 50 reps)',
        markersize=5, linewidth=2, alpha=0.8)

# No interrupt for reference
mu = d_noint['mu'].values
y  = d_noint['normalized_rate_unique_global'].values
v  = np.isfinite(y)
ax.plot(mu[v], y[v], 'x--', color='#1f77b4', label='No interrupt (reference)',
        markersize=5, linewidth=1.5, alpha=0.6)

ax.axvline(x=2.0, color='gray', linestyle=':', alpha=0.5, label='μ = 2')
ax.set_xlabel('μ (Lévy exponent)', fontsize=13)
ax.set_ylabel('λη (unique captures, pooled)', fontsize=13)
ax.set_title('Fine μ Resolution Around Optimum', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
p = outdir / 'fine_mu_step.png'
plt.savefig(p, dpi=150); plt.close()
print(f'  {p}')


# ══════════════════════════════════════════════════════════════════
# Plot 4: Summary — all experiments, unique-capture efficiency
# ══════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(12, 7))

configs = [
    (d_noint,  'No interrupt (ρ=0.02)',            '#1f77b4', '--', 'D', 5),
    (d_telep,  'Teleport lc=5 (ρ=0.02)',           '#ff7f0e', '--', 'D', 5),
    (d_destr,  'Viswanathan destr. (ρ=0.02)',      '#d62728', '-',  '^', 5),
    (d_base,   'Viswanathan non-destr. (ρ=0.02)',  '#2ca02c', '-',  'o', 5),
    (d_sparse, 'Viswanathan non-destr. (ρ=0.002)', '#1f77b4', '-',  's', 4),
    (d_dense,  'Viswanathan non-destr. (ρ=0.1)',   '#e377c2', '-',  'v', 4),
]

for df, label, color, ls, marker, ms in configs:
    mu = df['mu'].values
    y  = df['normalized_rate_unique_global'].values
    v  = np.isfinite(y)
    ax.plot(mu[v], y[v], marker=marker, linestyle=ls, color=color, label=label,
            markersize=ms, linewidth=2, alpha=0.85)

ax.axvline(x=2.0, color='gray', linestyle=':', alpha=0.5, linewidth=1)
ax.set_xlabel('μ (Lévy exponent)', fontsize=13)
ax.set_ylabel('λη (unique captures, pooled)', fontsize=13)
ax.set_title('Viswanathan Experiments — Summary', fontsize=14, fontweight='bold')
ax.legend(fontsize=9, loc='best')
ax.grid(True, alpha=0.3)
plt.tight_layout()
p = outdir / 'summary_all_experiments.png'
plt.savefig(p, dpi=150); plt.close()
print(f'  {p}')


# ══════════════════════════════════════════════════════════════════
# Print numeric summary
# ══════════════════════════════════════════════════════════════════
print('\n' + '='*75)
print('SUMMARY: Optimal μ* (unique captures, global pooled metric)')
print('='*75)

all_data = [
    ('No interrupt (ρ=0.02)',             d_noint),
    ('Teleport lc=5 (ρ=0.02)',           d_telep),
    ('Visw. destructive (ρ=0.02)',       d_destr),
    ('Visw. non-destr. (ρ=0.02, 100r)', d_base),
    ('Visw. non-destr. fine μ (80r)',    d_fine),
    ('Visw. non-destr. sparse (ρ=0.002)',d_sparse),
    ('Visw. non-destr. dense (ρ=0.1)',   d_dense),
]

for label, df in all_data:
    col = 'normalized_rate_unique_global'
    v = np.isfinite(df[col])
    if v.sum():
        idx = df.loc[v, col].idxmax()
        lam = get_lambda(df)
        rv  = df['rv'].iloc[0]
        print(f'  {label:42s}  μ* = {df.loc[idx,"mu"]:.2f}   '
              f'λη = {df.loc[idx,col]:.4f}   λ/rv = {lam/rv:.0f}')

print('='*75)
