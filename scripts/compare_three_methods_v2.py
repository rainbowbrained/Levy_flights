#!/usr/bin/env python3
"""
Comparison of three Lévy search strategies using pooled (global) metrics.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    csv_files = [
        'results_no_interrupt_large.csv',
        'results_teleport_large.csv',
        'results_viswanathan_100reps.csv',
    ]
    labels = [
        'No interrupt',
        'Teleport (lc=5)',
        'Viswanathan',
    ]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    outdir = Path('figures/comparison_v2')
    outdir.mkdir(parents=True, exist_ok=True)

    dfs = [pd.read_csv(f) for f in csv_files]

    # ---- Plot 1: λη unique (global pooled) ----
    fig, ax = plt.subplots(figsize=(10, 6))
    for df, label, color in zip(dfs, labels, colors):
        mu = df['mu'].values
        rate_global = df['normalized_rate_unique_global'].values
        valid = np.isfinite(rate_global)
        ax.plot(mu[valid], rate_global[valid], 'o-', color=color, label=label,
                markersize=5, linewidth=2.5, alpha=0.9)
    ax.set_xlabel('μ (Lévy exponent)', fontsize=13)
    ax.set_ylabel('λη (unique captures, pooled)', fontsize=13)
    ax.set_title('Search Efficiency — Unique Captures (pooled over all reps)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    p = outdir / 'comparison_unique_global.png'
    plt.savefig(p, dpi=150); plt.close()
    print(f'✅ {p}')

    # ---- Plot 2: λη unique — mean ± std per run ----
    fig, ax = plt.subplots(figsize=(10, 6))
    for df, label, color in zip(dfs, labels, colors):
        mu = df['mu'].values
        m = df['normalized_rate_unique_mean'].values
        s = df['normalized_rate_unique_std'].values
        v = np.isfinite(m) & np.isfinite(s)
        ax.plot(mu[v], m[v], 'o-', color=color, label=label,
                markersize=5, linewidth=2, alpha=0.85)
        ax.fill_between(mu[v], m[v]-s[v], m[v]+s[v], color=color, alpha=0.15)
    ax.set_xlabel('μ', fontsize=13)
    ax.set_ylabel('λη (unique, per-run mean ± σ)', fontsize=13)
    ax.set_title('Search Efficiency — Unique Captures (mean ± σ)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    p = outdir / 'comparison_unique_mean_std.png'
    plt.savefig(p, dpi=150); plt.close()
    print(f'✅ {p}')

    # ---- Plot 3: λη visits (global pooled) ----
    fig, ax = plt.subplots(figsize=(10, 6))
    for df, label, color in zip(dfs, labels, colors):
        mu = df['mu'].values
        rate_global = df['normalized_rate_visits_global'].values
        valid = np.isfinite(rate_global)
        ax.plot(mu[valid], rate_global[valid], 'o-', color=color, label=label,
                markersize=5, linewidth=2.5, alpha=0.9)
    ax.set_xlabel('μ', fontsize=13)
    ax.set_ylabel('λη (all visits, pooled)', fontsize=13)
    ax.set_title('Search Efficiency — All Visits (pooled)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    p = outdir / 'comparison_visits_global.png'
    plt.savefig(p, dpi=150); plt.close()
    print(f'✅ {p}')

    # ---- Plot 4: Two-panel — No interrupt & Teleport (left) vs Viswanathan (right) ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Left panel: No interrupt + Teleport
    for df, label, color in zip(dfs[:2], labels[:2], colors[:2]):
        mu = df['mu'].values
        m = df['normalized_rate_unique_mean'].values
        s = df['normalized_rate_unique_std'].values
        v = np.isfinite(m)
        ax1.plot(mu[v], m[v], 'o-', color=color, label=label,
                 markersize=5, linewidth=2, alpha=0.85)
        ax1.fill_between(mu[v], m[v]-s[v], m[v]+s[v], color=color, alpha=0.15)
    ax1.set_xlabel('μ', fontsize=12)
    ax1.set_ylabel('λη (unique captures)', fontsize=12)
    ax1.set_title('No Interrupt & Teleport', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)

    # Right panel: Viswanathan — global metric
    df = dfs[2]
    mu = df['mu'].values
    rg = df['normalized_rate_unique_global'].values
    v = np.isfinite(rg)
    ax2.plot(mu[v], rg[v], 'o-', color=colors[2], label='Viswanathan (pooled)',
             markersize=5, linewidth=2.5)
    # also show per-run mean ± std
    m = df['normalized_rate_unique_mean'].values
    s = df['normalized_rate_unique_std'].values
    v2 = np.isfinite(m)
    ax2.fill_between(mu[v2], m[v2]-s[v2], m[v2]+s[v2],
                     color=colors[2], alpha=0.12, label='± σ (per-run)')
    ax2.set_xlabel('μ', fontsize=12)
    ax2.set_ylabel('λη (unique captures)', fontsize=12)
    ax2.set_title('Viswanathan (100 reps)', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    p = outdir / 'comparison_two_panel.png'
    plt.savefig(p, dpi=150); plt.close()
    print(f'✅ {p}')

    # ---- Summary table ----
    print('\n' + '='*65)
    print('Optimal μ* (unique captures, global pooled metric)')
    print('='*65)
    for label, df in zip(labels, dfs):
        col = 'normalized_rate_unique_global'
        v = np.isfinite(df[col])
        if v.sum():
            idx = df.loc[v, col].idxmax()
            print(f'  {label:25s}  μ* = {df.loc[idx,"mu"]:.1f}   λη = {df.loc[idx,col]:.4f}')
    print('='*65)

if __name__ == '__main__':
    main()
