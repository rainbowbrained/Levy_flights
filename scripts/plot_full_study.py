#!/usr/bin/env python3
"""
Full study visualizations for the Lévy flight search diploma thesis.

Generates:
  1. Phase diagram: μ*(λ/rv) for all strategies
  2. Efficiency scaling: λη_max(λ/rv)
  3. Trapping heatmap: visits/unique ratio
  4. MFPT analysis
  5. lmin/lmax sensitivity
  6. Convergence analysis
  7. Analytical comparison
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
from glob import glob
import re

plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'legend.fontsize': 12,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'figure.dpi': 150,
})

import sys

# Allow switching between small, large and mega datasets via CLI
if '--mega' in sys.argv:
    _suffix = '_mega'
elif '--large' in sys.argv:
    _suffix = '_large'
else:
    _suffix = ''
_save_eps = '--eps' in sys.argv
PHASE_DIR = f'phase_diagram{_suffix}'
SENS_DIR  = f'sensitivity{_suffix}'
CONV_DIR  = f'convergence{_suffix}'
# No dedicated mfpt_mega; fall back to mfpt_large for mega runs
MFPT_DIR  = 'mfpt_large' if _suffix == '_mega' else f'mfpt{_suffix}'
OUTDIR = Path(f'figures/full_study{_suffix}')
OUTDIR.mkdir(parents=True, exist_ok=True)
print(f"Data dirs: {PHASE_DIR}/, {SENS_DIR}/, {CONV_DIR}/, {MFPT_DIR}/")
print(f"Output:    {OUTDIR}/")

if _save_eps:
    _orig_savefig = plt.savefig
    def _savefig_with_pdf(fname, *args, **kwargs):
        _orig_savefig(fname, *args, **kwargs)
        if isinstance(fname, (str, Path)):
            pdf_path = Path(fname).with_suffix('.pdf')
            _orig_savefig(str(pdf_path), *args, format='pdf', **kwargs)
            print(f'  + {pdf_path}')
    plt.savefig = _savefig_with_pdf

# ═══════════════════════════════════════════════════════════════
# Helper functions
# ═══════════════════════════════════════════════════════════════

def load_csv(path):
    try:
        df = pd.read_csv(path)
        if df.empty or len(df) == 0:
            print(f"  Warning: empty file {path}, skipping")
            return pd.DataFrame()
        return df
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        print(f"  Warning: cannot parse {path}: {e}, skipping")
        return pd.DataFrame()

def get_mu_star(df, metric='normalized_rate_unique_global'):
    """Return (mu*, value_at_mu*) from pooled global metric."""
    v = np.isfinite(df[metric])
    if v.sum() == 0:
        return np.nan, np.nan
    idx = df.loc[v, metric].idxmax()
    return df.loc[idx, 'mu'], df.loc[idx, metric]

def get_mu_star_band(df, metric='normalized_rate_unique_global',
                     std_col='normalized_rate_unique_std'):
    """Return (mu*, mu_low, mu_high): 1σ uncertainty band on peak location.

    The band spans all μ where η(μ) >= η(μ*) - σ(μ*).
    """
    v = np.isfinite(df[metric])
    if v.sum() == 0:
        return np.nan, np.nan, np.nan
    mu = df.loc[v, 'mu'].values
    y = df.loc[v, metric].values
    idx_star = np.argmax(y)
    mu_star = mu[idx_star]
    if std_col not in df.columns:
        return mu_star, mu_star, mu_star
    std = df.loc[v, std_col].values
    threshold = y[idx_star] - std[idx_star]
    within = mu[y >= threshold]
    return mu_star, within.min(), within.max()

def get_lambda_rv(df):
    rho = df['rho'].iloc[0]
    rv = df['rv'].iloc[0]
    lam = 1.0 / (2 * rho * rv)
    return lam / rv

# Strategy colors (consistent across all plots)
C_NOINT = '#1f77b4'
C_TELEP = '#ff7f0e'
C_VDESTR = '#d62728'
C_VNONDESTR = '#2ca02c'

_BAND_ALPHA = 0.3  # transparency for ±σ bands (requires PDF, not EPS)

STRATEGIES = {
    'no_interrupt': ('No interrupt', C_NOINT, 'o'),
    'teleport':     ('Teleport (lc=5rv)', C_TELEP, 's'),
    'visw_destr':   ('Visw. destructive', C_VDESTR, '^'),
    'visw_nondestr':('Visw. non-destr.', C_VNONDESTR, 'v'),
}

# ═══════════════════════════════════════════════════════════════
# 1. PHASE DIAGRAM: μ*(λ/rv)
# ═══════════════════════════════════════════════════════════════

def plot_phase_diagram():
    print("=== 1. Phase diagram ===")

    pdir = Path(PHASE_DIR)
    if not pdir.exists():
        print("  No phase_diagram/ directory found, skipping.")
        return

    # Collect μ*(λ/rv) for each strategy
    results = {k: {'lrv': [], 'mu_star': [], 'mu_lo': [], 'mu_hi': [],
                    'eta_max': []}
               for k in STRATEGIES}

    # Find all λ/rv values
    files = sorted(pdir.glob('lrv_*_no_interrupt.csv'))
    lrv_vals = []
    for f in files:
        m = re.search(r'lrv_(\d+)_', f.name)
        if m:
            lrv_vals.append(int(m.group(1)))

    for lrv in sorted(lrv_vals):
        for skey in STRATEGIES:
            csv_path = pdir / f'lrv_{lrv}_{skey}.csv'
            if not csv_path.exists():
                continue
            df = load_csv(csv_path)
            mu_s, eta_m = get_mu_star(df)
            _, mu_lo, mu_hi = get_mu_star_band(df)
            results[skey]['lrv'].append(lrv)
            results[skey]['mu_star'].append(mu_s)
            results[skey]['mu_lo'].append(mu_lo)
            results[skey]['mu_hi'].append(mu_hi)
            results[skey]['eta_max'].append(eta_m)

    # --- Plot 1a: μ*(λ/rv) --- two versions: clean and with σ bands
    for show_bands in [False, True]:
        fig, ax = plt.subplots(figsize=(10, 6))
        for skey, (label, color, marker) in STRATEGIES.items():
            x = np.array(results[skey]['lrv'])
            y = np.array(results[skey]['mu_star'])
            lo = np.array(results[skey]['mu_lo'])
            hi = np.array(results[skey]['mu_hi'])
            if len(x) == 0:
                continue
            valid = np.isfinite(y)
            if show_bands:
                ax.fill_between(x[valid], lo[valid], hi[valid],
                                color=color, alpha=_BAND_ALPHA, linewidth=0)
            ax.plot(x[valid], y[valid], marker=marker, linestyle='-', color=color,
                    label=label, markersize=8, linewidth=2.5)

        ax.axhline(y=2.0, color='gray', linestyle=':', alpha=0.6, label='μ = 2 (theory)')
        ax.set_xscale('log')
        ax.set_xlabel('λ / rv (sparsity parameter)')
        ax.set_ylabel('Optimal μ*')
        ax.set_title('Phase Diagram: Optimal Lévy Exponent vs Target Sparsity',
                     fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')
        ax.set_ylim(0.8, 3.2)
        plt.tight_layout()
        suffix = '' if show_bands else '_clean'
        p = OUTDIR / f'phase_diagram_mu_star{suffix}.png'
        plt.savefig(p); plt.close()
        print(f'  {p}')

    # --- Plot 1b: λη_max(λ/rv) — scaling ---
    fig, ax = plt.subplots(figsize=(10, 6))
    for skey, (label, color, marker) in STRATEGIES.items():
        x = np.array(results[skey]['lrv'])
        y = np.array(results[skey]['eta_max'])
        if len(x) == 0:
            continue
        valid = np.isfinite(y)
        ax.plot(x[valid], y[valid], marker=marker, linestyle='-', color=color,
                label=label, markersize=7, linewidth=2)

    ax.set_xscale('log')
    ax.set_xlabel('λ / rv', fontsize=16)
    ax.set_ylabel('λη_max (unique captures, pooled)', fontsize=16)
    ax.set_title('Efficiency Scaling: Maximum λη vs Target Sparsity',
                 fontsize=18, fontweight='bold')
    ax.legend(fontsize=13)
    ax.grid(True, alpha=0.3, which='both')
    plt.tight_layout()
    p = OUTDIR / 'efficiency_scaling.png'
    plt.savefig(p); plt.close()
    print(f'  {p}')

    # --- Plot 1c: Full η(μ) curves at selected λ/rv (2+3 layout) ---
    selected_lrv = [10, 50, 100, 500, 2000]
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(20, 13))
    gs = GridSpec(2, 6, figure=fig, hspace=0.32, wspace=0.45)
    # Top row: 2 wide panels (each spans 3 columns)
    ax0 = fig.add_subplot(gs[0, 0:3])
    ax1 = fig.add_subplot(gs[0, 3:6])
    # Bottom row: 3 panels (each spans 2 columns)
    ax2 = fig.add_subplot(gs[1, 0:2])
    ax3 = fig.add_subplot(gs[1, 2:4])
    ax4 = fig.add_subplot(gs[1, 4:6])
    axes = [ax0, ax1, ax2, ax3, ax4]

    for i, lrv in enumerate(selected_lrv):
        ax = axes[i]
        for skey, (label, color, _) in STRATEGIES.items():
            csv_path = pdir / f'lrv_{lrv}_{skey}.csv'
            if not csv_path.exists():
                continue
            df = load_csv(csv_path)
            mu = df['mu'].values
            y = df['normalized_rate_unique_global'].values
            v = np.isfinite(y)
            ax.plot(mu[v], y[v], '-', color=color, linewidth=2.5, alpha=0.85,
                    label=label if i == 0 else None)
            if 'normalized_rate_unique_std' in df.columns:
                std = df['normalized_rate_unique_std'].values
                ax.fill_between(mu[v], (y - std)[v], (y + std)[v],
                                color=color, alpha=_BAND_ALPHA, linewidth=0)
        ax.axvline(x=2.0, color='gray', linestyle=':', alpha=0.5)
        ax.set_xlabel('μ', fontsize=18)
        ax.set_title(f'λ/rv = {lrv}', fontsize=20, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.tick_params(axis='both', labelsize=16)
        if i == 0 or i == 2:
            ax.set_ylabel('λη (unique, pooled)', fontsize=18)

    axes[0].legend(fontsize=14, loc='best')
    plt.suptitle('λη(μ) Curves at Different Sparsity Levels',
                 fontsize=22, fontweight='bold', y=1.01)
    p = OUTDIR / 'eta_curves_by_sparsity.png'
    plt.savefig(p, bbox_inches='tight'); plt.close()
    print(f'  {p}')

    return results


# ═══════════════════════════════════════════════════════════════
# 2. TRAPPING HEATMAP
# ═══════════════════════════════════════════════════════════════

def plot_trapping_heatmap():
    print("=== 2. Trapping heatmap ===")

    pdir = Path(PHASE_DIR)
    if not pdir.exists():
        print("  Skipping (no data).")
        return

    files = sorted(pdir.glob('lrv_*_visw_nondestr.csv'))
    if not files:
        print("  No non-destructive files found.")
        return

    lrv_vals = []
    for f in files:
        m = re.search(r'lrv_(\d+)_', f.name)
        if m:
            lrv_vals.append(int(m.group(1)))
    lrv_vals = sorted(lrv_vals)

    # Build 2D array: ratio = visits_global / unique_global
    # First get the common mu grid
    df0 = load_csv(files[0])
    mu_vals = df0['mu'].values

    ratio_matrix = np.full((len(lrv_vals), len(mu_vals)), np.nan)

    for i, lrv in enumerate(lrv_vals):
        csv_path = pdir / f'lrv_{lrv}_visw_nondestr.csv'
        if not csv_path.exists():
            continue
        df = load_csv(csv_path)
        vis_g = df['normalized_rate_visits_global'].values
        uniq_g = df['normalized_rate_unique_global'].values
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(uniq_g > 0, vis_g / uniq_g, np.nan)
        ratio_matrix[i, :len(ratio)] = ratio

    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(ratio_matrix, aspect='auto',
                   extent=[mu_vals[0], mu_vals[-1], len(lrv_vals)-0.5, -0.5],
                   cmap='YlOrRd',
                   norm=mcolors.LogNorm(vmin=1, vmax=np.nanmax(ratio_matrix)))
    ax.set_yticks(range(len(lrv_vals)))
    ax.set_yticklabels([str(v) for v in lrv_vals])
    ax.set_xlabel('μ (Lévy exponent)', fontsize=16)
    ax.set_ylabel('λ / rv', fontsize=16)
    ax.set_title('Trapping Intensity: visits / unique (Viswanathan non-destructive)',
                 fontsize=18, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('visits / unique (log scale)')
    plt.tight_layout()
    p = OUTDIR / 'trapping_heatmap.png'
    plt.savefig(p); plt.close()
    print(f'  {p}')


# ═══════════════════════════════════════════════════════════════
# 3. MFPT ANALYSIS
# ═══════════════════════════════════════════════════════════════

def _mfpt_csv(lrv, strategy='visw_destr'):
    """Find best MFPT CSV: prefer phase_diagram (most flights), fall back to mfpt_large."""
    pdir = Path(PHASE_DIR)
    mdir = Path(MFPT_DIR)
    name = f'lrv_{lrv}_{strategy}.csv'
    if pdir.exists() and (pdir / name).exists():
        return pdir / name
    if mdir.exists() and (mdir / name).exists():
        return mdir / name
    return None

def plot_mfpt():
    print("=== 3. MFPT analysis ===")

    pdir = Path(PHASE_DIR)
    mdir = Path(MFPT_DIR)
    if not pdir.exists() and not mdir.exists():
        print("  Skipping (no data).")
        return

    # MFPT vs μ at different λ/rv for destructive Viswanathan
    selected_lrv = [20, 50, 100, 500, 1000]
    fig, ax = plt.subplots(figsize=(10, 6))

    cmap = plt.cm.viridis
    norm = mcolors.LogNorm(vmin=10, vmax=5000)

    for lrv in selected_lrv:
        csv_path = _mfpt_csv(lrv)
        if csv_path is None:
            continue
        df = load_csv(csv_path)
        mu = df['mu'].values
        mfpt = df['mfpt'].values
        lam = df['lambda'].iloc[0] if 'lambda' in df.columns else 1.0/(2*df['rho'].iloc[0]*df['rv'].iloc[0])
        # Normalize MFPT by λ
        mfpt_norm = mfpt / lam
        v = np.isfinite(mfpt_norm) & (mfpt_norm > 0)
        ax.plot(mu[v], mfpt_norm[v], 'o-', color=cmap(norm(lrv)),
                label=f'λ/rv={lrv}', markersize=5, linewidth=2)
        if 'mfpt_std' in df.columns:
            mfpt_std = df['mfpt_std'].values / lam
            ax.fill_between(mu[v], (mfpt_norm - mfpt_std)[v],
                            (mfpt_norm + mfpt_std)[v],
                            color=cmap(norm(lrv)), alpha=_BAND_ALPHA, linewidth=0)

    ax.set_xlabel('μ', fontsize=16)
    ax.set_ylabel('MFPT / λ', fontsize=16)
    ax.set_title('Mean First Passage Time (destructive Viswanathan)',
                 fontsize=18, fontweight='bold')
    ax.legend(fontsize=13)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    p = OUTDIR / 'mfpt_vs_mu.png'
    plt.savefig(p); plt.close()
    print(f'  {p}')

    # MFPT comparison: all strategies at baseline λ/rv=100
    fig, ax = plt.subplots(figsize=(10, 6))
    for skey, (label, color, marker) in STRATEGIES.items():
        csv_path = _mfpt_csv(100, skey)
        if csv_path is None:
            continue
        df = load_csv(csv_path)
        mu = df['mu'].values
        mfpt = df['mfpt'].values
        lam = df['lambda'].iloc[0] if 'lambda' in df.columns else 50.0
        mfpt_norm = mfpt / lam
        v = np.isfinite(mfpt_norm) & (mfpt_norm > 0)
        ax.plot(mu[v], mfpt_norm[v], marker=marker, linestyle='-', color=color,
                label=label, markersize=6, linewidth=2)
        if 'mfpt_std' in df.columns:
            mfpt_std = df['mfpt_std'].values / lam
            ax.fill_between(mu[v], (mfpt_norm - mfpt_std)[v],
                            (mfpt_norm + mfpt_std)[v],
                            color=color, alpha=_BAND_ALPHA, linewidth=0)

    ax.set_xlabel('μ', fontsize=16)
    ax.set_ylabel('MFPT / λ', fontsize=16)
    ax.set_title('MFPT Comparison at λ/rv = 100',
                 fontsize=18, fontweight='bold')
    ax.legend(fontsize=13)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    p = OUTDIR / 'mfpt_comparison_lrv100.png'
    plt.savefig(p); plt.close()
    print(f'  {p}')


# ═══════════════════════════════════════════════════════════════
# 4. SENSITIVITY TO lmin/lmax
# ═══════════════════════════════════════════════════════════════

def plot_sensitivity():
    print("=== 4. Sensitivity to lmin/lmax ===")

    sdir = Path(SENS_DIR)
    if not sdir.exists():
        print("  Skipping (no data).")
        return

    # --- lmax sensitivity ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    lmax_files = sorted(sdir.glob('lmax_ratio_*.csv'))
    lmax_mu_star = []
    lmax_eta_max = []
    lmax_ratios = []

    cmap = plt.cm.cool
    ratios_all = []
    for f in lmax_files:
        m = re.search(r'lmax_ratio_([\d.]+?)\.csv', f.name)
        if m:
            ratios_all.append(float(m.group(1)))
    if not ratios_all:
        print("  No lmax files found.")
        return
    norm = mcolors.LogNorm(vmin=min(ratios_all), vmax=max(ratios_all))

    for f in lmax_files:
        m = re.search(r'lmax_ratio_([\d.]+?)\.csv', f.name)
        if not m:
            continue
        ratio = float(m.group(1))
        df = load_csv(f)
        if df.empty:
            continue
        mu = df['mu'].values
        y = df['normalized_rate_unique_global'].values
        v = np.isfinite(y)
        c = cmap(norm(ratio))
        ax.plot(mu[v], y[v], 'o-', color=c,
                label=f'lmax/λ={ratio}', markersize=4, linewidth=1.8)
        if 'normalized_rate_unique_std' in df.columns:
            std = df['normalized_rate_unique_std'].values
            ax.fill_between(mu[v], (y - std)[v], (y + std)[v],
                            color=c, alpha=_BAND_ALPHA, linewidth=0)

        mu_s, eta_m = get_mu_star(df)
        lmax_ratios.append(ratio)
        lmax_mu_star.append(mu_s)
        lmax_eta_max.append(eta_m)

    ax.axvline(x=2.0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('μ')
    ax.set_ylabel('λη (unique, pooled)')
    ax.set_title('Effect of lmax / λ', fontweight='bold')
    ax.legend(fontsize=12, ncol=2)
    ax.grid(True, alpha=0.3)

    # --- lmin sensitivity ---
    ax = axes[1]
    lmin_files = sorted(sdir.glob('lmin_ratio_*.csv'))
    lmin_mu_star = []
    lmin_eta_max = []
    lmin_ratios = []

    ratios_all2 = []
    for f in lmin_files:
        m = re.search(r'lmin_ratio_([\d.]+?)\.csv', f.name)
        if m:
            ratios_all2.append(float(m.group(1)))
    if ratios_all2:
        norm2 = mcolors.Normalize(vmin=min(ratios_all2), vmax=max(ratios_all2))
    else:
        norm2 = mcolors.Normalize(0, 1)

    for f in lmin_files:
        m = re.search(r'lmin_ratio_([\d.]+?)\.csv', f.name)
        if not m:
            continue
        ratio = float(m.group(1))
        df = load_csv(f)
        if df.empty:
            continue
        mu = df['mu'].values
        y = df['normalized_rate_unique_global'].values
        v = np.isfinite(y)
        c2 = cmap(norm2(ratio))
        ax.plot(mu[v], y[v], 's-', color=c2,
                label=f'lmin/rv={ratio}', markersize=4, linewidth=1.8)
        if 'normalized_rate_unique_std' in df.columns:
            std = df['normalized_rate_unique_std'].values
            ax.fill_between(mu[v], (y - std)[v], (y + std)[v],
                            color=c2, alpha=_BAND_ALPHA, linewidth=0)

        mu_s, eta_m = get_mu_star(df)
        lmin_ratios.append(ratio)
        lmin_mu_star.append(mu_s)
        lmin_eta_max.append(eta_m)

    ax.axvline(x=2.0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('μ')
    ax.set_ylabel('λη (unique, pooled)')
    ax.set_title('Effect of lmin / rv', fontweight='bold')
    ax.legend(fontsize=12, ncol=2)
    ax.grid(True, alpha=0.3)

    plt.suptitle('Sensitivity to Power-Law Truncation (Destructive Viswanathan, λ/rv=100)',
                 fontsize=16, fontweight='bold', y=1.01)
    plt.tight_layout()
    p = OUTDIR / 'sensitivity_lmin_lmax_curves.png'
    plt.savefig(p, bbox_inches='tight'); plt.close()
    print(f'  {p}')

    # --- μ* and η_max vs lmax/λ and lmin/rv ---
    # Sort by ratio to fix point sequence bug
    if lmax_ratios:
        sort_idx = np.argsort(lmax_ratios)
        lmax_ratios = [lmax_ratios[i] for i in sort_idx]
        lmax_mu_star = [lmax_mu_star[i] for i in sort_idx]
        lmax_eta_max = [lmax_eta_max[i] for i in sort_idx]
    if lmin_ratios:
        sort_idx = np.argsort(lmin_ratios)
        lmin_ratios = [lmin_ratios[i] for i in sort_idx]
        lmin_mu_star = [lmin_mu_star[i] for i in sort_idx]
        lmin_eta_max = [lmin_eta_max[i] for i in sort_idx]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    if lmax_ratios:
        ax = axes[0, 0]
        ax.plot(lmax_ratios, lmax_mu_star, 'o-', color=C_VDESTR, markersize=7, linewidth=2)
        ax.axhline(y=2.0, color='gray', linestyle=':', alpha=0.5)
        ax.set_xlabel('lmax / λ')
        ax.set_ylabel('μ*')
        ax.set_title('μ* vs lmax/λ')
        ax.set_xscale('log')
        ax.grid(True, alpha=0.3)

        ax = axes[0, 1]
        ax.plot(lmax_ratios, lmax_eta_max, 'o-', color=C_VDESTR, markersize=7, linewidth=2)
        ax.set_xlabel('lmax / λ')
        ax.set_ylabel('λη_max')
        ax.set_title('Peak efficiency vs lmax/λ')
        ax.set_xscale('log')
        ax.grid(True, alpha=0.3)

    if lmin_ratios:
        ax = axes[1, 0]
        ax.plot(lmin_ratios, lmin_mu_star, 's-', color=C_VDESTR, markersize=7, linewidth=2)
        ax.axhline(y=2.0, color='gray', linestyle=':', alpha=0.5)
        ax.set_xlabel('lmin / rv')
        ax.set_ylabel('μ*')
        ax.set_title('μ* vs lmin/rv')
        ax.set_xscale('log')
        ax.grid(True, alpha=0.3)

        ax = axes[1, 1]
        ax.plot(lmin_ratios, lmin_eta_max, 's-', color=C_VDESTR, markersize=7, linewidth=2)
        ax.set_xlabel('lmin / rv')
        ax.set_ylabel('λη_max')
        ax.set_title('Peak efficiency vs lmin/rv')
        ax.set_xscale('log')
        ax.grid(True, alpha=0.3)

    plt.suptitle('Truncation Sensitivity Summary', fontsize=18, fontweight='bold')
    plt.tight_layout()
    p = OUTDIR / 'sensitivity_summary.png'
    plt.savefig(p); plt.close()
    print(f'  {p}')


# ═══════════════════════════════════════════════════════════════
# 5. CONVERGENCE ANALYSIS
# ═══════════════════════════════════════════════════════════════

def plot_convergence():
    print("=== 5. Convergence analysis ===")

    cdir = Path(CONV_DIR)
    if not cdir.exists():
        print("  Skipping (no data).")
        return

    # --- Convergence in flights (both lrv=100 and lrv=1000) ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    _lrv_colors = {100: C_VDESTR, 1000: '#9467bd'}
    _lrv_markers = {100: 'o', 1000: 's'}
    actual_reps = None

    for lrv in [100, 1000]:
        _flight_files = sorted(cdir.glob(f'lrv_{lrv}_flights_*.csv'),
                               key=lambda f: int(re.search(r'flights_(\d+)', f.name).group(1)))
        if not _flight_files:
            continue
        nf_vals = []
        mu_stars = []
        eta_maxes = []
        for f in _flight_files:
            nf = int(re.search(r'flights_(\d+)', f.name).group(1))
            df = load_csv(f)
            if df.empty:
                continue
            mu_s, eta_m = get_mu_star(df)
            nf_vals.append(nf)
            mu_stars.append(mu_s)
            eta_maxes.append(eta_m)
            if actual_reps is None:
                actual_reps = int(df['reps'].iloc[0])
        if nf_vals:
            c, m = _lrv_colors[lrv], _lrv_markers[lrv]
            axes[0].plot(nf_vals, mu_stars, f'{m}-', color=c, markersize=7, linewidth=2,
                         label=f'λ/rv={lrv}')
            axes[1].plot(nf_vals, eta_maxes, f'{m}-', color=c, markersize=7, linewidth=2,
                         label=f'λ/rv={lrv}')

    axes[0].axhline(y=2.0, color='gray', linestyle=':', alpha=0.5)
    for ax, ylabel in [(axes[0], 'μ*'), (axes[1], 'λη_max')]:
        ax.set_xscale('log')
        ax.set_xlabel('Number of flights')
        ax.set_ylabel(ylabel)
        ax.set_title(f'{ylabel} convergence (destructive)')
        ax.legend(fontsize=13)
        ax.grid(True, alpha=0.3)

    reps_str = f'reps={actual_reps}' if actual_reps else ''
    plt.suptitle(f'Convergence with Number of Flights ({reps_str})',
                 fontsize=18, fontweight='bold')
    plt.tight_layout()
    p = OUTDIR / 'convergence_flights.png'
    plt.savefig(p); plt.close()
    print(f'  {p}')

    # --- Convergence in reps (destructive + non-destructive, both lrv) ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    actual_flights = None

    for lrv in [100, 1000]:
        c_d = _lrv_colors[lrv]
        # Destructive reps
        _reps_files = sorted([f for f in cdir.glob(f'lrv_{lrv}_reps_*.csv') if 'nd_reps' not in f.name],
                             key=lambda f: int(re.search(r'reps_(\d+)', f.name).group(1)))
        nr_vals = []
        mu_stars_r = []
        eta_maxes_r = []
        for f in _reps_files:
            nr = int(re.search(r'reps_(\d+)', f.name).group(1))
            df = load_csv(f)
            if df.empty:
                continue
            mu_s, eta_m = get_mu_star(df)
            nr_vals.append(nr)
            mu_stars_r.append(mu_s)
            eta_maxes_r.append(eta_m)
            if actual_flights is None:
                actual_flights = int(df['flights'].iloc[0])
        if nr_vals:
            axes[0].plot(nr_vals, mu_stars_r, '^-', color=c_d, markersize=7, linewidth=2,
                         label=f'Destr. λ/rv={lrv}')
            axes[1].plot(nr_vals, eta_maxes_r, '^-', color=c_d, markersize=7, linewidth=2,
                         label=f'Destr. λ/rv={lrv}')

        # Non-destructive reps
        files_nd = sorted(cdir.glob(f'lrv_{lrv}_nd_reps_*.csv'),
                          key=lambda f: int(re.search(r'nd_reps_(\d+)', f.name).group(1)))
        nr_nd = []
        mu_nd = []
        eta_nd = []
        for f in files_nd:
            nr = int(re.search(r'nd_reps_(\d+)', f.name).group(1))
            df = load_csv(f)
            if df.empty:
                continue
            mu_s, eta_m = get_mu_star(df)
            nr_nd.append(nr)
            mu_nd.append(mu_s)
            eta_nd.append(eta_m)
        if nr_nd:
            c_nd = C_VNONDESTR if lrv == 100 else '#8c564b'
            axes[0].plot(nr_nd, mu_nd, 'v--', color=c_nd, markersize=7, linewidth=2,
                         label=f'Non-destr. λ/rv={lrv}')
            axes[1].plot(nr_nd, eta_nd, 'v--', color=c_nd, markersize=7, linewidth=2,
                         label=f'Non-destr. λ/rv={lrv}')

    axes[0].axhline(y=2.0, color='gray', linestyle=':', alpha=0.5)
    for ax, ylabel in [(axes[0], 'μ*'), (axes[1], 'λη_max')]:
        ax.set_xscale('log')
        ax.set_xlabel('Number of repetitions')
        ax.set_ylabel(ylabel)
        ax.set_title(f'{ylabel} convergence')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)

    flights_str = f'flights={actual_flights:,}' if actual_flights else ''
    plt.suptitle(f'Convergence with Number of Repetitions ({flights_str})',
                 fontsize=18, fontweight='bold')
    plt.tight_layout()
    p = OUTDIR / 'convergence_reps.png'
    plt.savefig(p); plt.close()
    print(f'  {p}')

    # --- Full curves at different reps (non-destructive, lrv_100) ---
    files_nd_100 = sorted(cdir.glob('lrv_100_nd_reps_*.csv'),
                          key=lambda f: int(re.search(r'nd_reps_(\d+)', f.name).group(1)))
    if not files_nd_100:
        files_nd_100 = sorted(cdir.glob('*nd_reps_*.csv'),
                              key=lambda f: int(re.search(r'nd_reps_(\d+)', f.name).group(1)))
    if files_nd_100:
        nr_nd_vals = [int(re.search(r'nd_reps_(\d+)', f.name).group(1)) for f in files_nd_100]
        fig, ax = plt.subplots(figsize=(10, 6))
        cmap = plt.cm.plasma
        norm = mcolors.LogNorm(vmin=max(5, min(nr_nd_vals)), vmax=max(nr_nd_vals))
        for f in files_nd_100:
            nr = int(re.search(r'nd_reps_(\d+)', f.name).group(1))
            df = load_csv(f)
            if df.empty:
                continue
            mu = df['mu'].values
            y = df['normalized_rate_unique_global'].values
            v = np.isfinite(y)
            cn = cmap(norm(nr))
            ax.plot(mu[v], y[v], 'o-', color=cn,
                    label=f'reps={nr}', markersize=4, linewidth=1.5)
            if 'normalized_rate_unique_std' in df.columns:
                std = df['normalized_rate_unique_std'].values
                ax.fill_between(mu[v], (y - std)[v], (y + std)[v],
                                color=cn, alpha=_BAND_ALPHA, linewidth=0)
        ax.axvline(x=2.0, color='gray', linestyle=':', alpha=0.5)
        ax.set_xlabel('μ')
        ax.set_ylabel('λη (unique, pooled)')
        ax.set_title('Non-destructive Viswanathan: Effect of Repetitions (λ/rv=100)',
                     fontweight='bold')
        ax.legend(fontsize=12)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        p = OUTDIR / 'convergence_nd_curves.png'
        plt.savefig(p); plt.close()
        print(f'  {p}')


# ═══════════════════════════════════════════════════════════════
# 6. ANALYTICAL COMPARISON
# ═══════════════════════════════════════════════════════════════

def _compute_mean_l(mu_arr, lmin, lmax):
    """Compute <l> for truncated power-law P(l) ∝ l^{-μ} on [lmin, lmax]."""
    lmin, lmax = float(lmin), float(lmax)
    ml = np.zeros_like(mu_arr, dtype=float)
    for i, m in enumerate(mu_arr):
        a = 1 - m
        a1 = a + 1  # = 2 - μ
        if abs(a) < 1e-10:       # μ = 1
            ml[i] = (lmax - lmin) / np.log(lmax / lmin)
        elif abs(a1) < 1e-10:    # μ = 2
            ml[i] = np.log(lmax / lmin) / (lmin**(-1) - lmax**(-1))
        else:
            num = lmax**a1 - lmin**a1
            den = lmax**a - lmin**a
            ml[i] = (a / a1) * num / den if abs(den) > 1e-30 else lmin
    return ml


def _compute_viswanathan_eta(mu_arr, rv, lam):
    """Viswanathan (1999) analytical prediction: η = 1/(⟨l⟩_MF · N).

    ⟨l⟩_MF: mean-field mean flight distance (eq 2) — uses untruncated P(l)~l^{-μ}
    on [rv, ∞) with flights longer than λ contributing λ (mean-field truncation).

    N = (λ/rv)^{(μ-1)/2}: mean number of flights between successive target
    encounters (eq 5, non-destructive scaling). This scaling is appropriate
    for the truncated simulation (l_max = λ) where the walker revisits areas.

    Returns λη (dimensionless efficiency).
    """
    rv, lam = float(rv), float(lam)
    R = lam / rv
    lam_eta = np.zeros_like(mu_arr, dtype=float)
    for i, m in enumerate(mu_arr):
        # Mean-field <l> from Viswanathan eq (2)
        if abs(m - 1.0) < 1e-8:
            ml = lam  # limit: all flights >> λ, truncated at λ
        elif abs(m - 2.0) < 1e-8:
            ml = rv * (np.log(R) + 1)
        else:
            Rp = R ** (2 - m)  # R^{2-μ}
            f = (m - 1) * (Rp - 1) / (2 - m) + Rp
            ml = rv * f
        # N from eq (5): number of flights between targets
        N = R ** ((m - 1) / 2)
        lam_eta[i] = lam / (ml * N)
    return lam_eta


def plot_analytical():
    """
    Compare simulation with Viswanathan (1999) analytical prediction.

    Uses η = 1/(⟨l⟩_MF · N) where:
      ⟨l⟩_MF — mean-field mean flight distance (eq 2)
      N = (λ/rv)^{(μ-1)/2} — flights between targets (eq 5)

    Normalization uses least-squares fit to the simulation data.
    """
    print("=== 6. Analytical comparison ===")

    pdir = Path(PHASE_DIR)
    if not pdir.exists():
        print("  Skipping (no data).")
        return

    # Use the cleanest data: destructive Viswanathan at moderate sparsity
    lrv_test = 100
    csv_path = pdir / f'lrv_{lrv_test}_visw_destr.csv'
    if not csv_path.exists():
        print(f"  No data for λ/rv={lrv_test}")
        return

    df = load_csv(csv_path)
    mu = df['mu'].values
    y = df['normalized_rate_unique_global'].values
    rv = df['rv'].iloc[0]
    lam = 1.0 / (2 * df['rho'].iloc[0] * rv)

    mu_theory = np.linspace(1.01, 3.0, 200)
    eta_theory = _compute_viswanathan_eta(mu_theory, rv, lam)

    v = np.isfinite(y)
    # Least-squares fit: minimize Σ(y - C·η_theory)²
    eta_at_sim = _compute_viswanathan_eta(mu[v], rv, lam)
    scale = np.dot(y[v], eta_at_sim) / np.dot(eta_at_sim, eta_at_sim)
    eta_theory_scaled = eta_theory * scale

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(mu[v], y[v], 'o-', color=C_VDESTR, label='Simulation (destructive)',
            markersize=6, linewidth=2)
    if 'normalized_rate_unique_std' in df.columns:
        std_vals = df['normalized_rate_unique_std'].values
        ax.fill_between(mu[v], (y - std_vals)[v], (y + std_vals)[v],
                        color=C_VDESTR, alpha=_BAND_ALPHA, linewidth=0,
                        label=r'$\pm 1\sigma$ (per-run)')
    ax.plot(mu_theory, eta_theory_scaled, '--', color='black',
            label=r'$\eta = 1/(\langle l \rangle \cdot N)$ (Viswanathan, 1999)',
            linewidth=2, alpha=0.7)
    ax.axvline(x=2.0, color='gray', linestyle=':', alpha=0.5, label='μ = 2')
    ax.set_xlabel('μ', fontsize=16)
    ax.set_ylabel('λη (unique, pooled)', fontsize=16)
    ax.set_title(f'Comparison with Analytical Prediction (λ/rv={lrv_test})',
                 fontsize=18, fontweight='bold')
    ax.legend(fontsize=14)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    p = OUTDIR / 'analytical_comparison.png'
    plt.savefig(p); plt.close()
    print(f'  {p}')

    # Multi-sparsity analytical comparison
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    test_lrv = [20, 100, 500, 2000]
    for idx, lrv in enumerate(test_lrv):
        ax = axes[idx]
        csv_path = pdir / f'lrv_{lrv}_visw_destr.csv'
        if not csv_path.exists():
            continue
        df = load_csv(csv_path)
        mu_d = df['mu'].values
        y_d = df['normalized_rate_unique_global'].values
        rv_d = df['rv'].iloc[0]
        lam_d = 1.0 / (2 * df['rho'].iloc[0] * rv_d)
        v = np.isfinite(y_d)

        et = _compute_viswanathan_eta(mu_theory, rv_d, lam_d)
        # Least-squares normalization
        et_sim = _compute_viswanathan_eta(mu_d[v], rv_d, lam_d)
        sc = np.dot(y_d[v], et_sim) / np.dot(et_sim, et_sim) if v.sum() > 0 else 1.0
        et_s = et * sc

        ax.plot(mu_d[v], y_d[v], 'o-', color=C_VDESTR, markersize=4, linewidth=2,
                label='Simulation')
        if 'normalized_rate_unique_std' in df.columns:
            std_d = df['normalized_rate_unique_std'].values
            ax.fill_between(mu_d[v], (y_d - std_d)[v], (y_d + std_d)[v],
                            color=C_VDESTR, alpha=_BAND_ALPHA, linewidth=0)
        ax.plot(mu_theory, et_s, '--', color='black', linewidth=1.5, alpha=0.7,
                label=r'$\eta = 1/(\langle l \rangle \cdot N)$')
        ax.axvline(x=2.0, color='gray', linestyle=':', alpha=0.4)
        ax.set_xlabel('μ')
        ax.set_ylabel('λη')
        ax.set_title(f'λ/rv = {lrv}')
        ax.legend(fontsize=12)
        ax.grid(True, alpha=0.3)

    plt.suptitle('Analytical vs Simulation at Different Sparsity Levels',
                 fontsize=18, fontweight='bold')
    plt.tight_layout()
    p = OUTDIR / 'analytical_multi_sparsity.png'
    plt.savefig(p); plt.close()
    print(f'  {p}')


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    results = plot_phase_diagram()
    plot_trapping_heatmap()
    plot_mfpt()
    plot_sensitivity()
    plot_convergence()
    plot_analytical()

    print(f'\nAll plots saved to {OUTDIR}/')
    print('Done.')
