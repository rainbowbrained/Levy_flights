#!/usr/bin/env python3
"""
Сравнение трёх методов Lévy поиска:
1. No interrupt - полёт продолжается после обнаружения
2. Teleport - телепорт на lc от таргета
3. Viswanathan - рестарт из позиции таргета
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def compare_three_methods(csv_files, labels, output_dir='figures/comparison'):
    """
    Сравнить три метода Lévy поиска.

    Args:
        csv_files: список из 3 путей к CSV файлам
        labels: список из 3 меток для легенды
        output_dir: директория для сохранения графиков
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Загрузить данные
    dfs = [pd.read_csv(f) for f in csv_files]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # blue, orange, green

    print("=" * 70)
    print("📊 СРАВНЕНИЕ ТРЁХ МЕТОДОВ LÉVY ПОИСКА")
    print("=" * 70)
    for i, (label, df) in enumerate(zip(labels, dfs)):
        print(f"{i+1}. {label}: {len(df)} точек, μ ∈ [{df['mu'].min():.1f}, {df['mu'].max():.1f}]")
    print()

    # ===== График 1: Normalized rate (unique) vs μ =====
    fig, ax = plt.subplots(figsize=(10, 6))

    for df, label, color in zip(dfs, labels, colors):
        mu = df['mu'].values
        rate_mean = df['normalized_rate_unique_mean'].values
        rate_std = df['normalized_rate_unique_std'].values

        # Фильтр NaN
        valid = np.isfinite(rate_mean) & np.isfinite(rate_std)
        mu = mu[valid]
        rate_mean = rate_mean[valid]
        rate_std = rate_std[valid]

        ax.plot(mu, rate_mean, 'o-', color=color, label=label,
               markersize=6, linewidth=2, alpha=0.8)
        ax.fill_between(mu, rate_mean - rate_std, rate_mean + rate_std,
                        color=color, alpha=0.2, linewidth=0)

    ax.set_xlabel('μ (Lévy exponent)', fontsize=12, fontweight='bold')
    ax.set_ylabel('λη (unique captures)', fontsize=12, fontweight='bold')
    ax.set_title('Normalized Search Efficiency vs μ (unique captures)',
                fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    output_path = Path(output_dir) / 'comparison_normalized_rate_unique.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✅ Сохранено: {output_path}")

    # ===== График 2: Normalized rate (visits) vs μ =====
    fig, ax = plt.subplots(figsize=(10, 6))

    for df, label, color in zip(dfs, labels, colors):
        mu = df['mu'].values
        rate_mean = df['normalized_rate_visits_mean'].values
        rate_std = df['normalized_rate_visits_std'].values

        valid = np.isfinite(rate_mean) & np.isfinite(rate_std)
        mu = mu[valid]
        rate_mean = rate_mean[valid]
        rate_std = rate_std[valid]

        ax.plot(mu, rate_mean, 'o-', color=color, label=label,
               markersize=6, linewidth=2, alpha=0.8)
        ax.fill_between(mu, rate_mean - rate_std, rate_mean + rate_std,
                        color=color, alpha=0.2, linewidth=0)

    ax.set_xlabel('μ (Lévy exponent)', fontsize=12, fontweight='bold')
    ax.set_ylabel('λη (all visits)', fontsize=12, fontweight='bold')
    ax.set_title('Normalized Search Efficiency vs μ (all visits)',
                fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    output_path = Path(output_dir) / 'comparison_normalized_rate_visits.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✅ Сохранено: {output_path}")

    # ===== График 3: Visits vs μ =====
    fig, ax = plt.subplots(figsize=(10, 6))

    for df, label, color in zip(dfs, labels, colors):
        mu = df['mu'].values
        visits_mean = df['visits_mean'].values
        visits_std = df['visits_std'].values

        valid = np.isfinite(visits_mean) & np.isfinite(visits_std)
        mu = mu[valid]
        visits_mean = visits_mean[valid]
        visits_std = visits_std[valid]

        ax.plot(mu, visits_mean, 'o-', color=color, label=label,
               markersize=6, linewidth=2, alpha=0.8)
        ax.fill_between(mu, visits_mean - visits_std, visits_mean + visits_std,
                        color=color, alpha=0.2, linewidth=0)

    ax.set_xlabel('μ (Lévy exponent)', fontsize=12, fontweight='bold')
    ax.set_ylabel('All visits (mean)', fontsize=12, fontweight='bold')
    ax.set_title('Total Visits vs μ', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    output_path = Path(output_dir) / 'comparison_visits.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✅ Сохранено: {output_path}")

    # ===== График 4: Unique captures vs μ =====
    fig, ax = plt.subplots(figsize=(10, 6))

    for df, label, color in zip(dfs, labels, colors):
        mu = df['mu'].values
        caps_mean = df['unique_caps_mean'].values
        caps_std = df['unique_caps_std'].values

        valid = np.isfinite(caps_mean) & np.isfinite(caps_std)
        mu = mu[valid]
        caps_mean = caps_mean[valid]
        caps_std = caps_std[valid]

        ax.plot(mu, caps_mean, 'o-', color=color, label=label,
               markersize=6, linewidth=2, alpha=0.8)
        ax.fill_between(mu, caps_mean - caps_std, caps_mean + caps_std,
                        color=color, alpha=0.2, linewidth=0)

    ax.set_xlabel('μ (Lévy exponent)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Unique captures (mean)', fontsize=12, fontweight='bold')
    ax.set_title('Unique Target Captures vs μ', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    output_path = Path(output_dir) / 'comparison_unique_captures.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✅ Сохранено: {output_path}")

    # ===== График 5: MFPT vs μ =====
    fig, ax = plt.subplots(figsize=(10, 6))

    for df, label, color in zip(dfs, labels, colors):
        mu = df['mu'].values
        mfpt = df['mfpt'].values

        valid = np.isfinite(mfpt)
        mu = mu[valid]
        mfpt = mfpt[valid]

        if len(mu) > 0:
            ax.plot(mu, mfpt, 'o-', color=color, label=label,
                   markersize=6, linewidth=2, alpha=0.8)

    ax.set_xlabel('μ (Lévy exponent)', fontsize=12, fontweight='bold')
    ax.set_ylabel('MFPT (Mean First Passage Time)', fontsize=12, fontweight='bold')
    ax.set_title('MFPT vs μ', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    output_path = Path(output_dir) / 'comparison_mfpt.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✅ Сохранено: {output_path}")

    # ===== Сводная таблица оптимальных μ =====
    print()
    print("=" * 70)
    print("📊 ОПТИМАЛЬНЫЕ ЗНАЧЕНИЯ μ*")
    print("=" * 70)

    for label, df in zip(labels, dfs):
        # Найти μ* для unique
        valid = np.isfinite(df['normalized_rate_unique_mean'])
        if valid.sum() > 0:
            idx_max = df.loc[valid, 'normalized_rate_unique_mean'].idxmax()
            mu_star = df.loc[idx_max, 'mu']
            rate_max = df.loc[idx_max, 'normalized_rate_unique_mean']
            print(f"{label}:")
            print(f"  μ* (unique) = {mu_star:.1f}, λη_max = {rate_max:.4f}")

        # Найти μ* для visits
        valid = np.isfinite(df['normalized_rate_visits_mean'])
        if valid.sum() > 0:
            idx_max = df.loc[valid, 'normalized_rate_visits_mean'].idxmax()
            mu_star = df.loc[idx_max, 'mu']
            rate_max = df.loc[idx_max, 'normalized_rate_visits_mean']
            print(f"  μ* (visits) = {mu_star:.1f}, λη_max = {rate_max:.4f}")
        print()

    print("=" * 70)
    print(f"✅ Все графики сохранены в {output_dir}/")
    print("=" * 70)


if __name__ == '__main__':
    csv_files = [
        'results_no_interrupt.csv',
        'results_teleport.csv',
        'results_viswanathan.csv'
    ]

    labels = [
        'No interrupt (продолжение полёта)',
        'Teleport (lc=5)',
        'Viswanathan (рестарт из таргета)'
    ]

    # Проверить что все файлы существуют
    for f in csv_files:
        if not Path(f).exists():
            print(f"❌ Файл не найден: {f}")
            sys.exit(1)

    compare_three_methods(csv_files, labels, output_dir='figures/comparison')
