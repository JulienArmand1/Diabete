#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Time spent in hypoglycaemia, before vs. after a cutoff date.

Runs at two thresholds and writes one figure per threshold:
  < 3.9 mmol/L — all lows
  < 3.0 mmol/L — severe lows only

The statistical unit is the DAY, not the individual reading: readings 5 min
apart are not independent, so treating them as such would overstate the
evidence (pseudoreplication).

Significance comes from a permutation test on the day labels, and the 95%
confidence interval from a bootstrap resampling whole days. Neither assumes
a normal distribution, which matters because daily time-below is skewed.

Usage: python3 time_in_hypo_stats.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ══════════════════════════════════════════════════════════════════
# Paths are resolved from this file, so they hold no matter which directory
# you run the script from. This script lives in analysis/, one level down.
PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR    = PROJECT_DIR / "data"
FIGURES_DIR = PROJECT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

GLUCOSE_CSV = DATA_DIR / "glucose_2026-07-20.csv"

# Anonymised subject identifier shown on the charts.
PATIENT_LABEL = "Patient A"

CUTOFF_DATE = pd.Timestamp("2026-07-01")
LOW_LIMITS  = [3.9, 3.0]   # one figure per threshold, in mmol/L
N_PERM      = 100_000      # permutations
N_BOOT      = 100_000      # bootstrap resamples
ONE_SIDED   = True         # H1: after < before
MAX_GAP_MIN = 20           # cap per-reading duration, to ignore sensor outages

COLOR_BEFORE = "#5C6BC0"
COLOR_AFTER  = "#26A69A"
# ══════════════════════════════════════════════════════════════════


def load_readings():
    """Load the export, keep automatic readings, and credit each with a duration."""
    df = pd.read_csv(GLUCOSE_CSV, skiprows=1, low_memory=False)
    # Record type 0 = automatic glucose reading; other types are notes, insulin, etc.
    df = df[df["Tietueen tyyppi"] == 0].copy()
    df["timestamp"] = pd.to_datetime(df["Laitteen aikaleima"], format="%d-%m-%Y %H:%M")
    # The export uses a comma as decimal separator.
    df["glucose"] = df["Aikaisempi glukoosi mmol/L"].str.replace(",", ".").astype(float)
    df = df.sort_values("timestamp").reset_index(drop=True)
    df["date"] = df["timestamp"].dt.date

    # Each reading is credited with the time until the next one, capped so that
    # a sensor outage does not count as hours spent at its last value.
    df["dt_min"] = df["timestamp"].diff().shift(-1).dt.total_seconds() / 60
    df["dt_min"] = df["dt_min"].clip(upper=MAX_GAP_MIN).fillna(0)
    return df


def daily_percent_below(df, limit):
    """Return, per day, the percentage of recorded time spent below `limit`."""
    total_min = df.groupby("date")["dt_min"].sum()
    low_min   = (df[df["glucose"] < limit].groupby("date")["dt_min"].sum()
                 .reindex(total_min.index).fillna(0))
    daily = (100 * low_min / total_min).replace([np.inf, -np.inf], np.nan).dropna()

    dates = pd.to_datetime(daily.index)
    return daily[dates < CUTOFF_DATE].values, daily[dates >= CUTOFF_DATE].values


def permutation_test(before, after):
    """Shuffle the day labels to build the null distribution of the difference."""
    observed = after.mean() - before.mean()      # negative if lows became rarer
    pool = np.concatenate([before, after])
    n_before = len(before)

    null_diffs = np.empty(N_PERM)
    for i in range(N_PERM):
        perm = np.random.permutation(pool)
        null_diffs[i] = perm[n_before:].mean() - perm[:n_before].mean()

    if ONE_SIDED:                                # H1: after < before
        p_value = (np.sum(null_diffs <= observed) + 1) / (N_PERM + 1)
    else:
        p_value = (np.sum(np.abs(null_diffs) >= abs(observed)) + 1) / (N_PERM + 1)
    return observed, null_diffs, p_value


def bootstrap_reduction(before, after):
    """Resample whole days to get the sampling distribution of the reduction."""
    boot_abs = np.empty(N_BOOT)
    boot_rel = np.empty(N_BOOT)
    for i in range(N_BOOT):
        a = np.random.choice(before, len(before), replace=True)
        b = np.random.choice(after,  len(after),  replace=True)
        boot_abs[i] = a.mean() - b.mean()
        boot_rel[i] = 100 * (a.mean() - b.mean()) / a.mean()
    return boot_abs, boot_rel


def ci95(x):
    return np.percentile(x, 2.5), np.percentile(x, 97.5)


def draw_group(ax, values, x, color):
    """Scatter one group's daily values, with its mean as a thick bar."""
    jitter = np.random.RandomState(1).uniform(-0.09, 0.09, len(values))
    ax.scatter(x + jitter, values, s=26, color=color, alpha=0.75,
               edgecolors="white", linewidths=0.5, zorder=3)
    ax.hlines(values.mean(), x - 0.26, x + 0.26, color=color, lw=2.6, zorder=4)
    ax.text(x + 0.30, values.mean(), f"{values.mean():.2f}%", va="center",
            ha="left", fontsize=10, fontweight="bold", color=color, zorder=5)


def make_figure(limit, before, after, observed, null_diffs, p_text,
                boot_abs, abs_drop, lo_abs, hi_abs, out_png):
    """Three panels: the daily values, the null distribution, the bootstrap CI."""
    fig, (ax_days, ax_perm, ax_boot) = plt.subplots(1, 3, figsize=(16, 5.5))
    fig.patch.set_facecolor("white")

    # ── Panel 1: the daily values themselves ──────────────────────
    draw_group(ax_days, before, 0, COLOR_BEFORE)
    draw_group(ax_days, after,  1, COLOR_AFTER)

    ax_days.set_xticks([0, 1])
    ax_days.set_xticklabels([f"Before {CUTOFF_DATE:%b %d}\nn = {len(before)} days",
                             f"{CUTOFF_DATE:%b %d} onwards\nn = {len(after)} days"],
                            fontsize=10, fontweight="bold")
    ax_days.set_xlim(-0.5, 1.5)
    ax_days.set_ylabel(f"% of day below {limit} mmol/L", fontsize=10)
    ax_days.set_title("Daily time spent low\none point = one day, bar = mean",
                      fontsize=11, fontweight="bold")
    ax_days.grid(axis="y", alpha=0.25)
    ax_days.set_facecolor("#FAFAFA")

    # ── Panel 2: permutation null distribution ────────────────────
    ax_perm.hist(null_diffs, bins=80, color="#B0BEC5", edgecolor="white", linewidth=0.4)
    ax_perm.axvline(observed, color="#C62828", lw=2.4, zorder=4)
    ax_perm.text(observed, ax_perm.get_ylim()[1] * 0.96, f" observed {observed:+.2f}",
                 color="#C62828", fontsize=9, fontweight="bold", va="top")

    ax_perm.set_xlabel("Difference in means (points), under H₀", fontsize=9)
    ax_perm.set_ylabel("Permutations", fontsize=9)
    ax_perm.set_title(f"Permutation test — {N_PERM:,} shuffles\n"
                      f"p = {p_text} ({'one-sided' if ONE_SIDED else 'two-sided'})",
                      fontsize=11, fontweight="bold")
    ax_perm.grid(axis="y", alpha=0.25)
    ax_perm.set_facecolor("#FAFAFA")

    # ── Panel 3: bootstrap distribution of the reduction ──────────
    ax_boot.hist(boot_abs, bins=80, color="#80CBC4", edgecolor="white", linewidth=0.4)
    ax_boot.axvspan(lo_abs, hi_abs, color="#26A69A", alpha=0.16, zorder=0)
    ax_boot.axvline(abs_drop, color="#00695C", lw=2.4, zorder=4)
    ax_boot.axvline(0, color="#757575", lw=1.2, ls="--", zorder=3)
    ax_boot.text(abs_drop, ax_boot.get_ylim()[1] * 0.96, f" {abs_drop:.2f} pts",
                 color="#00695C", fontsize=9, fontweight="bold", va="top")

    ax_boot.set_xlabel("Reduction in time spent low (percentage points)", fontsize=9)
    ax_boot.set_ylabel("Bootstrap resamples", fontsize=9)
    ax_boot.set_title(f"Bootstrap 95% CI — {N_BOOT:,} resamples\n"
                      f"{abs_drop:.2f} pts  [{lo_abs:.2f} ; {hi_abs:.2f}]",
                      fontsize=11, fontweight="bold")
    ax_boot.grid(axis="y", alpha=0.25)
    ax_boot.set_facecolor("#FAFAFA")

    severity = "all lows" if limit >= 3.9 else "severe lows"
    fig.suptitle(f"Time spent below {limit} mmol/L ({severity}) — before vs. after "
                 f"{CUTOFF_DATE:%d %b %Y} — {PATIENT_LABEL}",
                 fontsize=13, fontweight="bold")

    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight", facecolor="white")


def run(df, limit):
    """Full analysis at one threshold: stats to the console, one figure to disk."""
    np.random.seed(0)          # each threshold is independently reproducible

    before, after = daily_percent_below(df, limit)
    m_before, m_after = before.mean(), after.mean()
    abs_drop = m_before - m_after
    rel_drop = 100 * abs_drop / m_before

    observed, null_diffs, p_value = permutation_test(before, after)
    boot_abs, boot_rel = bootstrap_reduction(before, after)
    lo_abs, hi_abs = ci95(boot_abs)
    lo_rel, hi_rel = ci95(boot_rel)

    # With 100k permutations the smallest reachable p-value is 1/(N_PERM+1); a
    # result sitting on that floor means "smaller than this test can resolve".
    p_floor = 1 / (N_PERM + 1)
    p_text = f"< {p_floor:.5f}" if p_value <= p_floor else f"{p_value:.5f}"

    out_png = FIGURES_DIR / f"time_in_hypo_{str(limit).replace('.', '_')}.png"

    print("═" * 60)
    print(f"Subject               : {PATIENT_LABEL}")
    print(f"Low threshold         : < {limit} mmol/L")
    print(f"Days before / after   : {len(before)} / {len(after)}")
    print(f"Mean % low, before    : {m_before:.2f} %  (median {np.median(before):.2f})")
    print(f"Mean % low, after     : {m_after:.2f} %  (median {np.median(after):.2f})")
    print("─" * 60)
    print(f"Observed difference   : {observed:+.2f} percentage points")
    print(f"p-value ({'one-sided' if ONE_SIDED else 'two-sided'})    : {p_text}")
    print("─" * 60)
    print(f"Absolute reduction    : {abs_drop:.2f} pts   [95% CI {lo_abs:.2f} ; {hi_abs:.2f}]")
    print(f"Relative reduction    : {rel_drop:.1f} %     [95% CI {lo_rel:.1f} ; {hi_rel:.1f}]")
    print("═" * 60)

    make_figure(limit, before, after, observed, null_diffs, p_text,
                boot_abs, abs_drop, lo_abs, hi_abs, out_png)
    print(f"Saved: {out_png.relative_to(PROJECT_DIR)}\n")


df = load_readings()
for low_limit in LOW_LIMITS:
    run(df, low_limit)

plt.show()