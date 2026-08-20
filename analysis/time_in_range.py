#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Time-in-range comparison: before 2026-07-01 vs. on/after 2026-07-01.

Reads a LibreView CSV export, splits the readings at a cutoff date, and plots
two stacked bars showing the percentage of readings falling in each glucose
range for each period.

Usage: python3 time_in_range.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Cutoff separating the "before" period from the "on/after" period.
CUTOFF_DATE = pd.Timestamp("2026-07-01")

# Anonymised subject identifier shown on the chart — keep real names out of the repo.
PATIENT_LABEL = "Patient A"

# Paths are resolved from this file, so they hold no matter which directory
# you run the script from. This script lives in analysis/, one level down.
PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR    = PROJECT_DIR / "data"
FIGURES_DIR = PROJECT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

GLUCOSE_CSV = DATA_DIR / "glucose_2026-07-20.csv"
OUTPUT_PNG  = FIGURES_DIR / "time_in_range_comparison.png"

# LibreView ranges (bounds in mmol/L), ordered from lowest to highest.
RANGES = [
    ("Very low",  -np.inf, 3.0,   "#B71C1C"),
    ("Low",        3.0,    3.9,   "#E53935"),
    ("In target",  3.9,   10.0,   "#43A047"),
    ("High",      10.0,   13.9,   "#FDD835"),
    ("Very high", 13.9,   np.inf, "#FB8C00"),
]

# ── Loading (same parsing as the other scripts) ───────────────────
df = pd.read_csv(GLUCOSE_CSV, skiprows=1, low_memory=False)
# Record type 0 = automatic glucose reading; other types are notes, insulin, etc.
df = df[df["Tietueen tyyppi"] == 0].copy()
df["timestamp"] = pd.to_datetime(df["Laitteen aikaleima"], format="%d-%m-%Y %H:%M")
# The export uses a comma as decimal separator.
df["glucose"]   = df["Aikaisempi glukoosi mmol/L"].str.replace(",", ".").astype(float)
df = df.sort_values("timestamp").reset_index(drop=True)


# ── Percentage of readings per range ──────────────────────────────
def distribution(sub):
    """Return (percentage per range, number of readings) for a subset."""
    n = len(sub)
    if n == 0:
        return {name: 0.0 for name, *_ in RANGES}, 0
    out = {}
    for name, lo, hi, _ in RANGES:
        mask = (sub["glucose"] >= lo) & (sub["glucose"] < hi)
        out[name] = 100 * mask.sum() / n
    return out, n


def date_span(sub):
    """Return a 'first date → last date' label for the readings in a subset."""
    if len(sub) == 0:
        return "no data"
    return f"{sub['timestamp'].min():%Y-%m-%d} → {sub['timestamp'].max():%Y-%m-%d}"


def range_label(name, lo, hi):
    """Return a legend label pairing the range name with its bounds.

    Bounds are half-open [lo, hi), matching how distribution() classifies.
    """
    if lo == -np.inf:
        bounds = f"< {hi:.1f}"
    elif hi == np.inf:
        bounds = f"≥ {lo:.1f}"
    else:
        bounds = f"{lo:.1f} – {hi:.1f}"
    return f"{name}  ({bounds})"


before = df[df["timestamp"] < CUTOFF_DATE]
after  = df[df["timestamp"] >= CUTOFF_DATE]

dist_before, n_before = distribution(before)
dist_after,  n_after  = distribution(after)

# ── Chart: two stacked bars side by side ──────────────────────────
fig, ax = plt.subplots(figsize=(7, 8))
fig.patch.set_facecolor("white")

# Each tick states the period and the actual date span it covers.
labels = [
    f"Before {CUTOFF_DATE:%b %d}\n{date_span(before)}",
    f"{CUTOFF_DATE:%b %d} onwards\n{date_span(after)}",
]
x = [0, 1]
width = 0.55

for period, dist in [(0, dist_before), (1, dist_after)]:
    bottom = 0
    for name, lo, hi, color in RANGES:
        value = dist[name]
        ax.bar(period, value, width, bottom=bottom, color=color,
               edgecolor="white", linewidth=1.2)
        if value >= 3:  # label only segments big enough to stay readable
            ax.text(period, bottom + value / 2, f"{value:.0f}%",
                    ha="center", va="center", fontsize=10,
                    fontweight="bold",
                    color="white" if name != "High" else "#333")
        bottom += value

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=10, fontweight="bold")
ax.set_ylim(0, 100)
ax.set_ylabel("% of time (share of readings)", fontsize=10)
ax.set_title(f"Time in range — comparison\n{PATIENT_LABEL}",
             fontsize=13, fontweight="bold", pad=12)
ax.set_facecolor("#FAFAFA")
ax.grid(axis="y", alpha=0.25)

# Legend, ordered highest range first to match the stacking order.
legend_handles = [Patch(facecolor=c, label=range_label(n, lo, hi))
                  for n, lo, hi, c in reversed(RANGES)]
legend = ax.legend(handles=legend_handles, loc="center left",
                   bbox_to_anchor=(1.02, 0.5), fontsize=9, frameon=False,
                   title="Glucose (mmol/L)")
legend.get_title().set_fontweight("bold")
legend.get_title().set_fontsize(9)

# Sample-size caption
ax.text(0.5, -0.12,
        f"n = {n_before} readings (before)   |   n = {n_after} readings (after)",
        transform=ax.transAxes, ha="center", fontsize=8, color="#555")

plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=150, bbox_inches="tight", facecolor="white")
print(f"Saved: {OUTPUT_PNG.relative_to(PROJECT_DIR)}")
print(f"\nBefore {CUTOFF_DATE:%Y-%m-%d} (n={n_before}):",
      {k: round(v, 1) for k, v in dist_before.items()})
print(f"After  {CUTOFF_DATE:%Y-%m-%d} (n={n_after}):",
      {k: round(v, 1) for k, v in dist_after.items()})
plt.show()