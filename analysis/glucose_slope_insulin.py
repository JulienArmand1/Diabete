#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Daily glucose, slope and insulin panels, with stars marking rapid excursions.

For each day in the export, three stacked panels are drawn:
  1. glucose trace with the target band shaded,
  2. first derivative (slope), filled orange when rising and blue when falling,
  3. insulin — basal rate as a step area, bolus doses as stems.

A star marks the first reading of each episode where the slope crosses a
threshold, so a long excursion gets one marker rather than a cloud of points.

Usage: python3 glucose_slope_insulin.py
"""

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker

# ══════════════════════════════════════════════════════════════════
# PATHS — resolved from this file, so they hold no matter which directory you
# run the script from. This script lives in analysis/, one level down.
PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR    = PROJECT_DIR / "data"
FIGURES_DIR = PROJECT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

GLUCOSE_CSV = DATA_DIR / "glucose_2026-07-20.csv"   # CGM export (LibreView)
BOLUS_CSV   = DATA_DIR / "bolus_data_1.csv"         # pump export; optional
BASAL_CSV   = DATA_DIR / "basal_data_1.csv"         # pump export; optional

# Anonymised subject identifier shown on the chart — keep real names out of the repo.
PATIENT_LABEL = "Patient A"

# THRESHOLDS — adjust to taste.
HIGH_THRESHOLD =  6.0   # mmol/L/h → orange star (rapid rise)
LOW_THRESHOLD  = -4.5   # mmol/L/h → blue star   (rapid fall)

# PNG export quality. 300 = print quality; 600 gets very large when the figure
# spans many days.
PNG_DPI = 300

# Output basename (a .pdf and a .png are written into figures/).
OUTPUT_BASE = FIGURES_DIR / "glucose_slope_insulin"
# ══════════════════════════════════════════════════════════════════


# ── Loading ───────────────────────────────────────────────────────
def read_export(path):
    """Read a pump or CGM export from data/, or return None if absent.

    Accepts a plain .csv or the same name still inside a .zip sitting in
    data/, so a downloaded archive works without unpacking it.
    """
    if path.exists():
        print(f"  [ok] {path.name}")
        return pd.read_csv(path, skiprows=1, low_memory=False)
    for archive in sorted(DATA_DIR.glob("*.zip")):
        try:
            with zipfile.ZipFile(archive) as zf:
                member = next((n for n in zf.namelist()
                               if n.split("/")[-1] == path.name), None)
                if member:
                    print(f"  [ok] {path.name}  (in {archive.name})")
                    with zf.open(member) as f:
                        return pd.read_csv(f, skiprows=1)
        except zipfile.BadZipFile:
            pass
    print(f"  [missing] {path.name}")
    return None


def find_column(columns, needle):
    """Find the column containing `needle` (tolerates non-breaking spaces)."""
    for c in columns:
        if needle.lower() in str(c).replace("\xa0", " ").lower():
            return c
    return None


def to_float(s):
    """Parse a numeric column that uses a comma decimal separator."""
    return pd.to_numeric(s.astype(str).str.replace("\xa0", "", regex=False)
                          .str.replace(",", ".", regex=False), errors="coerce")


print("Loading data...")
raw = read_export(GLUCOSE_CSV)
if raw is None:
    raise SystemExit(f"Glucose export not found: {GLUCOSE_CSV}")

# Record type 0 = automatic glucose reading; other types are notes, insulin, etc.
df = raw[raw["Tietueen tyyppi"] == 0].copy()
df["timestamp"] = pd.to_datetime(df["Laitteen aikaleima"], format="%d-%m-%Y %H:%M")
# The export uses a comma as decimal separator.
df["glucose"]   = df["Aikaisempi glukoosi mmol/L"].str.replace(",", ".").astype(float)
df = df.sort_values("timestamp").reset_index(drop=True)

# ── Smoothed slope (centred 30 min window) ────────────────────────
# Readings are 5 min apart, so ×60/5 converts a step difference to mmol/L/h.
df["slope_raw"] = (df["glucose"].diff(1) / 5) * 60
df["slope"]     = df["slope_raw"].rolling(3, center=True, min_periods=3).mean()

# ── Insulin (bolus + basal), both optional ────────────────────────
bolus_df = basal_df = None

b = read_export(BOLUS_CSV)
if b is not None:
    units_col = find_column(b.columns, "Insuline administrée")
    bolus_df = pd.DataFrame({
        "timestamp": pd.to_datetime(b["Horodatage"], format="%Y-%m-%d %H:%M", errors="coerce"),
        "units":     to_float(b[units_col]),
        "type":      b[find_column(b.columns, "Type")].astype(str),
    }).dropna(subset=["timestamp", "units"])
    bolus_df = bolus_df[bolus_df["units"] > 0].sort_values("timestamp").reset_index(drop=True)
    print(f"       bolus: {len(bolus_df)} injections")

ba = read_export(BASAL_CSV)
if ba is not None:
    basal_df = pd.DataFrame({
        "timestamp": pd.to_datetime(ba["Horodatage"], format="%Y-%m-%d %H:%M", errors="coerce"),
        "rate":      to_float(ba["Taux"]),
    }).dropna(subset=["timestamp", "rate"]).sort_values("timestamp").reset_index(drop=True)
    print(f"       basal: {len(basal_df)} steps")


# ── Episode detection (first point of each consecutive group) ─────
def first_point_of_episodes(mask, gap_min=20):
    """Return a boolean mask keeping only the first point of each episode.

    Readings closer together than `gap_min` minutes count as the same episode,
    so one excursion yields one star instead of a cloud of points.
    """
    idx = mask[mask].index
    if len(idx) == 0:
        return pd.Series(False, index=mask.index)
    ts = df.loc[idx, "timestamp"]
    gaps = ts.diff().dt.total_seconds().fillna(9999) / 60
    result = pd.Series(False, index=mask.index)
    result.loc[idx[gaps > gap_min]] = True
    return result


episodes_high = first_point_of_episodes(df["slope"] >= HIGH_THRESHOLD)
episodes_low  = first_point_of_episodes(df["slope"] <= LOW_THRESHOLD)

df["date"] = df["timestamp"].dt.date
dates = sorted(df["date"].unique())

# ── Layout ────────────────────────────────────────────────────────
COLS = 4
ROWS = int(np.ceil(len(dates) / COLS))

fig = plt.figure(figsize=(COLS * 5, ROWS * 5))
fig.patch.set_facecolor("white")
outer = fig.add_gridspec(ROWS, COLS, wspace=0.3, hspace=0.45)

for idx, day in enumerate(dates):
    row, col = divmod(idx, COLS)
    inner = outer[row, col].subgridspec(3, 1, height_ratios=[2, 1, 1], hspace=0.08)
    ax_g = fig.add_subplot(inner[0])                  # glucose
    ax_s = fig.add_subplot(inner[1], sharex=ax_g)     # slope
    ax_i = fig.add_subplot(inner[2], sharex=ax_g)     # insulin

    day_df = df[df["date"] == day]
    ts     = day_df["timestamp"]
    glc    = day_df["glucose"]
    slope  = day_df["slope"]
    high_d = episodes_high[day_df.index]
    low_d  = episodes_low[day_df.index]

    # ── Glucose ───────────────────────────────────────────────────
    ax_g.axhspan(2.0,  3.9,  alpha=0.18, color="#FFCDD2", zorder=0)   # low
    ax_g.axhspan(3.9,  10.0, alpha=0.18, color="#C8E6C9", zorder=0)   # target
    ax_g.axhspan(10.0, 23.0, alpha=0.18, color="#FFF9C4", zorder=0)   # high
    ax_g.axhline(3.9,  color="#EF5350", lw=0.8, ls="--", alpha=0.7)
    ax_g.axhline(10.0, color="#FFA726", lw=0.8, ls="--", alpha=0.7)
    ax_g.plot(ts, glc, color="#1976D2", lw=1.3, zorder=3)

    # Stars on the glucose trace
    pts_high = day_df[high_d.values]
    pts_low  = day_df[low_d.values]
    if len(pts_high):
        ax_g.scatter(pts_high["timestamp"], pts_high["glucose"],
                     marker="*", s=120, color="#F57C00",
                     zorder=6, linewidths=0.5, edgecolors="white")
    if len(pts_low):
        ax_g.scatter(pts_low["timestamp"], pts_low["glucose"],
                     marker="*", s=120, color="#1565C0",
                     zorder=6, linewidths=0.5, edgecolors="white")

    ax_g.set_ylim(2.0, 23.0)
    ax_g.yaxis.set_major_locator(ticker.MultipleLocator(4))
    ax_g.tick_params(axis="y", labelsize=6.5, labelleft=True)
    ax_g.tick_params(axis="x", labelbottom=False)
    ax_g.set_ylabel("mmol/L", fontsize=7)
    ax_g.grid(True, alpha=0.2)
    ax_g.set_facecolor("#FAFAFA")

    # Daily stats and episode counts
    n_high = high_d.sum()
    n_low  = low_d.sum()
    ax_g.text(0.02, 0.97,
              f"min {glc.min():.1f}  mean {glc.mean():.1f}  max {glc.max():.1f}"
              f"   |  * {n_high} up  {n_low} down",
              transform=ax_g.transAxes, fontsize=6, va="top",
              bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                        alpha=0.85, edgecolor="none"))

    ax_g.set_title(pd.Timestamp(day).strftime("%d %b %Y"),
                   fontsize=8, fontweight="bold", pad=3)

    # ── Slope ─────────────────────────────────────────────────────
    ax_s.axhline(0, color="#333", lw=0.9, zorder=2)
    ax_s.axhspan(-1, 1, alpha=0.12, color="#E0E0E0", zorder=0)   # near-flat band
    ax_s.axhline(HIGH_THRESHOLD, color="#F57C00", lw=0.9, ls=":", alpha=0.8)
    ax_s.axhline(LOW_THRESHOLD,  color="#1565C0", lw=0.9, ls=":", alpha=0.8)

    ax_s.fill_between(ts, 0, slope, where=(slope >= 0), interpolate=True,
                      color="#F57C00", alpha=0.60, zorder=1)
    ax_s.fill_between(ts, 0, slope, where=(slope < 0), interpolate=True,
                      color="#1976D2", alpha=0.60, zorder=1)
    ax_s.plot(ts, slope, color="#333", lw=0.7, zorder=3)

    # Stars on the slope trace
    if len(pts_high):
        ax_s.scatter(pts_high["timestamp"], pts_high["slope"],
                     marker="*", s=120, color="#F57C00",
                     zorder=6, linewidths=0.5, edgecolors="white")
    if len(pts_low):
        ax_s.scatter(pts_low["timestamp"], pts_low["slope"],
                     marker="*", s=120, color="#1565C0",
                     zorder=6, linewidths=0.5, edgecolors="white")

    # Symmetric limits, always wide enough to show the thresholds.
    slope_max = max(np.nanmax(np.abs(slope.values)), abs(HIGH_THRESHOLD) + 1)
    slope_lim = np.ceil(slope_max / 2) * 2
    ax_s.set_ylim(-slope_lim, slope_lim)
    ax_s.yaxis.set_major_locator(ticker.MultipleLocator(slope_lim / 2))
    ax_s.tick_params(axis="y", labelsize=6.5)
    ax_s.tick_params(axis="x", labelbottom=False)   # hour labels live on ax_i
    ax_s.set_ylabel("mmol/L/h", fontsize=7)
    ax_s.grid(True, alpha=0.2)
    ax_s.set_facecolor("#FAFAFA")

    # ── Insulin: basal step area + bolus stems ────────────────────
    # Bounds come from the day's glucose so the x axis stays in sync.
    t0, t1 = ts.min(), ts.max()

    # A day counts as covered only if it overlaps the pump export's own span.
    # Without this test the basal branch below would carry the last known rate
    # forward and draw a flat band across days the export never described.
    def covers(frame):
        return (frame is not None and len(frame)
                and frame["timestamp"].min() <= t1
                and frame["timestamp"].max() >= t0)

    basal_covered = covers(basal_df)
    bolus_covered = covers(bolus_df)

    if basal_covered:
        # Carry the rate in force at t0 forward, then follow the day's changes.
        earlier = basal_df[basal_df["timestamp"] <= t0]
        rate0   = earlier["rate"].iloc[-1] if len(earlier) else basal_df["rate"].iloc[0]
        within  = basal_df[(basal_df["timestamp"] > t0) & (basal_df["timestamp"] <= t1)]
        xs = [t0, *within["timestamp"].tolist(), t1]
        ys = [rate0, *within["rate"].tolist(),
              (within["rate"].iloc[-1] if len(within) else rate0)]
        ax_i.fill_between(xs, 0, ys, step="post", color="#7E57C2", alpha=0.35, zorder=1)
        ax_i.step(xs, ys, where="post", color="#5E35B1", lw=1.0, zorder=2)
        rate_max = max(np.nanmax(ys), 0.1)
    else:
        rate_max = 1.0

    ax_i.set_ylim(0, np.ceil(rate_max / 0.5) * 0.5 * 1.15 + 0.05)
    ax_i.yaxis.set_major_locator(ticker.MaxNLocator(3))
    ax_i.tick_params(axis="y", labelsize=6.5)
    ax_i.set_ylabel("U/h basal", fontsize=7, color="#5E35B1")
    ax_i.grid(True, alpha=0.2)
    ax_i.set_facecolor("#FAFAFA")

    # Bolus on a twin axis — units, scaled independently of the basal rate.
    ax_b = ax_i.twinx()
    if bolus_covered:
        bol = bolus_df[(bolus_df["timestamp"] >= t0) & (bolus_df["timestamp"] <= t1)]
        if len(bol):
            ax_b.vlines(bol["timestamp"], 0, bol["units"],
                        color="#00897B", lw=1.4, zorder=4)
            ax_b.scatter(bol["timestamp"], bol["units"], marker="o", s=18,
                         color="#00897B", zorder=5, edgecolors="white", linewidths=0.5)
            for t, u in zip(bol["timestamp"], bol["units"]):
                ax_b.annotate(f"{u:.1f}", (t, u), textcoords="offset points",
                              xytext=(0, 2), ha="center", fontsize=5, color="#00695C")
            ax_b.set_ylim(0, max(bol["units"].max() * 1.30, 1.0))
    ax_b.set_ylabel("U bolus", fontsize=7, color="#00897B")
    ax_b.tick_params(axis="y", labelsize=6.5, colors="#00897B")

    # Say so outright when the pump export does not reach this day, rather than
    # leaving an empty panel that reads like "no insulin was delivered".
    if not basal_covered and not bolus_covered:
        ax_i.text(0.5, 0.5, "no insulin data", transform=ax_i.transAxes,
                  ha="center", va="center", fontsize=7.5, style="italic",
                  color="#9E9E9E", zorder=7)

    # ── Hour axis, on the bottom panel only ───────────────────────
    ax_i.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax_i.xaxis.set_major_locator(mdates.HourLocator(interval=4))
    plt.setp(ax_i.xaxis.get_majorticklabels(), rotation=45, fontsize=6)

# Hide unused cells in the last row
for idx in range(len(dates), ROWS * COLS):
    row, col = divmod(idx, COLS)
    fig.add_subplot(outer[row, col]).set_visible(False)

fig.suptitle(
    f"Daily glucose, slope and insulin — {PATIENT_LABEL}\n"
    f"Orange star = slope >= +{HIGH_THRESHOLD}  |  blue star = slope <= {LOW_THRESHOLD} mmol/L/h"
    f"   |  purple area = basal (U/h)  |  green stems = bolus (U)",
    fontsize=12, fontweight="bold", y=1.005
)

# ── Save ──────────────────────────────────────────────────────────
# PDF first: vector output stays sharp at any zoom.
plt.savefig(f"{OUTPUT_BASE}.pdf", bbox_inches="tight", facecolor="white")
plt.savefig(f"{OUTPUT_BASE}.png", dpi=PNG_DPI, bbox_inches="tight", facecolor="white")

shown = OUTPUT_BASE.relative_to(PROJECT_DIR)
print(f"Saved: {shown}.pdf (vector) + {shown}.png ({PNG_DPI} dpi)"
      f"   ({len(dates)} days)")
plt.show()
