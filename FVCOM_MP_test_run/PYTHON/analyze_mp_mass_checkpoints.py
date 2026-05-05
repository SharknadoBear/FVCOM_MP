#!/usr/bin/env python3
"""Analyze in-run mass-conservation checkpoints from mod_plast PLAST_MASS_CHECK output.

Reads the CSV files written by Write_Mass_Checkpoint in mod_plast.F and produces:

  1. Grand-total timeseries at CP0 (timestep entry) and CP12 (timestep exit)
     for cases b and c, with a comparison to confirm per-step net change.
  2. Per-stage mass-change timeseries: delta(grand_total) between consecutive
     checkpoints.  Stages that must be perfectly conservative (advection, vdif,
     MPI exchange) are highlighted so any non-zero delta is immediately visible.
  3. Per-stage mass-change heatmap: rows = timestep, cols = stage delta.
     Designed to show the onset day and which stage produces the first anomaly.
  4. NaN diagnostic table: which CP stages have NaN values, and how many rows.
     NaN in water-column sums indicates the underlying conc arrays were not
     halo-consistent at that stage -- itself diagnostic information.
  5. Summary CSV for each case with the pivoted per-stage deltas.

CSV format written by mod_plast.F:
  iint, T_days, checkpoint, iplast,
  wc_a_kg, wc_d_kg, wc_tot_kg,
  bed_a_kg, bed_d_kg, bed_tot_kg, grand_total_kg

Note: T_days in the CSV is the raw FVCOM model time (Julian day, e.g. 58485.677).
The script converts this to comparison_day = T_days - T_days.min() so the x-axis
starts at 0 and matches the "comparison day" convention used in memo_03.

Usage:
  python analyze_mp_mass_checkpoints.py               # auto-detect b/c CSVs
  python analyze_mp_mass_checkpoints.py --cases b c   # explicit
  python analyze_mp_mass_checkpoints.py --csv-dir /path/to/input_dir
  python analyze_mp_mass_checkpoints.py --day-zoom 2.5 4.0

Outputs go to ./output/mass_conservation/checkpoints/.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
TEST_ROOT = SCRIPT_DIR.parent
DEFAULT_CSV_DIR = TEST_ROOT / "INPUT"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "output" / "mass_conservation" / "checkpoints"

# Canonical checkpoint order as written by the Fortran code
CP_ORDER = [
    "CP0_start",
    "CP1_after_deposition",
    "CP2_after_erosion",
    "CP3_after_update_bottom",
    "CP3p5_after_adv_a",        # NEW: after adv_scal(conc_a), before adv_scal(conc_d)
    "CP4_after_advection",
    "CP5_after_vdif",
    "CP6_after_kinetics",
    "CP7_after_upper_clamp",
    "CP8_after_OBC",
    "CP9_after_PTsource",
    "CP10_after_neg_clamp",
    "CP11_after_MPI",
    "CP12_end",
]

# Stage delta labels: delta_N = grand_total(CP_N) - grand_total(CP_{N-1})
# (first entry is CP0 absolute, not a delta)
DELTA_LABELS = [f"d{CP_ORDER[i+1].split('_',1)[1]}_from_{CP_ORDER[i].split('_',1)[1]}"
                for i in range(len(CP_ORDER) - 1)]

# Stages where grand total MUST be exactly conserved (no source/sink expected).
# Any nonzero delta here is a mass-conservation bug.
CONSERVATIVE_STAGES = {
    "CP3→CP3p5",  # adv_a only: must be zero when conc_a=0
    "CP3p5→CP4",  # adv_d only: expected OBC flux only (same as case a)
    "CP4→CP5",    # vertical diffusion
    "CP5→CP6",    # kinetic exchange (conc_a+conc_d must be constant)
    "CP10→CP11",  # MPI halo exchange
    "CP11→CP12",  # recombination
}

# Short labels for the stage-delta axis (space-saving)
SHORT_DELTA = [
    "dep",       # CP0→CP1  deposition: water↓ bed↑
    "ero",       # CP1→CP2  erosion:    water↑ bed↓
    "bot",       # CP2→CP3  bed bookkeeping
    "adv_a*",    # CP3→CP3p5 advection of conc_a only ← must be zero if conc_a=0
    "adv_d*",    # CP3p5→CP4 advection of conc_d     ← must equal OBC flux only
    "vdif*",     # CP4→CP5  vertical diff  ← must be zero
    "kin*",      # CP5→CP6  kinetics       ← must be zero (a+d conserved)
    "clamp",     # CP6→CP7  upper clamp
    "OBC",       # CP7→CP8  OBC nudging
    "src",       # CP8→CP9  point source
    "neg",       # CP9→CP10 neg clamp
    "MPI*",      # CP10→CP11 MPI exchange  ← must be zero
    "recomb*",   # CP11→CP12 recombination ← must be zero
]

# For PLAST_FLOC=F (case a): CP3→CP3p5 captures the full horizontal advection
# AND vertical diffusion in one checkpoint window (adv_scal then vdif_scal before
# the next Write_MC_1Field call).  CP3p5→CP11 are reported from the same 'cnew'
# array, so they are identically zero and labelled accordingly.
SHORT_DELTA_A = [
    "dep",
    "ero",
    "bot",
    "transport*",      # CP3→CP3p5: adv_scal + vdif_scal combined
    "adv_d*\n(=0,NF)", # non-floc: same state as CP3p5
    "vdif*\n(=0,NF)",
    "kin*\n(=0,NF)",
    "clamp\n(=0,NF)",
    "OBC",
    "src",
    "neg",
    "MPI*",
    "recomb*\n(=0,NF)", # non-floc: same state as CP11
]

CASE_COLORS = {"a": "#2ca02c", "b": "#1f77b4", "c": "#ff7f0e"}


def _taxis(df: pd.DataFrame) -> tuple[pd.Series, str]:
    """Return (time_series, xlabel) using comparison_day when it has real spread,
    otherwise fall back to iint (for CSVs produced before the T_model/86400 fix)."""
    cd = df["comparison_day"]
    if cd.max() - cd.min() > 0.05:
        return cd, "Comparison day (days since run start)"
    # Fall back to iint index when all T_days are essentially the same value
    return pd.Series(df.index, index=df.index, name="iint"), "Timestep index (iint)"


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def find_csv(csv_dir: Path, case: str) -> Path:
    """Locate the mass-check CSV for a given case label."""
    patterns = [
        f"*{case}_mp_mass_check.csv",
        f"*{case.upper()}_mp_mass_check.csv",
    ]
    for pat in patterns:
        hits = sorted(csv_dir.glob(pat))
        if hits:
            return hits[0]
    raise FileNotFoundError(
        f"No mass-check CSV found for case '{case}' in {csv_dir}.\n"
        f"Expected a file matching *{case}_mp_mass_check.csv"
    )


def load_csv(path: Path) -> pd.DataFrame:
    """Read a mass-check CSV, strip whitespace from column names and checkpoint strings.
    CFLX_TRACE diagnostic lines written by the Fortran cflx-trace prints are
    excluded automatically (they don't match the 11-column CSV header)."""
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    # Drop any rows that don't have numeric iint (CFLX_TRACE lines land here)
    df = df[pd.to_numeric(df["iint"], errors="coerce").notna()].copy()
    df["iint"] = df["iint"].astype(int)
    df["checkpoint"] = df["checkpoint"].str.strip()
    # iplast is written with leading spaces by Fortran; cast to int after stripping
    df["iplast"] = pd.to_numeric(df["iplast"], errors="coerce").astype("Int64")
    # All numeric columns may be space-padded strings — coerce them all
    for col in ["T_days", "wc_a_kg", "wc_d_kg", "wc_tot_kg",
                "bed_a_kg", "bed_d_kg", "bed_tot_kg", "grand_total_kg"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------

def pivot_to_wide(df: pd.DataFrame, iplast: int = 1) -> pd.DataFrame:
    """
    For a single plastic class, pivot the long CSV so each row is one
    reporting timestep and each column is the grand_total at one checkpoint.

    Returns a DataFrame indexed by iint with columns named by CP_ORDER,
    plus a T_days column.
    """
    sub = df[df["iplast"] == iplast].copy()

    # Keep only checkpoints we know about (ignore any extras)
    sub = sub[sub["checkpoint"].isin(CP_ORDER)]

    # Pivot: index=iint, columns=checkpoint, values=grand_total_kg
    wide = sub.pivot_table(
        index="iint", columns="checkpoint", values="grand_total_kg", aggfunc="first"
    )

    # Reorder columns in canonical CP order (only those present)
    present = [cp for cp in CP_ORDER if cp in wide.columns]
    wide = wide[present]

    # Attach T_days (take the first T_days for each iint, they should all match)
    t_days = sub.groupby("iint")["T_days"].first()
    wide.insert(0, "T_days", t_days)

    wide.sort_index(inplace=True)

    # Compute comparison_day = days since run start (matches memo_03 convention)
    t_min = wide["T_days"].min()
    wide.insert(1, "comparison_day", wide["T_days"] - t_min)

    return wide


def compute_stage_deltas(wide: pd.DataFrame) -> pd.DataFrame:
    """
    Compute the per-step mass change for each consecutive CP pair.
    Columns: T_days, CP0_abs, then one delta column per consecutive stage.
    Positive = mass gained; negative = mass lost.
    """
    result = pd.DataFrame(index=wide.index)
    result["T_days"] = wide["T_days"]
    result["comparison_day"] = wide["comparison_day"]
    result["CP0_grand_total_kg"] = wide.get("CP0_start", np.nan)

    cp_present = [cp for cp in CP_ORDER if cp in wide.columns]
    for i in range(len(cp_present) - 1):
        cp_a = cp_present[i]
        cp_b = cp_present[i + 1]
        label = f"{cp_a}→{cp_b}"
        result[label] = wide[cp_b] - wide[cp_a]

    return result


def nan_diagnostic(df: pd.DataFrame) -> pd.DataFrame:
    """Count NaN rows per checkpoint in the raw long CSV."""
    rows = []
    for cp in CP_ORDER:
        sub = df[df["checkpoint"] == cp]
        n_total = len(sub)
        n_nan_wc = sub["wc_tot_kg"].isna().sum()
        n_nan_bed = sub["bed_tot_kg"].isna().sum()
        n_nan_grand = sub["grand_total_kg"].isna().sum()
        rows.append({
            "checkpoint": cp,
            "n_rows": n_total,
            "nan_wc_tot": n_nan_wc,
            "nan_bed_tot": n_nan_bed,
            "nan_grand_total": n_nan_grand,
            "frac_nan_grand": n_nan_grand / n_total if n_total else np.nan,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _save(fig: plt.Figure, path: Path, dpi: int = 150) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved -> {path.name}")


def plot_grand_total_timeseries(
    wide_dict: dict[str, pd.DataFrame],
    output_dir: Path,
    day_zoom: tuple[float, float] | None = None,
) -> None:
    """
    Panel 1: CP0 (start-of-step) and CP12 (end-of-step) grand totals for each case.
    Panel 2: CP12 - CP0 per step (should be near-zero if conservative).
    """
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    ax_tot, ax_diff = axes

    # Determine x-axis from first available wide DataFrame
    first_wide = next(iter(wide_dict.values()))
    _, xlabel = _taxis(first_wide)

    for case, wide in wide_dict.items():
        col = CASE_COLORS.get(case, "gray")
        t, _ = _taxis(wide)
        cp0 = wide.get("CP0_start")
        cp12 = wide.get("CP12_end")
        if cp0 is not None:
            ax_tot.plot(t, cp0, color=col, lw=1.4, label=f"Case {case} CP0 (start)")
        if cp12 is not None:
            ax_tot.plot(t, cp12, color=col, lw=1.4, ls="--", label=f"Case {case} CP12 (end)")
        if cp0 is not None and cp12 is not None:
            per_step_net = cp12 - cp0
            ax_diff.plot(t, per_step_net, color=col, lw=1.0, label=f"Case {case}")

    ax_tot.set_ylabel("Grand total mass [kg]")
    ax_tot.set_title("Domain-integrated MP mass at start (CP0) and end (CP12) of each timestep")
    ax_tot.legend(fontsize=8)
    ax_tot.grid(True, lw=0.4, alpha=0.5)

    ax_diff.axhline(0, color="k", lw=0.8, ls=":")
    ax_diff.set_ylabel("CP12 − CP0 per step [kg]")
    ax_diff.set_xlabel(xlabel)
    ax_diff.set_title("Per-step net mass change (should be ≈ 0 for conservative path)")
    ax_diff.legend(fontsize=8)
    ax_diff.grid(True, lw=0.4, alpha=0.5)

    if day_zoom:
        ax_tot.set_xlim(day_zoom)

    fig.tight_layout()
    _save(fig, output_dir / "cp_grand_total_timeseries.png")


def plot_stage_delta_timeseries(
    deltas_dict: dict[str, pd.DataFrame],
    output_dir: Path,
    day_zoom: tuple[float, float] | None = None,
) -> None:
    """
    One panel per consecutive CP stage (CP0→CP1 … CP11→CP12).
    Background color: red if the stage must be conservative but shows nonzero.
    """
    # Gather all stage columns present in any case
    all_stage_cols: list[str] = []
    for deltas in deltas_dict.values():
        for c in deltas.columns:
            if "→" in c and c not in all_stage_cols:
                all_stage_cols.append(c)

    n_stages = len(all_stage_cols)
    if n_stages == 0:
        return

    fig, axes = plt.subplots(n_stages, 1, figsize=(12, 2.0 * n_stages), sharex=True)
    if n_stages == 1:
        axes = [axes]

    for ax, col in zip(axes, all_stage_cols):
        # Extract "CPa→CPb" short name
        parts = col.split("→")
        cp_a_num = parts[0].split("_")[0].replace("CP", "")
        cp_b_num = parts[1].split("_")[0].replace("CP", "")
        stage_key = f"CP{cp_a_num}→CP{cp_b_num}"
        is_conservative = stage_key in CONSERVATIVE_STAGES

        if is_conservative:
            ax.set_facecolor("#fff0f0")

        any_plotted = False
        for case, deltas in deltas_dict.items():
            if col not in deltas.columns:
                continue
            c_color = CASE_COLORS.get(case, "gray")
            t, xlabel = _taxis(deltas)
            vals = deltas[col]
            ax.plot(t, vals, color=c_color, lw=0.9, label=f"Case {case}")
            any_plotted = True

        ax.axhline(0, color="k", lw=0.6, ls=":")

        # Short label
        idx = next((i for i, cp_pair in enumerate(
            [f"{CP_ORDER[i]}→{CP_ORDER[i+1]}" for i in range(len(CP_ORDER)-1)]
        ) if cp_pair == col), None)
        short = SHORT_DELTA[idx] if idx is not None and idx < len(SHORT_DELTA) else col
        suffix = "  ← MUST=0" if is_conservative else ""
        ax.set_ylabel(f"Δ[kg]\n{short}{suffix}", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, lw=0.3, alpha=0.4)
        if any_plotted:
            ax.legend(fontsize=7, loc="upper left")

    axes[-1].set_xlabel(xlabel)
    if day_zoom:
        axes[0].set_xlim(day_zoom)

    fig.suptitle("Per-stage mass change at each reporting step\n"
                 "(pink background = stage must be perfectly conservative)", y=1.002)
    fig.tight_layout()
    _save(fig, output_dir / "cp_stage_delta_timeseries.png")


def plot_stage_delta_heatmap(
    deltas_dict: dict[str, pd.DataFrame],
    output_dir: Path,
    day_zoom: tuple[float, float] | None = None,
    per_case_short_labels: dict[str, list[str]] | None = None,
) -> None:
    """
    2-panel heatmap (one per case): rows=timesteps, cols=CP stage deltas.
    Colormap: diverging around 0; colour saturates at ±1% of final total.
    """
    all_stage_cols = []
    for deltas in deltas_dict.values():
        for c in deltas.columns:
            if "→" in c and c not in all_stage_cols:
                all_stage_cols.append(c)
    if not all_stage_cols:
        return

    short_cols = []
    for col in all_stage_cols:
        idx = next((i for i, cp_pair in enumerate(
            [f"{CP_ORDER[i]}→{CP_ORDER[i+1]}" for i in range(len(CP_ORDER)-1)]
        ) if cp_pair == col), None)
        short_cols.append(SHORT_DELTA[idx] if idx is not None and idx < len(SHORT_DELTA) else col)

    # Font scale: 1.5× the current default for the heatmap figure
    _base_fs = plt.rcParams["font.size"]
    _fs = _base_fs * 1.5        # tick/label size  (10 → 15)
    _fs_title = _fs * 1.2       # panel titles     (15 → 18)
    _fs_suptitle = _fs * 1.3    # figure suptitle  (15 → 19.5)

    n_cases = len(deltas_dict)
    # Wider figure to accommodate larger fonts
    fig, axes = plt.subplots(1, n_cases, figsize=(7 * n_cases + 1, 11), sharey=False)
    if n_cases == 1:
        axes = [axes]

    for ax, (case, deltas) in zip(axes, deltas_dict.items()):
        t_ser, xlabel = _taxis(deltas)
        t = t_ser.values
        if day_zoom:
            mask = (t >= day_zoom[0]) & (t <= day_zoom[1])
            deltas = deltas[mask]
            t = t[mask]

        mat = deltas[[c for c in all_stage_cols if c in deltas.columns]].values
        # Determine colour scale: symmetric, saturate at 1% of max final total
        final_total = deltas["CP0_grand_total_kg"].dropna().max()
        vmax = max(abs(np.nanmax(mat)), abs(np.nanmin(mat)))
        vmax = min(vmax, 0.01 * final_total) if final_total and final_total > 0 else vmax
        vmax = vmax if vmax > 0 else 1e-6

        im = ax.imshow(
            mat.T,
            aspect="auto",
            cmap="RdBu_r",
            norm=mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax),
            origin="lower",
            extent=[t[0], t[-1], -0.5, mat.shape[1] - 0.5],
        )
        ax.set_yticks(range(len(short_cols)))
        _ylabels = (per_case_short_labels or {}).get(case, short_cols)
        ax.set_yticklabels(_ylabels, fontsize=_fs)
        ax.tick_params(axis="x", labelsize=_fs)
        ax.set_xlabel(xlabel, fontsize=_fs)
        ax.set_title(f"Case {case}: stage mass change [kg]", fontsize=_fs_title)

        # Mark conservative stages with dashed horizontal lines
        for j, col in enumerate(all_stage_cols):
            parts = col.split("→")
            cp_a_num = parts[0].split("_")[0].replace("CP", "")
            cp_b_num = parts[1].split("_")[0].replace("CP", "")
            if f"CP{cp_a_num}→CP{cp_b_num}" in CONSERVATIVE_STAGES:
                ax.axhline(j, color="k", lw=1.0, ls="--", alpha=0.6)

        cb = plt.colorbar(im, ax=ax, label="\u0394 mass [kg]", shrink=0.6)
        cb.ax.tick_params(labelsize=_fs)
        cb.set_label("\u0394 mass [kg]", fontsize=_fs)

    fig.suptitle("Stage mass-change heatmap (diverging: red = gain, blue = loss)\n"
                 "Dashed lines = stages that must be exactly conservative",
                 y=1.01, fontsize=_fs_suptitle)
    fig.tight_layout()
    _save(fig, output_dir / "cp_stage_heatmap.png")


def plot_day3_zoom(
    wide_dict: dict[str, pd.DataFrame],
    deltas_dict: dict[str, pd.DataFrame],
    output_dir: Path,
    onset_day: float = 2.957,
    window: float = 1.5,
) -> None:
    """
    Zoom into the onset-divergence window showing CP0 totals + conservative-stage
    deltas to pinpoint the first non-zero anomaly.
    """
    d_lo = max(0.0, onset_day - window / 2)
    d_hi = onset_day + window

    conservative_cols = []
    for deltas in deltas_dict.values():
        for col in deltas.columns:
            if "→" not in col:
                continue
            parts = col.split("→")
            cp_a_num = parts[0].split("_")[0].replace("CP", "")
            cp_b_num = parts[1].split("_")[0].replace("CP", "")
            stage_key = f"CP{cp_a_num}→CP{cp_b_num}"
            if stage_key in CONSERVATIVE_STAGES and col not in conservative_cols:
                conservative_cols.append(col)

    n_rows = 1 + len(conservative_cols)
    fig, axes = plt.subplots(n_rows, 1, figsize=(11, 2.5 * n_rows), sharex=True)
    if n_rows == 1:
        axes = [axes]

    # Determine x-axis type from first available wide DataFrame
    first_wide = next(iter(wide_dict.values()))
    _, xlabel = _taxis(first_wide)
    use_iint = (xlabel.startswith("Timestep"))

    ax_tot = axes[0]
    for case, wide in wide_dict.items():
        t_ser, _ = _taxis(wide)
        if use_iint:
            mask = pd.Series(True, index=wide.index)  # no zoom when using iint
        else:
            mask = (t_ser >= d_lo) & (t_ser <= d_hi)
        c_color = CASE_COLORS.get(case, "gray")
        t = t_ser[mask]
        if "CP0_start" in wide.columns:
            ax_tot.plot(t, wide.loc[mask, "CP0_start"], color=c_color, lw=1.4,
                        label=f"Case {case} CP0")
    if not use_iint:
        ax_tot.axvline(onset_day, color="red", lw=1.0, ls=":", label=f"onset day {onset_day}")
    ax_tot.set_ylabel("Grand total [kg]")
    ax_tot.set_title(f"Onset-window zoom: day {d_lo:.1f}–{d_hi:.1f}")
    ax_tot.legend(fontsize=8)
    ax_tot.grid(True, lw=0.4, alpha=0.5)

    for ax, col in zip(axes[1:], conservative_cols):
        parts = col.split("→")
        cp_a_num = parts[0].split("_")[0].replace("CP", "")
        cp_b_num = parts[1].split("_")[0].replace("CP", "")
        idx = next((i for i, cp_pair in enumerate(
            [f"{CP_ORDER[i]}→{CP_ORDER[i+1]}" for i in range(len(CP_ORDER)-1)]
        ) if cp_pair == col), None)
        short = SHORT_DELTA[idx] if idx is not None and idx < len(SHORT_DELTA) else col

        ax.set_facecolor("#fff0f0")
        for case, deltas in deltas_dict.items():
            if col not in deltas.columns:
                continue
            t_ser, _ = _taxis(deltas)
            if use_iint:
                mask = pd.Series(True, index=deltas.index)
            else:
                mask = (t_ser >= d_lo) & (t_ser <= d_hi)
            c_color = CASE_COLORS.get(case, "gray")
            ax.plot(t_ser[mask], deltas.loc[mask, col],
                    color=c_color, lw=0.9, label=f"Case {case}")
        ax.axhline(0, color="k", lw=0.6, ls=":")
        if not use_iint:
            ax.axvline(onset_day, color="red", lw=1.0, ls=":")
        ax.set_ylabel(f"Δ[kg] {short}\n(MUST=0)", fontsize=8)
        ax.legend(fontsize=7)
        ax.grid(True, lw=0.3, alpha=0.4)

    axes[-1].set_xlabel(xlabel)
    fig.tight_layout()
    _save(fig, output_dir / "cp_day3_zoom.png")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--cases", nargs="+", default=["b", "c"], metavar="CASE",
        help="Case labels to analyse (default: b c).",
    )
    parser.add_argument(
        "--csv-dir", type=Path, default=DEFAULT_CSV_DIR,
        help=f"Directory containing *_mp_mass_check.csv files (default: {DEFAULT_CSV_DIR}).",
    )
    parser.add_argument(
        "--iplast", type=int, default=1,
        help="Plastic class index to analyse (default: 1).",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for output figures and CSVs (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--day-zoom", type=float, nargs=2, default=None,
        metavar=("DAY_START", "DAY_END"),
        help="Restrict x-axis of timeseries plots to this day range.",
    )
    parser.add_argument(
        "--onset-day", type=float, default=2.957,
        help="Model day of first observed a-b divergence (default: 2.957 from memo_03).",
    )
    parser.add_argument(
        "--max-days", type=float, default=None, metavar="DAYS",
        help="Clip all cases to this many comparison days before plotting (e.g. 10.0).",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir: Path = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== FVCOM-MP Mass Checkpoint Analysis ===")
    print(f"CSV directory : {args.csv_dir}")
    print(f"Cases         : {args.cases}")
    print(f"Plastic class : {args.iplast}")
    print(f"Output dir    : {output_dir}\n")

    # ---- Load data ----
    raw_dict: dict[str, pd.DataFrame] = {}
    wide_dict: dict[str, pd.DataFrame] = {}
    deltas_dict: dict[str, pd.DataFrame] = {}

    for case in args.cases:
        try:
            csv_path = find_csv(args.csv_dir, case)
        except FileNotFoundError as exc:
            print(f"  WARNING: {exc}")
            continue

        print(f"  Loading case {case}: {csv_path.name}  ({csv_path.stat().st_size/1e6:.1f} MB)")
        raw = load_csv(csv_path)
        raw_dict[case] = raw

        wide = pivot_to_wide(raw, iplast=args.iplast)
        wide_dict[case] = wide

        deltas = compute_stage_deltas(wide)
        deltas_dict[case] = deltas

        # ---- NaN diagnostic ----
        nan_diag = nan_diagnostic(raw[raw["iplast"] == args.iplast])
        nan_path = output_dir / f"nan_diagnostic_{case}.csv"
        nan_diag.to_csv(nan_path, index=False)
        print(f"\n  Case {case}: NaN diagnostic")
        print(nan_diag[nan_diag["nan_grand_total"] > 0].to_string(index=False)
              or "  (no NaN rows in grand_total)")

        # ---- Summary statistics ----
        if "CP0_start" in wide.columns and "CP12_end" in wide.columns:
            cp0_vals = wide["CP0_start"].dropna()
            cp12_vals = wide["CP12_end"].dropna()
            if len(cp0_vals):
                print(f"\n  Case {case} summary:")
                print(f"    Reporting steps    : {len(wide)}")
                print(f"    Comparison day range: {wide['comparison_day'].min():.3f} – "
                      f"{wide['comparison_day'].max():.3f}")
                print(f"    Julian day range   : {wide['T_days'].min():.4f} – "
                      f"{wide['T_days'].max():.4f}")
                print(f"    CP0 final total    : {cp0_vals.iloc[-1]:.6g} kg")
                print(f"    CP12 final total   : {cp12_vals.iloc[-1]:.6g} kg")
                net = cp12_vals - cp0_vals
                n_nonzero = (net.abs() > 1e-10).sum()
                print(f"    Steps with CP12!=CP0: {n_nonzero} / {len(net)}")

        # ---- Per-stage summary ----
        stage_summary_path = output_dir / f"stage_deltas_{case}.csv"
        deltas.to_csv(stage_summary_path, index=True)
        print(f"\n  Stage delta CSV -> {stage_summary_path.name}")

    if not wide_dict:
        print("\nNo data loaded. Check --csv-dir path and case labels.")
        return 1

    # ---- Clip to max_days if requested ----
    if args.max_days is not None:
        for case in list(wide_dict.keys()):
            mask = wide_dict[case]["comparison_day"] <= args.max_days
            wide_dict[case] = wide_dict[case][mask]
            deltas_dict[case] = deltas_dict[case][mask]
        print(f"  (Data clipped to first {args.max_days} comparison days)")

    # ---- Plots ----
    print("\nGenerating plots ...")
    day_zoom = tuple(args.day_zoom) if args.day_zoom else None

    plot_grand_total_timeseries(wide_dict, output_dir, day_zoom=day_zoom)
    plot_stage_delta_timeseries(deltas_dict, output_dir, day_zoom=day_zoom)
    # Build per-case y-label overrides: case a uses non-floc labels
    _pcsl = {}
    if "a" in deltas_dict:
        _pcsl["a"] = SHORT_DELTA_A
    plot_stage_delta_heatmap(deltas_dict, output_dir, day_zoom=day_zoom,
                             per_case_short_labels=_pcsl if _pcsl else None)
    plot_day3_zoom(wide_dict, deltas_dict, output_dir, onset_day=args.onset_day)

    print(f"\nDone. All outputs in:\n  {output_dir}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
