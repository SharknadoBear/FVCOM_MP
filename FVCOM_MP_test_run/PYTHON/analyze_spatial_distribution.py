#!/usr/bin/env python3
"""Spatial distribution checking for FVCOM-MP cases a, b, and c.

Reads the spatial_snapshots_*.mat files produced by the MATLAB extractors
and generates plan-view maps, pointwise-max maps, cross-case difference maps,
and station time-series plots.

MAT file locations (relative to FVCOM_MP_test_run/):
    OUTPUT_a/spatial_snapshots_a.mat
    OUTPUT_b/spatial_snapshots_b.mat
    OUTPUT_c/spatial_snapshots_c.mat

Outputs written to:
    PYTHON/output/spatial_distribution/

Usage
-----
    cd FVCOM_source_repo_github/FVCOM_MP_test_run/PYTHON
    python analyze_spatial_distribution.py
    python analyze_spatial_distribution.py --cases a b c
    python analyze_spatial_distribution.py --no-diff-maps
    python analyze_spatial_distribution.py --show
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

try:
    import h5py
except ImportError as exc:
    raise SystemExit(
        "h5py is required. Install with: pip install -r requirements.txt"
    ) from exc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
TEST_ROOT   = SCRIPT_DIR.parent
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "output" / "spatial_distribution"

# ── Station names shared with MATLAB ──────────────────────────────────────
STATION_NAMES = [
    "Chesapeake Canal",
    "Main stream 1",
    "Main stream 2",
    "Main stream 3",
    "Main stream 4",
    "Main stream 5",
    "Upstream Delaware",
    "Upstream Schuylkill",
]

# ── Colormap configuration per variable group ──────────────────────────────
CMAPS = {
    "mp1":      "plasma",
    "bot_mass": "YlOrRd",
    "sed":      "Blues",
    "diff":     "RdBu_r",
}

# ── Domain extent for zoomed maps ─────────────────────────────────────────
# Lon stored in 0–360° convention (284–286° range for Delaware Bay).
EXTENT_LAT_MIN = 39.4    # southern boundary (°N)
EXTENT_LON_MAX = 285.5   # eastern boundary (°, 0–360 convention = −74.5° W)

# Floor concentration for LogNorm vmin (avoids log(0))
LOG_VMIN_FLOOR = 1.0e-9  # kg/m³ or kg/m² — below any physically meaningful value


# ============================================================================
#  Entry point
# ============================================================================

def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.show:
        matplotlib.use("TkAgg")

    output_dir: Path = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load all requested cases ──────────────────────────────────────────
    cases: dict[str, dict] = {}
    for case_id in args.cases:
        mat_path = resolve_mat(args.root, case_id)
        print(f"Loading case {case_id}: {mat_path}")
        cases[case_id] = load_snapshot_mat(mat_path)

    # ── Build a shared Triangulation from the first available case ─────────
    first = next(iter(cases.values()))
    triang = build_triangulation(first)
    print(f"Triangulation: {triang.x.size} nodes, {triang.triangles.shape[0]} elements")

    # ── Plots ─────────────────────────────────────────────────────────────
    # 1. Per-case day-mean plan-view maps
    for case_id, data in cases.items():
        print(f"  Plan-view maps: case {case_id}")
        plot_plan_view_maps(case_id, data, triang, output_dir, args.show)

    # 2. Per-case pointwise-max maps
    for case_id, data in cases.items():
        print(f"  Pointwise-max maps: case {case_id}")
        plot_max_maps(case_id, data, triang, output_dir, args.show)

    # 3. Cross-case difference maps
    if not args.no_diff_maps and len(cases) > 1:
        case_list = list(cases.keys())
        diff_pairs = []
        if "b" in cases and "a" in cases:
            diff_pairs.append(("b", "a"))
        if "c" in cases and "a" in cases:
            diff_pairs.append(("c", "a"))
        if "c" in cases and "b" in cases:
            diff_pairs.append(("c", "b"))
        for case_hi, case_lo in diff_pairs:
            print(f"  Difference maps: case {case_hi} − case {case_lo}")
            plot_diff_maps(case_hi, case_lo, cases[case_hi], cases[case_lo],
                           triang, output_dir, args.show)

    # 4. Station time-series plots
    print("  Station time-series plots")
    plot_station_timeseries(cases, output_dir, args.show)

    print(f"\nAll outputs written to: {output_dir}")
    return 0


# ============================================================================
#  Argument parsing
# ============================================================================

def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Spatial distribution checking for FVCOM-MP cases."
    )
    p.add_argument(
        "--root",
        type=Path,
        default=TEST_ROOT,
        help="FVCOM_MP_test_run root containing OUTPUT_a/OUTPUT_b/OUTPUT_c.",
    )
    p.add_argument(
        "--cases",
        nargs="+",
        default=["a", "b", "c"],
        help="Case IDs to process (default: a b c).",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for PNG figures.",
    )
    p.add_argument(
        "--no-diff-maps",
        action="store_true",
        help="Skip cross-case difference maps.",
    )
    p.add_argument(
        "--show",
        action="store_true",
        help="Show plots interactively in addition to saving PNG files.",
    )
    return p.parse_args(argv)


def resolve_mat(root: Path, case_id: str) -> Path:
    path = (root / f"OUTPUT_{case_id}" / f"spatial_snapshots_{case_id}.mat").resolve()
    if not path.is_file():
        raise FileNotFoundError(f"MAT file not found for case {case_id}: {path}")
    return path


# ============================================================================
#  MAT file loading  (v7.3 HDF5)
# ============================================================================

def load_snapshot_mat(mat_path: Path) -> dict:
    """Return a flat dict of arrays from spatial_snapshots_*.mat."""
    data: dict[str, Any] = {}
    with h5py.File(mat_path, "r") as f:
        root_group = f["out"] if "out" in f else f

        # Grid
        grid = root_group["grid"]
        data["lon"]    = _read_1d(grid["lon"])
        data["lat"]    = _read_1d(grid["lat"])
        data["lonc"]   = _read_1d(grid["lonc"])
        data["latc"]   = _read_1d(grid["latc"])
        data["h"]      = _read_1d(grid["h"])
        data["art1"]   = _read_1d(grid["art1"])
        data["nv"]     = _read_2d(grid["nv"])        # (nele, 3) or (3, nele)
        data["dzfrac"] = _read_2d(grid["dzfrac"])    # (node, siglay) or transposed
        data["n_nodes"]  = int(np.asarray(grid["n_nodes"]).flat[0])
        data["n_siglay"] = int(np.asarray(grid["n_siglay"]).flat[0])

        # Snapshot metadata
        snaps = root_group["snaps"]
        data["snap_time"] = _read_1d(snaps["snap_time"])
        data["n_snaps"]   = int(np.asarray(snaps["n_snaps"]).flat[0])

        # Snapshot fields  (n_nodes x n_snaps)
        for key in ["mp1_surf", "mp1_bot", "mp1_davg",
                    "mp1_agg_surf", "mp1_agg_bot", "mp1_agg_davg",
                    "mp1_dis_surf", "mp1_dis_bot", "mp1_dis_davg",
                    "bot_mass", "bot_mass_agg", "bot_mass_dis",
                    "sed_surf", "sed_bot", "sed_davg", "sed_bedfrac"]:
            if key in snaps:
                data[key] = _read_2d(snaps[key])  # (n_nodes, n_snaps) after transpose

        # Station metadata
        st = root_group["stations"]
        data["station_time"]    = _read_1d(st["time"])
        data["station_node_idx"]     = _read_1d(st["node_idx"]).astype(int) - 1  # 0-based
        data["station_node_lat"]     = _read_1d(st["node_lat"])
        data["station_node_lon"]     = _read_1d(st["node_lon"])
        data["station_node_dist_km"] = _read_1d(st["node_dist_km"])

        # Station mp1 profiles  (n_stations x n_siglay x n_total)
        for key in ["mp1", "mp1_agg", "mp1_dis", "sed"]:
            src_key = key
            if src_key in st:
                data[f"station_{key}"] = _read_nd(st[src_key])

    return data


def _read_1d(ds: h5py.Dataset) -> np.ndarray:
    arr = np.asarray(ds, dtype=np.float64).ravel()
    return arr


def _read_2d(ds: h5py.Dataset) -> np.ndarray:
    """Read a 2D dataset; always return (first_dim x second_dim) in C order."""
    arr = np.asarray(ds, dtype=np.float64)
    # MATLAB HDF5: stored column-major → h5py reads as (ncols, nrows)
    # We want (nrows, ncols) i.e. the natural MATLAB array orientation.
    if arr.ndim == 2:
        return arr.T
    return arr


def _read_nd(ds: h5py.Dataset) -> np.ndarray:
    """Read N-D dataset; return transposed to match MATLAB row-major layout."""
    arr = np.asarray(ds, dtype=np.float64)
    if arr.ndim >= 2:
        return arr.T
    return arr


# ============================================================================
#  Triangulation
# ============================================================================

def build_triangulation(data: dict) -> mtri.Triangulation:
    lon = data["lon"]
    lat = data["lat"]
    nv  = data["nv"].astype(int)   # (nele, 3)
    # Convert to 0-based if necessary
    if nv.min() >= 1:
        nv = nv - 1
    return mtri.Triangulation(lon, lat, nv)


# ============================================================================
#  Plan-view maps  (day-mean over the last N-hour window)
# ============================================================================

def plot_plan_view_maps(case_id: str, data: dict, triang: mtri.Triangulation,
                        output_dir: Path, show: bool) -> None:
    """Day-mean plan-view maps for all available fields for one case."""

    snap_time = data["snap_time"]
    lon = data["lon"]
    lat = data["lat"]
    plot_specs = _build_plot_specs(data)

    for spec in plot_specs:
        field_mean = np.nanmean(spec["field"], axis=1)   # (n_nodes,)
        t_label = "day-mean ({:.3f}\u2013{:.3f} days)".format(snap_time[0], snap_time[-1])

        fig, ax = plt.subplots(figsize=(10, 7))
        _tripcolor(ax, triang, field_mean, spec["cmap"], spec["vmin"], spec["vmax"],
                   log_scale=spec["log_scale"], label=spec["title"])
        _apply_map_extent(ax, lon, lat)
        _add_stations(ax, data)
        ax.set_title("Case {}  |  {}  |  {}".format(case_id, spec["title"], t_label), fontsize=10)
        ax.set_xlabel("Longitude (°, 0–360)")
        ax.set_ylabel("Latitude (°N)")

        fname = output_dir / "planview_dayavg_{}_case{}.png".format(spec["tag"], case_id)
        _savefig(fig, fname, show)


def plot_max_maps(case_id: str, data: dict, triang: mtri.Triangulation,
                  output_dir: Path, show: bool) -> None:
    """Pointwise-max plan-view maps (worst-case over last N hours)."""

    snap_time = data["snap_time"]
    lon = data["lon"]
    lat = data["lat"]
    plot_specs = _build_plot_specs(data)

    for spec in plot_specs:
        field_max = np.nanmax(spec["field"], axis=1)
        t_label = "max over {:.3f}\u2013{:.3f} days".format(snap_time[0], snap_time[-1])

        fig, ax = plt.subplots(figsize=(10, 7))
        _tripcolor(ax, triang, field_max, spec["cmap"], spec["vmin"], spec["vmax"],
                   log_scale=spec["log_scale"], label=spec["title"])
        _apply_map_extent(ax, lon, lat)
        _add_stations(ax, data)
        ax.set_title("Case {}  |  {}  |  {}".format(case_id, spec["title"], t_label), fontsize=10)
        ax.set_xlabel("Longitude (°, 0–360)")
        ax.set_ylabel("Latitude (°N)")

        fname = output_dir / "planview_max_{}_case{}.png".format(spec["tag"], case_id)
        _savefig(fig, fname, show)


def _build_plot_specs(data: dict) -> list[dict]:
    """Return list of {field, title, tag, cmap, vmin, vmax, log_scale} dicts."""
    specs = []

    def _add(key, title, cmap_key, use_log=True):
        if key not in data:
            return
        arr = data[key]
        if arr.ndim != 2:
            return
        pos = arr[np.isfinite(arr) & (arr > 0)]
        if use_log:
            vmin = max(float(np.percentile(pos, 1)) if pos.size else LOG_VMIN_FLOOR,
                       LOG_VMIN_FLOOR)
            vmax = max(float(np.nanpercentile(arr, 99)), vmin * 10)
        else:
            vmin = 0.0
            vmax = max(float(np.nanpercentile(arr, 99)), 1e-12)
        specs.append({"field": arr, "title": title, "tag": key,
                      "cmap": CMAPS[cmap_key], "vmin": vmin, "vmax": vmax,
                      "log_scale": use_log})

    _add("mp1_surf",     "mp1 surface [kg/m3]",            "mp1")
    _add("mp1_bot",      "mp1 bottom [kg/m3]",             "mp1")
    _add("mp1_davg",     "mp1 depth-avg [kg/m3]",          "mp1")
    _add("mp1_agg_surf", "mp1_agg surface [kg/m3]",        "mp1")
    _add("mp1_agg_bot",  "mp1_agg bottom [kg/m3]",         "mp1")
    _add("mp1_agg_davg", "mp1_agg depth-avg [kg/m3]",      "mp1")
    _add("mp1_dis_surf", "mp1_dis surface [kg/m3]",        "mp1")
    _add("mp1_dis_bot",  "mp1_dis bottom [kg/m3]",         "mp1")
    _add("mp1_dis_davg", "mp1_dis depth-avg [kg/m3]",      "mp1")
    _add("bot_mass",     "bot_mass_mp1 [kg/m2]",           "bot_mass")
    _add("bot_mass_agg", "bot_mass_agg_mp1 [kg/m2]",       "bot_mass")
    _add("bot_mass_dis", "bot_mass_dis_mp1 [kg/m2]",       "bot_mass")
    _add("sed_surf",     "coarse_sand_1 surface [g/L]",    "sed")
    _add("sed_bot",      "coarse_sand_1 bottom [g/L]",     "sed")
    _add("sed_davg",     "coarse_sand_1 depth-avg [g/L]",  "sed")
    _add("sed_bedfrac",  "sed bed fraction [-]",            "sed", use_log=False)

    return specs


# ============================================================================
#  Cross-case difference maps
# ============================================================================

def plot_diff_maps(case_hi: str, case_lo: str,
                   data_hi: dict, data_lo: dict,
                   triang: mtri.Triangulation,
                   output_dir: Path, show: bool) -> None:
    """Day-mean difference maps: case_hi minus case_lo (linear scale)."""

    diff_keys = [
        ("mp1_surf",  "mp1 surface",    "kg m\u207b\u00b3"),
        ("mp1_davg",  "mp1 depth-avg",  "kg m\u207b\u00b3"),
        ("bot_mass",  "bot_mass_mp1",   "kg m\u207b\u00b2"),
    ]

    snap_time = data_hi["snap_time"]
    t_label = "day-mean ({:.3f}\u2013{:.3f} days)".format(snap_time[0], snap_time[-1])
    lon = data_hi["lon"]
    lat = data_hi["lat"]

    for key, label, units in diff_keys:
        if key not in data_hi or key not in data_lo:
            continue

        mean_hi = np.nanmean(data_hi[key], axis=1)
        mean_lo = np.nanmean(data_lo[key], axis=1)
        diff    = mean_hi - mean_lo

        abs_max = np.nanpercentile(np.abs(diff), 99)
        abs_max = max(abs_max, 1e-12)

        fig, ax = plt.subplots(figsize=(10, 7))
        tcf = ax.tripcolor(
            triang, diff,
            cmap=CMAPS["diff"],
            vmin=-abs_max, vmax=abs_max,
            shading="gouraud",
        )
        plt.colorbar(tcf, ax=ax,
                     label="\u0394 {} [{}]".format(label, units),
                     shrink=0.8)
        _apply_map_extent(ax, lon, lat)
        _add_stations(ax, data_hi)
        ax.set_title(
            "Case {} \u2212 Case {}  |  {}  |  {}".format(case_hi, case_lo, label, t_label),
            fontsize=10,
        )
        ax.set_xlabel("Longitude (°, 0\u2013360)")
        ax.set_ylabel("Latitude (°N)")

        fname = output_dir / "diff_{}_case{}_minus_case{}.png".format(key, case_hi, case_lo)
        _savefig(fig, fname, show)


# ============================================================================
#  Station time series
# ============================================================================

def plot_station_timeseries(cases: dict[str, dict],
                             output_dir: Path, show: bool) -> None:
    """Depth-averaged station time series for mp1 and coarse_sand_1."""

    CASE_COLORS = {"a": "#1f77b4", "b": "#ff7f0e", "c": "#2ca02c"}
    n_stations = len(STATION_NAMES)

    # ── mp1 surface-layer time series ─────────────────────────────────────
    _plot_station_var(
        cases=cases,
        station_data_key="station_mp1",     # (n_stations, n_siglay, n_time)
        layer="surface",
        var_label="mp1 surface-layer [kg m⁻³]",
        tag="mp1_surf",
        case_colors=CASE_COLORS,
        output_dir=output_dir,
        show=show,
    )

    # ── mp1 depth-averaged ────────────────────────────────────────────────
    _plot_station_var(
        cases=cases,
        station_data_key="station_mp1",
        layer="davg",
        var_label="mp1 depth-averaged [kg m⁻³]",
        tag="mp1_davg",
        case_colors=CASE_COLORS,
        output_dir=output_dir,
        show=show,
    )

    # ── coarse_sand_1 surface ─────────────────────────────────────────────
    has_sed = any("station_sed" in d for d in cases.values())
    if has_sed:
        _plot_station_var(
            cases=cases,
            station_data_key="station_sed",
            layer="surface",
            var_label="coarse_sand_1 surface [g L⁻¹]",
            tag="sed_surf",
            case_colors=CASE_COLORS,
            output_dir=output_dir,
            show=show,
        )


def _plot_station_var(cases: dict[str, dict],
                      station_data_key: str,
                      layer: str,
                      var_label: str,
                      tag: str,
                      case_colors: dict,
                      output_dir: Path,
                      show: bool) -> None:
    """One figure with n_station subplots, all cases overlaid."""

    n_stations = len(STATION_NAMES)
    ncols = 2
    nrows = int(np.ceil(n_stations / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 3.5 * nrows), sharex=False)
    axes = axes.ravel()

    for is_idx in range(n_stations):
        ax = axes[is_idx]
        ax.set_title(STATION_NAMES[is_idx], fontsize=9)

        for case_id, data in cases.items():
            if station_data_key not in data:
                continue
            prof = data[station_data_key]   # (n_stations, n_siglay, n_time)
            time = data["station_time"]

            if is_idx >= prof.shape[0]:
                continue

            if layer == "surface":
                ts = prof[is_idx, 0, :]     # sigma layer k=0 (surface)
            elif layer == "bottom":
                ts = prof[is_idx, -1, :]    # sigma layer k=-1 (bottom)
            elif layer == "davg":
                # Weighted by dzfrac — use uniform weights here (no vol available)
                ts = np.nanmean(prof[is_idx, :, :], axis=0)
            else:
                ts = prof[is_idx, 0, :]

            # Convert FVCOM day to elapsed days from start
            t0 = time[0] if len(time) > 0 else 0.0
            t_elapsed = time - t0

            ax.plot(t_elapsed, ts, color=case_colors.get(case_id, "k"),
                    linewidth=1.2, label=f"case {case_id}")

        ax.set_xlabel("Elapsed time (days)", fontsize=8)
        ax.set_ylabel(var_label, fontsize=8)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, linestyle=":", alpha=0.5)

    # Hide unused subplots
    for idx in range(n_stations, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(f"Station time series — {var_label}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fname = output_dir / f"station_timeseries_{tag}.png"
    _savefig(fig, fname, show)


# ============================================================================
#  Plotting utilities
# ============================================================================

def _tripcolor(ax: plt.Axes, triang: mtri.Triangulation,
               values: np.ndarray, cmap: str,
               vmin: float, vmax: float,
               log_scale: bool = True,
               label: str = "") -> None:
    """Plot a scalar field on the unstructured mesh with optional LogNorm."""
    if log_scale:
        # Clamp values to [vmin, inf) so log(0) is avoided
        v = np.where(np.isfinite(values) & (values > 0), values, vmin)
        norm = LogNorm(vmin=vmin, vmax=vmax)
        tcf = ax.tripcolor(triang, v, cmap=cmap, norm=norm, shading="gouraud")
    else:
        tcf = ax.tripcolor(triang, values, cmap=cmap,
                           vmin=vmin, vmax=vmax, shading="gouraud")
    plt.colorbar(tcf, ax=ax, label=label, shrink=0.8, pad=0.02)


def _apply_map_extent(ax: plt.Axes, lon: np.ndarray, lat: np.ndarray) -> None:
    """Apply fixed zoom extent and equal-aspect ratio to a map axes."""
    lon_min = float(lon.min())
    lat_max = float(lat.max())
    ax.set_xlim(lon_min, EXTENT_LON_MAX)
    ax.set_ylim(EXTENT_LAT_MIN, lat_max)
    ax.set_aspect("equal", adjustable="box")


def _add_stations(ax: plt.Axes, data: dict) -> None:
    """Overlay station markers on a map axes."""
    node_lon = data.get("station_node_lon")
    node_lat = data.get("station_node_lat")
    if node_lon is None or node_lat is None:
        return
    ax.scatter(node_lon, node_lat, s=40, c="red", marker="^",
               zorder=5, linewidths=0.5, edgecolors="white",
               label="stations")
    for i, name in enumerate(STATION_NAMES):
        if i < len(node_lon):
            ax.annotate(
                name.split()[0],
                xy=(node_lon[i], node_lat[i]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=6,
                color="darkred",
            )


def _savefig(fig: plt.Figure, path: Path, show: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"    Saved: {path.name}")
    if show:
        plt.show()
    plt.close(fig)


# ============================================================================
#  Run
# ============================================================================

if __name__ == "__main__":
    sys.exit(main())
