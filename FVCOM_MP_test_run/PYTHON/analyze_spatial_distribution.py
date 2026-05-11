#!/usr/bin/env python3
"""Spatial distribution checking for FVCOM-MP cases a, b, and c.

Reads the spatial_snapshots_*.mat files produced by the MATLAB extractors
and generates plan-view maps, pointwise-max maps, cross-case difference maps,
station time-series plots, and a reproducible CSV summary.

MAT file locations (relative to FVCOM_MP_test_run/):
    Day 10: OUTPUT_a/spatial_snapshots_a.mat
            OUTPUT_b/spatial_snapshots_b.mat
            OUTPUT_c/spatial_snapshots_c.mat
    Day 30: OUTPUT_a/spatial_snapshots_a_day30.mat
            OUTPUT_b/spatial_snapshots_b_day30.mat
            OUTPUT_c/spatial_snapshots_c_day30.mat
    Day 90: OUTPUT_a/spatial_snapshots_a_day90.mat
            OUTPUT_b/spatial_snapshots_b_day90.mat
            OUTPUT_c/spatial_snapshots_c_day90.mat

Outputs written by default to:
    PYTHON/output/spatial_distribution/day10/
    PYTHON/output/spatial_distribution/day30/
    PYTHON/output/spatial_distribution/day90/
"""

from __future__ import annotations

import argparse
import csv
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
from matplotlib.colors import LogNorm
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
TEST_ROOT = SCRIPT_DIR.parent
DEFAULT_OUTPUT_BASE = SCRIPT_DIR / "output" / "spatial_distribution"
VALID_DAYS = (10, 30, 90)

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

CMAPS = {
    "mp1": "plasma",
    "bot_mass": "YlOrRd",
    "sed": "Blues",
    "diff": "RdBu_r",
}

# Lon is stored in 0-360 convention, around 284-286 deg for Delaware Bay.
EXTENT_LAT_MIN = 39.4
EXTENT_LON_MAX = 285.5

# Floor concentration for LogNorm vmin, avoiding log(0).
LOG_VMIN_FLOOR = 1.0e-9

PLOT_FIELDS = [
    ("mp1_surf", "mp1 surface [kg/m3]", "mp1", True),
    ("mp1_bot", "mp1 bottom [kg/m3]", "mp1", True),
    ("mp1_davg", "mp1 depth-avg [kg/m3]", "mp1", True),
    ("mp1_agg_surf", "mp1_agg surface [kg/m3]", "mp1", True),
    ("mp1_agg_bot", "mp1_agg bottom [kg/m3]", "mp1", True),
    ("mp1_agg_davg", "mp1_agg depth-avg [kg/m3]", "mp1", True),
    ("mp1_dis_surf", "mp1_dis surface [kg/m3]", "mp1", True),
    ("mp1_dis_bot", "mp1_dis bottom [kg/m3]", "mp1", True),
    ("mp1_dis_davg", "mp1_dis depth-avg [kg/m3]", "mp1", True),
    ("bot_mass", "bot_mass_mp1 [kg/m2]", "bot_mass", True),
    ("bot_mass_agg", "bot_mass_agg_mp1 [kg/m2]", "bot_mass", True),
    ("bot_mass_dis", "bot_mass_dis_mp1 [kg/m2]", "bot_mass", True),
    ("sed_surf", "coarse_sand_1 surface [g/L]", "sed", True),
    ("sed_bot", "coarse_sand_1 bottom [g/L]", "sed", True),
    ("sed_davg", "coarse_sand_1 depth-avg [g/L]", "sed", True),
    ("sed_bedfrac", "sed bed fraction [-]", "sed", False),
]

DIFF_FIELDS = [
    ("mp1_surf", "mp1 surface", "kg m-3"),
    ("mp1_davg", "mp1 depth-avg", "kg m-3"),
    ("mp1_agg_surf", "mp1_agg surface", "kg m-3"),
    ("mp1_agg_davg", "mp1_agg depth-avg", "kg m-3"),
    ("mp1_dis_surf", "mp1_dis surface", "kg m-3"),
    ("mp1_dis_davg", "mp1_dis depth-avg", "kg m-3"),
    ("bot_mass", "bot_mass_mp1", "kg m-2"),
    ("bot_mass_agg", "bot_mass_agg_mp1", "kg m-2"),
    ("bot_mass_dis", "bot_mass_dis_mp1", "kg m-2"),
]

SUMMARY_FIELDS = [
    "mp1_surf",
    "mp1_bot",
    "mp1_davg",
    "mp1_agg_davg",
    "mp1_dis_davg",
    "bot_mass",
    "bot_mass_agg",
    "bot_mass_dis",
    "sed_davg",
    "sed_bedfrac",
]


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    output_dir = resolve_output_dir(args.output_dir, args.day)
    output_dir.mkdir(parents=True, exist_ok=True)

    cases: dict[str, dict] = {}
    for case_id in args.cases:
        mat_path = resolve_mat(args.root, case_id, args.day)
        print(f"Loading case {case_id}: {mat_path}")
        data = load_snapshot_mat(mat_path)
        data["case_id"] = case_id
        data["day"] = args.day
        data["mat_path"] = str(mat_path)
        normalize_time_axes(data)
        cases[case_id] = data

    first = next(iter(cases.values()))
    triang = build_triangulation(first)
    print(f"Triangulation: {triang.x.size} nodes, {triang.triangles.shape[0]} elements")

    for case_id, data in cases.items():
        print(f"  Plan-view maps: case {case_id}")
        plot_plan_view_maps(case_id, data, triang, output_dir, args.show)

    for case_id, data in cases.items():
        print(f"  Pointwise-max maps: case {case_id}")
        plot_max_maps(case_id, data, triang, output_dir, args.show)

    if not args.no_diff_maps and len(cases) > 1:
        for case_hi, case_lo in diff_pairs(cases):
            print(f"  Difference maps: case {case_hi} - case {case_lo}")
            plot_diff_maps(
                case_hi,
                case_lo,
                cases[case_hi],
                cases[case_lo],
                triang,
                output_dir,
                args.show,
            )

    print("  Station time-series plots")
    plot_station_timeseries(cases, output_dir, args.show)

    summary_path = output_dir / "summary_stats.csv"
    write_summary_csv(cases, summary_path)
    print(f"  Summary CSV: {summary_path.name}")

    print(f"\nAll outputs written to: {output_dir}")
    return 0


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
        "--day",
        type=int,
        choices=VALID_DAYS,
        default=10,
        help="Simulation day snapshot set to process (default: 10).",
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
        default=None,
        help="Directory for PNG figures and summary CSV. Defaults to dayXX subfolder.",
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


def resolve_output_dir(output_dir_arg: Path | None, day: int) -> Path:
    if output_dir_arg is not None:
        return output_dir_arg.resolve()
    return (DEFAULT_OUTPUT_BASE / f"day{day}").resolve()


def resolve_mat(root: Path, case_id: str, day: int) -> Path:
    suffix = "" if day == 10 else f"_day{day}"
    path = (root / f"OUTPUT_{case_id}" / f"spatial_snapshots_{case_id}{suffix}.mat").resolve()
    if not path.is_file():
        raise FileNotFoundError(f"MAT file not found for case {case_id}, day {day}: {path}")
    return path


def load_snapshot_mat(mat_path: Path) -> dict:
    """Return a flat dict of arrays from spatial_snapshots_*.mat."""
    data: dict[str, Any] = {}
    with h5py.File(mat_path, "r") as f:
        root_group = f["out"] if "out" in f else f

        grid = root_group["grid"]
        data["lon"] = _read_1d(grid["lon"])
        data["lat"] = _read_1d(grid["lat"])
        data["lonc"] = _read_1d(grid["lonc"])
        data["latc"] = _read_1d(grid["latc"])
        data["h"] = _read_1d(grid["h"])
        data["art1"] = _read_1d(grid["art1"])
        data["nv"] = _read_2d(grid["nv"])
        data["dzfrac"] = _read_2d(grid["dzfrac"])
        data["n_nodes"] = int(np.asarray(grid["n_nodes"]).flat[0])
        data["n_siglay"] = int(np.asarray(grid["n_siglay"]).flat[0])

        snaps = root_group["snaps"]
        data["snap_time"] = _read_1d(snaps["snap_time"])
        data["n_snaps"] = int(np.asarray(snaps["n_snaps"]).flat[0])

        for key, _, _, _ in PLOT_FIELDS:
            if key in snaps:
                data[key] = _read_2d(snaps[key])

        st = root_group["stations"]
        data["station_time"] = _read_1d(st["time"])
        data["station_node_idx"] = _read_1d(st["node_idx"]).astype(int) - 1
        data["station_node_lat"] = _read_1d(st["node_lat"])
        data["station_node_lon"] = _read_1d(st["node_lon"])
        data["station_node_dist_km"] = _read_1d(st["node_dist_km"])

        for key in ["mp1", "mp1_agg", "mp1_dis", "sed"]:
            if key in st:
                data[f"station_{key}"] = _read_nd(st[key])

    return data


def _read_1d(ds: h5py.Dataset) -> np.ndarray:
    return np.asarray(ds, dtype=np.float64).ravel()


def _read_2d(ds: h5py.Dataset) -> np.ndarray:
    arr = np.asarray(ds, dtype=np.float64)
    if arr.ndim == 2:
        return arr.T
    return arr


def _read_nd(ds: h5py.Dataset) -> np.ndarray:
    arr = np.asarray(ds, dtype=np.float64)
    if arr.ndim >= 2:
        return arr.T
    return arr


def normalize_time_axes(data: dict) -> None:
    """Add simulation-day plotting axes without mutating the raw MAT times."""
    station_time = data.get("station_time")
    snap_time = data.get("snap_time")

    if station_time is None or len(station_time) == 0:
        station_ref = 0.0
    else:
        station_ref = float(station_time[0])

    if station_time is not None and np.nanmax(station_time) > 10000.0:
        data["station_time_plot"] = station_time - station_ref
        data["time_axis_label"] = "Simulation day"
    else:
        data["station_time_plot"] = station_time
        data["time_axis_label"] = "Simulation day"

    if snap_time is not None and np.nanmax(snap_time) > 10000.0:
        data["snap_time_plot"] = snap_time - station_ref
    else:
        data["snap_time_plot"] = snap_time


def build_triangulation(data: dict) -> mtri.Triangulation:
    lon = data["lon"]
    lat = data["lat"]
    nv = data["nv"].astype(int)
    if nv.min() >= 1:
        nv = nv - 1
    return mtri.Triangulation(lon, lat, nv)


def plot_plan_view_maps(
    case_id: str,
    data: dict,
    triang: mtri.Triangulation,
    output_dir: Path,
    show: bool,
) -> None:
    lon = data["lon"]
    lat = data["lat"]
    plot_specs = build_plot_specs(data)

    for spec in plot_specs:
        field_mean = np.nanmean(spec["field"], axis=1)
        t_label = time_window_label(data, prefix="mean")

        fig, ax = plt.subplots(figsize=(10, 7))
        tripcolor(
            ax,
            triang,
            field_mean,
            spec["cmap"],
            spec["vmin"],
            spec["vmax"],
            log_scale=spec["log_scale"],
            label=spec["title"],
        )
        apply_map_extent(ax, lon, lat)
        add_stations(ax, data)
        ax.set_title(f"Case {case_id} | {spec['title']} | {t_label}", fontsize=10)
        ax.set_xlabel("Longitude (deg, 0-360)")
        ax.set_ylabel("Latitude (deg N)")

        fname = output_dir / f"planview_dayavg_{spec['tag']}_case{case_id}.png"
        savefig(fig, fname, show)


def plot_max_maps(
    case_id: str,
    data: dict,
    triang: mtri.Triangulation,
    output_dir: Path,
    show: bool,
) -> None:
    lon = data["lon"]
    lat = data["lat"]
    plot_specs = build_plot_specs(data)

    for spec in plot_specs:
        field_max = np.nanmax(spec["field"], axis=1)
        t_label = time_window_label(data, prefix="max")

        fig, ax = plt.subplots(figsize=(10, 7))
        tripcolor(
            ax,
            triang,
            field_max,
            spec["cmap"],
            spec["vmin"],
            spec["vmax"],
            log_scale=spec["log_scale"],
            label=spec["title"],
        )
        apply_map_extent(ax, lon, lat)
        add_stations(ax, data)
        ax.set_title(f"Case {case_id} | {spec['title']} | {t_label}", fontsize=10)
        ax.set_xlabel("Longitude (deg, 0-360)")
        ax.set_ylabel("Latitude (deg N)")

        fname = output_dir / f"planview_max_{spec['tag']}_case{case_id}.png"
        savefig(fig, fname, show)


def build_plot_specs(data: dict) -> list[dict]:
    specs = []
    for key, title, cmap_key, use_log in PLOT_FIELDS:
        if key not in data:
            continue
        arr = data[key]
        if arr.ndim != 2:
            continue

        pos = arr[np.isfinite(arr) & (arr > 0)]
        if use_log:
            vmin = max(
                float(np.percentile(pos, 1)) if pos.size else LOG_VMIN_FLOOR,
                LOG_VMIN_FLOOR,
            )
            vmax = max(float(np.nanpercentile(arr, 99)), vmin * 10.0)
        else:
            vmin = 0.0
            vmax = max(float(np.nanpercentile(arr, 99)), 1.0e-12)

        specs.append(
            {
                "field": arr,
                "title": title,
                "tag": key,
                "cmap": CMAPS[cmap_key],
                "vmin": vmin,
                "vmax": vmax,
                "log_scale": use_log,
            }
        )
    return specs


def plot_diff_maps(
    case_hi: str,
    case_lo: str,
    data_hi: dict,
    data_lo: dict,
    triang: mtri.Triangulation,
    output_dir: Path,
    show: bool,
) -> None:
    lon = data_hi["lon"]
    lat = data_hi["lat"]
    t_label = time_window_label(data_hi, prefix="mean")

    for key, label, units in DIFF_FIELDS:
        if key not in data_hi or key not in data_lo:
            continue

        mean_hi = np.nanmean(data_hi[key], axis=1)
        mean_lo = np.nanmean(data_lo[key], axis=1)
        diff = mean_hi - mean_lo

        abs_max = float(np.nanpercentile(np.abs(diff), 99))
        abs_max = max(abs_max, 1.0e-12)

        fig, ax = plt.subplots(figsize=(10, 7))
        tcf = ax.tripcolor(
            triang,
            diff,
            cmap=CMAPS["diff"],
            vmin=-abs_max,
            vmax=abs_max,
            shading="gouraud",
        )
        plt.colorbar(tcf, ax=ax, label=f"Delta {label} [{units}]", shrink=0.8)
        apply_map_extent(ax, lon, lat)
        add_stations(ax, data_hi)
        ax.set_title(
            f"Case {case_hi} - Case {case_lo} | {label} | {t_label}",
            fontsize=10,
        )
        ax.set_xlabel("Longitude (deg, 0-360)")
        ax.set_ylabel("Latitude (deg N)")

        fname = output_dir / f"diff_{key}_case{case_hi}_minus_case{case_lo}.png"
        savefig(fig, fname, show)


def plot_station_timeseries(cases: dict[str, dict], output_dir: Path, show: bool) -> None:
    case_colors = {"a": "#1f77b4", "b": "#ff7f0e", "c": "#2ca02c"}

    station_specs = [
        ("station_mp1", "surface", "mp1 surface-layer [kg m-3]", "mp1_surf"),
        ("station_mp1", "davg", "mp1 depth-averaged [kg m-3]", "mp1_davg"),
        ("station_mp1_agg", "surface", "mp1_agg surface-layer [kg m-3]", "mp1_agg_surf"),
        ("station_mp1_agg", "davg", "mp1_agg depth-averaged [kg m-3]", "mp1_agg_davg"),
        ("station_mp1_dis", "surface", "mp1_dis surface-layer [kg m-3]", "mp1_dis_surf"),
        ("station_mp1_dis", "davg", "mp1_dis depth-averaged [kg m-3]", "mp1_dis_davg"),
        ("station_sed", "surface", "coarse_sand_1 surface [g L-1]", "sed_surf"),
    ]

    for station_data_key, layer, label, tag in station_specs:
        if any(station_data_key in d for d in cases.values()):
            plot_station_var(
                cases=cases,
                station_data_key=station_data_key,
                layer=layer,
                var_label=label,
                tag=tag,
                case_colors=case_colors,
                output_dir=output_dir,
                show=show,
            )


def plot_station_var(
    cases: dict[str, dict],
    station_data_key: str,
    layer: str,
    var_label: str,
    tag: str,
    case_colors: dict,
    output_dir: Path,
    show: bool,
) -> None:
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

            prof = data[station_data_key]
            time = data.get("station_time_plot", data["station_time"])

            if is_idx >= prof.shape[0]:
                continue

            if layer == "surface":
                ts = prof[is_idx, 0, :]
            elif layer == "bottom":
                ts = prof[is_idx, -1, :]
            elif layer == "davg":
                ts = np.nanmean(prof[is_idx, :, :], axis=0)
            else:
                ts = prof[is_idx, 0, :]

            ax.plot(
                time,
                ts,
                color=case_colors.get(case_id, "k"),
                linewidth=1.2,
                label=f"case {case_id}",
            )

        ax.set_xlabel("Simulation day", fontsize=8)
        ax.set_ylabel(var_label, fontsize=8)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, linestyle=":", alpha=0.5)

    for idx in range(n_stations, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(f"Station time series - {var_label}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fname = output_dir / f"station_timeseries_{tag}.png"
    savefig(fig, fname, show)


def write_summary_csv(cases: dict[str, dict], path: Path) -> None:
    rows = []
    day = next(iter(cases.values())).get("day", "")

    for case_id, data in cases.items():
        for key in SUMMARY_FIELDS:
            if key not in data:
                continue
            rows.extend(field_summary_rows(day, case_id, key, data[key]))

        rows.extend(split_summary_rows(day, case_id, data))

    for case_hi, case_lo in diff_pairs(cases):
        rows.extend(diff_summary_rows(day, case_hi, case_lo, cases[case_hi], cases[case_lo]))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "day",
                "section",
                "case",
                "case_hi",
                "case_lo",
                "field",
                "metric",
                "value",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def field_summary_rows(day: int, case_id: str, field: str, arr: np.ndarray) -> list[dict]:
    mean_field = np.nanmean(arr, axis=1)
    max_field = np.nanmax(arr, axis=1)
    n_nodes = mean_field.size
    positive_mean = int(np.count_nonzero(np.isfinite(mean_field) & (mean_field > 0.0)))
    positive_max = int(np.count_nonzero(np.isfinite(max_field) & (max_field > 0.0)))

    metrics = {
        "day_mean_max": np.nanmax(mean_field),
        "day_mean_mean": np.nanmean(mean_field),
        "day_mean_positive_nodes": positive_mean,
        "day_mean_positive_fraction": positive_mean / max(n_nodes, 1),
        "pointwise_max_max": np.nanmax(max_field),
        "pointwise_max_mean": np.nanmean(max_field),
        "pointwise_max_positive_nodes": positive_max,
        "pointwise_max_positive_fraction": positive_max / max(n_nodes, 1),
    }
    return [
        summary_row(day, "case_field", case_id, "", "", field, metric, value)
        for metric, value in metrics.items()
    ]


def split_summary_rows(day: int, case_id: str, data: dict) -> list[dict]:
    rows = []
    split_groups = [
        ("mp1_split_davg", "mp1_davg", "mp1_agg_davg", "mp1_dis_davg"),
        ("mp1_split_surf", "mp1_surf", "mp1_agg_surf", "mp1_dis_surf"),
        ("bot_mass_split", "bot_mass", "bot_mass_agg", "bot_mass_dis"),
    ]

    for label, total_key, agg_key, dis_key in split_groups:
        if total_key not in data or agg_key not in data or dis_key not in data:
            continue
        err = data[agg_key] + data[dis_key] - data[total_key]
        metrics = {
            "max_abs_error": np.nanmax(np.abs(err)),
            "mean_abs_error": np.nanmean(np.abs(err)),
        }
        for metric, value in metrics.items():
            rows.append(summary_row(day, "split_consistency", case_id, "", "", label, metric, value))

    return rows


def diff_summary_rows(
    day: int,
    case_hi: str,
    case_lo: str,
    data_hi: dict,
    data_lo: dict,
) -> list[dict]:
    rows = []
    for key, _, _ in DIFF_FIELDS:
        if key not in data_hi or key not in data_lo:
            continue
        diff = np.nanmean(data_hi[key], axis=1) - np.nanmean(data_lo[key], axis=1)
        metrics = {
            "abs_max": np.nanmax(np.abs(diff)),
            "abs_p99": np.nanpercentile(np.abs(diff), 99),
            "signed_mean": np.nanmean(diff),
            "positive_nodes": int(np.count_nonzero(np.isfinite(diff) & (diff > 0.0))),
            "negative_nodes": int(np.count_nonzero(np.isfinite(diff) & (diff < 0.0))),
        }
        for metric, value in metrics.items():
            rows.append(summary_row(day, "cross_case_diff", "", case_hi, case_lo, key, metric, value))
    return rows


def summary_row(
    day: int,
    section: str,
    case: str,
    case_hi: str,
    case_lo: str,
    field: str,
    metric: str,
    value: Any,
) -> dict:
    if isinstance(value, (np.floating, float)):
        value_out = f"{float(value):.10e}"
    elif isinstance(value, (np.integer, int)):
        value_out = str(int(value))
    else:
        value_out = str(value)
    return {
        "day": day,
        "section": section,
        "case": case,
        "case_hi": case_hi,
        "case_lo": case_lo,
        "field": field,
        "metric": metric,
        "value": value_out,
    }


def diff_pairs(cases: dict[str, dict]) -> list[tuple[str, str]]:
    pairs = []
    if "b" in cases and "a" in cases:
        pairs.append(("b", "a"))
    if "c" in cases and "a" in cases:
        pairs.append(("c", "a"))
    if "c" in cases and "b" in cases:
        pairs.append(("c", "b"))
    return pairs


def time_window_label(data: dict, prefix: str) -> str:
    snap_time = data.get("snap_time_plot", data["snap_time"])
    if snap_time is None or len(snap_time) == 0:
        return f"{prefix} over selected window"
    return f"{prefix} over days {snap_time[0]:.3f}-{snap_time[-1]:.3f}"


def tripcolor(
    ax: plt.Axes,
    triang: mtri.Triangulation,
    values: np.ndarray,
    cmap: str,
    vmin: float,
    vmax: float,
    log_scale: bool = True,
    label: str = "",
) -> None:
    if log_scale:
        v = np.where(np.isfinite(values) & (values > 0.0), values, vmin)
        norm = LogNorm(vmin=vmin, vmax=vmax)
        tcf = ax.tripcolor(triang, v, cmap=cmap, norm=norm, shading="gouraud")
    else:
        tcf = ax.tripcolor(triang, values, cmap=cmap, vmin=vmin, vmax=vmax, shading="gouraud")
    plt.colorbar(tcf, ax=ax, label=label, shrink=0.8, pad=0.02)


def apply_map_extent(ax: plt.Axes, lon: np.ndarray, lat: np.ndarray) -> None:
    lon_min = float(lon.min())
    lat_max = float(lat.max())
    ax.set_xlim(lon_min, EXTENT_LON_MAX)
    ax.set_ylim(EXTENT_LAT_MIN, lat_max)
    ax.set_aspect("equal", adjustable="box")


def add_stations(ax: plt.Axes, data: dict) -> None:
    node_lon = data.get("station_node_lon")
    node_lat = data.get("station_node_lat")
    if node_lon is None or node_lat is None:
        return

    ax.scatter(
        node_lon,
        node_lat,
        s=40,
        c="red",
        marker="^",
        zorder=5,
        linewidths=0.5,
        edgecolors="white",
        label="stations",
    )
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


def savefig(fig: plt.Figure, path: Path, show: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"    Saved: {path.name}")
    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
