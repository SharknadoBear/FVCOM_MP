#!/usr/bin/env python3
"""Spatial distribution comparison for FVCOM-MP case pairs.

This script reads the compact MAT files written by:

    MATLAB/extract_spatial_distribution_v2_day10.m
    MATLAB/extract_spatial_distribution_v2_day30.m
    MATLAB/extract_spatial_distribution_v2_day90.m

By default it processes cases c and d for days 10, 30, and 90, then writes
per-case maps, configurable difference maps, CSV summaries, and a short README
under:

    PYTHON/output/spatial_distribution_v2/day10
    PYTHON/output/spatial_distribution_v2/day30
    PYTHON/output/spatial_distribution_v2/day90
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
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
DEFAULT_OUTPUT_BASE = SCRIPT_DIR / "output" / "spatial_distribution_v2"
DEFAULT_DAYS = (10, 30, 90)
DEFAULT_CASES = ("c", "d")

# Lon is stored in 0-360 convention, around 284-286 deg for Delaware Bay.
EXTENT_LAT_MIN = 39.4
EXTENT_LON_MAX = 285.5
LOG_VMIN_FLOOR = 1.0e-12


@dataclass(frozen=True)
class FieldSpec:
    source: str
    tag: str
    label: str
    units: str
    cmap: str
    log_scale: bool


FIELD_SPECS = [
    FieldSpec("sed_davg", "sed_davg", "coarse_sand_1 depth-avg", "model units", "Blues", True),
    FieldSpec("sed_surf", "sed_surf", "coarse_sand_1 surface", "model units", "Blues", True),
    FieldSpec("sed_bot", "sed_bot", "coarse_sand_1 bottom", "model units", "Blues", True),
    FieldSpec("wset_floc_davg", "wset_floc_davg", "floc settling speed depth-avg", "m s-1", "viridis", True),
    FieldSpec("wset_floc_surf", "wset_floc_surf", "floc settling speed surface", "m s-1", "viridis", True),
    FieldSpec("wset_floc_bot", "wset_floc_bot", "floc settling speed bottom", "m s-1", "viridis", True),
    FieldSpec("d_floc_eff_um_davg", "dfloc_eff_um_davg", "Stokes-equivalent floc size depth-avg", "um", "magma", False),
    FieldSpec("d_floc_eff_um_surf", "dfloc_eff_um_surf", "Stokes-equivalent floc size surface", "um", "magma", False),
    FieldSpec("d_floc_eff_um_bot", "dfloc_eff_um_bot", "Stokes-equivalent floc size bottom", "um", "magma", False),
    FieldSpec("lambda_a_davg", "lambda_a_davg", "aggregation rate depth-avg", "s-1", "PuRd", True),
    FieldSpec("lambda_d_davg", "lambda_d_davg", "disaggregation rate depth-avg", "s-1", "PuBuGn", True),
    FieldSpec("mp1_davg", "mp1_davg", "mp1 total depth-avg", "kg m-3", "plasma", True),
    FieldSpec("mp1_agg_davg", "mp1_agg_davg", "mp1 aggregated depth-avg", "kg m-3", "plasma", True),
    FieldSpec("mp1_dis_davg", "mp1_dis_davg", "mp1 dispersed depth-avg", "kg m-3", "plasma", True),
    FieldSpec("agg_fraction_davg", "agg_fraction_davg", "aggregated mp1 fraction depth-avg", "-", "cividis", False),
    FieldSpec("bot_mass", "bot_mass", "bottom mp1 mass", "kg m-2", "YlOrRd", True),
    FieldSpec("bot_mass_agg", "bot_mass_agg", "bottom aggregated mp1 mass", "kg m-2", "YlOrRd", True),
    FieldSpec("bot_mass_dis", "bot_mass_dis", "bottom dispersed mp1 mass", "kg m-2", "YlOrRd", True),
]


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    output_base = args.output_base.resolve()

    processed_days = 0
    for day in args.days:
        print(f"\nDay {day}")
        output_dir = output_base / f"day{day}"
        output_dir.mkdir(parents=True, exist_ok=True)

        cases = load_cases(root, args.cases, day, strict=args.strict)
        if not cases:
            print(f"  No MAT files found for day {day}; skipping.")
            continue

        triang = build_triangulation(next(iter(cases.values())))
        for case_id, data in cases.items():
            print(f"  Per-case maps: case {case_id}")
            plot_case_maps(case_id, data, triang, output_dir, args.show)

        diff_hi = args.diff_hi.lower()
        diff_lo = args.diff_lo.lower()
        if diff_hi in cases and diff_lo in cases and not args.no_diff_maps:
            print(f"  Difference maps: case {diff_hi} - case {diff_lo}")
            plot_diff_maps(
                diff_hi,
                diff_lo,
                cases[diff_hi],
                cases[diff_lo],
                triang,
                output_dir,
                args.show,
            )
        elif not args.no_diff_maps:
            print(
                "  Difference maps skipped; both comparison cases "
                f"{diff_hi} and {diff_lo} are required."
            )

        write_summary_stats(cases, output_dir / "summary_stats.csv", day)
        write_diff_summary_stats(
            cases,
            output_dir / "diff_summary_stats.csv",
            day,
            diff_hi,
            diff_lo,
        )
        write_day_readme(output_dir / "README.md", day, cases, diff_hi, diff_lo)
        processed_days += 1

    if processed_days == 0:
        return 1

    print(f"\nOutputs written under: {output_base}")
    return 0


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate spatial c/d comparison maps from FVCOM-MP v2 MAT files."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=TEST_ROOT,
        help="FVCOM_MP_test_run root containing OUTPUT_<case> folders.",
    )
    parser.add_argument(
        "--days",
        nargs="+",
        type=int,
        default=list(DEFAULT_DAYS),
        help="Simulation days to process. Default: 10 30 90.",
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        default=list(DEFAULT_CASES),
        help="Case IDs to process. Default: c d.",
    )
    parser.add_argument(
        "--output-base",
        type=Path,
        default=DEFAULT_OUTPUT_BASE,
        help="Base directory for dayXX output folders.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Raise an error when an expected MAT file is missing.",
    )
    parser.add_argument(
        "--no-diff-maps",
        action="store_true",
        help="Skip case-difference maps.",
    )
    parser.add_argument(
        "--diff-hi",
        default="d",
        help="Case ID for the positive side of difference maps/stats. Default: d.",
    )
    parser.add_argument(
        "--diff-lo",
        default="c",
        help="Case ID for the negative side of difference maps/stats. Default: c.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show plots interactively in addition to saving PNG files.",
    )
    return parser.parse_args(argv)


def load_cases(root: Path, case_ids: list[str], day: int, strict: bool) -> dict[str, dict]:
    cases: dict[str, dict] = {}
    for case_id in case_ids:
        case_id = case_id.lower()
        mat_path = resolve_mat_path(root, case_id, day)
        if not mat_path.is_file():
            msg = f"  Missing case {case_id} MAT file: {mat_path}"
            if strict:
                raise FileNotFoundError(msg)
            print(msg)
            continue

        print(f"  Loading case {case_id}: {mat_path}")
        data = load_mat_file(mat_path)
        data["case_id"] = case_id
        data["day"] = day
        data["mat_path"] = str(mat_path)
        add_derived_fields(data)
        cases[case_id] = data
    return cases


def resolve_mat_path(root: Path, case_id: str, day: int) -> Path:
    return (
        root
        / f"OUTPUT_{case_id}"
        / f"spatial_distribution_v2_{case_id}_day{day}.mat"
    ).resolve()


def load_mat_file(mat_path: Path) -> dict[str, Any]:
    data: dict[str, Any] = {}
    with h5py.File(mat_path, "r") as handle:
        root_group = handle["out"] if "out" in handle else handle

        grid = root_group["grid"]
        data["lon"] = read_1d(grid["lon"])
        data["lat"] = read_1d(grid["lat"])
        data["lonc"] = read_1d(grid["lonc"]) if "lonc" in grid else None
        data["latc"] = read_1d(grid["latc"]) if "latc" in grid else None
        data["h"] = read_1d(grid["h"])
        data["nv"] = read_2d(grid["nv"]).astype(int)
        data["n_nodes"] = int(read_scalar(grid["n_nodes"])) if "n_nodes" in grid else data["lon"].size

        snaps = root_group["snaps"]
        data["n_records"] = int(read_scalar(snaps["n_records"])) if "n_records" in snaps else 0
        data["record_time"] = read_1d(snaps["record_time"]) if "record_time" in snaps else np.array([])
        data["units"] = read_units(snaps["units"]) if "units" in snaps else {}

        for key in snaps:
            obj = snaps[key]
            if isinstance(obj, h5py.Group) or key in {"vars_available"}:
                continue
            arr = try_read_numeric(obj)
            if arr is None:
                continue
            if arr.ndim == 0:
                data[key] = float(arr)
            else:
                data[key] = np.asarray(arr, dtype=np.float64).squeeze()
    return data


def read_units(group: h5py.Group) -> dict[str, str]:
    units: dict[str, str] = {}
    for key in group:
        value = read_matlab_string(group[key])
        if value:
            units[key] = value
    return units


def read_matlab_string(obj: h5py.Dataset) -> str:
    try:
        arr = np.asarray(obj)
    except Exception:
        return ""
    if arr.dtype.kind in {"S", "U"}:
        return "".join(np.asarray(arr).astype(str).ravel()).strip()
    if np.issubdtype(arr.dtype, np.integer):
        chars = [chr(int(v)) for v in arr.ravel() if int(v) != 0]
        return "".join(chars).strip()
    return ""


def try_read_numeric(obj: h5py.Dataset) -> np.ndarray | None:
    try:
        if not np.issubdtype(obj.dtype, np.number):
            return None
        return np.asarray(obj)
    except TypeError:
        return None


def read_1d(ds: h5py.Dataset) -> np.ndarray:
    return np.asarray(ds, dtype=np.float64).squeeze().ravel()


def read_2d(ds: h5py.Dataset) -> np.ndarray:
    arr = np.asarray(ds)
    if arr.ndim != 2:
        return np.asarray(arr)
    arr = arr.T
    if arr.shape[0] == 3 and arr.shape[1] != 3:
        arr = arr.T
    return arr


def read_scalar(ds: h5py.Dataset) -> float:
    return float(np.asarray(ds).squeeze())


def add_derived_fields(data: dict[str, Any]) -> None:
    add_fraction(data, "mp1_agg_surf", "mp1_surf", "agg_fraction_surf")
    add_fraction(data, "mp1_agg_bot", "mp1_bot", "agg_fraction_bot")
    add_fraction(data, "mp1_agg_davg", "mp1_davg", "agg_fraction_davg")

    # Accept either spelling if hand-edited MAT files are introduced later.
    for suffix in ("surf", "bot", "davg"):
        src = f"d_floc_eff_um_{suffix}"
        alt = f"dfloc_eff_um_{suffix}"
        if src not in data and alt in data:
            data[src] = data[alt]


def add_fraction(data: dict[str, Any], num_key: str, den_key: str, out_key: str) -> None:
    if num_key not in data or den_key not in data:
        return
    num = ensure_1d(data[num_key])
    den = ensure_1d(data[den_key])
    frac = np.full_like(num, np.nan, dtype=np.float64)
    valid = np.isfinite(num) & np.isfinite(den) & (den > 0.0)
    frac[valid] = num[valid] / den[valid]
    data[out_key] = frac


def build_triangulation(data: dict[str, Any]) -> mtri.Triangulation:
    lon = ensure_1d(data["lon"])
    lat = ensure_1d(data["lat"])
    nv = np.asarray(data["nv"], dtype=int)
    if nv.ndim != 2:
        raise ValueError("Grid nv array is not two-dimensional.")
    if nv.shape[1] != 3 and nv.shape[0] == 3:
        nv = nv.T
    if nv.shape[1] != 3:
        raise ValueError(f"Grid nv array has unexpected shape {nv.shape}.")
    if nv.min() >= 1:
        nv = nv - 1
    return mtri.Triangulation(lon, lat, nv)


def plot_case_maps(
    case_id: str,
    data: dict[str, Any],
    triang: mtri.Triangulation,
    output_dir: Path,
    show: bool,
) -> None:
    for spec in FIELD_SPECS:
        if spec.source not in data:
            print(f"    skip {spec.tag}: missing")
            continue
        values = ensure_1d(data[spec.source])
        if values.size != triang.x.size:
            print(f"    skip {spec.tag}: expected {triang.x.size} nodes, got {values.size}")
            continue
        fig, ax = plt.subplots(figsize=(9.5, 7.0))
        draw_scalar_map(ax, triang, values, spec, diff=False)
        apply_map_extent(ax, data["lon"], data["lat"])
        ax.set_title(f"Case {case_id} | day {data['day']} | {spec.label}", fontsize=10)
        ax.set_xlabel("Longitude (deg, 0-360)")
        ax.set_ylabel("Latitude (deg N)")
        fname = output_dir / f"case{case_id}_{spec.tag}.png"
        savefig(fig, fname, show)


def plot_diff_maps(
    case_hi: str,
    case_lo: str,
    data_hi: dict[str, Any],
    data_lo: dict[str, Any],
    triang: mtri.Triangulation,
    output_dir: Path,
    show: bool,
) -> None:
    for spec in FIELD_SPECS:
        if spec.source not in data_hi or spec.source not in data_lo:
            print(f"    skip diff {spec.tag}: missing")
            continue
        hi = ensure_1d(data_hi[spec.source])
        lo = ensure_1d(data_lo[spec.source])
        if hi.size != lo.size or hi.size != triang.x.size:
            print(f"    skip diff {spec.tag}: incompatible array sizes")
            continue
        diff = hi - lo

        fig, ax = plt.subplots(figsize=(9.5, 7.0))
        draw_scalar_map(ax, triang, diff, spec, diff=True)
        apply_map_extent(ax, data_hi["lon"], data_hi["lat"])
        ax.set_title(
            f"Case {case_hi} - case {case_lo} | day {data_hi['day']} | {spec.label}",
            fontsize=10,
        )
        ax.set_xlabel("Longitude (deg, 0-360)")
        ax.set_ylabel("Latitude (deg N)")
        fname = output_dir / f"diff_{spec.tag}_case{case_hi}_minus_case{case_lo}.png"
        savefig(fig, fname, show)


def draw_scalar_map(
    ax: plt.Axes,
    triang: mtri.Triangulation,
    values: np.ndarray,
    spec: FieldSpec,
    diff: bool,
) -> None:
    plot_values = np.asarray(values, dtype=np.float64).copy()
    plot_values[~np.isfinite(plot_values)] = np.nan

    if diff:
        finite = plot_values[np.isfinite(plot_values)]
        vmax = float(np.nanpercentile(np.abs(finite), 99.0)) if finite.size else 1.0
        vmax = max(vmax, LOG_VMIN_FLOOR)
        tpc = ax.tripcolor(
            triang,
            plot_values,
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            shading="gouraud",
        )
        plt.colorbar(tpc, ax=ax, label=f"delta {spec.label} [{spec.units}]", shrink=0.82)
        return

    norm = None
    vmin, vmax = scalar_limits(plot_values, spec)
    if spec.log_scale and has_positive(plot_values):
        plot_values[plot_values <= 0.0] = np.nan
        norm = LogNorm(vmin=vmin, vmax=vmax)
        vmin = None
        vmax = None

    tpc = ax.tripcolor(
        triang,
        plot_values,
        cmap=spec.cmap,
        norm=norm,
        vmin=vmin,
        vmax=vmax,
        shading="gouraud",
    )
    plt.colorbar(tpc, ax=ax, label=f"{spec.label} [{spec.units}]", shrink=0.82)


def scalar_limits(values: np.ndarray, spec: FieldSpec) -> tuple[float | None, float | None]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0.0, 1.0

    if spec.log_scale and has_positive(values):
        positive = finite[finite > 0.0]
        vmin = max(float(np.nanpercentile(positive, 1.0)), LOG_VMIN_FLOOR)
        vmax = max(float(np.nanpercentile(positive, 99.0)), vmin * 10.0)
        return vmin, vmax

    if spec.tag.startswith("agg_fraction"):
        return 0.0, 1.0

    vmin = float(np.nanpercentile(finite, 1.0))
    vmax = float(np.nanpercentile(finite, 99.0))
    if np.isclose(vmin, vmax):
        pad = max(abs(vmax) * 0.05, LOG_VMIN_FLOOR)
        vmin -= pad
        vmax += pad
    return vmin, vmax


def has_positive(values: np.ndarray) -> bool:
    return bool(np.any(np.isfinite(values) & (values > 0.0)))


def apply_map_extent(ax: plt.Axes, lon: np.ndarray, lat: np.ndarray) -> None:
    lon = ensure_1d(lon)
    lat = ensure_1d(lat)
    ax.set_xlim(float(np.nanmin(lon)), min(EXTENT_LON_MAX, float(np.nanmax(lon))))
    ax.set_ylim(max(EXTENT_LAT_MIN, float(np.nanmin(lat))), float(np.nanmax(lat)))
    ax.set_aspect("equal", adjustable="box")


def savefig(fig: plt.Figure, path: Path, show: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"    saved {path.name}")
    if show:
        plt.show()
    plt.close(fig)


def write_summary_stats(cases: dict[str, dict], path: Path, day: int) -> None:
    rows: list[dict[str, Any]] = []
    for case_id, data in cases.items():
        for spec in FIELD_SPECS:
            if spec.source not in data:
                continue
            rows.append(summary_row(day, case_id, spec, ensure_1d(data[spec.source])))

    write_csv(
        path,
        rows,
        [
            "day",
            "case",
            "field",
            "units",
            "count_finite",
            "min",
            "p05",
            "mean",
            "median",
            "p95",
            "max",
            "positive_fraction",
        ],
    )


def write_diff_summary_stats(
    cases: dict[str, dict],
    path: Path,
    day: int,
    case_hi: str,
    case_lo: str,
) -> None:
    rows: list[dict[str, Any]] = []
    if case_hi not in cases or case_lo not in cases:
        write_csv(
            path,
            rows,
            [
                "day",
                "case_hi",
                "case_lo",
                "field",
                "units",
                "count_finite",
                "signed_mean",
                "abs_p95",
                "abs_max",
                "positive_nodes",
                "negative_nodes",
            ],
        )
        return

    data_hi = cases[case_hi]
    data_lo = cases[case_lo]
    for spec in FIELD_SPECS:
        if spec.source not in data_hi or spec.source not in data_lo:
            continue
        hi = ensure_1d(data_hi[spec.source])
        lo = ensure_1d(data_lo[spec.source])
        if hi.size != lo.size:
            continue
        rows.append(diff_summary_row(day, case_hi, case_lo, spec, hi - lo))

    write_csv(
        path,
        rows,
        [
            "day",
            "case_hi",
            "case_lo",
            "field",
            "units",
            "count_finite",
            "signed_mean",
            "abs_p95",
            "abs_max",
            "positive_nodes",
            "negative_nodes",
        ],
    )


def summary_row(day: int, case_id: str, spec: FieldSpec, values: np.ndarray) -> dict[str, Any]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {
            "day": day,
            "case": case_id,
            "field": spec.tag,
            "units": spec.units,
            "count_finite": 0,
            "min": "",
            "p05": "",
            "mean": "",
            "median": "",
            "p95": "",
            "max": "",
            "positive_fraction": "",
        }
    return {
        "day": day,
        "case": case_id,
        "field": spec.tag,
        "units": spec.units,
        "count_finite": finite.size,
        "min": format_float(np.nanmin(finite)),
        "p05": format_float(np.nanpercentile(finite, 5.0)),
        "mean": format_float(np.nanmean(finite)),
        "median": format_float(np.nanmedian(finite)),
        "p95": format_float(np.nanpercentile(finite, 95.0)),
        "max": format_float(np.nanmax(finite)),
        "positive_fraction": format_float(np.count_nonzero(finite > 0.0) / finite.size),
    }


def diff_summary_row(
    day: int,
    case_hi: str,
    case_lo: str,
    spec: FieldSpec,
    diff: np.ndarray,
) -> dict[str, Any]:
    finite = diff[np.isfinite(diff)]
    if finite.size == 0:
        return {
            "day": day,
            "case_hi": case_hi,
            "case_lo": case_lo,
            "field": spec.tag,
            "units": spec.units,
            "count_finite": 0,
            "signed_mean": "",
            "abs_p95": "",
            "abs_max": "",
            "positive_nodes": 0,
            "negative_nodes": 0,
        }
    return {
        "day": day,
        "case_hi": case_hi,
        "case_lo": case_lo,
        "field": spec.tag,
        "units": spec.units,
        "count_finite": finite.size,
        "signed_mean": format_float(np.nanmean(finite)),
        "abs_p95": format_float(np.nanpercentile(np.abs(finite), 95.0)),
        "abs_max": format_float(np.nanmax(np.abs(finite))),
        "positive_nodes": int(np.count_nonzero(finite > 0.0)),
        "negative_nodes": int(np.count_nonzero(finite < 0.0)),
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {path.name}")


def write_day_readme(
    path: Path,
    day: int,
    cases: dict[str, dict],
    diff_hi: str,
    diff_lo: str,
) -> None:
    case_lines = []
    for case_id, data in cases.items():
        n_records = data.get("n_records", "")
        mat_path = data.get("mat_path", "")
        case_lines.append(f"- Case {case_id}: `{mat_path}` ({n_records} records averaged)")
    case_block = "\n".join(case_lines) if case_lines else "- No case files loaded"

    text = f"""# Spatial Distribution v2: Day {day}

This folder contains plan-view maps and statistics for the day-{day} final-24-hour
average extracted from the FVCOM-MP case outputs.

## Inputs
{case_block}

Difference figures and `diff_summary_stats.csv` use `{diff_hi} - {diff_lo}`.

## Fields
The maps include sediment concentration, reconstructed floc settling speed,
Stokes-equivalent floc size, MP aggregation/disaggregation rates, MP total and
split concentrations, aggregated MP fraction, and bottom MP mass fields.

The equivalent floc size saved in the MAT files is computed from:

```text
nu = 1.0e-6 * exp(-0.025 * (temp - 20))
mu = ambient_rho * nu
d_floc_eff = sqrt(18 * mu * settle_vel_floc_mp1 /
                  (9.81 * (1300.0 - ambient_rho)))
```

The MATLAB extractor bounds `nu` to `[0.5e-6, 2.0e-6]` and uses 20 degC if
`temp` is unavailable.
"""
    path.write_text(text, encoding="utf-8")
    print(f"  Wrote {path.name}")


def ensure_1d(value: Any) -> np.ndarray:
    return np.asarray(value, dtype=np.float64).squeeze().ravel()


def format_float(value: Any) -> str:
    if value is None or not np.isfinite(float(value)):
        return ""
    return f"{float(value):.10e}"


if __name__ == "__main__":
    sys.exit(main())
