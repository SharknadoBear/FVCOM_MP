#!/usr/bin/env python3
"""Compare FVCOM-MP mp1 mass-budget diagnostics for cases a, b, and c.

The MATLAB extractors save compact v7.3 MAT files in:

    ../OUTPUT_a/mp1_mass_budget_a.mat
    ../OUTPUT_b/mp1_mass_budget_b.mat
    ../OUTPUT_c/mp1_mass_budget_c.mat

This script reads those files, normalizes the different field names used by
the baseline and floc-enabled diagnostics, and writes comparison plots plus
CSV summaries under ``output/mass_conservation``.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

try:
    import h5py
except ImportError as exc:  # pragma: no cover - exercised only on missing envs
    raise SystemExit(
        "This script needs h5py to read MATLAB v7.3 MAT files. "
        "Install the Python requirements with: pip install -r requirements.txt"
    ) from exc

try:
    import netCDF4
    _HAS_NETCDF4 = True
except ImportError:
    _HAS_NETCDF4 = False

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
TEST_ROOT = SCRIPT_DIR.parent
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "output" / "mass_conservation"


@dataclass
class CaseBudget:
    """Normalized mass-budget fields for one experiment case."""

    case_id: str
    mat_file: Path
    raw_fields: dict[str, np.ndarray]
    frame: pd.DataFrame
    summary: dict[str, float | str | int]


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.mat_file is not None and len(args.cases) != 1:
        raise SystemExit("--mat-file can only be used when exactly one case is requested.")

    cases = []
    for case_id in args.cases:
        mat_file = resolve_case_mat_file(args.root, case_id, args.mat_file)
        print(f"Reading case {case_id}: {mat_file}")
        cases.append(load_case_budget(case_id, mat_file))

    time_origin_days = determine_time_origin(cases, args.time_origin_days)
    apply_common_time_origin(cases, time_origin_days)
    print(f"\nComparison time origin: FVCOM day {time_origin_days:.6f}")

    write_tables(cases, output_dir)
    river_nc = args.river_nc
    if river_nc is None:
        default_nc = SCRIPT_DIR.parent / "INPUT" / "waterPACT_riv_floc_MP.nc"
        if default_nc.is_file():
            river_nc = default_nc

    make_plots(cases, output_dir, time_origin_days=time_origin_days,
               river_nc=river_nc, iramp=args.iramp,
               extstep_seconds=args.extstep_seconds, isplit=args.isplit,
               show=args.show)

    print("\nMass-conservation comparison complete.")
    print(f"Outputs written to: {output_dir}")
    print("Key files:")
    print(f"  {output_dir / 'mass_summary.csv'}")
    print(f"  {output_dir / 'mass_timeseries_long.csv'}")
    for name in [
        "total_mass_by_case.png",
        "mass_change_by_case.png",
        "max_concentration_by_case.png",
        "cap_threshold_counts.png",
        "water_bed_partition.png",
        "floc_split_components_b_c.png",
        "split_consistency_b_c.png",
        "deposition_erosion_exchange.png",
    ]:
        path = output_dir / name
        if path.exists():
            print(f"  {path}")

    return 0


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot and summarize mp1 mass conservation diagnostics."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=TEST_ROOT,
        help="FVCOM_MP_test_run root containing OUTPUT_a/OUTPUT_b/OUTPUT_c.",
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        default=["a", "b", "c"],
        help="Case IDs to compare. Defaults to a b c.",
    )
    parser.add_argument(
        "--mat-file",
        type=Path,
        default=None,
        help=(
            "Optional exact MAT file. Only valid for one case; otherwise "
            "defaults to OUTPUT_<case>/mp1_mass_budget_<case>.mat."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for CSV tables and PNG figures.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display plots interactively in addition to saving PNG files.",
    )
    parser.add_argument(
        "--time-origin-days",
        type=float,
        default=None,
        help=(
            "FVCOM day used as x=0 in comparison plots. "
            "Default: earliest finite output time among loaded cases."
        ),
    )
    parser.add_argument(
        "--river-nc",
        type=Path,
        default=None,
        help=(
            "River NetCDF file (waterPACT_riv_floc_MP.nc) used to compute "
            "the analytic M0 + int(Q*C dt) river-source mass line. "
            "Defaults to INPUT/waterPACT_riv_floc_MP.nc next to this script."
        ),
    )
    parser.add_argument(
        "--iramp",
        type=int,
        default=216000,
        help=(
            "FVCOM IRAMP value (number of external barotropic steps over which "
            "the tanh ramp factor reaches ~1). Set 0 to disable ramp correction. "
            "Default: 216000 (matching RUN_a/b/c nml)."
        ),
    )
    parser.add_argument(
        "--extstep-seconds",
        type=float,
        default=0.4,
        dest="extstep_seconds",
        help="FVCOM EXTSTEP_SECONDS (s). Default: 0.4.",
    )
    parser.add_argument(
        "--isplit",
        type=int,
        default=5,
        help=(
            "FVCOM ISPLIT (external steps per internal step). "
            "Used in ramp: TMP = t_sec / (ISPLIT * EXTSTEP). Default: 5."
        ),
    )
    return parser.parse_args(argv)


def resolve_case_mat_file(root: Path, case_id: str, explicit_mat: Path | None) -> Path:
    if explicit_mat is not None:
        if len(case_id) != 1:
            raise ValueError("--mat-file is intended for a single designated case.")
        return explicit_mat.resolve()
    return (root / f"OUTPUT_{case_id}" / f"mp1_mass_budget_{case_id}.mat").resolve()


def load_case_budget(case_id: str, mat_file: Path) -> CaseBudget:
    if not mat_file.is_file():
        raise FileNotFoundError(f"Missing MAT file for case {case_id}: {mat_file}")

    with h5py.File(mat_file, "r") as handle:
        diagnostic = find_diagnostic_group(handle)
        budget_group = diagnostic["budget"]
        raw_fields = read_budget_fields(budget_group)

    frame = normalize_budget_frame(case_id, raw_fields)
    summary = summarize_case(case_id, mat_file, frame, raw_fields)
    return CaseBudget(case_id=case_id, mat_file=mat_file, raw_fields=raw_fields, frame=frame, summary=summary)


def determine_time_origin(cases: Iterable[CaseBudget], requested_origin: float | None) -> float:
    if requested_origin is not None:
        return float(requested_origin)

    starts = [
        first_finite(case.frame["time_days"].to_numpy(dtype=float))
        for case in cases
    ]
    starts = [value for value in starts if np.isfinite(value)]
    if not starts:
        raise ValueError("Could not determine a common time origin from the loaded cases.")
    return float(min(starts))


def apply_common_time_origin(cases: Iterable[CaseBudget], time_origin_days: float) -> None:
    for case in cases:
        frame = case.frame
        frame["comparison_days"] = frame["time_days"] - time_origin_days
        case.summary["comparison_time_origin_days"] = float(time_origin_days)
        case.summary["comparison_start_days"] = float(first_finite(frame["comparison_days"].to_numpy()))
        case.summary["comparison_end_days"] = float(last_finite(frame["comparison_days"].to_numpy()))


def find_diagnostic_group(handle: h5py.File) -> h5py.Group:
    """Return the diagnostic struct group from a MATLAB v7.3 MAT file."""

    if "diagnostic" in handle:
        return handle["diagnostic"]

    if "case_diagnostic" in handle:
        return handle["case_diagnostic"]

    if "combined" in handle and "diagnostics" in handle["combined"]:
        refs = np.asarray(handle["combined/diagnostics"])
        if refs.size == 0:
            raise KeyError("combined/diagnostics is empty.")
        return handle[refs.flat[0]]

    available = ", ".join(handle.keys())
    raise KeyError(f"Could not find diagnostic struct in MAT file. Root keys: {available}")


def read_budget_fields(budget_group: h5py.Group) -> dict[str, np.ndarray]:
    fields: dict[str, np.ndarray] = {}
    for name, obj in budget_group.items():
        if isinstance(obj, h5py.Dataset) and obj.dtype.kind != "O":
            fields[name] = read_numeric_dataset(obj)
    return fields


def read_numeric_dataset(dataset: h5py.Dataset) -> np.ndarray:
    data = np.asarray(dataset)
    matlab_class = decode_attr(dataset.attrs.get("MATLAB_class", b""))
    if matlab_class == "char":
        return np.asarray([decode_matlab_char(data)])
    return np.asarray(data, dtype=float).squeeze()


def decode_attr(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "decode"):
        return value.decode("utf-8", errors="replace")
    return str(value)


def decode_matlab_char(data: np.ndarray) -> str:
    chars = []
    for code in data.ravel(order="F"):
        code_int = int(code)
        if code_int:
            chars.append(chr(code_int))
    return "".join(chars)


def normalize_budget_frame(case_id: str, fields: dict[str, np.ndarray]) -> pd.DataFrame:
    time_days = vector(fields, "time_days")
    elapsed_days = time_days - first_finite(time_days)

    frame = pd.DataFrame(
        {
            "case": case_id,
            "time_days": time_days,
            "elapsed_days": elapsed_days,
            "record_index": optional_vector(fields, "record_index", len(time_days)),
            "iint": optional_vector(fields, "iint", len(time_days)),
        }
    )

    frame["water_mass_kg"] = pick_vector(
        fields, len(frame), "water_mass_kg", "water_mp1_mass_kg"
    )
    frame["bed_mass_kg"] = pick_vector(
        fields, len(frame), "bed_mass_all_nodes_kg", "bed_mp1_mass_kg"
    )
    frame["total_mass_kg"] = pick_vector(
        fields, len(frame), "total_mass_all_bed_kg", "total_combined_mass_kg"
    )
    frame["total_mass_change_kg"] = pick_vector(
        fields, len(frame), "total_mass_change_kg", "total_combined_change_kg"
    )
    frame["total_mass_rate_kg_per_day"] = pick_vector(
        fields, len(frame), "total_mass_rate_kg_per_day", "total_combined_rate_kg_per_day"
    )
    frame["water_negative_mass_kg"] = pick_vector(
        fields, len(frame), "water_negative_mass_kg", "water_mp1_negative_mass_kg"
    )
    frame["bed_negative_mass_kg"] = pick_vector(
        fields, len(frame), "bot_mass_negative_kg", "bed_mp1_negative_mass_kg"
    )
    frame["dep_flux_area_integral_kg"] = optional_vector(
        fields, "dep_flux_area_integral_kg", len(frame)
    )
    frame["ero_flux_area_integral_kg"] = optional_vector(
        fields, "ero_flux_area_integral_kg", len(frame)
    )
    frame["max_mp1_conc_kgm3"] = pick_vector(
        fields, len(frame), "max_mp1_conc_kgm3", "mp1_max_kgm3"
    )
    frame["max_mp1_agg_conc_kgm3"] = pick_vector(
        fields, len(frame), "max_mp1_agg_conc_kgm3", "mp1_agg_max_kgm3"
    )
    frame["max_mp1_dis_conc_kgm3"] = pick_vector(
        fields, len(frame), "max_mp1_dis_conc_kgm3", "mp1_dis_max_kgm3"
    )
    frame["mp1_count_ge_cap"] = optional_vector(fields, "mp1_count_ge_cap", len(frame))
    frame["mp1_count_near_or_ge_cap"] = optional_vector(
        fields, "mp1_count_near_or_ge_cap", len(frame)
    )
    frame["mp1_agg_count_ge_cap"] = optional_vector(
        fields, "mp1_agg_count_ge_cap", len(frame)
    )
    frame["mp1_agg_count_near_or_ge_cap"] = optional_vector(
        fields, "mp1_agg_count_near_or_ge_cap", len(frame)
    )
    frame["mp1_dis_count_ge_cap"] = optional_vector(
        fields, "mp1_dis_count_ge_cap", len(frame)
    )
    frame["mp1_dis_count_near_or_ge_cap"] = optional_vector(
        fields, "mp1_dis_count_near_or_ge_cap", len(frame)
    )
    frame["mp1_volume_near_or_ge_cap_m3"] = optional_vector(
        fields, "mp1_volume_near_or_ge_cap_m3", len(frame)
    )
    frame["mp1_agg_volume_near_or_ge_cap_m3"] = optional_vector(
        fields, "mp1_agg_volume_near_or_ge_cap_m3", len(frame)
    )
    frame["mp1_dis_volume_near_or_ge_cap_m3"] = optional_vector(
        fields, "mp1_dis_volume_near_or_ge_cap_m3", len(frame)
    )

    # Floc split fields are only present for cases b/c.
    frame["water_agg_mass_kg"] = optional_vector(fields, "water_agg_mass_kg", len(frame))
    frame["water_dis_mass_kg"] = optional_vector(fields, "water_dis_mass_kg", len(frame))
    frame["water_split_sum_mass_kg"] = optional_vector(
        fields, "water_split_sum_mass_kg", len(frame)
    )
    frame["water_split_minus_mp1_kg"] = optional_vector(
        fields, "water_split_minus_mp1_kg", len(frame)
    )
    frame["bed_agg_mass_kg"] = optional_vector(fields, "bed_agg_mass_kg", len(frame))
    frame["bed_dis_mass_kg"] = optional_vector(fields, "bed_dis_mass_kg", len(frame))
    frame["bed_split_sum_mass_kg"] = optional_vector(fields, "bed_split_sum_mass_kg", len(frame))
    frame["bed_split_minus_mp1_kg"] = optional_vector(fields, "bed_split_minus_mp1_kg", len(frame))
    frame["total_split_mass_kg"] = optional_vector(fields, "total_split_mass_kg", len(frame))
    frame["total_split_minus_combined_kg"] = optional_vector(
        fields, "total_split_minus_combined_kg", len(frame)
    )

    if frame["total_mass_change_kg"].isna().all():
        frame["total_mass_change_kg"] = frame["total_mass_kg"] - first_finite(
            frame["total_mass_kg"].to_numpy()
        )

    frame["dep_flux_sampled_cumsum_kg"] = np.nancumsum(
        frame["dep_flux_area_integral_kg"].to_numpy(dtype=float)
    )
    frame["ero_flux_sampled_cumsum_kg"] = np.nancumsum(
        frame["ero_flux_area_integral_kg"].to_numpy(dtype=float)
    )
    frame["bed_exchange_sampled_cumsum_kg"] = np.nancumsum(
        (
            frame["dep_flux_area_integral_kg"]
            - frame["ero_flux_area_integral_kg"]
        ).to_numpy(dtype=float)
    )
    return frame


def summarize_case(
    case_id: str,
    mat_file: Path,
    frame: pd.DataFrame,
    fields: dict[str, np.ndarray],
) -> dict[str, float | str | int]:
    total = frame["total_mass_kg"].to_numpy(dtype=float)
    water = frame["water_mass_kg"].to_numpy(dtype=float)
    bed = frame["bed_mass_kg"].to_numpy(dtype=float)

    summary: dict[str, float | str | int] = {
        "case": case_id,
        "mat_file": str(mat_file),
        "records": int(len(frame)),
        "time_start_days": float(first_finite(frame["time_days"].to_numpy())),
        "time_end_days": float(last_finite(frame["time_days"].to_numpy())),
        "elapsed_days": float(last_finite(frame["elapsed_days"].to_numpy())),
        "initial_total_mass_kg": float(first_finite(total)),
        "final_total_mass_kg": float(last_finite(total)),
        "final_total_change_kg": float(last_finite(frame["total_mass_change_kg"].to_numpy())),
        "max_abs_total_change_kg": float(np.nanmax(np.abs(frame["total_mass_change_kg"]))),
        "final_water_mass_kg": float(last_finite(water)),
        "final_bed_mass_kg": float(last_finite(bed)),
        "max_water_negative_mass_kg": float(np.nanmax(frame["water_negative_mass_kg"])),
        "max_bed_negative_mass_kg": float(np.nanmax(frame["bed_negative_mass_kg"])),
        "has_split_fields": bool("total_split_mass_kg" in fields),
        "max_mp1_conc_kgm3": nanmax_value(frame["max_mp1_conc_kgm3"]),
        "max_mp1_agg_conc_kgm3": nanmax_value(frame["max_mp1_agg_conc_kgm3"]),
        "max_mp1_dis_conc_kgm3": nanmax_value(frame["max_mp1_dis_conc_kgm3"]),
        "max_mp1_count_near_or_ge_cap": nanmax_value(frame["mp1_count_near_or_ge_cap"]),
        "max_mp1_agg_count_near_or_ge_cap": nanmax_value(frame["mp1_agg_count_near_or_ge_cap"]),
        "max_mp1_dis_count_near_or_ge_cap": nanmax_value(frame["mp1_dis_count_near_or_ge_cap"]),
    }

    for key in [
        "water_split_minus_mp1_kg",
        "bed_split_minus_mp1_kg",
        "total_split_minus_combined_kg",
    ]:
        values = frame[key].to_numpy(dtype=float) if key in frame else np.asarray([np.nan])
        summary[f"max_abs_{key}"] = float(nanmax_abs(values))

    return summary


def write_tables(cases: Iterable[CaseBudget], output_dir: Path) -> None:
    cases = list(cases)
    summary = pd.DataFrame([case.summary for case in cases])
    summary.to_csv(output_dir / "mass_summary.csv", index=False)

    timeseries = pd.concat([case.frame for case in cases], ignore_index=True)
    timeseries.to_csv(output_dir / "mass_timeseries_long.csv", index=False)

    print("\nSummary:")
    cols = [
        "case",
        "records",
        "elapsed_days",
        "comparison_start_days",
        "comparison_end_days",
        "initial_total_mass_kg",
        "final_total_mass_kg",
        "final_total_change_kg",
        "max_abs_total_split_minus_combined_kg",
    ]
    existing = [col for col in cols if col in summary.columns]
    print(summary[existing].to_string(index=False))


def compute_analytic_river_mass(
    river_nc: Path,
    time_origin_days: float,
    query_days: np.ndarray,
    iramp: int = 0,
    extstep_seconds: float = 1.0,
    isplit: int = 1,
) -> np.ndarray | None:
    """Integrate sum_j RAMP(t)*Q_j(t)*C_j(t) from the river NC file and return
    cumulative mass (kg) at each time in *query_days* (days since time_origin).

    iramp            FVCOM IRAMP (number of external steps).  0 = no ramp.
    extstep_seconds  FVCOM EXTSTEP_SECONDS (external barotropic step, s).
    isplit           FVCOM ISPLIT (external steps per internal step).

    RAMP(t) = tanh(TMP/IRAMP)  where  TMP = t_sec / (ISPLIT * EXTSTEP)
    (matches external_step.F: TMP = (IINT-1) + IEXT/ISPLIT)
    Returns None when netCDF4 is unavailable or the file cannot be read.
    """
    if not _HAS_NETCDF4:
        print("  [analytic line] netCDF4 not available, skipping.")
        return None
    try:
        ds = netCDF4.Dataset(river_nc)
        # time in 'days since 1858-11-17 00:00:00' (Modified Julian Day)
        riv_time_mjd = np.asarray(ds.variables["time"][:], dtype=float)
        q = np.asarray(ds.variables["river_flux"][:], dtype=float)   # (nt, nriv) m3/s
        c = np.asarray(ds.variables["mp1"][:], dtype=float)           # (nt, nriv) kg/m3
        ds.close()
    except Exception as exc:
        print(f"  [analytic line] could not read river NC: {exc}")
        return None

    # Total Q*C rate [kg/s] at each river-file timestep
    qc_rate = np.sum(q * c, axis=1)          # (nt,)

    # Convert river file time (MJD) to comparison days (days since time_origin)
    riv_comp_days = riv_time_mjd - time_origin_days

    # Apply FVCOM tanh ramp: RAMP(t) = tanh(TMP/IRAMP)
    # where TMP = t_sec / (ISPLIT * EXTSTEP)  [matches external_step.F]
    if iramp > 0 and extstep_seconds > 0.0 and isplit > 0:
        t_sec = np.maximum(riv_comp_days * 86400.0, 0.0)  # clamp pre-origin to 0
        tmp   = t_sec / (float(isplit) * extstep_seconds)
        ramp  = np.tanh(tmp / float(iramp))               # (nt,)
        qc_rate = qc_rate * ramp

    # Cumulative integral [kg] using trapezoidal rule
    dt_days = np.diff(riv_comp_days)               # (nt-1,) days
    dt_sec  = dt_days * 86400.0                    # seconds
    mid_rate = 0.5 * (qc_rate[:-1] + qc_rate[1:]) # (nt-1,) kg/s
    dM       = mid_rate * dt_sec                   # (nt-1,) kg per interval

    riv_cumulative = np.zeros(len(riv_comp_days))
    riv_cumulative[1:] = np.cumsum(dM)

    # Anchor so that cumulative mass added = 0 at comparison day 0 (model start)
    offset = float(np.interp(0.0, riv_comp_days, riv_cumulative))
    riv_cumulative -= offset

    # Interpolate onto query_days
    analytic = np.interp(
        query_days,
        riv_comp_days,
        riv_cumulative,
        left=np.nan,
        right=np.nan,
    )
    return analytic


def make_plots(
    cases: Iterable[CaseBudget],
    output_dir: Path,
    time_origin_days: float,
    river_nc: Path | None = None,
    iramp: int = 0,
    extstep_seconds: float = 1.0,
    isplit: int = 1,
    show: bool = False,
) -> None:
    cases = list(cases)
    plt.style.use("seaborn-v0_8-whitegrid")
    x_label = f"Days since FVCOM day {time_origin_days:.6f}"

    plot_total_mass(cases, output_dir, x_label,
                    time_origin_days=time_origin_days, river_nc=river_nc,
                    iramp=iramp, extstep_seconds=extstep_seconds, isplit=isplit)
    plot_mass_change(cases, output_dir, x_label)
    plot_max_concentration(cases, output_dir, x_label)
    plot_cap_counts(cases, output_dir, x_label)
    plot_water_bed_partition(cases, output_dir, x_label)
    plot_floc_split_components(cases, output_dir, x_label)
    plot_split_consistency(cases, output_dir, x_label)
    plot_deposition_erosion(cases, output_dir, x_label)

    if show:
        plt.show()
    else:
        plt.close("all")


def plot_total_mass(
    cases: list[CaseBudget],
    output_dir: Path,
    x_label: str,
    time_origin_days: float = 0.0,
    river_nc: Path | None = None,
    iramp: int = 0,
    extstep_seconds: float = 1.0,
    isplit: int = 1,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))

    # --- model output lines ---
    for case in cases:
        frame = case.frame
        ax.plot(
            frame["comparison_days"],
            frame["total_mass_kg"],
            label=f"case {case.case_id}",
            linewidth=2,
        )

    # --- analytic M^0 + int(Q*C dt) river-source line ---
    if river_nc is not None and river_nc.is_file():
        # Use a dense query grid spanning all loaded cases
        all_days = np.concatenate([
            case.frame["comparison_days"].to_numpy(dtype=float)
            for case in cases
        ])
        finite_days = all_days[np.isfinite(all_days)]
        if finite_days.size > 0:
            query_days = np.linspace(finite_days.min(), finite_days.max(), 2000)
            analytic_dM = compute_analytic_river_mass(
                river_nc, time_origin_days, query_days,
                iramp=iramp, extstep_seconds=extstep_seconds, isplit=isplit,
            )
            if analytic_dM is not None:
                # Anchor to the initial total mass of the first case
                M0 = first_finite(
                    cases[0].frame["total_mass_kg"].to_numpy(dtype=float)
                )
                analytic_mass = M0 + analytic_dM
                # Downsample to visible dots (every ~50th point)
                stride = max(1, len(query_days) // 120)
                ramp_label = (
                    r"$M^0+\int \mathrm{RAMP}(t)\,QC\,dt$ (analytic, IRAMP)"
                    if iramp > 0
                    else r"$M^0+\int QC\,dt$ (analytic)"
                )
                ax.plot(
                    query_days[::stride],
                    analytic_mass[::stride],
                    marker="o",
                    markersize=5,
                    markerfacecolor=(1.0, 0.0, 0.0, 0.35),
                    markeredgecolor=(0.7, 0.0, 0.0, 0.6),
                    linestyle="none",
                    label=ramp_label,
                    zorder=5,
                )

    ax.set_title("Domain-Integrated mp1 Total Mass")
    ax.set_xlabel(x_label)
    ax.set_ylabel("Water + bed mass (kg)")
    ax.legend()
    save_figure(fig, output_dir / "total_mass_by_case.png")


def plot_mass_change(cases: list[CaseBudget], output_dir: Path, x_label: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    for case in cases:
        frame = case.frame
        ax.plot(
            frame["comparison_days"],
            frame["total_mass_change_kg"],
            label=f"case {case.case_id}",
            linewidth=2,
        )
    ax.axhline(0.0, color="0.25", linewidth=1)
    ax.set_title("Total Mass Change From First Output Record")
    ax.set_xlabel(x_label)
    ax.set_ylabel("Mass change (kg)")
    ax.legend()
    save_figure(fig, output_dir / "mass_change_by_case.png")


def plot_max_concentration(cases: list[CaseBudget], output_dir: Path, x_label: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    any_values = False
    for case in cases:
        frame = case.frame
        if frame["max_mp1_conc_kgm3"].notna().any():
            any_values = True
            ax.plot(
                frame["comparison_days"],
                frame["max_mp1_conc_kgm3"],
                label=f"case {case.case_id} mp1",
                linewidth=2,
            )
        if frame["max_mp1_agg_conc_kgm3"].notna().any():
            any_values = True
            ax.plot(
                frame["comparison_days"],
                frame["max_mp1_agg_conc_kgm3"],
                label=f"case {case.case_id} agg",
                linewidth=1.8,
                linestyle="--",
            )
        if frame["max_mp1_dis_conc_kgm3"].notna().any():
            any_values = True
            ax.plot(
                frame["comparison_days"],
                frame["max_mp1_dis_conc_kgm3"],
                label=f"case {case.case_id} dis",
                linewidth=1.8,
                linestyle=":",
            )
    if not any_values:
        plt.close(fig)
        return
    ax.axhline(100.0, color="0.2", linewidth=1.2, linestyle="-.", label="100 kg/m3 cap")
    ax.set_title("Maximum Water-Column mp1 Concentration")
    ax.set_xlabel(x_label)
    ax.set_ylabel("Max concentration (kg/m3)")
    ax.legend()
    save_figure(fig, output_dir / "max_concentration_by_case.png")


def plot_cap_counts(cases: list[CaseBudget], output_dir: Path, x_label: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    any_values = False
    for case in cases:
        frame = case.frame
        if frame["mp1_count_near_or_ge_cap"].notna().any():
            any_values = True
            ax.plot(
                frame["comparison_days"],
                frame["mp1_count_near_or_ge_cap"],
                label=f"case {case.case_id} mp1",
                linewidth=2,
            )
        if frame["mp1_agg_count_near_or_ge_cap"].notna().any():
            any_values = True
            ax.plot(
                frame["comparison_days"],
                frame["mp1_agg_count_near_or_ge_cap"],
                label=f"case {case.case_id} agg",
                linewidth=1.8,
                linestyle="--",
            )
        if frame["mp1_dis_count_near_or_ge_cap"].notna().any():
            any_values = True
            ax.plot(
                frame["comparison_days"],
                frame["mp1_dis_count_near_or_ge_cap"],
                label=f"case {case.case_id} dis",
                linewidth=1.8,
                linestyle=":",
            )
    if not any_values:
        plt.close(fig)
        return
    ax.set_title("Grid Cells Near or Above 100 kg/m3 Cap")
    ax.set_xlabel(x_label)
    ax.set_ylabel("Cell count")
    ax.legend()
    save_figure(fig, output_dir / "cap_threshold_counts.png")


def plot_water_bed_partition(cases: list[CaseBudget], output_dir: Path, x_label: str) -> None:
    fig, axes = plt.subplots(len(cases), 1, figsize=(9, 3.2 * len(cases)), sharex=False)
    axes = np.atleast_1d(axes)
    for ax, case in zip(axes, cases):
        frame = case.frame
        ax.plot(frame["comparison_days"], frame["water_mass_kg"], label="water", linewidth=2)
        ax.plot(frame["comparison_days"], frame["bed_mass_kg"], label="bed", linewidth=2)
        ax.plot(frame["comparison_days"], frame["total_mass_kg"], label="total", linewidth=1.8, linestyle="--")
        ax.set_title(f"Case {case.case_id}: Water/Bed Partition")
        ax.set_ylabel("Mass (kg)")
        ax.legend(loc="best")
    axes[-1].set_xlabel(x_label)
    save_figure(fig, output_dir / "water_bed_partition.png")


def plot_floc_split_components(cases: list[CaseBudget], output_dir: Path, x_label: str) -> None:
    floc_cases = [case for case in cases if case.frame["water_agg_mass_kg"].notna().any()]
    if not floc_cases:
        return

    fig, axes = plt.subplots(len(floc_cases), 2, figsize=(12, 3.5 * len(floc_cases)), sharex=False)
    axes = np.atleast_2d(axes)
    for row, case in zip(axes, floc_cases):
        frame = case.frame
        ax = row[0]
        ax.plot(frame["comparison_days"], frame["water_agg_mass_kg"], label="water agg", linewidth=2)
        ax.plot(frame["comparison_days"], frame["water_dis_mass_kg"], label="water dis", linewidth=2)
        ax.plot(frame["comparison_days"], frame["water_mass_kg"], label="water mp1", linestyle="--", linewidth=1.8)
        ax.set_title(f"Case {case.case_id}: Water Split")
        ax.set_ylabel("Mass (kg)")
        ax.legend()

        ax = row[1]
        ax.plot(frame["comparison_days"], frame["bed_agg_mass_kg"], label="bed agg", linewidth=2)
        ax.plot(frame["comparison_days"], frame["bed_dis_mass_kg"], label="bed dis", linewidth=2)
        ax.plot(frame["comparison_days"], frame["bed_mass_kg"], label="bed mp1", linestyle="--", linewidth=1.8)
        ax.set_title(f"Case {case.case_id}: Bed Split")
        ax.legend()
    for ax in axes[-1, :]:
        ax.set_xlabel(x_label)
    save_figure(fig, output_dir / "floc_split_components_b_c.png")


def plot_split_consistency(cases: list[CaseBudget], output_dir: Path, x_label: str) -> None:
    floc_cases = [case for case in cases if case.frame["total_split_minus_combined_kg"].notna().any()]
    if not floc_cases:
        return

    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for case in floc_cases:
        frame = case.frame
        label = f"case {case.case_id}"
        axes[0].plot(frame["comparison_days"], frame["water_split_minus_mp1_kg"], label=label, linewidth=2)
        axes[1].plot(frame["comparison_days"], frame["bed_split_minus_mp1_kg"], label=label, linewidth=2)
        axes[2].plot(frame["comparison_days"], frame["total_split_minus_combined_kg"], label=label, linewidth=2)

    axes[0].set_title("Split Consistency: Water")
    axes[0].set_ylabel("agg + dis - mp1 (kg)")
    axes[1].set_title("Split Consistency: Bed")
    axes[1].set_ylabel("split - total (kg)")
    axes[2].set_title("Split Consistency: Total")
    axes[2].set_ylabel("split - combined (kg)")
    axes[2].set_xlabel(x_label)
    for ax in axes:
        ax.axhline(0.0, color="0.25", linewidth=1)
        ax.legend()
    save_figure(fig, output_dir / "split_consistency_b_c.png")


def plot_deposition_erosion(cases: list[CaseBudget], output_dir: Path, x_label: str) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=False)
    any_flux = False
    for case in cases:
        frame = case.frame
        if frame["dep_flux_area_integral_kg"].notna().any():
            any_flux = True
            axes[0].plot(
                frame["comparison_days"],
                frame["dep_flux_area_integral_kg"],
                label=f"case {case.case_id} dep",
                linewidth=2,
            )
        if frame["ero_flux_area_integral_kg"].notna().any():
            any_flux = True
            axes[0].plot(
                frame["comparison_days"],
                frame["ero_flux_area_integral_kg"],
                label=f"case {case.case_id} ero",
                linestyle="--",
                linewidth=2,
            )
        if frame["bed_exchange_sampled_cumsum_kg"].notna().any():
            axes[1].plot(
                frame["comparison_days"],
                frame["bed_exchange_sampled_cumsum_kg"],
                label=f"case {case.case_id}",
                linewidth=2,
            )

    if not any_flux:
        plt.close(fig)
        return

    axes[0].set_title("Area-Integrated Deposition/Erosion Flux Per Output Record")
    axes[0].set_ylabel("kg per output record")
    axes[0].legend()
    axes[1].set_title("Output-Sampled Cumulative Bed Exchange: cumsum(dep - ero)")
    axes[1].set_xlabel(x_label)
    axes[1].set_ylabel("kg")
    axes[1].axhline(0.0, color="0.25", linewidth=1)
    axes[1].legend()
    save_figure(fig, output_dir / "deposition_erosion_exchange.png")


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")


def vector(fields: dict[str, np.ndarray], name: str) -> np.ndarray:
    if name not in fields:
        raise KeyError(f"Missing required budget field: {name}")
    values = np.asarray(fields[name], dtype=float).squeeze()
    return np.atleast_1d(values)


def optional_vector(fields: dict[str, np.ndarray], name: str, length: int) -> np.ndarray:
    if name not in fields:
        return np.full(length, np.nan)
    values = np.asarray(fields[name], dtype=float).squeeze()
    values = np.atleast_1d(values)
    if values.size == 1 and length != 1:
        return np.full(length, float(values[0]))
    return values


def pick_vector(fields: dict[str, np.ndarray], length: int, *names: str) -> np.ndarray:
    for name in names:
        if name in fields:
            return optional_vector(fields, name, length)
    return np.full(length, np.nan)


def first_finite(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    return float(finite[0]) if finite.size else np.nan


def last_finite(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    return float(finite[-1]) if finite.size else np.nan


def nanmax_abs(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).any():
        return np.nan
    return float(np.nanmax(np.abs(values)))


def nanmax_value(values: pd.Series | np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).any():
        return np.nan
    return float(np.nanmax(values))


if __name__ == "__main__":
    sys.exit(main())
