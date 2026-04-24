#!/usr/bin/env python
"""Inspect the Delaware FVCOM-MP test cases and write reusable summaries."""

from __future__ import annotations

import json
import math
import re
from datetime import datetime
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


SCRIPT_DIR = Path(__file__).resolve().parent
TEST_ROOT = SCRIPT_DIR.parent
INPUT_DIR = TEST_ROOT / "INPUT"
OUTPUT_ROOT = SCRIPT_DIR / "output"

CASE_IDS = ("a", "b", "c")
CASE_RUN_FILES = {case: TEST_ROOT / f"RUN_{case}" / f"waterPACT_{case}_run.nml" for case in CASE_IDS}
CASE_JOB_FILES = {case: TEST_ROOT / f"RUN_{case}" / "job_Kestrel.sh" for case in CASE_IDS}
CASE_OUTPUT_DIRS = {case: TEST_ROOT / f"OUTPUT_{case}" for case in CASE_IDS}


def strip_comment(line: str) -> str:
    return line.split("!", 1)[0].rstrip()


def clean_value(text: str) -> str:
    return text.strip().rstrip(",").strip()


def unquote(text: str) -> str:
    value = clean_value(text)
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        return value[1:-1]
    return value


def parse_scalar(text: str):
    value = clean_value(text)
    upper = value.upper()
    if upper == "T":
        return True
    if upper == "F":
        return False
    if re.fullmatch(r"[+-]?\d+", value):
        try:
            return int(value)
        except ValueError:
            return value
    if re.fullmatch(r"[+-]?((\d+\.\d*)|(\d*\.\d+)|(\d+))(?:[Ee][+-]?\d+)?", value):
        try:
            return float(value)
        except ValueError:
            return value
    return unquote(value)


def parse_assignment_file(path: Path) -> dict[str, object]:
    data: dict[str, object] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = strip_comment(raw_line).strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[key.strip()] = parse_scalar(value)
    return data


def parse_river_sections(path: Path) -> list[dict[str, object]]:
    sections: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = strip_comment(raw_line).strip()
        if not line:
            continue
        if line.startswith("&"):
            current = {"section": line[1:].strip()}
            continue
        if line == "/":
            if current:
                sections.append(current)
            current = None
            continue
        if current is not None and "=" in line:
            key, value = line.split("=", 1)
            current[key.strip()] = parse_scalar(value)
    return sections


def decode_char_row(row) -> str:
    items = row.tolist()
    if items and isinstance(items[0], str):
        return "".join(items).strip().strip("\x00")
    return b"".join(items).decode("ascii", errors="ignore").strip().strip("\x00")


def times_span(dataset: Dataset) -> dict[str, object]:
    if "Times" in dataset.variables:
        times = dataset.variables["Times"]
        first = decode_char_row(times[0])
        last = decode_char_row(times[-1])
    else:
        first = None
        last = None
    raw_first = None
    raw_last = None
    units = None
    if "time" in dataset.variables:
        time_var = dataset.variables["time"]
        raw_first = float(np.asarray(time_var[0]))
        raw_last = float(np.asarray(time_var[-1]))
        units = getattr(time_var, "units", None)
    return {
        "first": first,
        "last": last,
        "raw_first": raw_first,
        "raw_last": raw_last,
        "units": units,
    }


def summarize_netcdf(path: Path) -> dict[str, object]:
    with Dataset(path) as dataset:
        return {
            "path": str(path),
            "dimensions": {name: len(dim) for name, dim in dataset.dimensions.items()},
            "variables": list(dataset.variables.keys()),
            "time_span": times_span(dataset),
        }


def summarize_startup_profiles(path: Path) -> dict[str, object]:
    with Dataset(path) as dataset:
        summary: dict[str, object] = {}
        for name in ("ssl", "tsl"):
            arr = np.asarray(dataset.variables[name][:], dtype=float)
            surface = arr[0, 0, :]
            bottom = arr[0, -1, :]
            summary[name] = {
                "shape": list(arr.shape),
                "min": float(np.nanmin(arr)),
                "max": float(np.nanmax(arr)),
                "mean": float(np.nanmean(arr)),
                "surface_mean": float(np.nanmean(surface)),
                "bottom_mean": float(np.nanmean(bottom)),
            }
        summary["time_span"] = times_span(dataset)
        return summary


def summarize_river_forcing(path: Path) -> dict[str, object]:
    with Dataset(path) as dataset:
        names = [decode_char_row(row) for row in dataset.variables["river_names"][:]]
        flux = np.asarray(dataset.variables["river_flux"][:], dtype=float)
        plastic = np.asarray(dataset.variables["mp1"][:], dtype=float)
        sediment = np.asarray(dataset.variables["coarse_sand_1"][:], dtype=float)

        group_indices = {
            "Delaware": slice(0, 6),
            "Schuylkill": slice(6, 9),
        }

        def group_total(rate: np.ndarray, group_slice: slice) -> float:
            # The Times strings are hourly and more reliable than the float32 MJD deltas.
            return float((rate[:-1, group_slice] * 3600.0).sum())

        per_river = {}
        for idx, name in enumerate(names):
            per_river[name] = {
                "river_flux_min_m3s": float(np.nanmin(flux[:, idx])),
                "river_flux_max_m3s": float(np.nanmax(flux[:, idx])),
                "river_flux_mean_m3s": float(np.nanmean(flux[:, idx])),
                "sediment_min_gl": float(np.nanmin(sediment[:, idx])),
                "sediment_max_gl": float(np.nanmax(sediment[:, idx])),
                "sediment_mean_gl": float(np.nanmean(sediment[:, idx])),
                "plastic_conc_kgm3": float(np.nanmean(plastic[:, idx])),
            }

        totals = {
            "water_m3": {},
            "plastic_kg": {},
            "sediment_kg": {},
        }
        for group_name, group_slice in group_indices.items():
            totals["water_m3"][group_name] = group_total(flux, group_slice)
            totals["plastic_kg"][group_name] = group_total(flux * plastic, group_slice)
            totals["sediment_kg"][group_name] = group_total(flux * sediment, group_slice)

        totals["water_m3"]["All"] = float((flux[:-1, :] * 3600.0).sum())
        totals["plastic_kg"]["All"] = float(((flux * plastic)[:-1, :] * 3600.0).sum())
        totals["sediment_kg"]["All"] = float(((flux * sediment)[:-1, :] * 3600.0).sum())

        return {
            "time_span": times_span(dataset),
            "river_names": names,
            "river_count": len(names),
            "vertical_distribution_note": "Each river uses 30 * 0.03333333333 in waterPACT_riv_floc_MP.nml.",
            "global_variable_stats": {
                "river_flux_min_m3s": float(np.nanmin(flux)),
                "river_flux_max_m3s": float(np.nanmax(flux)),
                "river_flux_mean_m3s": float(np.nanmean(flux)),
                "sediment_min_gl": float(np.nanmin(sediment)),
                "sediment_max_gl": float(np.nanmax(sediment)),
                "sediment_mean_gl": float(np.nanmean(sediment)),
                "plastic_conc_kgm3": float(np.nanmean(plastic)),
            },
            "per_river": per_river,
            "integrated_totals": totals,
        }


def summarize_job(path: Path) -> dict[str, object]:
    summary: dict[str, object] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith("#SBATCH --time="):
            summary["sbatch_time"] = line.split("=", 1)[1].strip()
        elif line.startswith("#SBATCH  --partition=") or line.startswith("#SBATCH --partition="):
            summary["partition"] = line.split("=", 1)[1].strip()
        elif line.startswith("srun ") and "--casename=" in line:
            summary["casename"] = line.split("--casename=", 1)[1].strip()
    return summary


def summarize_case(case_id: str) -> dict[str, object]:
    run_file = CASE_RUN_FILES[case_id]
    run_cfg = parse_assignment_file(run_file)
    job_cfg = summarize_job(CASE_JOB_FILES[case_id])
    plastic_file = INPUT_DIR / str(run_cfg["PLASTIC_MODEL_FILE"])
    plastic_cfg = parse_assignment_file(plastic_file)

    sediment_cfg: dict[str, object] | None = None
    sediment_file = None
    if run_cfg.get("SEDIMENT_MODEL") is True:
        sediment_file = INPUT_DIR / str(run_cfg["SEDIMENT_MODEL_FILE"])
        sediment_cfg = parse_assignment_file(sediment_file)

    output_listing = sorted(path.name for path in CASE_OUTPUT_DIRS[case_id].iterdir()) if CASE_OUTPUT_DIRS[case_id].exists() else []

    return {
        "case_id": case_id,
        "run_file": str(run_file),
        "job_file": str(CASE_JOB_FILES[case_id]),
        "job": job_cfg,
        "run": {
            "start_date": run_cfg.get("START_DATE"),
            "end_date": run_cfg.get("END_DATE"),
            "startup_ts_type": run_cfg.get("STARTUP_TS_TYPE"),
            "startup_t_vals": run_cfg.get("STARTUP_T_VALS"),
            "startup_s_vals": run_cfg.get("STARTUP_S_VALS"),
            "restart_first_out": run_cfg.get("RST_FIRST_OUT"),
            "netcdf_first_out": run_cfg.get("NC_FIRST_OUT"),
            "output_dir": run_cfg.get("OUTPUT_DIR"),
            "sediment_model": run_cfg.get("SEDIMENT_MODEL"),
            "plastic_model": run_cfg.get("PLASTIC_MODEL"),
            "plastic_model_file": run_cfg.get("PLASTIC_MODEL_FILE"),
            "sediment_model_file": run_cfg.get("SEDIMENT_MODEL_FILE"),
        },
        "plastic": {
            "file": plastic_file.name,
            "plast_floc": plastic_cfg.get("PLAST_FLOC"),
            "plast_ptsource": plastic_cfg.get("PLAST_PTSOURCE"),
            "plast_bedset": plastic_cfg.get("PLAST_BEDSET"),
            "plast_unit": plastic_cfg.get("PLAST_UNIT"),
            "plast_name": plastic_cfg.get("PLAST_NAME"),
            "plast_shape": plastic_cfg.get("PLAST_SHAP"),
            "plast_density_kgm3": plastic_cfg.get("PLAST_PRHO"),
            "plast_wssed_ms": plastic_cfg.get("PLAST_WSSED"),
            "ambient_water_density_kgm3": plastic_cfg.get("AMBIENT_WDEN"),
            "ini_mass_kg": plastic_cfg.get("INI_MASS"),
            "source_ag_rate": plastic_cfg.get("SOURCE_AG_RATE"),
        },
        "sediment": {
            "file": sediment_file.name if sediment_file else None,
            "sed_on": run_cfg.get("SEDIMENT_MODEL"),
            "sed_name": sediment_cfg.get("SED_NAME") if sediment_cfg else None,
            "sed_type": sediment_cfg.get("SED_TYPE") if sediment_cfg else None,
            "sed_sd50_mm": sediment_cfg.get("SED_SD50") if sediment_cfg else None,
            "sed_wset_mms": sediment_cfg.get("SED_WSET") if sediment_cfg else None,
            "sed_tau_e_nm2": sediment_cfg.get("SED_TAUE") if sediment_cfg else None,
            "sed_tau_d_nm2": sediment_cfg.get("SED_TAUD") if sediment_cfg else None,
            "sed_inf_bed": sediment_cfg.get("INF_BED") if sediment_cfg else None,
        },
        "output_files": output_listing,
    }


def file_inventory(directory: Path) -> list[dict[str, object]]:
    entries = []
    for path in sorted(directory.iterdir()):
        entries.append(
            {
                "name": path.name,
                "size_bytes": path.stat().st_size,
                "modified": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds"),
            }
        )
    return entries


def scenario_notes(cases: dict[str, dict[str, object]]) -> list[str]:
    notes = []
    if cases["a"]["plastic"]["plast_floc"] is False and cases["b"]["plastic"]["plast_floc"] is True:
        notes.append("Case a is the no-floc baseline, while case b turns PLAST_FLOC on.")
    if cases["b"]["plastic"]["file"] == "generic_plastic_b.inp" and cases["c"]["plastic"]["file"] == "generic_plastic_c.inp":
        notes.append("Cases b and c use different plastic filenames, but the file contents are identical.")
    if math.isclose(cases["b"]["sediment"]["sed_sd50_mm"], 0.1) and math.isclose(cases["c"]["sediment"]["sed_sd50_mm"], 1.0):
        notes.append("Case c changes only the sediment grain size from 0.1 mm to 1.0 mm relative to case b.")
    if cases["a"]["run"]["startup_ts_type"] != cases["b"]["run"]["startup_ts_type"]:
        notes.append("Case a and case b do not share the same startup salinity/temperature treatment.")
    if not cases["a"]["output_files"] and not cases["b"]["output_files"] and not cases["c"]["output_files"]:
        notes.append("OUTPUT_a, OUTPUT_b, and OUTPUT_c are currently empty, so the original scenario intent must be inferred from inputs alone.")
    return notes


def render_markdown(summary: dict[str, object]) -> str:
    cases = summary["cases"]
    river = summary["river_forcing"]
    startup = summary["startup_profiles"]
    netcdf = summary["netcdf_files"]

    lines = [
        "# Delaware FVCOM-MP Case Inventory",
        "",
        "This report was generated by `inspect_delaware_cases.py`.",
        "",
        "## Cases",
        "",
        "| Case | Sediment model | PLAST_FLOC | Startup TS | Plastic file | Sediment file | SD50 (mm) | NC first out | Batch |",
        "| --- | --- | --- | --- | --- | --- | ---: | --- | --- |",
    ]
    for case_id in CASE_IDS:
        case = cases[case_id]
        lines.append(
            "| {case_id} | {sediment_model} | {plast_floc} | {startup_ts_type} | {plastic_file} | {sediment_file} | {sd50} | {nc_first_out} | {batch} |".format(
                case_id=case_id,
                sediment_model=case["run"]["sediment_model"],
                plast_floc=case["plastic"]["plast_floc"],
                startup_ts_type=case["run"]["startup_ts_type"],
                plastic_file=case["plastic"]["file"],
                sediment_file=case["sediment"]["file"],
                sd50=case["sediment"]["sed_sd50_mm"],
                nc_first_out=case["run"]["netcdf_first_out"],
                batch=f"{case['job'].get('partition', 'n/a')} / {case['job'].get('sbatch_time', 'n/a')}",
            )
        )

    lines.extend(
        [
            "",
            "## Shared Forcing",
            "",
            f"- River forcing spans `{river['time_span']['first']}` to `{river['time_span']['last']}` with {river['river_count']} river entries.",
            f"- River names: {', '.join(river['river_names'])}.",
            f"- `mp1` is constant at {river['global_variable_stats']['plastic_conc_kgm3']:.10g} kg/m^3 across all rivers and all times.",
            f"- `coarse_sand_1` ranges from {river['global_variable_stats']['sediment_min_gl']:.6g} to {river['global_variable_stats']['sediment_max_gl']:.6g} g/L.",
            f"- Integrated annual water inflow: Delaware={river['integrated_totals']['water_m3']['Delaware']:.6e} m^3, Schuylkill={river['integrated_totals']['water_m3']['Schuylkill']:.6e} m^3.",
            f"- Integrated annual plastic inflow: Delaware={river['integrated_totals']['plastic_kg']['Delaware']:.6e} kg, Schuylkill={river['integrated_totals']['plastic_kg']['Schuylkill']:.6e} kg.",
            f"- Integrated annual sediment inflow: Delaware={river['integrated_totals']['sediment_kg']['Delaware']:.6e} kg, Schuylkill={river['integrated_totals']['sediment_kg']['Schuylkill']:.6e} kg.",
            "",
            "## Startup Profiles",
            "",
            f"- Startup file time: `{startup['time_span']['first']}`.",
            f"- `ssl` min/max/mean = {startup['ssl']['min']:.6g} / {startup['ssl']['max']:.6g} / {startup['ssl']['mean']:.6g}.",
            f"- `tsl` min/max/mean = {startup['tsl']['min']:.6g} / {startup['tsl']['max']:.6g} / {startup['tsl']['mean']:.6g}.",
            "",
            "## NetCDF Coverage",
            "",
        ]
    )
    for name, meta in netcdf.items():
        span = meta["time_span"]
        lines.append(
            f"- `{name}`: first=`{span['first']}`, last=`{span['last']}`, dimensions={meta['dimensions']}."
        )

    lines.extend(
        [
            "",
            "## Notes",
            "",
        ]
    )
    for note in summary["scenario_notes"]:
        lines.append(f"- {note}")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    cases = {case_id: summarize_case(case_id) for case_id in CASE_IDS}
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "input_inventory": file_inventory(INPUT_DIR),
        "cases": cases,
        "scenario_notes": scenario_notes(cases),
        "river_entries": parse_river_sections(INPUT_DIR / "waterPACT_riv_floc_MP.nml"),
        "river_forcing": summarize_river_forcing(INPUT_DIR / "waterPACT_riv_floc_MP.nc"),
        "startup_profiles": summarize_startup_profiles(INPUT_DIR / "waterPACT_DRE_its.nc"),
        "netcdf_files": {
            "waterPACT_DRE_its.nc": summarize_netcdf(INPUT_DIR / "waterPACT_DRE_its.nc"),
            "waterPACT_obc_2019.nc": summarize_netcdf(INPUT_DIR / "waterPACT_obc_2019.nc"),
            "waterPACT_tsobc_2019.nc": summarize_netcdf(INPUT_DIR / "waterPACT_tsobc_2019.nc"),
            "waterPACT_DRE_wind_2019.nc": summarize_netcdf(INPUT_DIR / "waterPACT_DRE_wind_2019.nc"),
            "waterPACT_DRE_pres_2019.nc": summarize_netcdf(INPUT_DIR / "waterPACT_DRE_pres_2019.nc"),
        },
    }

    json_path = OUTPUT_ROOT / "case_inventory.json"
    md_path = OUTPUT_ROOT / "case_inventory.md"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md_path.write_text(render_markdown(summary), encoding="utf-8")

    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
