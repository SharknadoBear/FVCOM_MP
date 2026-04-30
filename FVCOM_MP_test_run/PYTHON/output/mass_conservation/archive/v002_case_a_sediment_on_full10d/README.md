# Mass Conservation Output Archive v002

Archive label: `v002_case_a_sediment_on_full10d`

Created: 2026-04-29

Purpose: Preserve the mass-conservation analysis outputs after rerunning case `a` with `SEDIMENT_MODEL = T` for the full 10-day comparison window.

Context:
- Case `a`: `PLAST_FLOC = F`, `SEDIMENT_MODEL = T`.
- Case `b`: `PLAST_FLOC = T`, null aggregation/disaggregation setup, `SEDIMENT_MODEL = T`.
- Case `c`: `PLAST_FLOC = T`, active floc interaction setup, `SEDIMENT_MODEL = T`.
- This snapshot isolates the remaining no-floc versus null-floc difference after removing the previous sediment-runtime asymmetry.

Key summary from `mass_summary.csv`:
- Case `a` final total: 582.406247 kg.
- Case `b` final total: 371.850952 kg.
- Case `c` final total: 528.551468 kg.

Archived files:
- `mass_summary.csv`
- `mass_timeseries_long.csv`
- `total_mass_by_case.png`
- `mass_change_by_case.png`
- `cap_threshold_counts.png`
- `max_concentration_by_case.png`
- `water_bed_partition.png`
- `deposition_erosion_exchange.png`
- `floc_split_components_b_c.png`
- `split_consistency_b_c.png`
