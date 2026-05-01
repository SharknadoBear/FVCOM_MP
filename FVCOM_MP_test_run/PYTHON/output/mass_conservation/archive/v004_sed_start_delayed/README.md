# Mass Conservation Output Archive v004

Archive label: `v004_sed_start_delayed`

Created: 2026-04-30

Purpose: Preserve the completed 10-day mass-conservation analysis after rerunning
case `a` with sediment initialized but active `Advance_Sed` delayed beyond the
run length.

Context:
- Case `a`: `PLAST_FLOC = F`, `SEDIMENT_MODEL = T`, sediment input
  `generic_sediment_a_inactive.inp`.
- The sediment input sets `SED_START = 999999999`, so `Advance_Sed` returns at
  its startup guard and does not execute active sediment transport.
- Cases `b` and `c` are unchanged from the previous comparison.
- This snapshot separates sediment setup/compile/forcing effects from active
  runtime sediment transport.

Key summary from `mass_summary.csv`:
- Case `a` final total: 1001.251815 kg.
- Case `b` final total: 371.850952 kg.
- Case `c` final total: 528.551468 kg.

Interpretation:
- Case `a` returns to the old `v001`-like mass level when active
  `Advance_Sed` is delayed.
- This strongly localizes the sediment/plastic interference to code executed
  inside active `Advance_Sed`, rather than sediment setup, namelist presence, or
  the sediment-enabled river forcing path.

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
