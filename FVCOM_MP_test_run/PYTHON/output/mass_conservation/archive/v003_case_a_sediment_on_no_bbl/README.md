# Mass Conservation Output Archive v003

Archive label: `v003_case_a_sediment_on_no_bbl`

Created: 2026-04-30

Purpose: Preserve the mass-conservation analysis outputs after rerunning case `a`
with `SEDIMENT_MODEL = T` and a case-a sediment input that disables all BBL
roughness closures.

Context:
- Case `a`: `PLAST_FLOC = F`, `SEDIMENT_MODEL = T`, sediment input
  `generic_sediment_a_no_bbl.inp`.
- Case `a` BBL switches were disabled: `MB_BBL_USE = F`, `SG_BBL_USE = F`,
  `SSW_BBL_USE = F`, and `SSW_CALC_ZNOT = F`.
- Cases `b` and `c` are unchanged from the previous comparison.
- This snapshot tests whether the sediment-on case-a mass reduction is mainly
  caused by sediment BBL roughness feedback.

Key summary from `mass_summary.csv`:
- Case `a` final total: 573.277505 kg.
- Case `b` final total: 371.850952 kg.
- Case `c` final total: 528.551468 kg.

Interpretation:
- The case-a final mass differs from the `v002` sediment-on/BBL-on case-a value
  of 582.406247 kg by only 9.128742 kg over the 10-day window.
- This small change demotes BBL roughness as the leading explanation for the
  sediment/plastic interference.

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
