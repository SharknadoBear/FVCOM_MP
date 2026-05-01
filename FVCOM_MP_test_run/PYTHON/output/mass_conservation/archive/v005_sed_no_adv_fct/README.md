# Mass Conservation Output Archive v005

Archive label: `v005_sed_no_adv_fct`

Created: 2026-05-01

Purpose: Preserve the completed 10-day mass-conservation analysis after rerunning
case `a` with active sediment runtime enabled, but with the horizontal sediment
`adv_scal_backward` call and `fct_sed` limiter bypassed in `mod_sed.F`.

Context:
- Case `a`: `PLAST_FLOC = F`, `SEDIMENT_MODEL = T`, sediment input
  `mp_sediment_a_no_adv.inp`.
- The sediment input starts active sediment immediately with `SED_START = 1`.
- The code build keeps the rest of `Advance_Sed` active but bypasses the
  sediment horizontal scalar advection/FCT block.
- Cases `b` and `c` are unchanged from the previous comparison.

Key summary from `mass_summary.csv`:
- Case `a` final total: 994.200315 kg.
- Case `b` final total: 371.850952 kg.
- Case `c` final total: 528.551468 kg.

Interpretation:
- Case `a` returns to the `v001`/`v004`-like mass level when sediment runtime
  is active but the sediment horizontal advection/FCT block is bypassed.
- This is strong evidence that the unwanted sediment/plastic interaction is
  caused by the active sediment `adv_scal_backward`/`fct_sed` stage or by shared
  state modified there before `Advance_Plast` runs.

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
- `source_mass_budget_mat/mp1_mass_budget_a.mat`
- `source_mass_budget_mat/mp1_mass_budget_b.mat`
- `source_mass_budget_mat/mp1_mass_budget_c.mat`
