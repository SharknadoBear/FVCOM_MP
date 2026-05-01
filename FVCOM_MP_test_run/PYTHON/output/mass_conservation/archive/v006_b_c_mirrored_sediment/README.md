# Mass Conservation Output Archive v006

Archive label: `v006_b_c_mirrored_sediment`

Created: 2026-05-01

Purpose: Preserve the mass-conservation analysis after rerunning cases `b` and
`c` with the same sediment input used by case `a`, while keeping the present
source build with sediment `adv_scal_backward` and `fct_sed` bypassed in
`mod_sed.F`.

Context:
- Cases `a`, `b`, and `c` all use `SEDIMENT_MODEL = T`.
- Cases `a`, `b`, and `c` all point to `mp_sediment_a_no_adv.inp`.
- Case `a` keeps `PLAST_FLOC = F`.
- Cases `b` and `c` keep `PLAST_FLOC = T`.
- The sediment input starts active sediment immediately with `SED_START = 1`,
  but the code build bypasses sediment horizontal advection/FCT.

Key summary from `mass_summary.csv`:
- Case `a` final total: 994.200315 kg.
- Case `b` final total: 347.563342 kg.
- Case `c` final total: 348.522527 kg.

Interpretation:
- Cases `b` and `c` now evolve almost identically, confirming that the old
  b/c difference was controlled by sediment-input differences rather than by
  separate plastic behavior.
- Cases `b` and `c` track case `a` for roughly the first three days, then both
  enter an oscillatory lower-mass branch while case `a` continues to rise
  smoothly.
- The divergence is now localized to the plastic floc/split branch activated in
  b/c, not to different sediment input settings.

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
- `run_config/waterPACT_a_run.nml`
- `run_config/waterPACT_b_run.nml`
- `run_config/waterPACT_c_run.nml`
- `run_config/mp_sediment_a_no_adv.inp`
