# Spatial Distribution v2: Day 10

This folder contains plan-view maps and statistics for the day-10 final-24-hour
average extracted from the FVCOM-MP case c/d NetCDF outputs.

## Inputs
- Case c: `C:\Users\huan111\OneDrive - PNNL\Desktop\WaterPACT_Local\1_Model_Build\Model_develop\FVCOM_source_repo_github\FVCOM_MP_test_run\OUTPUT_c\spatial_distribution_v2_c_day10.mat` (24.0 records averaged)
- Case d: `C:\Users\huan111\OneDrive - PNNL\Desktop\WaterPACT_Local\1_Model_Build\Model_develop\FVCOM_source_repo_github\FVCOM_MP_test_run\OUTPUT_d\spatial_distribution_v2_d_day10.mat` (24.0 records averaged)

Case `c` is the comparison run. Case `d` is the ORIG_SED effective-floc coupling
run. Difference figures and `diff_summary_stats.csv` use `d - c`.

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
