# FVCOM-MP Production Run Project Handoff

This handoff is for a future agent that will create clean production-run
workspaces away from the current code-development/debugging repository.

Current repo root used as the template source:

```text
c:\Users\huan111\OneDrive - PNNL\Desktop\WaterPACT_Local\1_Model_Build\Model_develop\FVCOM_source_repo_github
```

Primary template folder inside that repo:

```text
FVCOM_MP_test_run/
```

Do not treat this current folder as the production workspace. It contains
debugging history, memo/report artifacts, archive outputs, and code-development
state. The new production folders should be created as separate sibling or
external project folders.

## Kestrel And Credential Rules

Use Huan's Kestrel skill/context when connecting to NREL Kestrel.

- Kestrel username: `yhuang168`
- Host: `kestrel.nrel.gov`
- Required SSH form:

```bash
ssh -m hmac-sha2-256 yhuang168@kestrel.nrel.gov
```

- Use the same MAC setting for `scp` and `rsync`:

```bash
scp -o MACs=hmac-sha2-256 local_file yhuang168@kestrel.nrel.gov:/remote/path/
rsync -av --progress -e "ssh -m hmac-sha2-256" local_dir/ yhuang168@kestrel.nrel.gov:/remote/path/
```

- Never write, store, echo, commit, or log Huan's password or OTP.
- Ask for the current OTP only when ready to start a new authenticated session.
- Before any `sbatch`, inspect the Slurm script for account, partition, time,
  node count, module stack, executable path, run directory, and output/error
  paths.

## Current Model State The New Agent Should Know

The current Delaware FVCOM-MP debugging work has converged on the following
configuration state.

### Case Meanings

Current test-run templates live in:

```text
FVCOM_MP_test_run/RUN_a/waterPACT_a_run.nml
FVCOM_MP_test_run/RUN_b/waterPACT_b_run.nml
FVCOM_MP_test_run/RUN_c/waterPACT_c_run.nml
```

The present useful interpretation is:

| Case | Role | Run namelist sediment model | Plastic file | Floc behavior |
|---|---|---|---|---|
| `a` | MP-only/no-floc control | `SEDIMENT_MODEL = T`, `SEDIMENT_MODEL_FILE = 'mp_sediment_a_no_adv.inp'` | `generic_plastic_a.inp` | `PLAST_FLOC = F` |
| `b` | Coupled fine-sediment case | `SEDIMENT_MODEL = T`, `SEDIMENT_MODEL_FILE = 'mp_sediment_a_no_adv.inp'` | `generic_plastic_b.inp` | `PLAST_FLOC = T`, but aggregation gate excludes the 100 um MP |
| `c` | Coupled coarse-sediment case | `SEDIMENT_MODEL = T`, `SEDIMENT_MODEL_FILE = 'generic_sediment_c.inp'` | `generic_plastic_c.inp` | `PLAST_FLOC = T`, aggregation active |

Important nuance: the published-style `a` case is conceptually "MP only /
no floc", but the current debug template still initializes the sediment model
through `mp_sediment_a_no_adv.inp`. Confirm whether the production MP-only run
should keep this exact tested path or truly set `SEDIMENT_MODEL = F`.

### Sediment And Plastic Inputs

Template input files live in:

```text
FVCOM_MP_test_run/INPUT/
```

Key files:

- `generic_plastic_a.inp`
- `generic_plastic_b.inp`
- `generic_plastic_c.inp`
- `mp_sediment_a_no_adv.inp`
- `generic_sediment.inp`
- `generic_sediment_c.inp`

Current MP setup:

- MP species: `mp1`
- Shape: sphere
- Axes: `PLAST_AAXI = PLAST_BAXI = PLAST_CAXI = 0.100` mm
- Density: `PLAST_PRHO = 950` kg/m3
- `SOURCE_AG_RATE = 0.0`
- `PLAST_FLOC = F` in case `a`; `T` in cases `b` and `c`

Current sediment setup:

- Sediment species: `coarse_sand_1`
- Fine sediment case: `SED_SD50 = 0.1` mm
- Coarse sediment case: `SED_SD50 = 1.0` mm
- `SED_WSET = 3` mm/s in both fine and coarse sediment files

Model behavior established by memos:

- For case `b`, D50 = 0.1 mm gives an aggregation gate of about 34 um, so the
  100 um MP is excluded from aggregation. Cases `a` and `b` are numerically
  identical in the current diagnostics.
- For case `c`, D50 = 1.0 mm gives a gate of about 162 um, so 100 um MP can
  aggregate. Case `c` shows active aggregation and bed deposition.
- Free/disaggregated MP is buoyant in seawater. A code-formula estimate at
  `rho_w ~= 1025 kg/m3` is about `-0.39 mm/s`.
- Aggregated MP uses the sediment settling branch, currently `+3.0 mm/s`.
  The branch difference is therefore about `3.4 mm/s`.

### Known Fixes And Configuration Cautions

Carry forward these lessons from the debug/memo phase:

- Fix B+ (`isnan()` guards on erosion/deposition fluxes) resolved the major
  mass-conservation bug from NaN erosion/deposition fluxes.
- `adv_scal` is safe in the tested path, but `adv_scal_backward` / FCT sediment
  paths previously corrupted shared flux arrays. Do not re-enable backward/FCT
  without a dedicated verification run.
- Case `c` was previously misconfigured to use the same D50 = 0.1 mm sediment
  as case `b`. Production coupled cases must explicitly point each run namelist
  to the intended sediment file.
- With `SOURCE_AG_RATE = 0.0`, current split consistency is tight. If production
  changes `SOURCE_AG_RATE > 0`, revisit the multi-species flux-limiting issue
  before trusting long runs.
- Confirm `STARTUP_TS_TYPE = 'observed'` and `STARTUP_FILE = 'waterPACT_DRE_its.nc'`
  remain intentional.
- Confirm source ramp behavior. The v012 archive used a ramp-corrected analytic
  source treatment for the first 10-day mass check.
- Confirm OBC export diagnostics before interpreting total domain mass over
  multi-month production windows.

## Project 1: Delaware Production Run Workspace

Create a new local project folder outside this repo. Recommended local parent:

```text
c:\Users\huan111\OneDrive - PNNL\Desktop\WaterPACT_Local\2_Production_Runs\
```

Recommended project root:

```text
Delaware_FVCOM_MP_production/
```

Create this directory structure:

```text
Delaware_FVCOM_MP_production/
  README.md
  MANIFEST.md
  INPUT/
    common/
    case_mp_only/
    case_coupled_fine/
    case_coupled_coarse/
    case_coupled_sensitivity/
  RUN/
    mp_only/
    coupled_fine/
    coupled_coarse/
    coupled_sensitivity/
  OUTPUT/
    mp_only/
    coupled_fine/
    coupled_coarse/
    coupled_sensitivity/
  ANALYSIS/
    scripts/
    notebooks/
    figures/
    summaries/
```

Planned Delaware cases:

| Production case | Template to start from | Purpose |
|---|---|---|
| `mp_only` | `RUN_a`, `generic_plastic_a.inp`, `mp_sediment_a_no_adv.inp` | no-floc/no-aggregation baseline |
| `coupled_fine` | `RUN_b`, `generic_plastic_b.inp`, `mp_sediment_a_no_adv.inp` or `generic_sediment.inp` | coupled path with fine sediment; expected near-MP-only behavior unless physics changes |
| `coupled_coarse` | `RUN_c`, `generic_plastic_c.inp`, `generic_sediment_c.inp` | active floc/deposition case |
| `coupled_sensitivity` | duplicate `RUN_c` then vary a controlled parameter | sensitivity/calibration case, e.g. D50, `SED_WSET`, source strength, or `SOURCE_AG_RATE` after safety review |

Essential files to copy from `FVCOM_MP_test_run/INPUT/` into
`Delaware_FVCOM_MP_production/INPUT/common/`:

```text
waterPACT_grd.dat
waterPACT_dep.dat
waterPACT_cor.dat
waterPACT_sig.dat
waterPACT_spg.dat
waterPACT_obc.dat
waterPACT_obc_2019.nc
waterPACT_tsobc_2019.nc
waterPACT_DRE_wind_2019.nc
waterPACT_DRE_pres_2019.nc
waterPACT_DRE_its.nc
waterPACT_riv_floc_MP.nml
waterPACT_riv_floc_MP.nc
```

Essential model-configuration files to copy and then case-specialize:

```text
generic_plastic_a.inp
generic_plastic_b.inp
generic_plastic_c.inp
mp_sediment_a_no_adv.inp
generic_sediment.inp
generic_sediment_c.inp
```

Essential run templates to copy:

```text
FVCOM_MP_test_run/RUN_a/waterPACT_a_run.nml
FVCOM_MP_test_run/RUN_a/job_Kestrel.sh
FVCOM_MP_test_run/RUN_b/waterPACT_b_run.nml
FVCOM_MP_test_run/RUN_b/job_Kestrel.sh
FVCOM_MP_test_run/RUN_c/waterPACT_c_run.nml
FVCOM_MP_test_run/RUN_c/job_Kestrel.sh
```

After copying, rename the case namelists and casenames consistently. Example:

```text
RUN/mp_only/waterPACT_mp_only_run.nml
RUN/coupled_fine/waterPACT_coupled_fine_run.nml
RUN/coupled_coarse/waterPACT_coupled_coarse_run.nml
RUN/coupled_sensitivity/waterPACT_coupled_sensitivity_run.nml
```

Update every namelist's:

- `CASE_TITLE`
- `START_DATE`, `END_DATE`
- `INPUT_DIR`
- `OUTPUT_DIR`
- `RST_FIRST_OUT`, `RST_OUT_INTERVAL`, `RST_OUTPUT_STACK`
- `NC_FIRST_OUT`, `NC_OUT_INTERVAL`, `NC_OUTPUT_STACK`
- `SEDIMENT_MODEL`, `SEDIMENT_MODEL_FILE`
- `PLASTIC_MODEL`, `PLASTIC_MODEL_FILE`
- river, wind, pressure, open-boundary, startup T/S file names if paths change

Update every Slurm script's:

- `#SBATCH --job-name`
- `#SBATCH --account` and `--partition` if Kestrel policy changed
- `#SBATCH --time`
- `#SBATCH --nodes`
- `#SBATCH --ntasks-per-node`
- `#SBATCH -o` and `#SBATCH -e`
- module stack
- `srun ./fvcom_TLMF --casename=<case_name>`

Do not copy heavy debug output archives into the new production workspace.
Create empty `OUTPUT/<case>/` folders or placeholder `.gitkeep` files only.

Recommended analysis templates to copy into `ANALYSIS/scripts/`:

```text
FVCOM_MP_test_run/PYTHON/analyze_spatial_distribution.py
FVCOM_MP_test_run/MATLAB/extract_mp1_mass_budget_a.m
FVCOM_MP_test_run/MATLAB/extract_mp1_mass_budget_b_c.m
FVCOM_MP_test_run/MATLAB/extract_spatial_snapshots_a.m
FVCOM_MP_test_run/MATLAB/extract_spatial_snapshots_b_c.m
```

For production, adapt the analysis scripts to use the new project root and
case names instead of hard-coded `OUTPUT_a`, `OUTPUT_b`, `OUTPUT_c`.

## Project 2: Lab Flume Floc Calibration Workspace

Create a second new local project folder outside this repo. Recommended root:

```text
Lab_Flume_FVCOM_MP_floc_calibration/
```

Use the same high-level structure:

```text
Lab_Flume_FVCOM_MP_floc_calibration/
  README.md
  MANIFEST.md
  INPUT/
    common/
    case_mp_only/
    case_fine_sediment/
    case_coarse_sediment/
    case_calibration_sweep/
  RUN/
    mp_only/
    fine_sediment/
    coarse_sediment/
    calibration_sweep/
  OUTPUT/
    mp_only/
    fine_sediment/
    coarse_sediment/
    calibration_sweep/
  ANALYSIS/
    scripts/
    lab_data/
    figures/
    summaries/
```

No dedicated lab-flume template was found in the current repo. Start by copying
the Delaware run/config templates, then replace geometry, forcing, and
observation-specific files with flume-specific versions when they are available.

Initial lab-flume case plan:

| Flume case | Template to start from | Purpose |
|---|---|---|
| `mp_only` | Delaware `RUN_a` + `generic_plastic_a.inp` | control run for MP transport without floc interaction |
| `fine_sediment` | Delaware `RUN_b` + fine sediment config | check expected no/weak aggregation for D50 = 0.1 mm |
| `coarse_sediment` | Delaware `RUN_c` + coarse sediment config | active aggregation/deposition reference |
| `calibration_sweep` | duplicate `coarse_sediment` | tune against lab data; vary one parameter at a time |

Lab-flume files that must be supplied or generated later:

- flume grid, bathymetry, sigma, Coriolis handling, and sponge/open-boundary
  files appropriate for a laboratory domain
- inflow or source forcing for MP and sediment injection
- water-level/velocity boundary forcing if needed
- lab observation CSVs for suspended MP, aggregated MP, sediment concentration,
  bed deposition, or image-derived floc size/time series
- station/probe definitions matching the lab measurement locations

Configuration details to revisit for flume calibration:

- `START_DATE` / `END_DATE` should use the lab experiment time window, not the
  Delaware 2019 calendar forcing by default.
- Surface wind, air pressure, tide, and salinity/temperature open-boundary files
  may be irrelevant for the flume and should be turned off or replaced.
- The model may need much shorter output intervals than Delaware production.
- Calibration sweeps should record exact values of D50, `SED_WSET`, MP density,
  MP size, source concentration, source timing, and any aggregation/disaggregation
  controls changed from the template.
- Keep `SOURCE_AG_RATE = 0.0` until the agent confirms the intended physics and
  verifies mass conservation with the nonzero source aggregation path.

## Minimum Verification Before Submitting Production Jobs

For each new case in either project:

1. Run a file-presence check for all files named in the run namelist.
2. Diff the new namelist against its template and write a short case-specific
   note explaining every intentional change.
3. Confirm the output directory is empty or intentionally resumable.
4. Confirm the Slurm script points at the correct executable and casename.
5. On Kestrel, run a short smoke test first, not the full production walltime.
6. After smoke test, check:
   - FVCOM started and wrote NetCDF/restart output.
   - no immediate NaN/Inf messages.
   - mass budget is plausible.
   - split consistency is plausible for coupled cases.
   - no concentration cap artifacts.
7. Only then submit full production jobs.

## Suggested Handoff Deliverables For The New Agent

The future agent should produce:

- the two clean local project directories
- a `MANIFEST.md` in each project listing every copied template and source path
- a `README.md` in each project describing the case matrix
- case-specialized `RUN/*/*.nml` files
- case-specialized Slurm scripts
- copied/shared `INPUT/common` forcing/config files
- empty `OUTPUT/<case>/` folders
- copied and lightly adapted `ANALYSIS/scripts/`
- a short `SETUP_NOTES.md` documenting unresolved decisions before Kestrel upload

Do not submit Kestrel jobs until Huan explicitly confirms the final production
case matrix, run duration, Kestrel remote path, and resource settings.
