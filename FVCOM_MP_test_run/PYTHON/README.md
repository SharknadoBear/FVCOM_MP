# Delaware Test Diagnostics

This folder holds the reusable Python tooling for the Delaware River FVCOM-MP
test cases in `../INPUT`, `../RUN_a`, `../RUN_b`, and `../RUN_c`.

## Files

- `inspect_delaware_cases.py`: inspects the run namelists, plastic/sediment
  input files, and shared NetCDF forcing files; writes machine-readable and
  Markdown summaries under `output/`.
- `analyze_mp1_mass_conservation.py`: reads the MATLAB mass-budget MAT files
  in `../OUTPUT_a`, `../OUTPUT_b`, and `../OUTPUT_c`, then writes comparison
  plots and CSV summaries under `output/mass_conservation/`.
- `setup_env.ps1`: creates a local virtual environment in `.venv` while still
  exposing the host scientific packages through `--system-site-packages`.
- `requirements.txt`: minimal package list if a clean environment needs to be
  recreated elsewhere.

## Suggested Usage

```powershell
cd FVCOM_MP_test_run\PYTHON
.\setup_env.ps1
.\.venv\Scripts\python.exe .\inspect_delaware_cases.py
```

If the host machine already has the needed packages, the script can also be run
directly with:

```powershell
python .\inspect_delaware_cases.py
```

Mass-budget comparison plots can be regenerated with:

```powershell
python .\analyze_mp1_mass_conservation.py
```

The generated summaries are:

- `output/case_inventory.json`
- `output/case_inventory.md`
- `output/mass_conservation/mass_summary.csv`
- `output/mass_conservation/mass_timeseries_long.csv`
- `output/mass_conservation/*.png`
