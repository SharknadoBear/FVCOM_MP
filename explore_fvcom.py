import re, os

base = r"C:\Users\huan111\OneDrive - PNNL\Desktop\WaterPACT_Local\1_Model_Build\Model_develop\FVCOM_source_repo_github"

def read_around(fname, pattern, before=5, after=35):
    path = os.path.join(base, fname)
    if not os.path.exists(path):
        print(f"File not found: {path}")
        return
    with open(path, encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    p = re.compile(pattern, re.IGNORECASE)
    found = False
    for i, l in enumerate(lines):
        if p.search(l):
            s = max(0, i - before)
            e = min(len(lines), i + after)
            print(f"\n=== {fname} (lines {s+1}-{e}) ===")
            for j in range(s, e):
                print(f"{j+1:5d}  {lines[j]}", end='')
            found = True
            # For some, we might want to continue searching, but here we break to match requested logic
            break
    if not found:
        print(f"Pattern '{pattern}' not found in {fname}")

# 1. external_step.F - Look for any number starting with 1013 (often 101325 or 1013.25)
read_around("external_step.F", r"1013", before=10, after=20)

# 2. mod_force.F - SURFACE_AIRPRESSURE or UPDATE_AIRPRESSURE
read_around("mod_force.F", r"SUBROUTINE SURFACE_AIRPRESSURE", before=2, after=70)

# 3. mod_input.F - Search for NetCDF reading logic
read_around("mod_input.F", r"air_pressure", before=5, after=40)
