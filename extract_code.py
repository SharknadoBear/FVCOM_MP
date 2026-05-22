import os
import re

base = r"C:\Users\huan111\OneDrive - PNNL\Desktop\WaterPACT_Local\1_Model_Build\Model_develop\FVCOM_source_repo_github"

def show_lines(fpath, start, end):
    if not os.path.exists(fpath):
        print(f"File not found: {fpath}")
        return
    with open(fpath, encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    print(f"\n--- {os.path.basename(fpath)}: lines {start} to {end} ---")
    for i in range(start-1, min(end, len(lines))):
        print(f"{i+1:5d}: {lines[i]}", end='')

def find_patterns(fpath, pattern, context_before=5, context_after=30, max_matches=1):
    if not os.path.exists(fpath):
        print(f"File not found: {fpath}")
        return
    import re
    with open(fpath, encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    pat = re.compile(pattern, re.IGNORECASE)
    matches = 0
    for i, l in enumerate(lines):
        if pat.search(l):
            start = max(0, i - context_before)
            end = min(len(lines), i + context_after + 1)
            print(f"\n--- match pattern '{pattern}' at {os.path.basename(fpath)} line {i+1} ---")
            for j in range(start, end):
                print(f"{j+1:5d}: {lines[j]}", end='')
            matches += 1
            if matches >= max_matches:
                break

print("=== START EXTRACTION ===")
show_lines(os.path.join(base, "external_step.F"), 545, 590)

print("\n=== mod_force.F: AIR_PRESSURE_STORE ===")
find_patterns(os.path.join(base, "mod_force.F"), r"SUBROUTINE\s+AIR_PRESSURE_STORE", context_before=2, context_after=45)

print("\n=== mod_force.F: PA declaration and assignment ===")
find_patterns(os.path.join(base, "mod_force.F"), r"\bPA\s*=", context_before=2, context_after=5, max_matches=3)
find_patterns(os.path.join(base, "mod_force.F"), r"REAL.*::.*PA", context_before=2, context_after=5, max_matches=2)

print("\n=== mod_input.F: reading air_pressure from NC ===")
find_patterns(os.path.join(base, "mod_input.F"), r"AIRPRESSURE_FILE|air_pressure", context_before=10, context_after=30, max_matches=3)
