import re
import sys

files = sys.argv[1:]
patterns = re.compile(r'air_pressure|AIRPRESS|PATM|patm|baro|101325|P_ATM|AIRPRESSURE_ON', re.IGNORECASE)

for fpath in files:
    try:
        with open(fpath, encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
        hits = [(i+1, l.rstrip()) for i,l in enumerate(lines) if patterns.search(l)]
        if hits:
            print(f"\n=== {fpath} ===")
            for ln, txt in hits:
                print(f"  {ln}: {txt}")
        else:
            print(f"\n=== {fpath} === NO MATCHES")
    except Exception as e:
        print(f"ERROR reading {fpath}: {e}")
