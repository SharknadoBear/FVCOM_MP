"""Quick stats for the case-a mass-check CSV."""
import re

fn = 'waterPACT_a_mp_mass_check.csv'
rows = open(fn).readlines()
data = [l for l in rows if re.match(r'^\s+\d', l)]

vals = {}; tdays = {}
nan_count = 0
for l in data:
    p = l.split(',')
    iint = p[0].strip(); cp = p[2].strip()
    grand = p[10].strip(); t = p[1].strip()
    try:
        g = float(grand)
        if g != g:  # isnan
            nan_count += 1
        vals[(iint, cp)] = g
        tdays[iint] = float(t)
    except Exception:
        pass

all_iints = sorted(set(k[0] for k in vals), key=lambda x: int(x))
t0 = min(tdays.values())
tmax = max(tdays.values())
n_steps = len(all_iints)

print(f"Reporting steps: {n_steps}")
print(f"Julian day range: {t0:.4f} – {tmax:.4f}")
print(f"Duration: {tmax - t0:.4f} days")
print(f"NaN events (grand_total): {nan_count}")

# Check all checkpoints present
CPS = ['CP0_start','CP1_after_deposition','CP2_after_erosion',
       'CP3_after_update_bottom','CP3p5_after_adv_a','CP4_after_advection',
       'CP5_after_vdif','CP6_after_kinetics','CP7_after_upper_clamp',
       'CP8_after_OBC','CP9_after_PTsource','CP10_after_neg_clamp',
       'CP11_after_MPI','CP12_end']
cp_counts = {cp: sum(1 for k in vals if k[1]==cp) for cp in CPS}
print("\nCheckpoint row counts:")
for cp, n in cp_counts.items():
    print(f"  {cp}: {n}")

# CP0->CP12 drift
deltas = []
for iint in all_iints:
    v0 = vals.get((iint,'CP0_start'))
    v12 = vals.get((iint,'CP12_end'))
    if v0 is not None and v12 is not None:
        deltas.append((tdays[iint]-t0, v12-v0, v0, iint))

if deltas:
    nets = [d[1] for d in deltas]
    cumulative = sum(nets)
    final_mass = deltas[-1][2]
    print(f"\nCP0->CP12 per-step: mean={sum(nets)/len(nets):+.4e} kg, "
          f"min={min(nets):+.4e}, max={max(nets):+.4e}")
    print(f"CP0->CP12 cumulative drift: {cumulative:+.4f} kg "
          f"({cumulative/final_mass*100:+.4f}%)")
    print(f"Final domain mass (CP0): {final_mass:.2f} kg")

# Per-stage mean absolute delta
print("\nPer-stage mean |delta| (kg):")
CP_PAIRS = list(zip(CPS[:-1], CPS[1:]))
for c0, c1 in CP_PAIRS:
    ds = []
    for iint in all_iints:
        v0 = vals.get((iint,c0)); v1 = vals.get((iint,c1))
        if v0 is not None and v1 is not None:
            ds.append(abs(v1-v0))
    if ds:
        nz = sum(1 for d in ds if d > 1e-12)
        print(f"  {c0[:8]}→{c1[:8]}: mean={sum(ds)/len(ds):.4e}  nonzero={nz}/{len(ds)}")
