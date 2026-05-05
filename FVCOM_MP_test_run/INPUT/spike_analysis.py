import re

fn = 'waterPACT_b_mp_mass_check.csv'
rows = open(fn).readlines()
data = [l for l in rows if re.match(r'^\s+\d', l)]
vals = {}; tdays = {}
for l in data:
    p = l.split(',')
    iint = p[0].strip(); cp = p[2].strip()
    grand = p[10].strip(); t = p[1].strip()
    try:
        vals[(iint, cp)] = float(grand)
        tdays[iint] = float(t)
    except Exception:
        pass

CP_NEED = ['CP0_start', 'CP3p5_after_adv_a', 'CP4_after_advection', 'CP12_end']
all_iints = sorted(set(k[0] for k in vals), key=lambda x: int(x))
t0 = min(tdays.values())
deltas = []
for iint in all_iints:
    row = {cp: vals.get((iint, cp)) for cp in CP_NEED}
    if any(v is None for v in row.values()):
        continue
    net = row['CP12_end'] - row['CP0_start']
    adv_d = row['CP4_after_advection'] - row['CP3p5_after_adv_a']
    deltas.append((tdays[iint] - t0, net, adv_d, iint))

dbn = sorted(deltas, key=lambda x: x[1])
print('5 largest per-step losses:')
for t, net, adv_d, i in dbn[:5]:
    print(f'  day={t:.4f}  net={net:+.4e}  adv_d={adv_d:+.4e}  iint={i}')
print()
print('5 largest per-step gains:')
for t, net, adv_d, i in reversed(dbn[-5:]):
    print(f'  day={t:.4f}  net={net:+.4e}  adv_d={adv_d:+.4e}  iint={i}')
print()

early = [abs(d[2]) for d in deltas if d[0] < 5]
late = [abs(d[2]) for d in deltas if d[0] >= 15]
print(f'Mean |adv_d| day 0-5:   {sum(early)/len(early):.4e} kg')
print(f'Mean |adv_d| day 15-20: {sum(late)/len(late):.4e} kg')
print(f'Amplitude growth factor: {sum(late)/len(late)/(sum(early)/len(early)):.1f}x')
print()

CP0_KEY = 'CP0_start'
for day in [0, 2, 5, 10, 15, 19]:
    c = min(deltas, key=lambda x: abs(x[0] - day))
    mass = vals.get((c[3], CP0_KEY), 0.0)
    print(f'  CP0 at day {day:2d}: {mass:.2f} kg')
