import re

for case in ['b', 'c']:
    fn = f'waterPACT_{case}_mp_mass_check.csv'
    rows = open(fn).readlines()

    data = [l for l in rows if re.match(r'^\s+\d', l)]
    cflx = [l for l in rows if l.startswith('CFLX')]

    nan_by_cp = {}
    tot_by_cp = {}
    grand_cp0 = []
    grand_cp12 = []
    grand_cp3 = []
    grand_cp3p5 = []
    grand_cp4 = []

    for l in data:
        parts = l.split(',')
        cp = parts[2].strip()
        has_nan = 'NaN' in l
        nan_by_cp[cp] = nan_by_cp.get(cp, 0) + (1 if has_nan else 0)
        tot_by_cp[cp] = tot_by_cp.get(cp, 0) + 1
        if not has_nan:
            val = float(parts[10].strip())
            if cp == 'CP0_start':               grand_cp0.append(val)
            if cp == 'CP12_end':                grand_cp12.append(val)
            if cp == 'CP3_after_update_bottom': grand_cp3.append(val)
            if cp == 'CP3p5_after_adv_a':       grand_cp3p5.append(val)
            if cp == 'CP4_after_advection':     grand_cp4.append(val)

    steps = tot_by_cp.get('CP0_start', 0)
    print(f'\n=== Case {case} ===')
    print(f'  Steps: {steps}')
    print(f'  NaN by checkpoint:')
    for cp in ['CP2_after_erosion', 'CP3_after_update_bottom', 'CP3p5_after_adv_a',
               'CP4_after_advection', 'CP5_after_vdif', 'CP10_after_neg_clamp']:
        n = nan_by_cp.get(cp, 0)
        t = tot_by_cp.get(cp, 0)
        marker = '  <-- FIXED' if cp == 'CP2_after_erosion' and n == 0 else ''
        print(f'    {cp:42s}  {n:4d}/{t}  ({100*n/max(t,1):.1f}%){marker}')

    if grand_cp0:
        print(f'  CP0  grand_total: min={min(grand_cp0):.4e}  max={max(grand_cp0):.4e}')
    if grand_cp12:
        print(f'  CP12 grand_total: min={min(grand_cp12):.4e}  max={max(grand_cp12):.4e}')

    # adv_a delta: CP3 -> CP3p5 (must be 0 when conc_a=0)
    if grand_cp3 and grand_cp3p5 and len(grand_cp3) == len(grand_cp3p5):
        deltas_a = [abs(a - b) for a, b in zip(grand_cp3p5, grand_cp3)]
        print(f'  adv_a delta (CP3->CP3p5): '
              f'max={max(deltas_a):.4e}  mean={sum(deltas_a)/len(deltas_a):.4e}'
              + ('  <-- ZERO: adv_a clean' if max(deltas_a) == 0.0 else ''))
    else:
        print(f'  adv_a delta: cannot compute (CP3 steps={len(grand_cp3)}, CP3p5 steps={len(grand_cp3p5)})')

    # adv_d delta: CP3p5 -> CP4 (tidal OBC signal)
    if grand_cp3p5 and grand_cp4 and len(grand_cp3p5) == len(grand_cp4):
        deltas_d = [a - b for a, b in zip(grand_cp4, grand_cp3p5)]
        print(f'  adv_d delta (CP3p5->CP4): '
              f'min={min(deltas_d):.4e}  max={max(deltas_d):.4e}')

    nz_cflx = sum(1 for l in cflx if '0.000000E+00' not in l.split('=')[-1])
    print(f'  CFLX_TRACE non-zero: {nz_cflx}/{len(cflx)}')

    # NaN temporal trend at CP2 (10 bins)
    cp2rows = [l for l in data if 'CP2_after_erosion' in l]
    n_bins = 5
    bin_size = max(1, len(cp2rows) // n_bins)
    print(f'  CP2 NaN trend ({n_bins} bins):')
    for i in range(n_bins):
        chunk = cp2rows[i*bin_size:(i+1)*bin_size]
        nn = sum(1 for l in chunk if 'NaN' in l)
        print(f'    bin {i+1}: {nn}/{len(chunk)} ({100*nn/max(len(chunk),1):.1f}%)')
