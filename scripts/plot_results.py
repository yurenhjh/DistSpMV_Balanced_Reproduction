#!/usr/bin/env python3
"""Plot cant ablation medians: Balanced vs Naive Comm Time.

Reads `data/cant_ablation_medians.csv` (preferred) or falls back to
`data/results_perf.csv` to compute medians, then saves
`data/ablation_trend.png`.
"""
import csv
import os
import re
import statistics
from collections import defaultdict


def read_medians_file(path, np_target=4):
    mapping = defaultdict(dict)  # mapping[round(frac,3)][mode] = comm
    if not os.path.exists(path):
        return mapping
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        # expected header: frac,np,mode,...,comm_s_median,...
        for r in reader:
            frac_s = r.get('frac','').strip()
            if frac_s == '':
                continue
            try:
                frac = float(frac_s.replace('p','.'))
            except Exception:
                try:
                    frac = float(frac_s)
                except Exception:
                    continue
            try:
                npv = int(r.get('np','0') or 0)
            except Exception:
                npv = 0
            mode = (r.get('mode') or '').lower()
            # comm column name in medians file
            comm_val = r.get('comm_s_median') or r.get('comm_s') or r.get('comm')
            if comm_val is None or comm_val == '':
                continue
            try:
                comm = float(comm_val)
            except Exception:
                continue
            if npv == np_target and mode in ('balanced','naive'):
                mapping[round(frac,3)][mode] = comm
    return mapping


def compute_medians_from_results(path, np_target=4):
    temp = defaultdict(lambda: defaultdict(list))
    if not os.path.exists(path):
        return {}
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for r in reader:
            mname = r.get('Matrix','')
            frac = None
            frac_s = (r.get('lower_bound_frac') or '').strip()
            if frac_s:
                try:
                    frac = float(frac_s)
                except Exception:
                    frac = None
            if frac is None:
                mo = re.search(r'frac([0-9p\.]+)', mname)
                if mo:
                    frac = float(mo.group(1).replace('p','.'))
            if frac is None:
                continue
            try:
                npv = int(r.get('np','0') or 0)
            except Exception:
                npv = 0
            mode = (r.get('mode') or '').lower()
            try:
                comm = float(r.get('comm_s') or r.get('comm') or '')
            except Exception:
                continue
            if npv == np_target and mode in ('balanced','naive'):
                temp[round(frac,3)][mode].append(comm)
    out = defaultdict(dict)
    for frac_k, d in temp.items():
        for mode, arr in d.items():
            if arr:
                out[frac_k][mode] = statistics.median(sorted(arr))
    return out


def main():
    medians_path = os.path.join('data','cant_ablation_medians.csv')
    results_path = os.path.join('data','results_perf.csv')
    desired_fracs = [0.5, 0.8, 1.0]
    np_target = 4

    data_map = read_medians_file(medians_path, np_target=np_target)
    # If any desired frac is missing, fallback to recomputing from results_perf
    missing = [f for f in desired_fracs if round(f,3) not in data_map]
    if missing:
        recomputed = compute_medians_from_results(results_path, np_target=np_target)
        # merge
        for k,v in recomputed.items():
            data_map.setdefault(k, {}).update(v)

    # build plotting arrays
    xs = []
    ys_bal = []
    ys_naive = []
    for f in desired_fracs:
        key = round(f,3)
        xs.append(f)
        entry = data_map.get(key, {})
        ys_bal.append(entry.get('balanced', float('nan')))
        ys_naive.append(entry.get('naive', float('nan')))

    # plotting
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception as e:
        print('matplotlib is required to run this script:', e)
        return 2

    plt.figure(figsize=(6,4))
    plt.plot(xs, ys_bal, marker='o', linestyle='-', color='#1f77b4', label='Balanced')
    plt.plot(xs, ys_naive, marker='s', linestyle='--', color='#ff7f0e', label='Naive')
    plt.xticks(xs)
    plt.xlabel('lower_bound_frac')
    plt.ylabel('Comm Time (s)')
    plt.title('cant: Threshold Ablation (Comm Time median, np=4)')
    plt.grid(alpha=0.3)
    plt.legend()
    out_dir = 'data'
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir,'ablation_trend.png')
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    print('Wrote', out_path)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
