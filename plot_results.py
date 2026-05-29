#!/usr/bin/env python3
"""
plot_results.py
读取 CSV（results_main.csv、results_ablation.csv、results_balance.csv），生成 figures/ 下 PNG。
"""
import os, csv, math
import matplotlib.pyplot as plt
from collections import defaultdict

ROOT = os.path.abspath(os.path.dirname(__file__))
MAIN_CSV = os.path.join(ROOT, 'results_main.csv')
ABL_CSV  = os.path.join(ROOT, 'results_ablation.csv')
BAL_CSV  = os.path.join(ROOT, 'results_balance.csv')
FIG_DIR  = os.path.join(ROOT, 'figures')

os.makedirs(FIG_DIR, exist_ok=True)

def read_csv_rows(path):
    if not os.path.exists(path):
        print(f"Missing {path}, skipping.")
        return []
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        return list(reader)

main_rows = read_csv_rows(MAIN_CSV)
abl_rows  = read_csv_rows(ABL_CSV)
bal_rows  = read_csv_rows(BAL_CSV)

# Normalize key for matrix (strip path)
def basename(mat):
    return os.path.basename(mat)

# Build lookup from main
main = {}
for r in main_rows:
    key = (basename(r['Matrix']), r['Algorithm'], r['np'])
    main[key] = r

# Figure 1: GFlops comparison (np=1 and np=2)
mat_order = ['cant.mtx', 'ecology1.mtx', 'bcsstk30.mtx']
alg_order = ['naive', 'metis_naive', 'balanced']
np_list = ['1','2']

for idx, npv in enumerate(np_list):
    fig, ax = plt.subplots(figsize=(8,4))
    x = range(len(mat_order))
    width = 0.2
    for i, alg in enumerate(alg_order):
        vals = []
        for m in mat_order:
            key = (m, alg, npv)
            v = main.get(key)
            if v and v.get('GFlops') and v.get('GFlops')!='NA':
                vals.append(float(v['GFlops']))
            else:
                vals.append(0.0)
        positions = [p + (i-1)*width for p in x]
        ax.bar(positions, vals, width=width, label=alg)
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace('.mtx','') for m in mat_order])
    ax.set_ylabel('GFlops')
    ax.set_title(f'Performance Comparison (np={npv})')
    ax.legend()
    plt.tight_layout()
    figpath = os.path.join(FIG_DIR, f'fig1_gflops_comparison_np{npv}.png')
    fig.savefig(figpath)
    print('Saved', figpath)
    plt.close(fig)

# Combined two-subplot version
fig, axes = plt.subplots(1,2, figsize=(12,4), sharey=True)
for ax, npv in zip(axes, np_list):
    x = range(len(mat_order))
    width = 0.2
    for i, alg in enumerate(alg_order):
        vals = []
        for m in mat_order:
            key = (m, alg, npv)
            v = main.get(key)
            vals.append(float(v['GFlops']) if v and v.get('GFlops') and v.get('GFlops')!='NA' else 0.0)
        positions = [p + (i-1)*width for p in x]
        ax.bar(positions, vals, width=width, label=alg)
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace('.mtx','') for m in mat_order])
    ax.set_title(f'np={npv}')
axes[0].set_ylabel('GFlops')
fig.suptitle('Performance Comparison (np=1 and np=2)')
fig.legend(alg_order, loc='upper right')
fig.tight_layout(rect=[0,0,1,0.95])
figpath = os.path.join(FIG_DIR, 'fig1_gflops_comparison.png')
fig.savefig(figpath)
print('Saved', figpath)
plt.close(fig)

# Figure 2: Communication time comparison (two subplots np=1/2)
fig, axes = plt.subplots(1,2, figsize=(12,4), sharey=True)
for ax, npv in zip(axes, np_list):
    x = range(len(mat_order))
    width = 0.2
    for i, alg in enumerate(alg_order):
        vals = []
        for m in mat_order:
            key = (m, alg, npv)
            v = main.get(key)
            vals.append(float(v['CommTime']) if v and v.get('CommTime') and v.get('CommTime')!='NA' else 0.0)
        positions = [p + (i-1)*width for p in x]
        ax.bar(positions, vals, width=width, label=alg)
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace('.mtx','') for m in mat_order])
    ax.set_title(f'np={npv}')
axes[0].set_ylabel('CommTime (s)')
fig.suptitle('Communication Time Comparison')
fig.tight_layout(rect=[0,0,1,0.95])
figpath = os.path.join(FIG_DIR, 'fig2_comm_time.png')
fig.savefig(figpath)
print('Saved', figpath)
plt.close(fig)

# Figure 3: Speedup of Balanced over Baselines (np=2)
# Speedup = baseline_total_time / balanced_total_time
pairs = []
mat_labels = [m.replace('.mtx','') for m in mat_order]
speedup_vs_naive = []
speedup_vs_metis = []
for m in mat_order:
    t_naive = main.get((m,'naive','2'))
    t_metis  = main.get((m,'metis_naive','2'))
    t_bal    = main.get((m,'balanced','2'))
    def get_total(r):
        if not r or r.get('TotalTime') in (None,'NA',''):
            return None
        try:
            return float(r['TotalTime'])
        except:
            return None
    tn = get_total(t_naive)
    tm = get_total(t_metis)
    tb = get_total(t_bal)
    sn = (tn/tb) if (tn and tb and tb>0) else 0.0
    sm = (tm/tb) if (tm and tb and tb>0) else 0.0
    speedup_vs_naive.append(sn)
    speedup_vs_metis.append(sm)

x = range(len(mat_order))
width = 0.35
fig, ax = plt.subplots(figsize=(8,4))
ax.bar([p-width/2 for p in x], speedup_vs_naive, width=width, label='Balanced/Naive')
ax.bar([p+width/2 for p in x], speedup_vs_metis, width=width, label='Balanced/METIS+Naive')
ax.set_xticks(x)
ax.set_xticklabels(mat_labels)
ax.set_ylabel('Speedup (time)')
ax.set_title('Speedup of Balanced over Baselines (np=2)')
ax.axhline(1.0, color='k', linestyle='--')
ax.legend()
fig.tight_layout()
figpath = os.path.join(FIG_DIR, 'fig3_speedup.png')
fig.savefig(figpath)
print('Saved', figpath)
plt.close(fig)

# Figure 4: Ablation study on lower_bound_frac (cant, np=2)
if abl_rows:
    fracs = []
    gflops = []
    comms  = []
    for r in abl_rows:
        if r.get('np') != '2': continue
        try:
            fr = float(r['frac'])
        except:
            continue
        fracs.append(fr)
        gflops.append(float(r['GFlops']) if r.get('GFlops') and r.get('GFlops')!='NA' else 0.0)
        comms.append(float(r['CommTime']) if r.get('CommTime') and r.get('CommTime')!='NA' else 0.0)
    if fracs:
        # sort by frac
        idxs = sorted(range(len(fracs)), key=lambda i: fracs[i])
        fracs_s = [fracs[i] for i in idxs]
        gfl_s = [gflops[i] for i in idxs]
        comm_s = [comms[i] for i in idxs]
        fig, ax1 = plt.subplots(figsize=(8,4))
        ax2 = ax1.twinx()
        ax1.plot(fracs_s, gfl_s, '-o', color='blue', label='GFlops')
        ax2.plot(fracs_s, comm_s, '-s', color='red', label='CommTime')
        ax1.set_xlabel('lower_bound_frac')
        ax1.set_ylabel('GFlops', color='blue')
        ax2.set_ylabel('CommTime (s)', color='red')
        ax1.set_title('Ablation Study on lower_bound_frac (cant, np=2)')
        fig.tight_layout()
        figpath = os.path.join(FIG_DIR, 'fig4_ablation.png')
        fig.savefig(figpath)
        print('Saved', figpath)
        plt.close(fig)
else:
    print('No ablation data; skipping fig4.')

# Figure 5: Load Imbalance Ratio (LIR)
if bal_rows:
    # organize by matrix -> alg -> lir
    data = defaultdict(dict)
    for r in bal_rows:
        if r.get('np') != '2': continue
        m = basename(r['Matrix'])
        alg = r['Algorithm']
        lir = r.get('LIR','NA')
        try:
            lir_f = float(lir) if lir not in (None,'NA','') else 0.0
        except:
            lir_f = 0.0
        data[m][alg] = lir_f
    mats = [m for m in mat_order]
    fig, ax = plt.subplots(figsize=(8,4))
    x = range(len(mats))
    width = 0.2
    for i, alg in enumerate(alg_order):
        vals = [ data.get(m,{}).get(alg,0.0) for m in mats]
        pos = [p + (i-1)*width for p in x]
        ax.bar(pos, vals, width=width, label=alg)
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace('.mtx','') for m in mats])
    ax.set_ylabel('Load Imbalance Ratio (LIR)')
    ax.set_title('Load Imbalance Ratio Comparison (np=2)')
    ax.legend()
    fig.tight_layout()
    figpath = os.path.join(FIG_DIR, 'fig5_load_balance.png')
    fig.savefig(figpath)
    print('Saved', figpath)
    plt.close(fig)
else:
    print('No balance data; skipping fig5.')

# Figure 6: Per-process remote_nnz comparison (cant, np=2)
# Try to get balanced remote from ablation (frac=1.0) and metis total nnz from balance CSV
remote_bal = None
for r in abl_rows:
    try:
        if r.get('np')=='2' and float(r.get('frac',-1))==1.0:
            remote_bal = (int(r.get('RemoteNNZ_p0',0)), int(r.get('RemoteNNZ_p1',0)))
            break
    except:
        continue
metis_tot = None
for r in bal_rows:
    if basename(r['Matrix'])=='cant.mtx' and r['Algorithm']=='metis_naive' and r['np']=='2':
        try:
            metis_tot = (int(r['NNZ_p0']), int(r['NNZ_p1']))
        except:
            metis_tot = None
        break
if remote_bal is not None and metis_tot is not None:
    labels = ['p0','p1']
    x = range(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(6,4))
    ax.bar([p-width/2 for p in x], metis_tot, width=width, label='METIS+Naive (NNZ total proxy)')
    ax.bar([p+width/2 for p in x], remote_bal, width=width, label='Balanced (remote_nnz)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('NNZ / remote nnz')
    ax.set_title('Per-process Remote NNZ (cant, np=2)')
    ax.legend()
    fig.tight_layout()
    figpath = os.path.join(FIG_DIR, 'fig6_remote_nnz.png')
    fig.savefig(figpath)
    print('Saved', figpath)
    plt.close(fig)
else:
    print('Insufficient data for fig6; skipping.')

print('Plotting finished.')
