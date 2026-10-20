"""Run time of scoring cells against one reference with the Python package rfindsc (Supplementary Table 3, sheet
"Run time and memory"; the Python counterpart of fig5_runtime_score.R).

Data: the counts written by fig5_runtime_export.R into $RFIND_RUNTIME_EXPORT (the 75,572 cells of the data of Hammond et al. after
quality control), the first n cells of the same random order as the R script; reference: DAM, the top L genes per direction of the
published list (L = 500 or all), the genes present in the counts. Scoring: rfindsc.score_panel_pairs() with its own gene filter
and baseline from the cells given. The time covers scoring only, not reading the counts.
Output: one line appended to figs/Fig5_benchmark/output/results/fig5_runtime_raw.tsv
Run (from the project root; rfindsc importable, e.g. PYTHONPATH=tools/rfindsc-py/src):
  RFIND_RUNTIME_EXPORT=<folder> python3 figs/Fig5_benchmark/analysis/fig5_runtime_score.py n=1000 L=500 cores=1
"""
import os
import math
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
import rfindsc

kv = dict(a.split('=', 1) for a in sys.argv[1:])
N, LARG, CORES = kv.get('n', '1000'), kv.get('L', '500'), int(kv.get('cores', '1'))
D = os.environ['RFIND_RUNTIME_EXPORT']
OUT = 'figs/Fig5_benchmark/output/results/fig5_runtime_raw.tsv'

ng, nc = (int(v) for v in open(os.path.join(D, 'dims.txt')).read().split())
X = sp.csc_matrix((np.fromfile(os.path.join(D, 'x.bin'), dtype='<f8'), np.fromfile(os.path.join(D, 'i.bin'), dtype='<i4'),
                   np.fromfile(os.path.join(D, 'p.bin'), dtype='<i4')), shape=(ng, nc))
genes = open(os.path.join(D, 'genes.txt')).read().split('\n')[:ng]
cells = np.array(open(os.path.join(D, 'cells.txt')).read().split('\n')[:nc])
order = np.array([int(v) for v in open(os.path.join(D, 'order.txt')).read().split()]) - 1     # R's 1-based order
n = nc if N == 'all' else int(N)
keep = np.sort(order[:n])
X = X[:, keep]
libsize = np.asarray(X.sum(axis=0)).ravel()
tab = pd.read_csv('data/Hammond2019/output/references/dam_KerenShaul.csv')
LV = None if LARG == 'all' else int(LARG)
gset = set(genes)
upl = tab[tab['Fold Change'] > 0].sort_values('Fold Change', ascending=False, kind='stable')['Symbol'].tolist()
dnl = tab[tab['Fold Change'] < 0].sort_values('Fold Change', ascending=True, kind='stable')['Symbol'].tolist()
up = [g for g in (upl if LV is None else upl[:LV]) if g in gset]
dn = [g for g in (dnl if LV is None else dnl[:LV]) if g in gset]
print(f'rfindsc_py: {X.shape[1]} cells x {X.shape[0]} genes; DAM L = {LARG}: up {len(up)}, down {len(dn)}; {CORES} core(s)', flush=True)

t0 = time.perf_counter()
out = rfindsc.score_panel_pairs(X, {'DAM': {'up': up, 'down': dn}}, config=rfindsc.rfindsc_config(n_cores=CORES, top_n=math.inf if LV is None else LV), genes=genes,
                                cells=cells[keep].tolist(), libsize=libsize)
secs = time.perf_counter() - t0
s = out['scores']
s = np.asarray(s['DAM'] if isinstance(s, pd.DataFrame) else s, dtype=float).ravel()
assert len(s) == X.shape[1] and np.isfinite(s).any()
print(f'scoring time {secs:.1f} s', flush=True)
new = not os.path.exists(OUT)
with open(OUT, 'a') as g:
    if new:
        g.write('method\tcells\tL\tcores\tseconds\tup\tdown\n')
    g.write(f'rfindsc_py\t{X.shape[1]}\t{LARG}\t{CORES}\t{secs:.2f}\t{len(up)}\t{len(dn)}\n')
