"""RFind-sc demo (Python): the DAM-index of 1,000 microglia (Hammond et al. 2019; white matter after LPC vs saline injection,
500 each) with the DAM reference of Keren-Shaul et al. 2017 (the top 500 genes per direction, 408 up and 427 down in the data), as in Figs. 2 and 5 of the paper.
Uses the Python copy of the Running Fisher kernel in this repository (analysis/running_fisher.py with analysis/_hyper.py; the kernel
of rfindsc-py 0.3.5) and numpy, scipy, pandas. The same steps as demo/run_demo.R:
  1. Delta = (x - e) / sqrt(e + 1), x = UMI count, e = baseline share of the gene x total UMI of the cell (all genes);
  2. UP = the 500 genes with the largest positive Delta, DOWN = the 500 most negative;
  3. S = RF(UP, up) - RF(UP, down) - RF(DOWN, up) + RF(DOWN, down) (Running Fisher; universe P = the genes of the demo).
Then the AUC of LPC vs saline (Mann-Whitney) and the check against the paper's stored DAM-index (demo/expected_output.csv).
Run from the project root:  python demo/run_demo.py
Output: demo/output_py.csv (cell_id, condition, DAM_index)
"""
import os
import sys
import time

t_start = time.perf_counter()
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.io import mmread  # noqa: E402
from scipy.stats import rankdata  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'analysis'))
from running_fisher import running_fisher  # noqa: E402

D = os.path.join(ROOT, 'demo', 'data')
NTOP, CLIP = 500, 300.0

# ---- data ----
try:
    X = mmread(os.path.join(D, 'counts.mtx.gz'), spmatrix=False)          # genes x cells, raw UMI counts
except TypeError:                                                         # scipy < 1.15
    X = mmread(os.path.join(D, 'counts.mtx.gz'))
X = X.toarray()                                                           # dense: 5,369 x 1,000 (43 MB)
genes = pd.read_csv(os.path.join(D, 'genes.tsv.gz'), sep='\t', dtype={'gene': str, 'baseline_share': float})
cells = pd.read_csv(os.path.join(D, 'labels.tsv'), sep='\t', dtype={'cell_id': str, 'condition': str, 'total_umi': float})
barcodes = pd.read_csv(os.path.join(D, 'barcodes.tsv.gz'), header=None)[0].tolist()
assert barcodes == cells.cell_id.tolist() and X.shape == (len(genes), len(cells))
ref = pd.read_csv(os.path.join(D, 'dam_reference.tsv'), sep='\t', dtype={'direction': str, 'rank': int, 'gene': str})
ref_up = ref.loc[ref.direction == 'up', 'gene'].tolist()                  # ranked by fold change, lower case
ref_dn = ref.loc[ref.direction == 'down', 'gene'].tolist()
glow = np.array([g.lower() for g in genes.gene])
P = len(set(glow))                                                        # universe of the pair cell x reference (pair_universe)
share = genes.baseline_share.to_numpy()
print(f'demo: {X.shape[0]} genes (P = {P}) x {X.shape[1]} cells; reference {len(ref_up)} up / {len(ref_dn)} down genes')


def rf(a, b):
    return running_fisher(a, b, P, clip=CLIP)


# ---- per-cell Delta, UP / DOWN lists, signed Running Fisher score ----
score = np.empty(X.shape[1])
for j in range(X.shape[1]):
    x = X[:, j]
    e = share * cells.total_umi.iat[j]                                    # expected count
    d = (x - e) / np.sqrt(e + 1)                                          # Delta
    pos = np.flatnonzero(d > 0); pos = pos[np.argsort(-d[pos], kind='stable')]
    neg = np.flatnonzero(d < 0); neg = neg[np.argsort(d[neg], kind='stable')]
    up, dn = list(glow[pos[:NTOP]]), list(glow[neg[:NTOP]])
    s = rf(up, ref_up) - rf(up, ref_dn) - rf(dn, ref_up) + rf(dn, ref_dn)
    score[j] = max(-CLIP, min(s, CLIP))

lpc = (cells.condition == 'LPC').to_numpy()
r = rankdata(score); n1, n0 = lpc.sum(), (~lpc).sum()
auc = (r[lpc].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
out = pd.DataFrame({'cell_id': cells.cell_id, 'condition': cells.condition, 'DAM_index': [repr(float(v)) for v in score]})
out.to_csv(os.path.join(ROOT, 'demo', 'output_py.csv'), index=False)
t_run = time.perf_counter() - t_start

print('\nfirst cells:')
print(pd.DataFrame({'cell_id': cells.cell_id, 'condition': cells.condition, 'DAM_index': score.round(4)}).head(5).to_string(index=False))
print(f'\nmedian DAM-index: LPC {np.median(score[lpc]):.3f}, saline {np.median(score[~lpc]):.3f}')
print(f'AUC (LPC vs saline): {auc:.4f}')
fx = os.path.join(ROOT, 'demo', 'expected_output.csv')
if os.path.exists(fx):
    ex = pd.read_csv(fx, dtype={'cell_id': str, 'condition': str, 'DAM_index': float})
    assert ex.cell_id.tolist() == cells.cell_id.tolist()
    dmax = float(np.max(np.abs(ex.DAM_index.to_numpy() - score)))
    print(f"max |difference| from the paper's stored DAM-index: {dmax:.3g} -> {'OK' if dmax < 1e-8 else 'DIFFERS'}")
print(f'written demo/output_py.csv; run time {t_run:.1f} s (Python {sys.version.split()[0]}, numpy {np.__version__})')
