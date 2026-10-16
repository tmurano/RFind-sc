# -*- coding: utf-8 -*-
"""De Donno 2026 — Stress vs No Stress pseudo-bulk axis for one cell type (control line), to be matched
against the Perturb-seq KO effect vectors.

  <python> analysis/perturbseq/02_stress_axis.py

Cells: our QC (analysis/dedonno/01_qc_embed.R -> data/DeDonno2026/output/cell_meta_qc.csv), the authors'
cell types with number suffixes removed (21 types), line == 'Ctrl' (flox/flox without Cre = receptors intact; the
only line kept by 01_qc_embed.R, so this filter keeps every cell).
For each cell type, raw counts of all No Stress cells and of all Stress cells are summed (two pseudo-bulk
libraries) and turned into an axis in three ways:
  log2FC        = log2(CPM_Stress + 1) - log2(CPM_NoStress + 1)
  P2delta       = P2 residual of Stress minus that of No Stress (baseline = the two libraries pooled)
  P2_vsNoStress = (Stress count - expected) / sqrt(expected + 1), expected = No Stress profile x Stress UMI
                  -> the axis used for the KO matching
Ambient genes (haemoglobin + choroid plexus; AMBIENT below) are only FLAGGED in the `ambient` column; they are
kept in the axis (the reference is not hand-edited).

Input : data/DeDonno2026/raw/dedonno2026.h5ad (raw/X), data/DeDonno2026/output/cell_meta_qc.csv
Output: data/Perturbseq/axis_dedonno_<DG|CA3>_stress.csv
        (Symbol, log2FC, P2delta, P2_vsNoStress, CPM_NoStress, CPM_Stress, ambient)
"""
import os, csv, time, numpy as np, h5py, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))
H5 = os.path.join(BASE, 'data', 'DeDonno2026', 'raw', 'dedonno2026.h5ad')
META = os.path.join(BASE, 'data', 'DeDonno2026', 'output', 'cell_meta_qc.csv')
OUT = os.path.join(BASE, 'data', 'Perturbseq'); os.makedirs(OUT, exist_ok=True)
AXES = {'DG': 'DG Glut', 'CA3': 'CA3 Glut', 'CA1': 'CA1 Glut', 'GABA': 'GABA'}   # tag -> author cell type (numbers removed)
GR_TARGETS = ['Fkbp5', 'Tsc22d3', 'Sgk1', 'Ddit4', 'Klf9', 'Zbtb16', 'Plin4', 'Sult1a1', 'Per1']
## ambient list: haemoglobin + choroid plexus genes
AMBIENT = {'Hba-a1', 'Hba-a2', 'Hba-x', 'Hbb-bs', 'Hbb-bt', 'Hbb-b1', 'Hbb-b2', 'Hbb-y', 'Alas2', 'Hbq1a', 'Hbq1b',
           'Ttr', 'Folr1', 'Prlr', 'Kl', 'Clic6', 'Aqp1', 'Enpp2', 'Kcnj13', 'Otx2os1', '1500015O10Rik', 'Htr2c', 'Sostdc1'}
t0 = time.time(); say = lambda *a: print(f'[{time.time()-t0:6.1f}s]', *a, flush=True)
dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])

f = h5py.File(H5, 'r')
def col(path):
    o = f[path]
    if isinstance(o, h5py.Group):                        # categorical column
        return dec(o['categories'][:])[o['codes'][:]]
    return dec(o[:])
ids = col('obs/_index') if '_index' in f['obs'] else col('obs/index')
genes = col('raw/var/feature_name')                      # columns of raw/X
X = f['raw/X']; indptr = X['indptr'][:]
assert len(indptr) - 1 == len(ids)
say(f'h5ad: {len(ids):,} cells x {len(genes):,} genes')

meta = pd.read_csv(META, usecols=['cell_id', 'cell_type', 'condition', 'line'])
meta = meta[meta.line == 'Ctrl'].set_index('cell_id')
pos = pd.Series(np.arange(len(ids)), index=ids)

for tag, ct in AXES.items():
    m = meta[meta.cell_type == ct]
    rows = {c: np.sort(pos.loc[m.index[m.condition == c]].values) for c in ['No Stress', 'Stress']}
    say(f'{tag} ({ct}, control line): No Stress {len(rows["No Stress"]):,} / Stress {len(rows["Stress"]):,} cells')
    sums = {}
    for c, rr in rows.items():
        s = np.zeros(len(genes), np.float64)
        for r in rr:                                         # CSR rows of the selected cells
            lo, hi = indptr[r], indptr[r + 1]
            s[X['indices'][lo:hi]] += X['data'][lo:hi]
        sums[c] = s
    A, B = sums['No Stress'], sums['Stress']
    cpmA, cpmB = A / A.sum() * 1e6, B / B.sum() * 1e6
    log2fc = np.log2(cpmB + 1) - np.log2(cpmA + 1)
    M2 = np.c_[A, B]; base = M2.sum(1) / M2.sum(); ex = np.outer(base, M2.sum(0))
    V2 = (M2 - ex) / np.sqrt(ex + 1); p2 = V2[:, 1] - V2[:, 0]
    exS = cpmA * B.sum() / 1e6                                 # counts expected in Stress from the No Stress profile
    p2ns = (B - exS) / np.sqrt(exS + 1)
    amb = np.isin(genes, list(AMBIENT))
    say(f'  total UMI No Stress {A.sum():,.0f} / Stress {B.sum():,.0f}')
    say('  GR targets (log2FC): ' + ', '.join(
        f'{g} {log2fc[np.flatnonzero(genes == g)[0]]:+.3f}' for g in GR_TARGETS if (genes == g).any()))
    say(f'  ambient genes flagged: {int(amb.sum())} (e.g. ' + ', '.join(
        f'{g} {log2fc[np.flatnonzero(genes == g)[0]]:+.2f}' for g in ['Hba-a1', 'Hbb-bs', 'Ttr'] if (genes == g).any()) + ')')
    o = np.argsort(-p2ns)
    say('  top 12 up (P2_vsNoStress): ' + ', '.join(genes[o[:12]]))
    say('  top 12 down (P2_vsNoStress): ' + ', '.join(genes[o[::-1][:12]]))
    out = os.path.join(OUT, f'axis_dedonno_{tag}_stress.csv')
    with open(out, 'w', newline='') as g:
        w = csv.writer(g)
        w.writerow(['Symbol', 'log2FC', 'P2delta', 'P2_vsNoStress', 'CPM_NoStress', 'CPM_Stress', 'ambient'])
        for j in o:
            w.writerow([genes[j], f'{log2fc[j]:.5f}', f'{p2[j]:.4f}', f'{p2ns[j]:.4f}',
                        f'{cpmA[j]:.3f}', f'{cpmB[j]:.3f}', bool(amb[j])])
    say(f'  wrote {os.path.relpath(out, BASE)}')
f.close()
