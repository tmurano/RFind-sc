# -*- coding: utf-8 -*-
"""Pool the three GABAergic subclass pseudo-bulks (Pvalb, Sst, Lamp5) into one GABA context.

  <python> analysis/perturbseq/01c_combine_gaba.py

Why pooled: the De Donno GABA type is not split into interneuron classes, so the Perturb-seq
side is pooled the same way. Caveat: after pooling, KO columns can differ in their subclass make-up, which biases
the mixed group 017 (shown by a homogeneity check not included in this repository). The share of each subclass
per knock-out is therefore printed.

Input : data/Perturbseq/pb_{Pvalb,Sst,Lamp5}.npz   (01b_scan_subclasses.py)
Output: data/Perturbseq/pb_GABA3.npz (same layout)
"""
import os, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data', 'Perturbseq')
PARTS = ['Pvalb', 'Sst', 'Lamp5']

D = {p: np.load(os.path.join(DATA, f'pb_{p}.npz'), allow_pickle=True) for p in PARTS}
genes = D[PARTS[0]]['genes']
assert all(np.array_equal(D[p]['genes'], genes) for p in PARTS), 'gene axes differ'
targets = sorted(set().union(*[set(D[p]['targets'].astype(str)) for p in PARTS]))
tix = {t: j for j, t in enumerate(targets)}
M = np.zeros((len(genes), len(targets))); n = np.zeros((len(PARTS), len(targets)), np.int64)
u = np.zeros(len(targets)); det = np.zeros(len(genes), np.int64)
for k, p in enumerate(PARTS):
    d = D[p]; j = np.array([tix[t] for t in d['targets'].astype(str)])
    M[:, j] += d['counts']; n[k, j] = d['n_cells']; u[j] += d['tot_umi']; det += d['det_cells']
ntot = n.sum(0)
order = np.lexsort((np.array(targets), -ntot))                        # most cells first, as in the other pseudo-bulks
np.savez_compressed(os.path.join(DATA, 'pb_GABA3.npz'), counts=M[:, order].astype(np.float32), genes=genes,
                    targets=np.array(targets)[order], n_cells=ntot[order], tot_umi=u[order], det_cells=det)
ko = np.array([not (t == 'Non_target' or t.startswith('Safe_target')) for t in np.array(targets)[order]])
fr = n[:, order][:, ko] / np.maximum(ntot[order][ko], 1)
print(f'pb_GABA3.npz: {len(targets):,} targets, {int(ntot.sum()):,} cells, KO median {int(np.median(ntot[order][ko]))} cells/KO')
for k, p in enumerate(PARTS):
    print(f'  {p:<6} fraction per KO: mean {fr[k].mean():.2f}, range {fr[k].min():.2f}-{fr[k].max():.2f}')
