# -*- coding: utf-8 -*-
"""Supplementary Fig. S3 — RFind-sc between all pairs of columns (knock-outs and safe-targeting groups) of one pseudo-bulk.

Same columns, genes and Delta as Fig. 4 (analysis/perturbseq/03_axis_vs_ko.py): knock-outs and safe-targeting groups with >= 10
nuclei; genes with CPM >= 1 in the non-targeting column; Delta = (x - e)/sqrt(e + 1), e = the share of each gene in the summed
counts of all these columns and the non-targeting column x column UMI. Each column is a signed list: its top 500 genes (UP) and
bottom 500 (DOWN) by Delta. RFind score of a pair = the
4-term Running Fisher of the two signed lists (UPxUP - UPxDOWN - DOWNxUP + DOWNxDOWN; each term -log10 of the smallest
hypergeometric p over every pair of cut-offs, universe = the genes compared; analysis/running_fisher.py, clip 3000), the same score
as Fig. 4. The score is symmetric, so only i < j is computed. ~64 min on one core; run in parallel.

  python3 figs/Supp_S3_Perturbseq/analysis/figS3_rfind_allpairs.py [n_workers] [delta|ratio] [g151|DG|CA1|CA3|GABA|HIP]
  (project root)

The second argument chooses the per-column value the genes are ranked by: 'delta' (default, as above) or 'ratio', the simple
fold change log2((CPM + 0.1)/(baseline CPM + 0.1)) against the same baseline, so that
the two differ only in the transform.

The third argument chooses the pseudo-bulk (default g151). DG, CA1 and CA3 are the single subclasses of Fig. 4
(pb_g017_DG.npz, pb_CA1.npz, pb_g017_CA3.npz) and GABA the GABAergic neurons of Fig. 4 (pb_GABA3.npz); the paper uses g151, HIP and
GABA (Supplementary Fig. S3b-e). HIP = dentate gyrus,
CA1 and CA3 combined (load_hip; Delta only: the second argument must be delta for HIP). g151 = the thalamic glutamatergic neurons '151 TH Prkcd Grin2c Glut'
(pb_g151.npz), the one cell type that Shi et al. name in their main-text comparison of knock-outs (Fig. 2I SWI/SNF,
highlighted in Fig. 2H); the cell type of Supplementary Fig. S3b-d.

Output: figs/Supp_S3_Perturbseq/output/data/rfind_allpairs_<context>[_ratio].npz (S = n x n float32, targets, kind, nuclei) and
figs/Supp_S3_Perturbseq/output/figS3_rfind_numbers_<context>[_ratio].txt
"""
import os, sys, time
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(BASE, 'analysis'))
from running_fisher import running_fisher

OUT = os.path.join(BASE, 'figs', 'Supp_S3_Perturbseq', 'output')
os.makedirs(os.path.join(OUT, 'data'), exist_ok=True)
CUT, MINC, N_SIG, CLIP = 1.0, 10, 500, 3000.0          # as Fig. 4 (analysis/perturbseq/03_axis_vs_ko.py)
PB = {'DG': 'pb_g017_DG.npz', 'CA1': 'pb_CA1.npz', 'CA3': 'pb_g017_CA3.npz', 'g151': 'pb_g151.npz', 'GABA': 'pb_GABA3.npz'}   # GABA = the GABAergic neurons of Fig. 4
PSEUDO = 0.1                                           # CPM pseudo-count of the simple fold change

def load_hip():
    """HIP: dentate gyrus, CA1 and CA3 combined without mixing their nuclei.
    Each cell type keeps its own genes, baseline and Delta (as load()); the Delta of a column is summed over the three cell types
    (genes kept in all three; knock-outs and safe-targeting groups kept, with >= MINC nuclei, in all three). n = nuclei summed."""
    parts = []
    for c in ('DG', 'CA1', 'CA3'):
        d = np.load(os.path.join(BASE, 'data', 'Perturbseq', PB[c]), allow_pickle=True)
        M, t, n, g = d['counts'].astype(np.float64), d['targets'].astype(str), d['n_cells'], d['genes'].astype(str)
        jb = int(np.where(t == 'Non_target')[0][0])
        keep = (M / M.sum(0, keepdims=True) * 1e6)[:, jb] >= CUT
        inb = (n >= MINC) | (np.arange(len(t)) == jb)
        base = M[keep][:, inb].sum(1) / M[keep][:, inb].sum()
        sel = (n >= MINC) & (t != 'Non_target')
        X = M[keep][:, sel]; ex = np.outer(base, X.sum(0))
        parts.append(((X - ex) / np.sqrt(ex + 1), g[keep], t[sel], n[sel]))
    genes = sorted(set.intersection(*[set(p[1]) for p in parts]))
    cols = sorted(set.intersection(*[set(p[2]) for p in parts]))
    V = np.zeros((len(genes), len(cols))); nn = np.zeros(len(cols), int)
    for Vx, gx, tx, nx in parts:
        gi = {x: i for i, x in enumerate(gx)}; ci = {x: i for i, x in enumerate(tx)}
        V += Vx[np.ix_([gi[x] for x in genes], [ci[x] for x in cols])]
        nn += np.array([nx[ci[x]] for x in cols])
    t = np.array(cols)
    kind = np.array(['safe' if x.startswith('Safe_target') else 'KO' for x in t])
    O = np.argsort(-V, axis=0, kind='stable')
    return O[:N_SIG].T.copy(), O[::-1][:N_SIG].T.copy(), V.shape[0], t, kind, nn

def load(transform='delta', ctx='g151'):
    d = np.load(os.path.join(BASE, 'data', 'Perturbseq', PB[ctx]), allow_pickle=True)
    M, t, n = d['counts'].astype(np.float64), d['targets'].astype(str), d['n_cells']
    jb = int(np.where(t == 'Non_target')[0][0])
    keep = (M / M.sum(0, keepdims=True) * 1e6)[:, jb] >= CUT
    ## baseline = the share of each gene in the summed counts of ALL columns kept for Fig. 4 (knock-outs and safe-targeting
    ## with >= MINC nuclei, and the non-targeting column), exactly as analysis/perturbseq/03_axis_vs_ko.py. The non-targeting
    ## column enters the baseline but is not one of the columns compared.
    inb = (n >= MINC) | (np.arange(len(t)) == jb)
    base = M[keep][:, inb].sum(1) / M[keep][:, inb].sum()
    sel = (n >= MINC) & (t != 'Non_target')
    M, t, n = M[keep][:, sel], t[sel], n[sel]
    kind = np.array(['safe' if x.startswith('Safe_target') else 'KO' for x in t])
    if transform == 'ratio':
        cpm = M / M.sum(0, keepdims=True) * 1e6
        V = np.log2((cpm + PSEUDO) / (base[:, None] * 1e6 + PSEUDO))
    else:
        ex = np.outer(base, M.sum(0))
        V = (M - ex) / np.sqrt(ex + 1)
    O = np.argsort(-V, axis=0, kind='stable')
    return O[:N_SIG].T.copy(), O[::-1][:N_SIG].T.copy(), V.shape[0], t, kind, n

G = {}
def init(up, dn, P):
    G['UP'] = [u.tolist() for u in up]; G['DN'] = [x.tolist() for x in dn]; G['P'] = P

def row(i):
    UP, DN, P = G['UP'], G['DN'], G['P']
    out = np.empty(len(UP) - i - 1, np.float32)
    for k, j in enumerate(range(i + 1, len(UP))):
        s = (running_fisher(UP[i], UP[j], universe_size=P, clip=CLIP) - running_fisher(UP[i], DN[j], universe_size=P, clip=CLIP)
             - running_fisher(DN[i], UP[j], universe_size=P, clip=CLIP) + running_fisher(DN[i], DN[j], universe_size=P, clip=CLIP))
        out[k] = min(max(s, -CLIP), CLIP)
    return i, out

if __name__ == '__main__':
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    tf = sys.argv[2] if len(sys.argv) > 2 else 'delta'
    assert tf in ('delta', 'ratio')
    ctx = sys.argv[3] if len(sys.argv) > 3 else 'g151'
    assert ctx in PB or ctx == 'HIP'
    sfx = '' if tf == 'delta' else '_ratio'
    tag = f'_{ctx}'
    t0 = time.time()
    UP, DN, P, tgt, kind, nuc = load_hip() if ctx == 'HIP' else load(tf, ctx)
    n = len(tgt)
    print(f'{n:,} columns ({(kind == "KO").sum():,} KO, {(kind == "safe").sum()} safe), {P:,} genes, {n * (n - 1) // 2:,} pairs, '
          f'{nw} workers', flush=True)
    S = np.zeros((n, n), np.float32)
    mpc = mp.get_context('fork')                        # (not 'ctx': that name is the pseudo-bulk context)
    order = list(range(n - 1))                          # long rows first would balance better; imap_unordered is enough
    done = 0
    with mpc.Pool(nw, initializer=init, initargs=(UP, DN, P)) as pool:
        for i, out in pool.imap_unordered(row, order, chunksize=4):
            S[i, i + 1:] = out; S[i + 1:, i] = out
            done += 1
            if done % 100 == 0:
                print(f'  {done:,}/{n - 1:,} rows, {time.time() - t0:6.0f} s', flush=True)
    np.fill_diagonal(S, np.nan)
    np.savez_compressed(os.path.join(OUT, 'data', f'rfind_allpairs_{ctx}{sfx}.npz'), S=S, targets=tgt, kind=kind, nuclei=nuc)
    ko, sf = kind == 'KO', kind == 'safe'
    tri = lambda m: S[np.ix_(m, m)][np.triu_indices(int(m.sum()), 1)]
    kk, ss, ks = tri(ko), tri(sf), S[np.ix_(ko, sf)].ravel()
    lines = [f'RFind-sc between all pairs, {ctx} ({n:,} columns, {P:,} genes; UP/DOWN {N_SIG} by {"Delta" if tf == "delta" else "simple ratio"}; 4-term Running Fisher)',
             f'KOxKO   mean {kk.mean():+.3f}  median {np.median(kk):+.3f}  sd {kk.std():.3f}  positive {(kk > 0).mean():.3f}  n {len(kk):,}',
             f'safexsafe mean {ss.mean():+.3f}  median {np.median(ss):+.3f}  sd {ss.std():.3f}  positive {(ss > 0).mean():.3f}  n {len(ss):,}',
             f'KOxsafe mean {ks.mean():+.3f}  median {np.median(ks):+.3f}  sd {ks.std():.3f}  positive {(ks > 0).mean():.3f}  n {len(ks):,}']
    iu = np.triu_indices(int(ko.sum()), 1); Sk = S[np.ix_(ko, ko)][iu]; tk = tgt[ko]
    top = np.argsort(-Sk)[:15]
    lines.append('top KO pairs: ' + ', '.join(f'{tk[iu[0][k]]}x{tk[iu[1][k]]} {Sk[k]:+.1f}' for k in top))
    lines.append(f'time {time.time() - t0:.0f} s')
    open(os.path.join(OUT, f'figS3_rfind_numbers{tag}{sfx}.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
