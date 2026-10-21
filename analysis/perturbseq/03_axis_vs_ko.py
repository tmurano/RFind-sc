# -*- coding: utf-8 -*-
"""Match a De Donno stress axis against every Perturb-seq KO effect vector of one context (RFind).

  <python> analysis/perturbseq/03_axis_vs_ko.py

Settings:
  KO side   targets with >= 10 cells (+ non); genes kept if the non column has CPM >= 1;
            P2 residual (x - e)/sqrt(e + 1) with e = the share of each gene in the summed counts of ALL columns
            (knock-outs, safe-targeting, non-targeting) x column UMI, in every context. A baseline of the safe-targeting
            groups alone would carry its own estimation error into every knock-out, so that knock-outs would share it and
            look alike.
  axis side P2_vsNoStress of 02_stress_axis.py, ALL genes kept (the reference is not hand-edited)
  RFind     4-term Running Fisher on the top/bottom 500 genes: each pair of lists scored by the smallest p over every
            top-G x top-K cut-off pair; clip 3000; kernel = running_fisher of analysis/running_fisher.py, a copy of the
            kernel of rfindsc-py 0.3.5, unchanged in 0.3.6 (same computation as analysis/running_fisher.R)
Runs (axis x context), RUNS below (Fig. 4c = CA1; Fig. 4d, e and Supplementary Tables 1 and 2 = all four):
  DG stress  x g017_DG  (DG granule cells only)
  CA3 stress x g017_CA3 (CA3 only)
  CA1 stress x CA1 (016 CA1-ProS Glut); GABA stress x GABA3 (Pvalb + Sst + Lamp5 nuclei, summed)

Output (analysis/perturbseq/results/; PERTURB_RES overrides the folder, e.g. for checks): axis_<axis>_vs_<context>.csv — one
row per KO and safe guide: RFind, rank among the KOs, z (against KOs), n_cells, target CPM in the non column (is the target
expressed in this context); summary printed.
"""
import os, sys, csv, time, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(BASE, 'analysis'))   # analysis/running_fisher.py: the Running Fisher kernel of rfindsc-py 0.3.5, copied
from running_fisher import running_fisher
DATA = os.path.join(BASE, 'data', 'Perturbseq')
RES = os.environ.get('PERTURB_RES', os.path.join(BASE, 'analysis', 'perturbseq', 'results')); os.makedirs(RES, exist_ok=True)
CUT, MINC, N_SIG, CLIP = 1.0, 10, 500, 3000.0
MINC = int(os.environ.get('PERTURB_MINC', MINC))    # override of the nuclei cut-off, for checks only (e.g. 5, 20)
AXCOL = 'P2_vsNoStress'
PB = {'g017_DG': 'pb_g017_DG.npz', 'g017_CA3': 'pb_g017_CA3.npz', 'g151': 'pb_g151.npz',
      'CA1': 'pb_CA1.npz', 'GABA3': 'pb_GABA3.npz'}        # CA1 / GABA3 from 01b_scan_subclasses.py + 01c_combine_gaba.py
RUNS = [('DG', 'g017_DG'), ('CA3', 'g017_CA3'), ('CA1', 'CA1'), ('GABA', 'GABA3')]
## PERTURB_RUNS overrides the list for one-off work without touching the pipeline default, e.g.
##   PERTURB_RUNS='DG:CA1,DG:g017_CA3' python analysis/perturbseq/03_axis_vs_ko.py
if os.environ.get('PERTURB_RUNS'):
    RUNS = [tuple(x.split(':')) for x in os.environ['PERTURB_RUNS'].split(',') if x.strip()]
WATCH = ['Nr3c1', 'Nr3c2', 'Ncoa1', 'Crebrf', 'Fkbp5', 'Sgk1', 'Creb1', 'Crhr1', 'Ntrk2', 'Bdnf', 'Grin2b', 'Camk2a']
t0 = time.time(); say = lambda *a: print(f'[{time.time()-t0:6.1f}s]', *a, flush=True)

def load_axis(tag):
    sym, val = [], []
    with open(os.path.join(DATA, f'axis_dedonno_{tag}_stress.csv')) as g:
        for r in csv.DictReader(g):
            sym.append(r['Symbol']); val.append(float(r[AXCOL]))           # ambient genes kept
    return np.array(sym), np.array(val)

def load_ko(ctx):
    d = np.load(os.path.join(DATA, PB[ctx]), allow_pickle=True)
    M0, gn, tg0, nc0 = d['counts'].astype(np.float64), d['genes'].astype(str), d['targets'].astype(str), d['n_cells']
    jb0 = int(np.where(tg0 == 'Non_target')[0][0])
    sel = (nc0 >= MINC) | (np.arange(len(tg0)) == jb0)
    M, tgt, nc = M0[:, sel], tg0[sel], nc0[sel]
    jb = int(np.where(tgt == 'Non_target')[0][0])
    cpm_non_all = M0[:, jb0] / M0[:, jb0].sum() * 1e6                    # target expression, before the gene cut
    keep = (M / np.maximum(M.sum(0, keepdims=True), 1) * 1e6)[:, jb] >= CUT
    M, gnk = M[keep], gn[keep]
    kind = np.array(['non' if t == 'Non_target' else 'safe' if t.startswith('Safe_target') else 'KO' for t in tgt])
    SAFE = kind == 'safe'
    base = M.sum(1) / M.sum()                                            # baseline = the share of each gene in the summed
                                                                         # counts of ALL columns
    ex = np.outer(base, M.sum(0)); V = (M - ex) / np.sqrt(ex + 1)
    n_ko_all = int(sum(1 for t in tg0 if t != 'Non_target' and not t.startswith('Safe_target')))
    return dict(V=V, genes=gnk, tgt=tgt, nc=nc, kind=kind, cpm_target=dict(zip(gn, cpm_non_all)), n_ko_all=n_ko_all)

def rf_pair(a, b, P):
    """Running Fisher score of two ranked gene lists (gene indices, most extreme first): -log10 of the smallest
    hypergeometric upper-tail p, P(X >= M), over every pair of cut-offs (top G of a, top K of b), universe P, clipped at
    CLIP — the same definition as analysis/running_fisher.R. Computed by the kernel running_fisher
    of analysis/running_fisher.py (log scale; only the corners at ranks of shared genes, identical to the full G x K scan)."""
    return running_fisher(a.tolist(), b.tolist(), universe_size=P, clip=CLIP)

def rfind(Vc, ax):
    """4-term RFind score of every column of Vc against the axis: top/bottom N_SIG genes of each side."""
    P = Vc.shape[0]
    oA = np.argsort(-ax); aUP, aDN = oA[:N_SIG], oA[::-1][:N_SIG]
    S = np.empty(Vc.shape[1])
    for j in range(Vc.shape[1]):
        o = np.argsort(-Vc[:, j]); UP, DN = o[:N_SIG], o[::-1][:N_SIG]
        S[j] = rf_pair(aUP, UP, P) - rf_pair(aUP, DN, P) - rf_pair(aDN, UP, P) + rf_pair(aDN, DN, P)
    return np.clip(S, -CLIP, CLIP)

cache = {}
for axtag, ctx in RUNS:
    if not os.path.exists(os.path.join(DATA, PB[ctx])):
        say(f'skip {axtag} x {ctx}: {PB[ctx]} not built yet'); continue
    sym, val = load_axis(axtag)
    K = cache.setdefault(ctx, load_ko(ctx))
    common, iK, iA = np.intersect1d(K['genes'], sym, return_indices=True)
    use = K['kind'] != 'non'
    Vc, tgt, knd, nc = K['V'][iK][:, use], K['tgt'][use], K['kind'][use], K['nc'][use]
    ax = val[iA]
    S = rfind(Vc, ax)
    KO, SF = knd == 'KO', knd == 'safe'
    sK = S[KO]
    mu, sd = sK.mean(), sK.std()
    order = np.argsort(-sK); rank = np.empty(KO.sum(), int); rank[order] = np.arange(1, KO.sum() + 1)
    rank_all = np.full(len(S), -1); rank_all[np.flatnonzero(KO)] = rank
    out = os.path.join(RES, f'axis_{axtag}_vs_{ctx}.csv')
    with open(out, 'w', newline='') as g:
        w = csv.writer(g)
        w.writerow(['target', 'kind', 'n_cells', 'rfind', 'rank_among_KO', 'z_vs_KO', 'target_cpm_in_non'])
        for j in np.argsort(-S):
            t = tgt[j]; gsym = t if knd[j] == 'KO' else ''
            w.writerow([t, knd[j], int(nc[j]), f'{S[j]:.3f}', rank_all[j] if knd[j] == 'KO' else '',
                        f'{(S[j]-mu)/sd:+.3f}',
                        f'{K["cpm_target"].get(gsym, np.nan):.2f}' if gsym else ''])
    tk = tgt[KO]
    print()
    say(f'=== axis {axtag} stress x {ctx}: {len(common):,} common genes; KO {KO.sum():,} of {K["n_ko_all"]:,} '
        f'with >= {MINC} cells (median {int(np.median(nc[KO]))} cells/KO); safe {SF.sum()} ===')
    say(f'  RFind  KO mean {mu:+.2f} sd {sd:.2f}')
    say('  aligned top 20 : ' + ', '.join(f'{tk[j]} {sK[j]:+.1f}' for j in order[:20]))
    say('  opposite top 20: ' + ', '.join(f'{tk[j]} {sK[j]:+.1f}' for j in order[::-1][:20]))
    for gname in WATCH:
        jj = np.flatnonzero(tk == gname)
        if len(jj):
            j = int(jj[0])
            say(f'    {gname:<7} {sK[j]:+7.2f}  rank {rank[j]:>5}/{KO.sum():,} (from the opposite end {KO.sum()-rank[j]+1:>5})'
                f'  z {(sK[j]-mu)/sd:+.2f}  target CPM in non {K["cpm_target"].get(gname, np.nan):.1f}')
    say(f'  -> {os.path.relpath(out, BASE)}')
