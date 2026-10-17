# -*- coding: utf-8 -*-
"""Fig. 4 with two other ways of matching: decoupleR's ULM and the scFOCAL-type score (Supplementary Table 2).

Question: does the Fig. 4 result depend on matching by RFind-sc, or do the matching methods of Fig. 5 give it too? Everything
is computed on the inputs of Fig. 4: the same columns, genes and Delta of the knock-outs (analysis/perturbseq/03_axis_vs_ko.py)
and the same stress response (P2_vsNoStress of 02_stress_axis.py). As in the paper, each knock-out is a reference (its top 500
and bottom 500 genes by Delta) and the stress response is the profile that is scored.
  RFind    the Fig. 4 score, read from analysis/perturbseq/results/axis_<cell type>_vs_<context>.csv (not recomputed)
  ULM      decoupleR run_ulm: the Delta of the knock-out's 1,000 genes as weights (mor), 0 for its other genes; the genes of
           the regression = the genes that are a target of at least one knock-out (decoupleR keeps the genes of the network);
           score = t of the slope (= t of the Pearson correlation). Same as decoupleR::run_ulm (decoupleR 2.16.0, R 4.5.2;
           minsize 5) on the DG run: largest difference of t 5.3e-14
  ULM_all  the same with every gene a target of every knock-out (the whole Delta profile as weights)
  scFOCAL  Spearman correlation between the knock-out's Delta and the stress response over the knock-out's 1,000 genes (the
           scFOCAL-type score of Fig. 5: Spearman between the reference values and the profile over the reference genes)
Measures, as in the paper:
  - positions of the five knock-outs of Fig. 4d (Nr3c1, Nr3c2, Hsp90ab1, Ncoa1, Ncoa3) and of Crebrf (Fig. 4c), counted from
    the opposed end
  - Supplementary Table 2: AUC of the five against the other knock-outs and P against 100,000 random sets of five (any, and
    matched by quartiles of target expression and nuclei), with the draws of fig4_random_set_null.py (seed 1; the draws do not
    depend on the scores, so every method gets the same draws and RFind must reproduce geneset_random_null.csv)
  - Spearman correlation of the knock-out scores of each method with those of RFind
In the paper: Supplementary Table 2, sheet "Other matching methods" (made by supp_table_random_sets.py from
fig4_other_methods_summary.csv): the positions of the five knock-outs of Fig. 4d and the comparison with random sets for RFind,
ULM and scFOCAL. ULM_all, Crebrf and the Spearman correlations with RFind are not in the paper.
Output: figs/Fig4_Perturbseq/output/fig4_other_methods.txt (summary), fig4_other_methods_summary.csv (cell type x method)
Run   : python figs/Fig4_Perturbseq/analysis/fig4_other_methods.py   (from the project root; numpy, scipy, pandas;
        ~2 min)
"""
import os, time
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.getcwd()                                     # the project root (inputs)
assert os.path.isdir(os.path.join(BASE, 'figs', 'Fig4_Perturbseq')), 'run from the project root'
DATA = os.path.join(BASE, 'data', 'Perturbseq')
RES = os.path.join(BASE, 'analysis', 'perturbseq', 'results')
FIG4 = os.path.join(BASE, 'figs', 'Fig4_Perturbseq', 'output')                 # Fig. 4 tables (inputs)
OUT = os.path.join(os.path.dirname(HERE), 'output'); os.makedirs(OUT, exist_ok=True)   # outputs of this script

CUT, MINC, N_SIG = 1.0, 10, 500                        # as 03_axis_vs_ko.py and figS3_rfind_allpairs.py
PB = {'g017_DG': 'pb_g017_DG.npz', 'g017_CA3': 'pb_g017_CA3.npz', 'CA1': 'pb_CA1.npz', 'GABA3': 'pb_GABA3.npz',
      'g151': 'pb_g151.npz'}
RUNS = [('DG', 'g017_DG'), ('CA1', 'CA1'), ('CA3', 'g017_CA3'), ('GABA', 'GABA3')]   # the order of fig4_random_set_null.py
SHEET = {'DG': 'DG Glut', 'CA1': 'CA1 Glut', 'CA3': 'CA3 Glut', 'GABA': 'GABA'}   # sheets of SuppTable_KO_rankings.xlsx
FIG4D = ['Nr3c1', 'Nr3c2', 'Hsp90ab1', 'Ncoa1', 'Ncoa3']
WATCH = FIG4D + ['Crebrf']
METHODS = ['RFind', 'ULM', 'ULM_all', 'scFOCAL']
B = 100_000
t0 = time.time()
LOG = []
def say(*a):
    s = ' '.join(str(x) for x in a); LOG.append(s); print(f'[{time.time() - t0:6.1f}s] {s}', flush=True)


## ---------- inputs (as 03_axis_vs_ko.py) ----------
def load_axis(tag):
    d = pd.read_csv(os.path.join(DATA, f'axis_dedonno_{tag}_stress.csv'), usecols=['Symbol', 'P2_vsNoStress'])
    return d['Symbol'].astype(str).to_numpy(), d['P2_vsNoStress'].to_numpy(float)

def load_delta(ctx):
    """Delta of every column with >= MINC nuclei (knock-outs and safe-targeting), genes with CPM >= CUT in the non-targeting
    column, baseline = the share of each gene in the summed counts of these columns and the non-targeting column
    (03_axis_vs_ko.py, figS3_rfind_allpairs.py)."""
    d = np.load(os.path.join(DATA, PB[ctx]), allow_pickle=True)
    M, gn, t, n = d['counts'].astype(np.float64), d['genes'].astype(str), d['targets'].astype(str), d['n_cells']
    jb = int(np.where(t == 'Non_target')[0][0])
    cpm_non = dict(zip(gn, M[:, jb] / M[:, jb].sum() * 1e6))
    keep = (M / np.maximum(M.sum(0, keepdims=True), 1) * 1e6)[:, jb] >= CUT
    inb = (n >= MINC) | (np.arange(len(t)) == jb)
    base = M[keep][:, inb].sum(1) / M[keep][:, inb].sum()
    sel = (n >= MINC) & (t != 'Non_target')
    X = M[keep][:, sel]; ex = np.outer(base, X.sum(0))
    kind = np.array(['safe' if x.startswith('Safe_target') else 'KO' for x in t[sel]])
    return (X - ex) / np.sqrt(ex + 1), gn[keep], t[sel], n[sel], kind, cpm_non

## ---------- the matching methods: references = columns of R (genes x references); profile(s) = columns of Q ----------
def ref_genes(R):
    """Top N_SIG and bottom N_SIG genes of every reference column (the argsort of 03_axis_vs_ko.py)."""
    idx = np.empty((R.shape[1], 2 * N_SIG), int)
    for j in range(R.shape[1]):
        o = np.argsort(-R[:, j]); idx[j] = np.r_[o[:N_SIG], o[::-1][:N_SIG]]
    return idx

def ulm_t(Q, R, idx=None):
    """t of decoupleR's ULM for every profile (column of Q) x reference (column of R): Pearson correlation over the genes of
    the network between the profile and the weights (Delta of the reference's genes idx, 0 elsewhere; idx None = all genes),
    t = r sqrt(df / (1 - r^2)), df = genes - 2. Returns profiles x references."""
    if idx is None:
        W, F = R, np.ones(R.shape[0], bool)
    else:
        W = np.zeros_like(R)
        for j in range(R.shape[1]):
            W[idx[j], j] = R[idx[j], j]
        F = np.zeros(R.shape[0], bool); F[np.unique(idx)] = True
    Qf, Wf = Q[F] - Q[F].mean(0), W[F] - W[F].mean(0)
    r = (Qf.T @ Wf) / np.outer(np.sqrt((Qf ** 2).sum(0)), np.sqrt((Wf ** 2).sum(0)))
    df = int(F.sum()) - 2
    return r * np.sqrt(df / ((1 - r) * (1 + r))), int(F.sum())

def scfocal(Q, R, idx):
    """Spearman correlation between each reference (Delta over its idx genes) and each profile over the same genes;
    ordinal ranks (Delta has no ties in practice). Returns profiles x references."""
    out = np.empty((Q.shape[1], R.shape[1]))
    m = (2 * N_SIG - 1) / 2
    for j in range(R.shape[1]):
        rq = np.argsort(np.argsort(Q[idx[j]], axis=0), axis=0) - m
        rr = np.argsort(np.argsort(R[idx[j], j])) - m
        out[:, j] = (rr @ rq) / (rr @ rr)
    return out


## ---------- Part 1: Fig. 4 (stress response x knock-outs) ----------
def auc_of(ranks, idx, n):
    k = idx.shape[-1]; r = ranks[idx].sum(axis=-1)
    return (r - k * (k + 1) / 2) / (k * (n - k))

def no_repeats(draw, redraw):
    while True:
        s = np.sort(draw, axis=1)
        bad = np.flatnonzero((s[:, 1:] == s[:, :-1]).any(axis=1))
        if not len(bad):
            return draw
        draw[bad] = redraw(len(bad))

def random_sets(tables, scores):
    """fig4_random_set_null.py for one method: tables = the sheets of SuppTable_KO_rankings.xlsx (row order = the draws), scores = {ct: {KO:
    score}}. One generator (seed 1) over the cell types in the order of RUNS, as there."""
    rng = np.random.default_rng(1); rows = {}
    for ct, _ in RUNS:
        d = tables[ct]; n = len(d)
        ranks = pd.Series([scores[ct][k] for k in d['Knock-out']]).rank(method='average').to_numpy()
        member = d['Knock-out'].isin(FIG4D).to_numpy(); m_idx = np.flatnonzero(member); k = len(FIG4D)
        obs = float(auc_of(ranks, m_idx, n))
        draw = no_repeats(rng.integers(0, n, (B, k)), lambda m: rng.integers(0, n, (m, k)))
        null1 = auc_of(ranks, draw, n)
        q_cpm = pd.qcut(d['Target CPM (non-targeting nuclei)'].rank(method='first'), 4, labels=False).to_numpy()
        q_nuc = pd.qcut(d['Nuclei'].rank(method='first'), 4, labels=False).to_numpy()
        cells = q_cpm * 4 + q_nuc
        pools = [np.flatnonzero((cells == cells[i]) & ~member) for i in m_idx]
        matched = lambda m: np.stack([p[rng.integers(0, len(p), m)] for p in pools], axis=1)
        null2 = auc_of(ranks, no_repeats(matched(B), matched), n)
        rows[ct] = dict(auc_fig4d=obs, p_any=(1 + int((null1 <= obs).sum())) / (1 + B),
                        p_matched=(1 + int((null2 <= obs).sum())) / (1 + B))
    return rows

def part1():
    tables = {ct: pd.read_excel(os.path.join(FIG4, 'SuppTable_KO_rankings.xlsx'), sheet_name=SHEET[ct])
              for ct, _ in RUNS}
    scores = {m: {} for m in METHODS}; nfeat = {}; rows = []
    for ct, ctx in RUNS:
        sym, val = load_axis(ct)
        V, gn, tgt, nc, kind, _ = load_delta(ctx)
        common, iK, iA = np.intersect1d(gn, sym, return_indices=True)
        ko = kind == 'KO'
        R, names, ax = V[iK][:, ko], tgt[ko], val[iA]
        rr = pd.read_csv(os.path.join(RES, f'axis_{ct}_vs_{ctx}.csv')).query("kind == 'KO'").set_index('target')
        rf = rr['rfind']
        assert set(rf.index) == set(names), ct
        idx = ref_genes(R)
        t_ulm, nfeat[ct] = ulm_t(ax[:, None], R, idx)
        t_all, _ = ulm_t(ax[:, None], R)
        sc = scfocal(ax[:, None], R, idx)
        S = {'RFind': rf.reindex(names).to_numpy(), 'ULM': t_ulm[0], 'ULM_all': t_all[0], 'scFOCAL': sc[0]}
        n = len(names)
        say(f'{ct} x {ctx}: {len(common):,} genes, {n:,} knock-outs; ULM regression over {nfeat[ct]:,} genes')
        P = {}
        for m in METHODS:
            scores[m][ct] = dict(zip(names, S[m]))
            opp = np.empty(n, int); opp[np.argsort(S[m], kind='stable')] = np.arange(1, n + 1)   # 1 = most opposed
            if m == 'RFind':                                 # the ranks of 03_axis_vs_ko.py (its file rounds the scores)
                opp = n - rr['rank_among_KO'].reindex(names).to_numpy().astype(int) + 1
            P[m] = opp
            pos = {g: int(opp[np.flatnonzero(names == g)[0]]) for g in WATCH}
            rho = stats.spearmanr(S[m], S['RFind']).statistic
            top20 = len(set(names[np.argsort(S[m])[:20]]) & set(names[np.argsort(S['RFind'])[:20]]))
            rows.append(dict(cell_type=ct, method=m, n_ko=n, **{f'pos_{g}': pos[g] for g in WATCH},
                             median_pos_fig4d=float(np.median([pos[g] for g in FIG4D])),
                             spearman_vs_RFind=round(float(rho), 4), opposed_top20_shared_with_RFind=top20))
            say(f'  {m:8s} from the opposed end: ' + ', '.join(f'{g} {pos[g]}' for g in WATCH)
                + f'  | Spearman with RFind {rho:+.3f}; opposed top 20 shared {top20}'
                + '  | opposed top 10: ' + ', '.join(names[np.argsort(S[m])[:10]]))
        rk = pd.DataFrame({'cell_type': ct, 'knock_out': names, **{m: S[m] for m in METHODS}})
        for m in METHODS:
            rk[f'pos_{m}'] = P[m]                            # position from the opposed end (1 = most opposed)
        part1.ranks.append(rk)
    ## Supplementary Table 2 for every method
    say('Supplementary Table 2 (five knock-outs of Fig. 4d against random sets of five; AUC < 0.5 = towards the opposed end):')
    ref = pd.read_csv(os.path.join(FIG4, 'geneset_random_null.csv')).set_index('cell_type')
    summ = pd.DataFrame(rows).set_index(['cell_type', 'method'])
    for m in METHODS:
        rs = random_sets(tables, scores[m])
        for ct, _ in RUNS:
            for k, v in rs[ct].items():
                summ.loc[(ct, m), k] = v
            if m == 'RFind':
                assert abs(rs[ct]['auc_fig4d'] - ref.loc[ct, 'auc_fig4d']) < 1e-4 and \
                    np.isclose(rs[ct]['p_any'], ref.loc[ct, 'p_any'], rtol=1e-9, atol=0) and \
                    np.isclose(rs[ct]['p_matched'], ref.loc[ct, 'p_matched'], rtol=1e-9, atol=0), ct   # CSV read: last digits
        say(f'  {m:8s} ' + '; '.join(f'{ct} AUC {rs[ct]["auc_fig4d"]:.3f} P {rs[ct]["p_any"]:.2g} / matched {rs[ct]["p_matched"]:.2g}'
                                     for ct, _ in RUNS))
    say('  RFind reproduces geneset_random_null.csv (AUC and both P, every cell type)')
    ## RFind positions must be those of fig4d_GR-complex_positions.csv
    gp = pd.read_csv(os.path.join(FIG4, 'fig4d_GR-complex_positions.csv'))
    for _, r in gp.iterrows():
        assert summ.loc[(r.cell_type, 'RFind'), f'pos_{r.KO}'] == r.rank_from_opposed_end, (r.cell_type, r.KO)
    say('  RFind reproduces fig4d_GR-complex_positions.csv (the 20 positions)')
    return summ.reset_index()
part1.ranks = []


if __name__ == '__main__':
    s1 = part1()
    s1.to_csv(os.path.join(OUT, 'fig4_other_methods_summary.csv'), index=False, float_format='%.6g')
    say(f'done in {time.time() - t0:.0f} s')
    open(os.path.join(OUT, 'fig4_other_methods.txt'), 'w').write('\n'.join(LOG) + '\n')
