"""Supplementary Table 2 (redo) — is the placement of the Fig. 4d gene group special among knock-outs of the same library?

The Perturb-seq library (Shi et al.) is not a random sample of genes: it favours genes expressed widely in the brain and
strongly linked to disease (transcription factors and receptors first). Public gene sets compared within this library
inherit that bias, so the position of the Fig. 4d group among public sets does not show how particular its placement is.
Here the null distribution is built from the library itself.

For each cell type (the Fig. 4 ranking of that cell type):
  - statistic : AUC = P(RFind score of a member > that of a non-member) for the five knock-outs of the Fig. 4d group
                (Nr3c1, Nr3c2, Hsp90ab1, Ncoa1, Ncoa3) against all other knock-outs; < 0.5 = towards the opposed end
                (P(member > non-member), as for a Wilcoxon rank-sum test)
  - null 1    : 100,000 sets of five knock-outs drawn at random from all knock-outs of the cell type
  - null 2    : 100,000 matched sets: for each of the five genes, one knock-out drawn from the knock-outs (other than the
                five) in the same quartile of target expression (CPM in the non-targeting nuclei) and the same quartile of
                the number of nuclei, so that the random sets have knock-outs of comparable strength
  - P         : (1 + number of random sets with AUC <= the observed AUC) / (1 + 100,000), one-sided towards the opposed end
Input : figs/Fig4_Perturbseq/output/SuppTable_KO_rankings.xlsx, one sheet per cell type (supp_table_ko_rankings.py)
Output: figs/Fig4_Perturbseq/output/geneset_random_null.csv
Run   : python figs/Fig4_Perturbseq/analysis/fig4_random_set_null.py   (from the project root; numpy, pandas)
"""
import numpy as np
import pandas as pd

OUT = 'figs/Fig4_Perturbseq/output'
CTS = ['DG', 'CA1', 'CA3', 'GABA']
SHEET = {'DG': 'DG Glut', 'CA1': 'CA1 Glut', 'CA3': 'CA3 Glut', 'GABA': 'GABA'}   # sheets of SuppTable_KO_rankings.xlsx
FIG4D = ['Nr3c1', 'Nr3c2', 'Hsp90ab1', 'Ncoa1', 'Ncoa3']
B = 100_000
rng = np.random.default_rng(1)


def auc_of(ranks, idx, n):
    """AUC of the members idx (rows of 'ranks', average ranks of the scores, 1 = lowest) against the other knock-outs."""
    k = idx.shape[-1]
    r = ranks[idx].sum(axis=-1)
    return (r - k * (k + 1) / 2) / (k * (n - k))


def no_repeats(draw, redraw):
    """Redraw the rows that contain the same knock-out twice (a set has five different knock-outs)."""
    while True:
        s = np.sort(draw, axis=1)
        bad = np.flatnonzero((s[:, 1:] == s[:, :-1]).any(axis=1))
        if not len(bad):
            return draw
        draw[bad] = redraw(len(bad))


rows = []
for ct in CTS:
    d = pd.read_excel(f'{OUT}/SuppTable_KO_rankings.xlsx', sheet_name=SHEET[ct])
    n = len(d)
    ranks = d['RFind score'].rank(method='average').to_numpy()
    member = d['Knock-out'].isin(FIG4D).to_numpy()
    m_idx = np.flatnonzero(member)
    assert len(m_idx) == len(FIG4D), (ct, d.loc[member, 'Knock-out'].tolist())
    obs = float(auc_of(ranks, m_idx, n))
    # null 1: any five different knock-outs
    k = len(FIG4D)
    draw = no_repeats(rng.integers(0, n, (B, k)), lambda m: rng.integers(0, n, (m, k)))
    null1 = auc_of(ranks, draw, n)
    # null 2: matched by quartiles of target expression and of nuclei (five different knock-outs, none of the five genes)
    q_cpm = pd.qcut(d['Target CPM (non-targeting nuclei)'].rank(method='first'), 4, labels=False).to_numpy()
    q_nuc = pd.qcut(d['Nuclei'].rank(method='first'), 4, labels=False).to_numpy()
    cells = q_cpm * 4 + q_nuc
    pools = [np.flatnonzero((cells == cells[i]) & ~member) for i in m_idx]
    pool_sizes = [len(p) for p in pools]
    matched = lambda m: np.stack([p[rng.integers(0, len(p), m)] for p in pools], axis=1)
    draw2 = no_repeats(matched(B), matched)
    null2 = auc_of(ranks, draw2, n)
    rows.append(dict(cell_type=ct, n_knockouts=n, auc_fig4d=round(obs, 4),
                     null_any_mean=round(float(null1.mean()), 4), null_any_sd=round(float(null1.std()), 4),
                     p_any=(1 + int((null1 <= obs).sum())) / (1 + B),
                     pct_any=round(100 * float((null1 <= obs).mean()), 3),
                     null_matched_mean=round(float(null2.mean()), 4), null_matched_sd=round(float(null2.std()), 4),
                     p_matched=(1 + int((null2 <= obs).sum())) / (1 + B),
                     pct_matched=round(100 * float((null2 <= obs).mean()), 3),
                     matched_pool_sizes=' '.join(map(str, pool_sizes)), draws=B, seed=1))
res = pd.DataFrame(rows)
res.to_csv(f'{OUT}/geneset_random_null.csv', index=False)
print(res.to_string(index=False))
