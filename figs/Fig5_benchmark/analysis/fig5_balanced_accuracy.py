"""Supplementary Table 3, sheet Balanced accuracy (Results 5, Methods) — how well the sign of each method's score separates the two
groups of cells of each task (the result is not in a figure).
Balanced accuracy at a threshold of 0: a cell is called positive when its score is above 0 and negative otherwise (0 or below);
  balanced accuracy = the mean of the fraction of positive cells with a score above 0 and the fraction of negative cells with a score
  of 0 or below (0.5 = chance; also 0.5 for a score that is above 0 in every cell, or never above 0).
Tasks: the eight type A tasks at L = 500 (the DAM and Glucocorticoid tasks H1, D1-D3; the PS2APP and Forced swim tasks DA1, AS1-AS3) and the three type B tasks
  (H2, H3, H5); the Maturity task (H4, ceiling) is left out, as in the comparisons of the AUC.
Methods: those of the AUC table at L = 500 (the main-figure methods; UCell with the up genes only, WMEAN and ORA of decoupleR only in
  Supplementary Table 3). Scores: the per-cell scores of the benchmark, the same cells and scores as the AUC; Beyondcell cells without
  a score = 0, as in the benchmark (fig5_bench_redo.R, BEYONDCELL_NAN = 'as_zero'); no other method has a missing score.
95% CIs: bootstrap over mice (fig5_mice.py: all cells of a drawn mouse; within the positive and the negative group, or from all mice
  in H2 and H5, whose mice have cells in both groups), 1,000 resamples, the same resamples for every method; default_rng(1), one
  generator, the tasks asking in the order H1, DA1, D1, D2, D3, AS1, AS2, AS3, H2, H3, H5, and the tasks with the same mice sharing the
  draws (H1, DA1; D1-D3, AS1-AS3; H2, H5). RFind-sc minus the method with the 95% CI of the paired
  difference; RFind-sc higher / lower when the CI excludes 0 and |difference| > 0.01, otherwise not distinguished (the rule of the AUC
  comparisons, fig5_redo_tables.py).
Input : figs/Fig5_benchmark/output/results/fig5_redo_percell_L500.csv.gz (H1, D1-D3, H2, H3, H5),
        extra_refs/fig5_extra_percell_L500.csv.gz (DA1, AS1-AS3)
Output: figs/Fig5_benchmark/output/results/fig5_balanced_accuracy_L500.csv (read by fig5_redo_tables.py); one row per task and method:
        n_pos, n_neg, pos_above0 / neg_above0 (fraction of positive / negative cells with a score above 0), balanced_accuracy, lo, hi
        (its 95% CI), diff_vs_RFind, d_lo, d_hi (RFind-sc minus the method and its 95% CI), result
Run   : <python with numpy and pandas> figs/Fig5_benchmark/analysis/fig5_balanced_accuracy.py   (from the project root; about 10 s)
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig5_mice import MouseDraws, mouse_of  # noqa: E402

O = 'figs/Fig5_benchmark/output/'
TYPE_A = ['H1', 'DA1', 'D1', 'D2', 'D3', 'AS1', 'AS2', 'AS3']   # drawing order (DA1 next to H1, AS1-AS3 after D1-D3)
TYPE_B = ['H2', 'H3', 'H5']
METH_A = ['RFind-sc', 'AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM (log2FC)', 'scFOCAL-type',
          'UCell (up genes)', 'WMEAN (log2FC)', 'ORA']               # RFind-sc first: column 0 is the reference of the differences
METH_B = ['RFind-sc', 'AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM', 'WMEAN', 'ORA']
NB, SEED = 1000, 1


def verdict(diff, lo, hi):
    return 'RFind-sc higher' if (lo > 0 and diff > 0.01) else ('RFind-sc lower' if (hi < 0 and diff < -0.01) else 'not distinguished')


pc = pd.concat([pd.read_csv(O + 'results/fig5_redo_percell_L500.csv.gz', low_memory=False),
                pd.read_csv(O + 'extra_refs/fig5_extra_percell_L500.csv.gz', low_memory=False)], ignore_index=True, sort=False)
assert pc.positive.dtype == bool and set(TYPE_A + TYPE_B) <= set(pc.id), sorted(set(pc.id))

rows = []
DRAWS = MouseDraws(np.random.default_rng(SEED), NB)
for t in TYPE_A + TYPE_B:
    g = pc[pc.id == t]
    pos = g.positive.to_numpy()
    ms = METH_A if t in TYPE_A else METH_B
    S = g[ms].to_numpy(float)
    S[:, ms.index('Beyondcell')] = np.nan_to_num(S[:, ms.index('Beyondcell')], nan=0.0)   # cells without a score = 0
    assert not np.isnan(S).any(), (t, [m for m, k in zip(ms, np.isnan(S).any(0)) if k])
    C = np.where(pos[:, None], S > 0, S <= 0).astype(float)          # 1 = the sign calls the cell's group correctly
    ba = 0.5 * (C[pos].mean(0) + C[~pos].mean(0))
    boot = np.empty((NB, len(ms)))
    for b, i in enumerate(DRAWS.cells(mouse_of(g.cell_id), pos)):
        Ci, pi = C[i], pos[i]
        boot[b] = 0.5 * (Ci[pi].mean(0) + Ci[~pi].mean(0))
    for k, m in enumerate(ms):
        lo, hi = np.percentile(boot[:, k], [2.5, 97.5])
        d_lo, d_hi = np.percentile(boot[:, 0] - boot[:, k], [2.5, 97.5])
        diff = ba[0] - ba[k]
        rows.append({'id': t, 'type': 'A' if t in TYPE_A else 'B', 'method': m, 'n_pos': int(pos.sum()), 'n_neg': int((~pos).sum()),
                     'pos_above0': (S[pos, k] > 0).mean(), 'neg_above0': (S[~pos, k] > 0).mean(),
                     'balanced_accuracy': ba[k], 'lo': lo, 'hi': hi,
                     'diff_vs_RFind': np.nan if k == 0 else diff, 'd_lo': np.nan if k == 0 else d_lo, 'd_hi': np.nan if k == 0 else d_hi,
                     'result': '' if k == 0 else verdict(diff, d_lo, d_hi)})
res = pd.DataFrame(rows)
assert np.allclose(res.balanced_accuracy, 0.5 * (res.pos_above0 + 1 - res.neg_above0), rtol=0, atol=1e-12)
os.makedirs(O + 'results', exist_ok=True)
out = O + 'results/fig5_balanced_accuracy_L500.csv'
res.to_csv(out, index=False)

## summary on the console: the mean over the eight type A tasks, and the results against RFind-sc by type
main = METH_A[:10] + ['ULM']                                         # type B: ULM with weight 1; WMEAN and ORA are table-only
x = res[res.method.isin(main)]
print('mean balanced accuracy, eight type A tasks (%):')
print((100 * x[x.type == 'A'].groupby('method', sort=False).balanced_accuracy.mean()).round(1).to_string())
print('\nresults against RFind-sc (main-figure methods):')
print(x[x.method != 'RFind-sc'].groupby(['type', 'result']).size().unstack(fill_value=0).to_string())
print('\nULM against RFind-sc, type A:', x[(x.type == 'A') & (x.method == 'ULM (log2FC)')].result.value_counts().to_dict())
print('wrote', out, f'({len(res)} rows)')
