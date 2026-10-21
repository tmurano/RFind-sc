"""Supplementary Table 2 — the capability benchmark of Fig. 5 and Supplementary Fig. S4 (the eight type A tasks and the three type B
tasks of these figures, the Maturity task at the ceiling, the ten methods of the figure, the methods reported only here).
Tasks (12): type A (differential expression lists with direction and size) = the DAM and Glucocorticoid tasks, H1 DAM and D1-D3
  Glucocorticoid (oligodendrocytes, astrocytes, endothelial cells), and the PS2APP and Forced swim tasks (fig5_extra_score.R), DA1, the
  cells of H1 with the PS2APP reference, and AS1-AS3, the cells of D1-D3 with the Forced swim reference; Fig. 5c draws H1 and D2,
  Supplementary Fig. S4b the other six; Figs. 5e and 5f average the eight type A tasks;
  type B (MSigDB Hallmark gene sets) = H2, H3 Interferon and H5 G2M (Fig. 5d); H4 Maturity (type A) is left out of the main figure
  because the median AUC of the methods is 0.99 or more at L = 500 (ceiling) and is shown here only, without win / tie / loss.
Methods in the main figure: RFind-sc, AddModuleScore, AUCell, UCell (its signed signature: up genes and down genes), singscore, ssGSEA,
  VISION, Beyondcell, ULM of decoupleR (log2 fold changes as weights for type A; weight 1 for gene sets) and scFOCAL-type (type A only).
Reported here only: UCell with the up genes only and decoupleR's WMEAN (log2FC weights; weight 1 for gene sets) and ORA.
Comparison (Methods): AUC over all task cells; 95% CI of the AUC and of the difference from RFind-sc by bootstrap over mice (all cells of
  each drawn mouse; 1,000 resamples, the same resamples for every method); RFind-sc higher / lower when the CI of the difference excludes 0 and |difference| > 0.01,
  otherwise not distinguished (the difference is 0.01 or less or cannot be told apart; this does not show that the two are at the same level).
Beyondcell: cells without a score are kept and set to 0, as in the benchmark.
Sheets Best length (block 1: fig5_best_length_8tasks.py, every method at its best reference length in the eight type A tasks, per
  task / leave one task out over the other seven, with the length re-selected in each resample; block 2: fig5_best_length.py, leave one
  task out within the four DAM and Glucocorticoid tasks) and Mean AUC differences (block 1: the mean over the eight type A tasks;
  block 2: the mean over the four DAM and Glucocorticoid tasks; RFind-sc minus each method); "comparable" = the whole 95% CI of the difference lies within
  -0.01..+0.01.
Sheet Run time and memory (fig5_runtime.sh; written only when its output exists): the run time and peak memory of scoring the cells of
  the data of Hammond et al. against the DAM reference and of the knock-out x knock-out matching of Supplementary Fig. S3, with the
  computer and the versions of the run.
The task ids (H1, D1-D3, DA1, AS1-AS3, H2, H3, H5, H4) are the keys of the input files and are not written to the sheets.
Input : figs/Fig5_benchmark/output/results/fig5_redo_auc_L<L>.csv, extra_refs/fig5_extra_auc_L<L>.csv, extra_refs/fig5_extra_percell_L500.csv.gz,
        results/fig5_redo_beyondcell_L500.csv,
        results/fig5_best_length.csv, results/fig5_mean_auc_diff.csv, results/fig5_best_length_8tasks.csv, results/fig5_mean_auc_diff_8tasks.csv;
        if present, results/fig5_runtime.csv and results/fig5_runtime_env.tsv
Output: figs/Fig5_benchmark/output/SuppTable2.xlsx
Run   : <python with pandas and openpyxl> figs/Fig5_benchmark/analysis/fig5_redo_tables.py   (from the project root)
"""
import os
import numpy as np
import pandas as pd

O = 'figs/Fig5_benchmark/output/'
LS = ['100', '250', '500', '1000', '2000', 'all']
TYPE_A_MAIN, TYPE_A_ADDED = ['H1', 'D1', 'D2', 'D3'], ['DA1', 'AS1', 'AS2', 'AS3']   # DAM and Glucocorticoid; PS2APP and Forced swim
TYPE_A = TYPE_A_MAIN + TYPE_A_ADDED                                # the eight type A tasks (Figs. 5e, f)
TYPE_A_C = ['H1', 'D2']                                            # the type A tasks drawn in Fig. 5c; the other six in Supplementary Fig. S4b
TYPE_B, CEIL = ['H2', 'H3', 'H5'], ['H4']
TASKS = TYPE_A + TYPE_B + CEIL
PAIR = {'DA1': 'H1', 'AS1': 'D1', 'AS2': 'D2', 'AS3': 'D3'}         # added task -> the DAM or Glucocorticoid task on the same cells
TSETNAME = {t: 'Fig. 5c' if t in TYPE_A_C else 'Supplementary Fig. S4b' for t in TYPE_A}   # where each type A task is drawn
MAIN_A = ['AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM (log2FC)', 'scFOCAL-type']
MAIN_B = ['AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM']
EXTRA_A = ['UCell (up genes)', 'WMEAN (log2FC)', 'ORA']
EXTRA_B = ['WMEAN', 'ORA']
DISPLAY = {'UCell (up genes)': 'UCell (up genes only)', 'ULM (log2FC)': 'ULM (decoupleR; log2FC weights)', 'WMEAN (log2FC)': 'WMEAN (decoupleR; log2FC weights)',
           'ORA': 'ORA (decoupleR)'}
FORM = {'RFind-sc': 'A: ranked up and down lists', 'scFOCAL-type': 'A: ranked list (log2 fold changes)',
        'ULM (log2FC)': 'A: log2 fold changes as weights', 'WMEAN (log2FC)': 'A: log2 fold changes as weights',
        'singscore': "A': up and down sets", 'ssGSEA': "A': up and down sets",
        'VISION': "A': up +1, down -1", 'Beyondcell': "A': up and down sets", 'UCell': "A': signed signature (up and down)", 'AddModuleScore': "A'': up set", 'UCell (up genes)': "A'': up set",
        'AUCell': "A'': up set", 'ORA': "A'': up set"}
TASKNAME = {'H1': 'DAM: LPC-injected vs saline-injected white matter (microglia)', 'H2': 'Interferon: cluster 9 vs other microglia',
            'H3': 'Interferon: P540 vs P100 microglia', 'H4': 'Maturity: P100 vs E14.5 and P4/P5 microglia',   # E14.5 = age 'E14' of the metadata
            'H5': 'G2M: cluster 2b vs other microglia', 'D1': 'Glucocorticoid: stress vs no stress, oligodendrocytes',
            'D2': 'Glucocorticoid: stress vs no stress, astrocytes', 'D3': 'Glucocorticoid: stress vs no stress, endothelial cells',
            'DA1': 'PS2APP: LPC-injected vs saline-injected white matter (microglia)',
            'AS1': 'Forced swim: stress vs no stress, oligodendrocytes', 'AS2': 'Forced swim: stress vs no stress, astrocytes',
            'AS3': 'Forced swim: stress vs no stress, endothelial cells'}


def methods_of(t):
    if t in TYPE_B:
        return ['RFind-sc'] + MAIN_B, EXTRA_B
    return ['RFind-sc'] + MAIN_A, EXTRA_A


def load(L):                                                       # the benchmark, then the added type A tasks
    d = pd.concat([pd.read_csv(f'{O}results/fig5_redo_auc_L{L}.csv', dtype={'L': str}),
                   pd.read_csv(f'{O}extra_refs/fig5_extra_auc_L{L}.csv', dtype={'L': str})], ignore_index=True)
    return d[d.id.isin(TASKS)]


def in_main(t):
    if t in CEIL:
        return 'no (ceiling)'
    return 'no (Supplementary Fig. S4b; in the means of Figs. 5e and 5f)' if t in TYPE_A and t not in TYPE_A_C else 'yes'


def judge(a, b):
    """RFind-sc (a) against a method (b): the CI of the difference (same bootstrap draws; columns diff_vs_RFind, d_lo, d_hi of b)"""
    if pd.isna(b.AUC):
        return ''
    if b.d_lo > 0 and b.diff_vs_RFind > 0.01:
        return 'RFind-sc higher'
    if b.d_hi < 0 and b.diff_vs_RFind < -0.01:
        return 'RFind-sc lower'
    return 'not distinguished'


def table(d, L):
    rows = []
    for t in TASKS:
        if L != '500' and t in TYPE_B:
            continue                                                   # gene sets have no length to vary
        main, extra = methods_of(t)
        x = d[d.id == t].set_index('method')
        if 'RFind-sc' not in x.index:
            continue
        rf = x.loc['RFind-sc']
        for m in main + extra:
            if m not in x.index:
                continue
            r = x.loc[m]
            res = '' if m == 'RFind-sc' else ('not judged (ceiling task)' if t in CEIL else judge(rf, r))
            rows.append({'task id': t, 'task': TASKNAME[t], 'type': 'B (gene set)' if t in TYPE_B else 'A (list with direction and size)',
                         'in the main figure': in_main(t),
                         'method': DISPLAY.get(m, m), 'method in the main figure': 'yes' if m in main else 'no (this table only)',
                         'form of the reference given': 'gene set (weight 1 for decoupleR)' if t in TYPE_B else FORM.get(m, ''),
                         'reference length per direction': 'gene set' if t in TYPE_B else L,
                         'reference genes up': r.ref_up, 'reference genes down': r.ref_down, 'positive cells': r.n_pos, 'negative cells': r.n_neg,
                         'AUC': r.AUC, 'AUC 95% CI low': r.lo, 'AUC 95% CI high': r.hi,
                         'RFind-sc minus method': np.nan if m == 'RFind-sc' else r.diff_vs_RFind,
                         'difference 95% CI low': np.nan if m == 'RFind-sc' else r.d_lo, 'difference 95% CI high': np.nan if m == 'RFind-sc' else r.d_hi,
                         'result against RFind-sc': res,
                         'note': r.note if isinstance(r.note, str) else ''})
    return pd.DataFrame(rows)


A = {L: load(L) for L in LS if os.path.exists(f'{O}results/fig5_redo_auc_L{L}.csv')}
a500 = table(A['500'], '500')
along = pd.concat([table(A[L], L) for L in LS if L != '500' and L in A], ignore_index=True)

## mean AUC over the type A tasks by reference length: the eight tasks (Fig. 5f) and the four DAM and Glucocorticoid tasks
def mean_length(tasks):
    rows = []
    for L, d in A.items():
        d = d[d.id.isin(tasks)]
        for m in ['RFind-sc'] + MAIN_A + EXTRA_A:
            g = d[(d.method == m) & d.AUC.notna()]
            rows.append({'method': DISPLAY.get(m, m), 'in the main figure': 'yes' if m in ['RFind-sc'] + MAIN_A else 'no', 'L': L,
                         'mean AUC': g.AUC.mean() if len(g) == len(tasks) else np.nan})
    x = pd.DataFrame(rows).pivot_table(index=['method', 'in the main figure'], columns='L', values='mean AUC', aggfunc='first')
    x = x.reindex(columns=[L for L in LS if L in A]); x.columns = [f'L = {c}' for c in x.columns]
    order = [DISPLAY.get(m, m) for m in ['RFind-sc'] + MAIN_A + EXTRA_A]
    x = x.reset_index(); x['o'] = x.method.map({m: i for i, m in enumerate(order)})
    return x.sort_values('o').drop(columns='o').reset_index(drop=True)


length8, length4 = mean_length(TYPE_A), mean_length(TYPE_A_MAIN)
length_notes = [['Mean AUC over the eight type A tasks by the reference length per direction (Fig. 5f): the four DAM and Glucocorticoid '
                 'tasks (DAM, and Glucocorticoid in the oligodendrocytes, astrocytes and endothelial cells) and the four PS2APP and Forced '
                 'swim tasks (the same cells with the PS2APP and the Forced swim references).'],
                ['Mean AUC over the four DAM and Glucocorticoid tasks by the reference length per direction.']]

## counts of the comparisons (main-figure methods; L = 500 and all genes; type A by task set)
cnt = []
SETS = [('A', 'the four DAM and Glucocorticoid tasks', TYPE_A_MAIN), ('A', 'the four PS2APP and Forced swim tasks', TYPE_A_ADDED),
        ('A', 'all eight type A tasks', TYPE_A), ('B', 'the three type B tasks (Fig. 5d)', TYPE_B)]
for L, tab in [('500', a500), ('all', along[along['reference length per direction'] == 'all'])]:
    for typ, tsname, ts in SETS:
        if L == 'all' and typ == 'B':
            continue
        y = tab[tab['task id'].isin(ts) & (tab['method in the main figure'] == 'yes') & (tab.method != 'RFind-sc')]
        for m, g in list(y.groupby('method', sort=False)) + [('all methods', y)]:
            cnt.append({'reference length per direction': L if typ == 'A' else 'gene set', 'type': typ, 'tasks': tsname, 'method': m,
                        'comparisons': len(g),
                        'RFind-sc higher': (g['result against RFind-sc'] == 'RFind-sc higher').sum(),
                        'not distinguished': (g['result against RFind-sc'] == 'not distinguished').sum(),
                        'RFind-sc lower': (g['result against RFind-sc'] == 'RFind-sc lower').sum()})
cnt = pd.DataFrame(cnt)

bc = pd.read_csv(f'{O}results/fig5_redo_beyondcell_L500.csv') if os.path.exists(f'{O}results/fig5_redo_beyondcell_L500.csv') else pd.DataFrame()
if len(bc):                                                              # cells without a Beyondcell score, set to 0 as in the benchmark
    bc = bc[bc.id.isin(TASKS) & (bc.rule == 'as_zero')]
    bc = pd.DataFrame({'task id': bc.id, 'task': bc.id.map(TASKNAME), 'cells': bc.n_cells, 'positive cells': bc.n_pos,
                       'cells without a Beyondcell score': bc.n_nan, 'fraction of cells without a score': bc.frac_nan,
                       'Beyondcell AUC (these cells set to 0)': bc.AUC, 'AUC 95% CI low': bc.lo, 'AUC 95% CI high': bc.hi,
                       'RFind-sc AUC': bc.RFind_AUC, 'RFind-sc minus Beyondcell': bc.diff_vs_RFind,
                       'difference 95% CI low': bc.d_lo, 'difference 95% CI high': bc.d_hi})
    ## the added type A tasks: the cells without a score from the per-cell file (L = 500), the AUC and the difference from the
    ## Beyondcell row of their AUC table (the same bootstrap draws as the two-column bootstrap of fig5_bench_redo.R)
    pe = pd.read_csv(f'{O}extra_refs/fig5_extra_percell_L500.csv.gz', usecols=['id', 'positive', 'Beyondcell'])
    xa = A['500'].set_index(['id', 'method'])
    ext = []
    for t in TYPE_A_ADDED:
        q = pe[pe.id == t]; r = xa.loc[(t, 'Beyondcell')]; nn = int(q.Beyondcell.isna().sum())
        assert nn == r.n_na and len(q) == r.n_pos + r.n_neg, t
        ext.append({'task id': t, 'task': TASKNAME[t], 'cells': len(q), 'positive cells': int(q.positive.sum()),
                    'cells without a Beyondcell score': nn, 'fraction of cells without a score': nn / len(q),
                    'Beyondcell AUC (these cells set to 0)': r.AUC, 'AUC 95% CI low': r.lo, 'AUC 95% CI high': r.hi,
                    'RFind-sc AUC': xa.loc[(t, 'RFind-sc'), 'AUC'], 'RFind-sc minus Beyondcell': r.diff_vs_RFind,
                    'difference 95% CI low': r.d_lo, 'difference 95% CI high': r.d_hi})
    bc = pd.concat([bc, pd.DataFrame(ext)], ignore_index=True)    # the added tasks after the tasks of the benchmark file

tasks = a500[a500.method == 'RFind-sc'][['task id', 'task', 'type', 'in the main figure', 'positive cells', 'negative cells',
                                         'reference genes up', 'reference genes down']].copy()
tasks.insert(2, 'reference', tasks['task id'].map({'H1': 'DAM (Keren-Shaul et al., 2017)',
                                                    'DA1': 'PS2APP (PS2APP vs non-transgenic cortical microglia, 13 months, bulk; Srinivasan et al., 2016)',
                                                    'AS1': 'Forced swim (6-min forced swim, 3 h, ventral hippocampus, bulk; von Ziegler et al., 2022)',
                                                    'AS2': 'Forced swim (6-min forced swim, 3 h, ventral hippocampus, bulk; von Ziegler et al., 2022)',
                                                    'AS3': 'Forced swim (6-min forced swim, 3 h, ventral hippocampus, bulk; von Ziegler et al., 2022)', 'H2': 'Interferon (MSigDB Hallmark)', 'H3': 'Interferon (MSigDB Hallmark)',
                                                    'H4': 'Maturity (Matcovitch-Natan et al., 2016)', 'H5': 'G2M (MSigDB Hallmark)',
                                                    'D1': 'Glucocorticoid (corticosterone 3 h, hippocampus, bulk)', 'D2': 'Glucocorticoid (corticosterone 3 h, hippocampus, bulk)',
                                                    'D3': 'Glucocorticoid (corticosterone 3 h, hippocampus, bulk)'}))
tasks.insert(4, 'how the positive cells were defined', tasks['task id'].map({'H1': 'experimental condition (injection)', 'H2': "the original authors' cluster",
                                                                            'H3': 'experimental condition (age)', 'H4': 'experimental condition (age)',
                                                                            'H5': "the original authors' cluster", 'D1': 'experimental condition (stress)',
                                                                            'D2': 'experimental condition (stress)', 'D3': 'experimental condition (stress)'}))
tasks.insert(5, 'basis in the original study', tasks['task id'].map({
    'H1': 'injury-responsive microglia share the DAM genes (Hammond et al., 2019)', 'H2': 'cluster 9 of the original study, high in interferon-response genes',
    'H3': 'the aged microglia up-regulate interferon-response genes (Hammond et al., 2019)',
    'H4': 'embryonic and early postnatal microglia differ from adult microglia (Hammond et al., 2019)',
    'H5': 'cluster 2b of the original study, a proliferating cluster whose top markers include Ube2c and Birc5 (Hammond et al., 2019)',
    'D1': 'among the non-neuronal cells, one of the three cell types with the largest stress-induced changes in expression (De Donno et al., 2026, Fig. 2b, c)',
    'D2': 'among the non-neuronal cells, one of the three cell types with the largest stress-induced changes in expression (De Donno et al., 2026, Fig. 2b, c)',
    'D3': 'among the non-neuronal cells, one of the three cell types with the largest stress-induced changes in expression (De Donno et al., 2026, Fig. 2b, c)'}))
for t, p0 in PAIR.items():                                         # the added tasks: the cells and their definition as the paired task
    tasks.loc[tasks['task id'] == t, 'how the positive cells were defined'] = tasks.loc[tasks['task id'] == p0, 'how the positive cells were defined'].iloc[0]
    tasks.loc[tasks['task id'] == t, 'basis in the original study'] = (f'the cells of the {TASKNAME[p0].split(":")[0]} task: '
                                                                       + tasks.loc[tasks['task id'] == p0, 'basis in the original study'].iloc[0])
med = a500[(a500['method in the main figure'] == 'yes')].groupby('task id').AUC.median()
tasks['median AUC of the main-figure methods (L = 500)'] = tasks['task id'].map(med)

## ============ sheets Best length and Mean AUC differences ============
YN = {True: 'yes', False: 'no'}
RULE = ('Result against RFind-sc: RFind-sc higher or lower when the 95% CI of the difference excludes 0 and the difference exceeds 0.01; '
        'otherwise not distinguished. Comparable: the whole 95% CI of the difference lies within -0.01 to +0.01 (a subset of not distinguished).')


## Best length: the eight type A tasks (fig5_best_length_8tasks.py), and leave one task out within the four DAM and Glucocorticoid tasks
bl = pd.read_csv(O + 'results/fig5_best_length_8tasks.csv', dtype={'RFind-sc length': str, 'method length': str})
bl4 = pd.read_csv(O + 'results/fig5_best_length.csv', dtype={'RFind-sc length': str, 'method length': str})


def best_table(x):
    return pd.DataFrame({
        'version': x['version'], 'task id': x['task id'], 'task': x['task id'].map(TASKNAME), 'task in': x['task id'].map(TSETNAME),
        'method': x.method.map(lambda m: DISPLAY.get(m, m)),
        'method in the main figure': x['method set'].map({'main': 'yes', 'table only': 'no (this table only)'}),
        'form of the reference given': x.method.map(lambda m: FORM.get(m, '')),
        'RFind-sc length per direction': x['RFind-sc length'], 'RFind-sc AUC': x['RFind-sc AUC'],
        'method length per direction': x['method length'], 'method AUC': x['method AUC'],
        'RFind-sc minus method': x['RFind-sc minus method'],
        'difference 95% CI low': x['difference 95% CI low'], 'difference 95% CI high': x['difference 95% CI high'],
        'result against RFind-sc': x['result'], 'comparable (CI within -0.01 to +0.01)': x['comparable'].map(YN),
        'length re-selected in each resample: difference 95% CI low': x['reselected difference 95% CI low'],
        'length re-selected in each resample: difference 95% CI high': x['reselected difference 95% CI high'],
        'length re-selected in each resample: result against RFind-sc': x['reselected result'],
        'length re-selected in each resample: comparable': x['reselected comparable'].map(YN),
        'share of resamples in which RFind-sc chose its length': x['RFind-sc share of resamples at its point-estimate length']})


def best_sorted(x):                                                 # main methods first, then table-only; version, then task order
    x = x.copy()
    x['o'] = x['method set'].map({'main': 0, 'table only': 1}); x['v'] = x.version.map({'per task': 0, 'leave one task out': 1})
    x = x.sort_values(['v', 'o'], kind='stable').drop(columns=['o', 'v'])
    return best_table(x).reset_index(drop=True)


assert set(bl.version) == {'per task', 'leave one task out'} and set(bl['task id']) == set(TYPE_A)
pt8, pt4 = bl[bl.version == 'per task'], bl4[bl4.version == 'per task']   # per task: the rows of the DAM and Glucocorticoid tasks are those of the 4-task file
for c in pt4.columns:
    a_, b_ = pt8[pt8['task id'].isin(TYPE_A_MAIN)][c].reset_index(drop=True), pt4[c].reset_index(drop=True)
    assert (np.allclose(a_, b_, rtol=0, atol=1e-12) if pd.api.types.is_float_dtype(b_) else (a_.astype(str) == b_.astype(str)).all()), c
best = best_sorted(bl)
best4 = best_sorted(bl4[(bl4.version == 'leave one task out') & (bl4['task set'] == 'main')])
best_note = [
    'Each method at its best reference length in the eight type A tasks (the four DAM and Glucocorticoid tasks: DAM, and Glucocorticoid in '
    'the oligodendrocytes, astrocytes and endothelial cells; the four PS2APP and Forced swim tasks: the same cells with the PS2APP and the '
    'Forced swim references; fig5_best_length_8tasks.py). The length per direction (100, 250, 500, 1,000, 2,000 or all genes) is chosen '
    'for each method separately, RFind-sc included. Per task: the length with the highest AUC in the task. Leave one task out: for each '
    'task, the length with the highest mean AUC over the other seven tasks.',
    'The length is chosen on the same cells on which the AUC is measured (per task), so these AUCs are upper bounds, alike for '
    'every method; leave one task out does not use the task itself for the choice. 95% CIs: bootstrap over mice (all cells of each drawn '
    'mouse), 1,000 resamples, the same resamples for every method and length. Columns "length re-selected in each resample": the length is chosen again in every '
    'resample (for every method, RFind-sc included); the difference shown is still the point-estimate difference, only the CI changes. ' + RULE]
best4_note = ['Leave one task out within the four DAM and Glucocorticoid tasks (fig5_best_length.py): for each task, the length with the '
              'highest mean AUC over the other three of these tasks. The per-task rows of these tasks are those of the table above.']

## Mean AUC differences: the mean over the eight type A tasks (results/fig5_mean_auc_diff_8tasks.csv) and over the four DAM and
## Glucocorticoid tasks (results/fig5_mean_auc_diff.csv)
def meandiff_table(md):
    return pd.DataFrame({
    'reference length': md.comparison, 'method': md.method.map(lambda m: DISPLAY.get(m, m)),
    'method in the main figure': md['method set'].map({'main': 'yes', 'table only': 'no (this table only)'}),
    'form of the reference given': md.method.map(lambda m: FORM.get(m, '')),
    'mean AUC, RFind-sc': md['mean AUC RFind-sc'], 'mean AUC, method': md['mean AUC method'],
    'RFind-sc minus method (mean AUC)': md.difference, 'difference 95% CI low': md['CI low'], 'difference 95% CI high': md['CI high'],
    'result against RFind-sc': md.result, 'comparable (CI within -0.01 to +0.01)': md.comparable.map(YN)})


meandiff, meandiff4 = meandiff_table(pd.read_csv(O + 'results/fig5_mean_auc_diff_8tasks.csv')), meandiff_table(pd.read_csv(O + 'results/fig5_mean_auc_diff.csv'))
meandiff_note = [
    'Mean AUC over the eight type A tasks (Figs. 5e, f: DAM and Glucocorticoid in the oligodendrocytes, astrocytes and endothelial cells, and '
    'the same cells with the PS2APP and the Forced swim references), RFind-sc minus each method (fig5_best_length_8tasks.py). Fixed '
    'length: every method at the same length per direction. Best length per task: each method (RFind-sc included) at its own best length '
    'in each task. Leave one task out: each task at '
    'the length that is best over the other seven tasks. The best lengths are chosen on the same cells (upper bounds, alike for every method).',
    '95% CIs: bootstrap over mice (all cells of each drawn mouse), 1,000 resamples, the same draws for the tasks with the same mice; in each '
    'resample the mean over the eight tasks of the resampled AUCs. ' + RULE]
meandiff4_note = ['Mean AUC over the four DAM and Glucocorticoid tasks, RFind-sc minus each method (fig5_best_length.py; the same comparisons, '
                  'leave one task out over the other three tasks).']

## Mean over the tasks of each type (Fig. 5e; results/fig5_task_means.csv of fig5_task_means.py)
tm = pd.read_csv(O + 'results/fig5_task_means.csv')
tm['result'] = tm.result.fillna('')
tm_cols = {'method': 'method', 'mean AUC': 'mean AUC', 'mean AUC 95% CI low': 'mean AUC 95% CI low', 'mean AUC 95% CI high': 'mean AUC 95% CI high',
           'RFind-sc minus method': 'RFind-sc minus method (mean AUC)', 'difference 95% CI low': 'difference 95% CI low',
           'difference 95% CI high': 'difference 95% CI high', 'result': 'result against RFind-sc'}
tmA, tmB = [tm[tm.type == k][list(tm_cols)].rename(columns=tm_cols) for k in ('A', 'B')]
tmA_note = ['Mean AUC over the eight type A tasks with references of 500 genes per direction (Fig. 5e, left: the four DAM and '
            'Glucocorticoid tasks and the four PS2APP and Forced swim tasks), and RFind-sc minus each method (fig5_task_means.py).',
            '95% CIs: bootstrap over mice (all cells of each drawn mouse), 1,000 resamples, the same draws for the tasks with the same mice and '
            'for every method; in each resample the mean over the tasks of the resampled AUCs. ' + RULE]
tmB_note = ['Mean AUC over the three type B tasks (Fig. 5e, right: the Hallmark gene sets of Fig. 5d), and RFind-sc minus each method; '
            'the same bootstrap (H2 and H5, defined by clusters, draw from all mice).']


def write_blocks(w, name, blocks):
    """blocks: list of (lines of text above the table, table); one blank row between blocks; widths from the tables only"""
    row = 0
    widths = {}
    for lines, df in blocks:
        for ln in lines:
            pd.DataFrame([[ln]]).to_excel(w, sheet_name=name, index=False, header=False, startrow=row)
            row += 1
        df.to_excel(w, sheet_name=name, index=False, startrow=row)
        for j, c in enumerate(df.columns):
            vals = [len(str(c))] + [len(str(v)) for v in df[c] if not (isinstance(v, float) and np.isnan(v))]
            widths[j] = max(widths.get(j, 0), max(vals))
        row += len(df) + 2
    ws = w.sheets[name]
    from openpyxl.utils import get_column_letter
    for j, wd in widths.items():
        ws.column_dimensions[get_column_letter(j + 1)].width = min(60, max(10, wd + 2))


readme = [
    ('Supplementary Table 2', 'Capability benchmark against existing scoring methods (Fig. 5 and Supplementary Fig. S4).'),
    ('Tasks', 'Twelve tasks: seven with the data and references of Fig. 2 and Fig. 3 (type A, the four DAM and Glucocorticoid tasks: DAM, '
              'and Glucocorticoid in the oligodendrocytes, astrocytes and endothelial cells; type B, Fig. 5d: two Interferon tasks and G2M); '
              'the four PS2APP and Forced swim tasks (type A) on the same cells as the DAM and Glucocorticoid tasks, with a second reference '
              'of each state (PS2APP for the DAM task; Forced swim for the three Glucocorticoid tasks); and the Maturity task, left out of the '
              'main figure because the median AUC of the methods was 0.99 or more (ceiling). Fig. 5c draws the DAM task and the '
              'Glucocorticoid task in the astrocytes, and Supplementary Fig. S4b the other six type A tasks; Figs. 5e and 5f average the '
              'eight type A tasks. Positive and negative cells are defined by the '
              'experimental condition or the cluster of the original study, not by any score of this study.'),
    ('Methods', 'Main figure: RFind-sc, AddModuleScore, AUCell, UCell (signed signature), singscore, ssGSEA, VISION, Beyondcell, ULM of decoupleR '
                'and scFOCAL-type (implemented in this study from the published formula; type A only). This table only: UCell with the up genes only, '
                'WMEAN and ORA of decoupleR. '
                "Column \"form of the reference given\": A, ranked lists or weights; A', names of up- and "
                "down-regulated genes; A'', names of up-regulated genes."),
    ('References given', 'Gene sets as they are (weight 1 for decoupleR). Type A references: the top L genes per direction (L = 500 in Fig. 5b-d; '
                         '100, 250, 500, 1,000, 2,000 and all genes in Fig. 5f), the same genes for every method.'),
    ('AUC and result', 'Mann-Whitney AUC of the score between the positive and negative cells, over all task cells; 95% CIs of the AUC and of the '
                       'difference from RFind-sc by a bootstrap over mice (all cells of each drawn mouse; mice drawn within the positive and the negative '
                       'group, or from all mice in the tasks defined by a cluster; 1,000 resamples, the same resamples for every method). Result against '
                       'RFind-sc, task by task: RFind-sc higher or lower when the CI of the difference excludes 0 and the difference exceeds 0.01; '
                       'otherwise not distinguished (the difference is 0.01 or less or cannot be told apart; this does not show that the two are at the same level). The '
                       'ceiling task is not judged.'),
    ('Sheets', 'Tasks; AUC, L = 500 (Fig. 5b-d and Supplementary Fig. S4b); Mean over tasks (Fig. 5e); AUC, other lengths (Fig. 5f; type A tasks); '
               'Reference length (mean AUC over the eight type A tasks, Fig. 5f, and over the four DAM and Glucocorticoid tasks); Comparisons (numbers of '
               'comparisons in which RFind-sc was higher, not distinguished or lower, by task set); Beyondcell missing '
               '(the number and fraction of cells to which Beyondcell gives no score, and the AUC with these cells set to 0, as in the main '
               'analysis).'),
    ('Sheet: Mean over tasks', 'Mean over tasks (Fig. 5e: the mean AUC of each method over the eight type A tasks, 500 genes per '
                               'direction, and over the three type B tasks, and RFind-sc minus each method, with 95% CIs and the result).'),
    ('Sheet: Best length', 'Best length (each method, RFind-sc included, at its own best reference length in the eight type A tasks: per task '
                           'and leave one task out over the other seven tasks; below, leave one task out within the four DAM and '
                           'Glucocorticoid tasks; the '
                           'length is chosen on the same cells, so the AUCs are upper '
                           'bounds alike for every method; with the CI and result when the length is chosen again in each resample; comparable = '
                           'the whole CI of the difference within -0.01 to +0.01).'),
    ('Sheet: Mean AUC differences', 'Mean AUC differences (mean AUC over the eight type A tasks, and below over the four DAM and '
                                    'Glucocorticoid tasks, RFind-sc minus each method, with 95% CIs, at '
                                    'each fixed length, at the best length per task and leave one task out; comparable as above).'),
]
## ============ sheet Run time and memory (fig5_runtime.sh; only when its output exists) ============
RT, RTENV = f'{O}results/fig5_runtime.csv', f'{O}results/fig5_runtime_env.tsv'
rt_blocks = []
if os.path.exists(RT):
    rt = pd.read_csv(RT).drop(columns=['Exit code'], errors='ignore')
    env = pd.read_csv(RTENV, sep='\t', dtype=str) if os.path.exists(RTENV) else pd.DataFrame(columns=['item', 'value'])
    rt_note = ['Run time (seconds of scoring, not of reading the data) and peak memory (maximum resident set size of the process; '
               'single-core runs only, because the workers of multi-core runs are separate processes) of scoring the 75,572 cells of the '
               'data of Hammond et al. (or the first 1,000 or 10,000 cells of a fixed random order) against the DAM reference (the top 500 '
               'genes per direction or the whole list, among the genes in the counts), and of the knock-out x knock-out matching of '
               'Supplementary Fig. S3 (2,048 profiles of TH Prkcd Grin2c Glut, 7 worker processes).',
               'Computer and versions: ' + '; '.join(f'{a} {b}' for a, b in env[['item', 'value']].itertuples(index=False)) + '.']
    rt_blocks = [(rt_note, rt)]
    readme.append(('Sheet: Run time and memory', 'Run time and peak memory of RFind-sc (R and Python packages), ULM and UCell for 1,000 to '
                   '75,572 cells, and of the knock-out x knock-out matching of Supplementary Fig. S3, with the computer and the versions used.'))
sheets = {'README': pd.DataFrame(readme, columns=['item', 'description']), 'Tasks': tasks, 'AUC, L = 500': a500, 'AUC, other lengths': along,
          'Comparisons': cnt, 'Beyondcell missing': bc}
added = {'Mean over tasks': [(tmA_note, tmA), (tmB_note, tmB)],
         'Reference length': [(length_notes[0], length8), (length_notes[1], length4)],
         'Best length': [(best_note, best), (best4_note, best4)],
         'Mean AUC differences': [(meandiff_note, meandiff), (meandiff4_note, meandiff4)]}
ORDER = ['README', 'Tasks', 'AUC, L = 500', 'Mean over tasks', 'AUC, other lengths', 'Reference length', 'Comparisons',
         'Beyondcell missing', 'Best length', 'Mean AUC differences']
if rt_blocks:
    added['Run time and memory'] = rt_blocks
    ORDER.append('Run time and memory')
if len(bc):                                                              # the AUC with the unscored cells set to 0 is that of the main table
    chk = a500[a500.method == 'Beyondcell'].set_index('task id')
    for c, k in [('Beyondcell AUC (these cells set to 0)', 'AUC'), ('AUC 95% CI low', 'AUC 95% CI low'), ('AUC 95% CI high', 'AUC 95% CI high')]:
        assert np.allclose(bc.set_index('task id')[c].reindex(chk.index), chk[k]), c
NOKEY = lambda df: df.drop(columns=['task id'], errors='ignore')        # task ids are internal keys
sheets = {k: NOKEY(v) for k, v in sheets.items()}
added = {k: [(lines, NOKEY(df)) for lines, df in b] for k, b in added.items()}
out = f'{O}SuppTable2.xlsx'
with pd.ExcelWriter(out, engine='openpyxl') as w:
    for name in ORDER:
        if name in added:
            write_blocks(w, name, added[name]); continue
        df = sheets[name]
        df.to_excel(w, sheet_name=name, index=False)
        ws = w.sheets[name]
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = min(60, max(10, max(len(str(c.value)) if c.value is not None else 0 for c in col) + 2))
print('wrote', out, '|', ', '.join(f'{k}: {len(v)} rows' for k, v in sheets.items()), '|',
      ', '.join(f'{k}: ' + ' + '.join(str(len(df)) for _, df in b) + ' rows' for k, b in added.items()))
print(cnt[cnt.method == 'all methods'].to_string(index=False))
