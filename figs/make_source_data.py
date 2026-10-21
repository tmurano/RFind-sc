"""Source Data for Figs. 2-5 and Supplementary Fig. S4: one Excel file per figure, one sheet per panel, with every statistic shown
in the figure or quoted in the text (effect sizes, P values, correlations, AUC, n of cells and mice) and the per-mouse values drawn as points.

Reads only the outputs of the figure scripts (figs/*/output/) and the per-cell scores of the pipelines; computes nothing new
except the per-mouse means drawn in Figs. 2c, 2d, 3c and 3d, and, for Fig. 5, what the figure scripts compute from the same files
while drawing (the ROC curve of Fig. 5b, the means of Fig. 5f, the higher / not distinguished / lower rule). Fig. 5
(figs/Fig5_benchmark/analysis/fig5_rebuild_figure_draft.py; a, table of methods; b, the DAM example; c, type A, the DAM task and the
Glucocorticoid task in the astrocytes; d, type B; e, the mean over the tasks of each type, from results/fig5_task_means.csv; f, reference
length over the eight type A tasks) and Supplementary Fig. S4 (figS4_capability_map.py; a, table of all tasks; b, the six type A tasks
not drawn in Fig. 5c) follow those scripts: their lists of methods and tasks are
copied below and must be changed with them. Only the tasks and methods drawn in these two figures are written; the other results
of the benchmark that the paper reports (the Maturity task, the methods not drawn in these figures) are in Supplementary Table 3
(figs/Fig5_benchmark/analysis/fig5_redo_tables.py).

Run from the project root:  python figs/make_source_data.py      (needs pandas, numpy, openpyxl)
Output: submission/Source_Data_Fig2.xlsx, Source_Data_Fig3.xlsx, Source_Data_Fig4.xlsx, Source_Data_Fig5.xlsx,
        Source_Data_Supplementary_FigS4.xlsx (the names used for the upload; submission/ is not tracked)
"""
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path('submission')                                     # the files uploaded to the journal (not tracked)
F2, F3, F4 = Path('figs/Fig2_Hammond/output'), Path('figs/Fig3_DeDonno/output'), Path('figs/Fig4_Perturbseq/output')
CLUSTERS = ['1', '2a', '2b', '2c', '3', '4', '5', '6', '7a', '7b', '7c', '8', '9']
IDX2 = {'Maturity': 'Maturity-index', 'DAM': 'DAM-index', 'Interferon': 'Interferon-index', 'G2M': 'G2M-index'}
IDX3 = {'CORT': 'Glucocorticoid-index', 'KA6h': 'Seizure-index', 'NFkB': 'NF-κB-index', 'UPR': 'UPR-index'}
CT3 = ['DG Glut', 'CA1 Glut', 'CA3 Glut', 'GABA', 'Astro', 'Microglia', 'Oligodendrocytes', 'Endothelial']


def cohend(a, b):
    """Cohen's d of b - a with the pooled standard deviation (as in the figure scripts)."""
    n1, n2 = len(a), len(b)
    sp = np.sqrt(((n1 - 1) * np.var(a, ddof=1) + (n2 - 1) * np.var(b, ddof=1)) / (n1 + n2 - 2))
    return (np.mean(b) - np.mean(a)) / sp


def write(path, readme, sheets):
    """README sheet (one line per row) first, then the panels (a list of tables = side by side in one sheet)."""
    with pd.ExcelWriter(path, engine='openpyxl') as xw:
        pd.DataFrame({'README': readme}).to_excel(xw, sheet_name='README', index=False)
        for name, df in sheets.items():
            if isinstance(df, pd.DataFrame):
                df.to_excel(xw, sheet_name=name, index=False)
            else:                                                  # several tables in one sheet, side by side, one empty column apart
                col = 0
                for d in df:
                    d.to_excel(xw, sheet_name=name, index=False, startcol=col); col += d.shape[1] + 1
    print('wrote', path)


def fig2():
    s = pd.read_csv('data/Hammond2019/output/scores_hammond_9axes.csv', dtype={'cluster': str})
    meta = pd.read_csv('data/Hammond2019/raw/cell_meta_hammond.csv', usecols=['cell_id', 'age'])
    s = s[s.cluster.isin(CLUSTERS)].merge(meta, on='cell_id', how='left')
    s['mouse'] = s.cell_id.str.replace(r'_[^_]+$', '', regex=True)

    st = pd.read_csv(F2 / 'Fig2_stats_table.csv', dtype={'cluster': str})
    vio = pd.DataFrame({
        'Panel': np.where(st.axis.isin(['Maturity', 'DAM']), '2c', '2d'),
        'Index': st.axis.map(IDX2), 'Cluster': st.cluster,
        'Cells (n)': st.n, 'Mice with >= 20 cells (n)': st.n_animals,
        'Median index': st['median'], 'Mean index': st['mean'], 'SD of index': st['sd'],
        'd0 (mean/SD; number above the violin)': st.d0,
        "Cohen's d, cluster vs all other clusters": st.d_vs_rest, 'AUC, cluster vs all other clusters': st.auc_vs_rest})
    vio['Index'] = pd.Categorical(vio.Index, list(IDX2.values()))
    vio = vio.sort_values(['Panel', 'Index', 'Cluster'], key=lambda c: c.map(CLUSTERS.index) if c.name == 'Cluster' else c)

    pm = (s.melt(id_vars=['cluster', 'mouse'], value_vars=list(IDX2), var_name='axis', value_name='v')
            .groupby(['axis', 'cluster', 'mouse']).v.agg(['mean', 'size']).reset_index())
    pm = pm[pm['size'] >= 20]
    pts = pd.DataFrame({'Panel': np.where(pm.axis.isin(['Maturity', 'DAM']), '2c', '2d'), 'Index': pm.axis.map(IDX2),
                        'Cluster': pm.cluster, 'Mouse (sample)': pm.mouse, 'Cells (n)': pm['size'],
                        'Mean index of the mouse (point)': pm['mean'].round(4)})

    # Figs. 2e, f: Pearson r over cells (written in the panels); group names as in Results 2
    # (age keys of the metadata: 'E14' = embryonic day 14.5 of Hammond et al., written E14.5 in the labels)
    s['grp'] = s.age.replace({'P100_Injured_WM': 'LPC', 'P100_Control_WM': 'Saline'})
    early, late, inj = ['E14', 'P4/P5'], ['P30', 'P100', 'P540'], ['Saline', 'LPC']
    groups = [('e/f left (development)', 'All development (E14.5–P540)', early + late),
              ('e/f left (development)', 'Early development (E14.5–P4/P5)', early),
              ('e/f left (development)', 'Maturation/ageing (P30–P540)', late),
              ('e/f right (white matter, P100)', 'Saline and LPC', inj),
              ('e/f right (white matter, P100)', 'Saline', ['Saline']), ('e/f right (white matter, P100)', 'LPC', ['LPC'])]
    rows = []
    for pan, (xa, ya) in [('2e', ('Maturity', 'DAM')), ('2f', ('DAM', 'Interferon'))]:
        for side, lab, gs in groups:
            x = s[s.grp.isin(gs)]
            rows.append({'Panel': pan + ' ' + side.split(' ', 1)[1], 'x': IDX2[xa], 'y': IDX2[ya], 'Group': lab,
                         'Cells (n)': len(x), 'Mice (n)': x.mouse.nunique(),
                         'Pearson r over cells': round(np.corrcoef(x[xa], x[ya])[0, 1], 3)})
    cor = pd.DataFrame(rows)


    write(OUT / 'Source_Data_Fig2.xlsx', [
        'Source Data for Fig. 2 (microglia of Hammond et al., 2019; 13 microglial clusters, Mono/Mac excluded).',
        'Index = the RFind-sc index of each cell (signed similarity to the reference; not a test of the cell).',
        'Fig2c-d_stats: per cluster; d0 = mean/SD of the index over the cells of the cluster (the number above each violin).',
        'Fig2c-d_points: per-mouse means drawn as white points (mice contributing >= 20 cells to the cluster).',
        'Fig2e-f_correlations: Pearson r between two indices over the cells of each group (the r written in the panels).',
        'Source tables: figs/Fig2_Hammond/output/Fig2_stats_table.csv, fig2_scatter_correlations.csv; '
        'per-cell scores: data/Hammond2019/output/scores_hammond_9axes.csv.'],
        {'Fig2c-d_stats': vio, 'Fig2c-d_points': pts, 'Fig2e-f_correlations': cor})


def fig3():
    st = pd.read_csv(F3 / 'fig3_violin_stats.csv')
    nc = pd.read_csv(F3 / 'Fig3_celltype_ncells.csv').set_index('cell_type')
    st['Index'] = st.axis.map(IDX3)
    vio = pd.DataFrame({
        'Panel': np.where(st.axis.isin(['CORT', 'KA6h']), '3c', '3d'), 'Index': st.Index, 'Cell type': st.cell_type,
        'No Stress cells (n)': st.cell_type.map(nc.n_NoStress), 'Stress cells (n)': st.cell_type.map(nc.n_Stress),
        'No Stress mice (n)': st.n_mice_NoStress, 'Stress mice (n)': st.n_mice_Stress,
        "Cohen's d, Stress - No Stress (per-mouse means; number above the violins)": st.d.round(3),
        'P, two-sided Wilcoxon rank-sum on per-mouse means (nominal)': st.p.round(4)})
    vio['Cell type'] = pd.Categorical(vio['Cell type'], CT3)
    vio = vio.sort_values(['Panel', 'Index', 'Cell type'])

    s = pd.read_csv('data/DeDonno2026/output/scores_dedonno.csv')
    s = s[(s.line == 'Ctrl') & s.cell_type.isin(CT3)]
    pm = (s.melt(id_vars=['cell_type', 'condition', 'donor_id'], value_vars=list(IDX3), var_name='axis', value_name='v')
            .groupby(['axis', 'cell_type', 'condition', 'donor_id']).v.agg(['mean', 'size']).reset_index())
    pm = pm[pm['size'] >= 10]
    pts = pd.DataFrame({'Panel': np.where(pm.axis.isin(['CORT', 'KA6h']), '3c', '3d'), 'Index': pm.axis.map(IDX3),
                        'Cell type': pm.cell_type, 'Condition': pm.condition, 'Mouse': pm.donor_id, 'Cells (n)': pm['size'],
                        'Mean index of the mouse (point)': pm['mean'].round(4)})

    cr = pd.read_csv(F3 / 'fig3_axis_correlations.csv')
    cor = pd.DataFrame({'Panel': np.where(cr.pair.str.startswith('Glucocorticoid'), '3e', '3f'), 'Pair': cr.pair.str.replace('_vs_', ' vs '),
                        'Cell type': cr.cell_type, 'Cells (n)': cr.n, 'Pearson r, No Stress cells': cr.r_NoStress,
                        'Pearson r, Stress cells': cr.r_Stress, 'Pearson r, all cells': cr.r_pooled,
                        'Spearman rho, all cells': cr.spearman_pooled})

    write(OUT / 'Source_Data_Fig3.xlsx', [
        'Source Data for Fig. 3 (hippocampus of De Donno et al., 2026; control line, No Stress vs Stress; 8 representative cell types).',
        "Fig3c-d_stats: Cohen's d (pooled SD) of Stress minus No Stress over per-mouse means (mice with >= 10 cells), as written above the "
        'violins; the summary panels of Figs. 3e and 3f plot the same d. P = nominal two-sided Wilcoxon rank-sum on the per-mouse means '
        '(for reference; smallest possible value for 6 vs 6 mice = 0.0022).',
        'Fig3c-d_points: per-mouse means drawn as points.',
        'Fig3e-f_correlations: correlation of the two indices over the cells of each cell type (not written in the figure).',
        'Source tables: figs/Fig3_DeDonno/output/fig3_violin_stats.csv, Fig3_celltype_ncells.csv, fig3_axis_correlations.csv; '
        'per-cell scores: data/DeDonno2026/output/scores_dedonno.csv.'],
        {'Fig3c-d_stats': vio, 'Fig3c-d_points': pts, 'Fig3e-f_correlations': cor})


def fig4():
    ## Fig. 4c = CA1 Glut (fig4c_rank_tables.py CA1 CA1); the 40 knock-outs are checked against the input of that script
    xl = pd.ExcelFile(F4 / 'SuppTable_KO_rankings.xlsx')                 # Supplementary Table 1 (README sheet first)
    ca1 = pd.read_excel(xl, sheet_name=next(n for n in ('CA1 Glut', 'CA1') if n in xl.sheet_names))
    c = pd.concat([ca1.head(20).assign(**{'List in Fig. 4c': 'opposed (20 most opposed)'}),
                   ca1.tail(20).iloc[::-1].assign(**{'List in Fig. 4c': 'aligned (20 most aligned)'})])
    ax = pd.read_csv('analysis/perturbseq/results/axis_CA1_vs_CA1.csv')
    ax = ax[ax.kind == 'KO'].sort_values('rank_among_KO')             # rank_among_KO 1 = most aligned
    assert list(ca1.head(20)['Knock-out']) == list(ax.target.iloc[::-1].head(20)) and \
        list(ca1.tail(20).iloc[::-1]['Knock-out']) == list(ax.target.head(20)), 'Fig. 4c lists differ from the CA1 rankings'
    pos = pd.read_csv(F4 / 'fig4d_GR-complex_positions.csv')
    de = pd.DataFrame({'Cell type': pos.cell_type, 'Knock-out': pos.KO, 'RFind score': pos.rfind,
                       'z among knock-outs': pos.z_vs_KO, 'Rank from the opposed end (1 = most opposed)': pos.rank_from_opposed_end,
                       'Knock-outs ranked (n)': pos.n_ko,
                       'Nuclei of the knock-out (n)': pos.nuclei})
    rn = pd.read_csv(F4 / 'geneset_random_null.csv')
    scan = pd.DataFrame({'Cell type': rn.cell_type, 'Knock-outs ranked (n)': rn.n_knockouts,
                         'AUC of the five knock-outs (< 0.5 = towards the opposed end)': rn.auc_fig4d,
                         'Random sets: mean AUC': rn.null_any_mean, 'Random sets: SD of AUC': rn.null_any_sd, 'P, random sets': rn.p_any,
                         'Matched sets: mean AUC': rn.null_matched_mean, 'Matched sets: SD of AUC': rn.null_matched_sd,
                         'P, matched sets': rn.p_matched})
    gr = pd.read_csv(F4 / 'fig4e_expression_vs_score.csv')
    NM = {'DG': 'DG Glut', 'CA1': 'CA1 Glut', 'CA3': 'CA3 Glut', 'GABA': 'GABA'}
    f = pd.DataFrame({'Gene': gr.gene, 'Cell type': gr.cell_type.map(NM), 'Expression in the No Stress cells of De Donno et al. (CPM)': gr.cpm_no_stress,
                      'RFind score of the knock-out': gr.rfind_score, 'Nuclei of the knock-out (n)': gr.ko_nuclei})
    grs = pd.read_csv(F4 / 'fig4e_expression_vs_score_stats.csv')
    fs = pd.DataFrame({'Gene': grs.gene, 'x': grs.x, 'y': grs.y, 'Statistic': grs.statistic + ' (descriptive; no test)', 'Value': grs.value,
                       'Cell types (n)': grs.n_cell_types})      # no P: no test with four cell types (Methods)
    write(OUT / 'Source_Data_Fig4.xlsx', [
        'Source Data for Fig. 4 (stress response of De Donno et al. matched with the knock-outs of the in vivo Perturb-seq of Shi et al., 2026).',
        'RFind score = the signed Running Fisher score between the stress response (Delta computed from the summed counts of Stress cells with those of No Stress cells as the baseline, per cell type) and the knock-out effect; not a test.',
        'Fig4c_CA1: the 40 knock-outs listed in Fig. 4c (CA1 Glut); all knock-outs = Supplementary Table 1.',
        'Fig4d_positions: the five genes of Fig. 4d in the four cell types.',
        'Fig4d_random_sets: the five genes as a set compared with 100,000 random (and 100,000 matched) sets of five knock-outs from the library (Supplementary Table 2).',
        'Fig4e_expression: Nr3c1 and Ncoa1 in the four cell types of Fig. 4e (CPM in the summed counts of the No Stress cells of De Donno et al.; RFind score of the knock-out of the same gene); in the table to the right, the Spearman rho of the two over the four cell types for each gene (descriptive; no test, as there are only four cell types).',
        'Source tables: figs/Fig4_Perturbseq/output/SuppTable_KO_rankings.xlsx (sheet CA1 Glut), fig4d_GR-complex_positions.csv, geneset_random_null.csv, fig4e_expression_vs_score.csv, fig4e_expression_vs_score_stats.csv.'],
        {'Fig4c_CA1': c, 'Fig4d_positions': de, 'Fig4d_random_sets': scan, 'Fig4e_expression': [f, fs]})


## ---- Fig. 5 and Supplementary Fig. S4 (the capability benchmark) ----
F5 = Path('figs/Fig5_benchmark/output')
RF = 'RFind-sc'
## as in fig5_rebuild_figure_draft.py: the methods of c and of the strip of b (CMETH), of d (SETM), the tasks and the lengths
M5C = [RF, 'AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM (log2FC)', 'scFOCAL-type']
M5D = [RF, 'AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM']
TA5, TB5 = ['H1', 'D1', 'D2', 'D3'], ['H2', 'H3', 'H5']             # the DAM and Glucocorticoid tasks; type B (d)
TC5 = ['H1', 'D2']                                                  # the type A tasks drawn in c
XA5 = ['DA1', 'AS1', 'AS2', 'AS3']                                  # the PS2APP and Forced swim tasks (in the means of Figs. 5e, f)
S4B = [t for t in TA5 if t not in TC5] + XA5                        # Supplementary Fig. S4b, in the order drawn: D1, D3, DA1, AS1-AS3
TNAME5 = {'H1': 'DAM, LPC vs saline', 'D1': 'Glucocorticoid, oligodendrocytes', 'D2': 'Glucocorticoid, astrocytes',
          'D3': 'Glucocorticoid, endothelial cells', 'H2': 'Interferon, cluster 9', 'H3': 'Interferon, P540 vs P100', 'H5': 'G2M, cluster 2b',
          'DA1': 'PS2APP, LPC vs saline', 'AS1': 'Forced swim, oligodendrocytes', 'AS2': 'Forced swim, astrocytes',
          'AS3': 'Forced swim, endothelial cells'}
LS5 = ['100', '250', '500', '1000', '2000', 'all']
LAB5 = {'UCell': 'UCell (signed)', 'ULM (log2FC)': 'ULM (decoupleR), log2FC weights', 'ULM': 'ULM (decoupleR), weights of 1'}
FORM5 = {RF: 'Ranked lists', 'AddModuleScore': 'Up names only', 'AUCell': 'Up names only', 'UCell': 'Up/down names only', 'singscore': 'Up/down names only',
         'ssGSEA': 'Up/down names only', 'VISION': 'Up/down names only', 'Beyondcell': 'Up/down names only', 'ULM (log2FC)': 'log2FC weights',
         'scFOCAL-type': 'Ranked lists'}                                                   # the column 'Type A tasks' of the table a
## the table a as drawn (rows of one text = the group lines of the table)
HEAD5A = ['Method', 'As designed: built for', 'As designed: direction', 'As designed: magnitude', 'As designed: gene sets', 'Given in this study: type A tasks', 'Given in this study: type B tasks', 'What is measured']
ROWS5A = [
    ('RFind-sc', 'External ranked DE lists; gene sets', '✓', 'Rank', '✓', 'Ranked lists', 'Gene set', 'Best overlap of cell and reference lists'),
    ('Names only (membership)',),
    ('AddModuleScore', 'Gene sets', '–', '–', '✓', 'Up names only', 'Gene set', 'Set expression above control genes'),
    ('AUCell', 'Gene sets', '–', '–', '✓', 'Up names only', 'Gene set', "Set genes among the cell's top 5%"),
    ('Names and direction (up and down sets)',),
    ('UCell', 'Gene sets (up/down allowed)', '✓', '–', '✓', 'Up/down names only', 'Gene set', 'Set-gene ranks; up − down'),
    ('singscore', 'Up/down gene sets', '✓', '–', '✓', 'Up/down names only', 'Gene set', 'Mean up / down set ranks'),
    ('ssGSEA', 'Gene sets', '△', '–', '✓', 'Up/down names only', 'Gene set', 'Rank enrichment (up − down)'),
    ('VISION', 'Signatures (±1)', '✓', '–', '✓', 'Up/down names only', 'Gene set', 'Signed mean, standardized'),
    ('Beyondcell', 'Drug signatures (up/down)', '✓', '–', '✓', 'Up/down names only', 'Gene set', 'Up − down expression, scaled'),
    ('Names, direction and magnitude (ranks or weights)',),
    ('ULM (decoupleR)', 'Weighted networks', '✓', 'Weight', '✓', 'log2FC weights', 'Gene set', 't value of a linear model'),
    ('scFOCAL-type', 'Drug signatures (log2FC)', '✓', 'Rank', '–', 'Ranked lists', 'Not run', 'Rank correlation with log2FC')]
RES ={'higher': 'RFind-sc higher', 'nd': 'not distinguished', 'lower': 'RFind-sc lower'}


def verdict5(diff, lo, hi):
    """The rule of the benchmark: higher / lower when the 95% CI of the difference excludes 0 and the difference exceeds 0.01."""
    if lo > 0 and diff > 0.01:
        return 'higher'
    if hi < 0 and diff < -0.01:
        return 'lower'
    return 'nd'


def length5(L):
    return int(L) if str(L).isdigit() else 'all genes'


def auc5(A, ids, meths, lowered=False):
    """Rows of one AUC strip (Figs. 5b-d, Supplementary Fig. S4b): AUC, 95% CI, RFind-sc minus the method with its 95% CI,
    result, symbol. lowered: the lower edge goes down just enough to show every CI of a point inside the range (S4, as drawn)."""
    rows = []
    for t in ids:
        v = A[A.id == t].set_index('method')
        ok = [m for m in meths if m in v.index and np.isfinite(v.loc[m, 'AUC'])]
        ylo = min(1.012, max(v.loc[m, 'hi'] if np.isfinite(v.loc[m, 'hi']) else v.loc[m, 'AUC'] for m in ok) + 0.012) - 0.2   # axis as drawn
        if lowered:
            vv = v.loc[ok]
            while ((vv.AUC >= ylo) & (vv.lo < ylo)).any():
                ylo = vv[(vv.AUC >= ylo) & (vv.lo < ylo)].lo.min() - 0.005
        for m in meths:
            r = v.loc[m]; rf = m == RF
            res = verdict5(r.diff_vs_RFind, r.d_lo, r.d_hi)
            rows.append({'Task id': t, 'Task': TNAME5[t], 'Positive cells (n)': r.n_pos, 'Negative cells (n)': r.n_neg,
                         'Reference genes up (n)': r.ref_up, 'Reference genes down (n)': r.ref_down,
                         'Method': LAB5.get(m, m), 'Reference given': 'gene set' if t in TB5 else FORM5[m], 'AUC': r.AUC, 'AUC 95% CI low': r.lo, 'AUC 95% CI high': r.hi,
                         'AUC of RFind-sc minus AUC of the method': np.nan if rf else r.diff_vs_RFind,
                         'Difference 95% CI low': np.nan if rf else r.d_lo, 'Difference 95% CI high': np.nan if rf else r.d_hi,
                         'Result': 'RFind-sc (reference of the comparison)' if rf else RES[res],
                         'Below the axis range (triangle at the bottom edge)': 'yes' if r.AUC < ylo else 'no'})
    return pd.DataFrame(rows)


def counts5(df):
    c = df[df.Method != RF].Result.value_counts()
    return ', '.join(f'{k} {c.get(k, 0)}' for k in RES.values()) + f' (of {int(c.sum())})'


GCCELL5 = {'D1': 'oligodendrocytes', 'D2': 'astrocytes', 'D3': 'endothelial cells'}


def gc_zero5(ids):
    """The note on the Glucocorticoid tasks `ids` of a sheet: the share of cells to which signed UCell gives 0 and Beyondcell gives no
    score (these cells set to 0, as in the benchmark)."""
    BY = pd.read_csv(F5 / 'results/fig5_redo_beyondcell_L500.csv')
    BY = BY[(BY.rule == 'as_zero') & BY.id.isin(ids)].set_index('id')
    UZ = pd.read_csv(F5 / 'results/fig5_ucell_zero.csv')
    UZ = UZ.assign(z=UZ.n_cells * UZ.share_score0).groupby('id')[['z', 'n_cells']].sum()
    u0 = ' and '.join(f'{100 * UZ.loc[t, "z"] / UZ.loc[t, "n_cells"]:.2f}%' for t in ids)
    b0 = ' and '.join(f'{100 * BY.frac_nan[t]:.1f}%' for t in ids)
    cells = ' and the '.join(GCCELL5[t] for t in ids)
    many = len(ids) > 1
    return (f'Signed UCell gave 0 to {u0} of the cells of the Glucocorticoid task{"s" if many else ""} in the {cells}'
            f'{", respectively" if many else ""}; Beyondcell gave no score to {b0} (these cells set to 0).')


def fig5():
    A = pd.read_csv(F5 / 'results/fig5_redo_auc_L500.csv')
    assert {'ULM (log2FC)'} <= set(A.method)                       # otherwise the figure draws ULM with weights of 1 in its place

    a = pd.DataFrame([r if len(r) > 1 else (r[0],) + ('',) * 7 for r in ROWS5A], columns=HEAD5A)

    ## b: the DAM task (H1): the cells on the t-SNE of Fig. 2, their DAM-index, the ROC curve; the AUC strip = the H1 rows of c
    pc = pd.read_csv(F5 / 'results/fig5_redo_percell_L500.csv.gz', usecols=['id', 'cell_id', 'positive', 'group', RF])
    pc = pc[pc.id == 'H1'].reset_index(drop=True)
    xy = pd.read_csv('data/Hammond2019/raw/cell_meta_hammond.csv', usecols=['cell_id', 'tsne_x', 'tsne_y'])
    e = pc.merge(xy, on='cell_id', how='left', validate='one_to_one')
    assert e[['tsne_x', 'tsne_y']].notna().all().all()
    v = e[RF].to_numpy(dtype=float); pos = e.positive.to_numpy(dtype=bool)
    o = np.argsort(-v, kind='mergesort'); pp = pos[o]                # the ROC curve as drawn (cells from the highest DAM-index)
    e = e.iloc[o].reset_index(drop=True)
    cells = pd.DataFrame({'Cell id': e.cell_id, 'Group': np.where(pp, 'LPC (positive)', 'saline (negative)'),
                          'Sample group (Hammond et al.)': e.group, 't-SNE 1 (Fig. 2)': e.tsne_x, 't-SNE 2 (Fig. 2)': e.tsne_y,
                          'DAM-index (RFind-sc)': e[RF], 'Rank by DAM-index (1 = highest)': np.arange(1, len(e) + 1),
                          'ROC false positive rate (cells down to this one)': np.cumsum(~pp) / (~pp).sum(),
                          'ROC true positive rate (cells down to this one)': np.cumsum(pp) / pp.sum()})
    b_auc = auc5(A, ['H1'], M5C)
    vl = np.percentile(np.abs(v), 99); lo, hi = np.nanpercentile(v, [0.5, 99.5])
    c = auc5(A, TC5, M5C)
    d = auc5(A, TB5, M5D)
    assert np.allclose(b_auc.AUC.values, c[c['Task id'] == 'H1'].AUC.values)

    ALL = pd.concat([pd.read_csv(F5 / f'results/fig5_redo_auc_L{L}.csv').assign(L=L) for L in LS5]
                    + [pd.read_csv(F5 / 'extra_refs' / f'fig5_extra_auc_L{L}.csv').assign(L=L) for L in LS5], ignore_index=True)
    TA8 = TA5 + XA5                                                # e: the eight type A tasks
    rows = []
    for m in M5C:
        for L in LS5:
            s = ALL[(ALL.method == m) & (ALL.L == L) & ALL.id.isin(TA8)].set_index('id')
            r = {'Method': LAB5.get(m, m), 'Reference length (genes per direction)': length5(L)}
            r.update({f'AUC, {TNAME5[t]}': s.AUC.get(t, np.nan) for t in TA8})
            r['Mean AUC of the 8 type A tasks (point drawn)'] = s.AUC.mean() if len(s) == len(TA8) else np.nan
            rows.append(r)
    el = pd.DataFrame(rows)

    TMS = pd.read_csv(F5 / 'results/fig5_task_means.csv')            # e: the mean over the tasks of each type (fig5_task_means.py)
    me = pd.concat([TMS[(TMS.type == 'A') & TMS.method.isin(M5C)].assign(type='A, eight tasks, 500 genes per direction'),
                    TMS[(TMS.type == 'B') & TMS.method.isin(M5D)].assign(type='B, three tasks (gene sets)')], ignore_index=True)
    me = me[['type', 'tasks', 'method', 'mean AUC', 'mean AUC 95% CI low', 'mean AUC 95% CI high', 'RFind-sc minus method',
             'difference 95% CI low', 'difference 95% CI high', 'result']]
    write(OUT / 'Source_Data_Fig5.xlsx', [
        'Source Data for Fig. 5 (capability benchmark against existing scoring methods).',
        'AUC = probability that a positive cell has a higher score than a negative cell. 95% CI of the AUC and of the difference in AUC: '
        'bootstrap over the mice of the task (all cells of each drawn mouse; 1,000 resamples, the same draws for every method).',
        'Result: RFind-sc higher / lower when the 95% CI of (AUC of RFind-sc minus AUC of the method) excludes 0 and the difference exceeds '
        '0.01, otherwise not distinguished (not marked in Fig. 5; given in Supplementary Table 3).',
        'Fig5a: the table as drawn (✓ uses as it is; △ only partly, ssGSEA up − down by convention; – cannot; scFOCAL-type, implemented '
        'from the published formula of scFOCAL). '
        'As designed: what each method uses (✓ as it is; △ only partly; – cannot). Given in this study: the form of the reference in the type A tasks (up/down '
        'names only and up names only: what was left of the lists for methods that cannot use the rank or the direction) and in the type B tasks (Not run: no gene sets accepted).',
        'Fig5b (left table): the 4,891 microglia of the DAM task (LPC-injected = positive, saline-injected = negative white matter, P100) with '
        'their t-SNE coordinates (those of Fig. 2, from the metadata of Hammond et al.; the other cells of the atlas are drawn in grey) and '
        'DAM-index (RFind-sc, top 500 genes per direction), sorted from the highest DAM-index; the ROC curve joins (0, 0) and the two '
        'rates of each row in this order.',
        f'Fig5b, drawing: colour scale of the DAM-index on the t-SNE ±{vl:.4f} (99th percentile of |DAM-index|, values beyond in the end '
        f'colours); histograms: 39 equal bins from {lo:.4f} to {hi:.4f} (0.5th and 99.5th percentiles, values beyond in the end bins), '
        'density per group.',
        'Fig5b (right table): the AUC of each method in the DAM task (the same rows as the DAM task of Fig5c).',
        "Reference given (Figs. 5b-d): the form of reference given to the method (as in the table a; gene set = the Hallmark set).",
        f'Fig5c: the type A tasks drawn (DAM; Glucocorticoid in the astrocytes), top 500 genes per direction (nominal; '
        f'{int(A[A.id == "D2"].ref_up.iloc[0])} up and {int(A[A.id == "D2"].ref_down.iloc[0])} down Glucocorticoid genes remained among '
        f'the genes analysed). The other six type A tasks are in Supplementary Fig. S4b. Comparisons: {counts5(c)}.',
        gc_zero5(['D2']),
        f'Fig5d: type B tasks (MSigDB Hallmark gene sets); the methods that take gene sets (ULM with weights of 1). Comparisons: {counts5(d)}. '
        'With a single gene set, ULM and VISION rank the cells in the same order (identical AUCs).',
        'Fig5e: the mean AUC of each method over the eight type A tasks (500 genes per direction) and over the three type B tasks, and '
        'RFind-sc minus each method, with 95% CIs (bootstrap over mice, the same draws for the tasks with the same mice; drawn: the '
        'differences; result by the rule above, not marked in the figure).',
        'Fig5f: mean AUC over the eight type A tasks by reference length (genes per direction; all genes = the whole list): the four '
        'DAM and Glucocorticoid tasks and the four PS2APP and Forced swim tasks (the same cells with the PS2APP and the Forced swim '
        'references); the AUC strips of the tasks not drawn in Fig5c are in the Source Data of Supplementary Fig. S4.',
        'Source tables: figs/Fig5_benchmark/output/results/fig5_redo_auc_L<L>.csv and extra_refs/fig5_extra_auc_L<L>.csv (L = 100, 250, 500, '
        '1000, 2000, all), results/fig5_task_means.csv, results/fig5_ucell_zero.csv, '
        'results/fig5_redo_beyondcell_L500.csv; per-cell scores: results/fig5_redo_percell_L500.csv.gz; t-SNE: data/Hammond2019/raw/cell_meta_hammond.csv.'],
        {'Fig5a': a, 'Fig5b': [cells, b_auc], 'Fig5c': c, 'Fig5d': d, 'Fig5e': me, 'Fig5f': el})


def figS4():
    """Supplementary Fig. S4 as drawn by figS4_capability_map.py: a, the table of all tasks; b, the AUC strips of the six type A tasks
    not drawn in Fig. 5c. (The comparisons of type A at each method's best length and of type B are in Supplementary Table 3.)"""
    A = pd.read_csv(F5 / 'results/fig5_redo_auc_L500.csv')
    X5 = pd.read_csv(F5 / 'extra_refs' / 'fig5_extra_auc_L500.csv')

    ## b: the AUC strips of the six tasks (the lower edge of a panel goes below the range of 0.2 only to show a whole CI)
    GB = [t for t in S4B if t in TA5]                               # D1, D3: the Glucocorticoid tasks not drawn in Fig. 5c
    top = auc5(pd.concat([A[A.id.isin(GB)], X5], ignore_index=True), S4B, M5C, lowered=True)

    ## a: every task of the comparison (Supplementary Table 3, sheet Tasks; Maturity, the ceiling task, is only in the table)
    tasks = pd.read_excel(F5 / 'SuppTable3.xlsx', sheet_name='Tasks')
    tasks = tasks[~tasks.task.str.startswith('Maturity')][['task', 'reference', 'type', 'how the positive cells were defined',
                                                            'positive cells', 'negative cells', 'in the main figure']].reset_index(drop=True)
    write(OUT / 'Source_Data_Supplementary_FigS4.xlsx', [
        'Source Data for Supplementary Fig. S4 (Fig. 5): a, the table of all tasks; b, the type A tasks not drawn in Fig. 5c. The comparisons of each '
        'method with RFind-sc in the eight type A tasks (at 500 genes per direction and at each method\'s best reference length) and in the '
        'type B tasks are in Supplementary Table 3 (sheets AUC, L = 500, Comparisons and Best length).',
        'FigS4a_tasks: every task of the comparison (from the sheet Tasks of Supplementary Table 3; Maturity, the ceiling task, is only '
        'in that table).',
        'FigS4b_typeA_tasks: the AUC of every method of Fig. 5c in the six type A tasks not drawn in Fig. 5c, top 500 genes per '
        'direction (Glucocorticoid in the oligodendrocytes and the endothelial cells: '
        f'{int(A[A.id == "D1"].ref_up.iloc[0])} up and {int(A[A.id == "D1"].ref_down.iloc[0])} down genes remained among the analysed genes; '
        'PS2APP, the cells of the DAM task with the PS2APP reference; Forced swim, the cells of the three Glucocorticoid tasks with the '
        f'Forced swim reference: {int(X5[X5.id == "AS1"].ref_up.iloc[0])} up and {int(X5[X5.id == "AS1"].ref_down.iloc[0])} down genes '
        'remained among the analysed genes); symbols and colours as in Fig. 5 (all filled; no comparison marked). '
        f'Comparisons: {counts5(top)}. The lower edge of a panel goes below the AUC range of 0.2 only so that no 95% CI is cut.',
        gc_zero5(GB),
        'Source tables: figs/Fig5_benchmark/output/SuppTable3.xlsx (Tasks), results/fig5_redo_auc_L500.csv, '
        'extra_refs/fig5_extra_auc_L500.csv, results/fig5_ucell_zero.csv and results/fig5_redo_beyondcell_L500.csv.'],
        {'FigS4a_tasks': tasks, 'FigS4b_typeA_tasks': top})


if __name__ == '__main__':
    OUT.mkdir(exist_ok=True)
    fig2(); fig3(); fig4(); fig5(); figS4()
