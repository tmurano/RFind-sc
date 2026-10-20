"""Source Data for Supplementary Figs. S1-S3 and S5: one Excel file per figure, one sheet (or a few) per panel, with the numbers
drawn in each panel and the statistics quoted in the legend or in Supplementary Notes 1-3 and 5. Same style as make_source_data.py
(README sheet first, then the panels; one header row per sheet; notes, when needed, as rows under the table).

Reads only the figure scripts' inputs and outputs (figs/Supp_S*/ and, for S5, the per-cell scores of analysis/kang); the constants
of the figures (toy lists of S1, complexes of S3d) are read from the figure scripts themselves (ast; the scripts are not run), so
the workbooks follow any change made there.
Computes nothing new except what the figure scripts compute on the fly from these inputs (the Running Fisher grid of S1, the
summaries of the pair distributions of S3c, the ranks of S3e, the per-donor means of S5), with the same code.
Not included: schematics without data (S2a, S2d, S3a), the text table of S2b (definitions; no numbers) and S5a, b (the authors'
t-SNE and the table of cell types and markers; as Fig. 3a, b).

Run from the project root:  python figs/make_source_data_supp.py      (needs pandas, numpy, scipy, openpyxl)
Output: submission/Source_Data_Supplementary_FigS1.xlsx, Source_Data_Supplementary_FigS2.xlsx,
        Source_Data_Supplementary_FigS3.xlsx, Source_Data_Supplementary_FigS5.xlsx (the names used for the upload; submission/
        is not tracked)
"""
import ast
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from scipy.stats import hypergeom

OUT = Path('submission')                                     # the files uploaded to the journal (not tracked)
S1, S2, S3 = Path('figs/Supp_S1_RunningFisher'), Path('figs/Supp_S2_Delta'), Path('figs/Supp_S3_Perturbseq')
sys.path.insert(0, str(S3 / 'analysis'))
from figS3_shi_pairs import SHI                      # pairs named as similar by Shi et al. (the list used by the figure)


def constants(script, names):
    """Literal assignments of a figure script (read with ast; the script is not run)."""
    tree = ast.parse(Path(script).read_text(encoding='utf-8'))
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Tuple) and not isinstance(node.value, ast.Tuple):
                    continue                          # unpacking of a computed value: not a constant
                elts = t.elts if isinstance(t, ast.Tuple) else [t]
                vals = node.value.elts if isinstance(t, ast.Tuple) else [node.value]
                for e, v in zip(elts, vals):
                    if isinstance(e, ast.Name) and e.id in names:
                        out[e.id] = ast.literal_eval(v)
    missing = set(names) - set(out)
    assert not missing, f'{script}: {missing} not found'
    return out


def write(path, readme, sheets):
    """README sheet (one line per row) first, then the panels; a sheet is a DataFrame or (DataFrame, [note rows under it])."""
    with pd.ExcelWriter(path, engine='openpyxl') as xw:
        pd.DataFrame({'README': readme}).to_excel(xw, sheet_name='README', index=False)
        for name, df in sheets.items():
            notes = []
            if isinstance(df, tuple):
                df, notes = df
            df.to_excel(xw, sheet_name=name, index=False)
            ws = xw.sheets[name]
            for k, s in enumerate(notes):
                ws.cell(row=len(df) + 3 + k, column=1, value=s)
    print('wrote', path)


# ------------------------------------------------------------------ Supplementary Fig. S1 (toy example of the Running Fisher test)
def figS1():
    c = constants(S1 / 'analysis/figS1_running_fisher.py', ['LISTS', 'EX_G', 'EX_K', 'P_UNIVERSE', 'TABLE'])
    L, U, EG, EK = c['LISTS'], c['P_UNIVERSE'], c['EX_G'], c['EX_K']

    def fisher(l1, l2, G, K):                         # as fisher_p() of the figure script
        M = len(set(l1[:G]) & set(l2[:K]))
        return M, hypergeom.sf(M - 1, U, K, G)

    def grid(l1, l2):
        return np.array([[-np.log10(fisher(l1, l2, g, k)[1]) for k in range(1, len(l2) + 1)] for g in range(1, len(l1) + 1)])

    # a: the four lists (each ranks all 26 genes: 13 up, 13 down) and the illustrative spreadsheet rows drawn in a
    a = pd.DataFrame({'Rank (1 = strongest change)': range(1, 14), **{k: L[k] for k in ['b1_up', 'b1_down', 'b2_up', 'b2_down']}})
    drawn = []
    for lst, rows in c['TABLE'].items():
        for i, r in enumerate(rows):
            if r is not None:
                drawn.append(f'{lst}: {r[0]} {r[1].replace(chr(8722), "-")}')
    a_notes = [f'Toy example; background (U) = the {U} genes GeneA-GeneZ. b1_up / b2_up start from the most increased gene, '
               'b1_down / b2_down from the most decreased gene.',
               'The log2FC values drawn in the tables of panel a are illustrative (not used in b-d): ' + '; '.join(drawn) + '.']

    # b: one pair (b1_up x b2_up) at one pair of cut-offs
    M, P = fisher(L['b1_up'], L['b2_up'], EG, EK)
    shared = [g for g in L['b1_up'][:EG] if g in L['b2_up'][:EK]]
    b = pd.DataFrame([
        ('Pair', 'b1_up × b2_up'), ('G (top genes of b1_up)', EG), ('K (top genes of b2_up)', EK),
        ('Genes in the top G of b1_up', ', '.join(L['b1_up'][:EG])), ('Genes in the top K of b2_up', ', '.join(L['b2_up'][:EK])),
        ('Shared genes', ', '.join(shared)), ('U (all genes)', U),
        ('2 × 2 table: in top K and in top G (M)', M), ('2 × 2 table: in top K, not in top G (K − M)', EK - M),
        ('2 × 2 table: not in top K, in top G (G − M)', EG - M), ('2 × 2 table: in neither (U − G − K + M)', U - EG - EK + M),
        ("P, one-sided Fisher's exact test (hypergeometric upper tail)", P), ('−log10 P', -np.log10(P))],
        columns=['Quantity', 'Value'])

    # c: the 13 x 13 grid for up x up
    H = grid(L['b1_up'], L['b2_up'])
    bg, bk = [int(i) + 1 for i in np.unravel_index(np.argmax(H), H.shape)]
    cg = pd.DataFrame(np.round(H, 4), columns=[f'K = {k}' for k in range(1, 14)])
    cg.insert(0, 'G (top genes of b1_up)', range(1, 14))
    ties = int(np.isclose(H, H.max()).sum())
    c_notes = ['Cell = −log10 P of the one-sided Fisher test of the top G genes of b1_up against the top K genes of b2_up (U = 26); '
               'colour of panel c.',
               f'Dashed square = the test of panel b (G = {EG}, K = {EK}): −log10 P = {H[EG - 1, EK - 1]:.4f}.',
               f'Red square = the smallest P (G = {bg}, K = {bk}): P = {10 ** -H.max():.4g}, −log10 P = {H.max():.4f} '
               f'(= Score of up × up in d; {ties} cell(s) reach this value, the first in row order is framed).']

    # d: the four pairs and S_total
    pairs = [('b1_up', 'b2_up', 'up × up', 'same', +1), ('b1_down', 'b2_down', 'down × down', 'same', +1),
             ('b1_up', 'b2_down', 'up × down', 'opposite', -1), ('b1_down', 'b2_up', 'down × up', 'opposite', -1)]
    rows = []
    for l1, l2, lab, dr, sg in pairs:
        h = grid(L[l1], L[l2]); g, k = [int(i) + 1 for i in np.unravel_index(np.argmax(h), h.shape)]
        rows.append({'Pair': lab, 'Lists': f'{l1} × {l2}', 'Direction': dr, 'Sign in S_total': '+' if sg > 0 else '−',
                     'G at the smallest P': g, 'K at the smallest P': k, 'M at the smallest P': fisher(L[l1], L[l2], g, k)[0],
                     'Smallest P': 10 ** -h.max(), 'Score = −log10 (smallest P) (bar)': round(h.max(), 4), '_s': sg * h.max()})
    d = pd.DataFrame(rows)
    st = d.pop('_s').sum()
    d = pd.concat([d, pd.DataFrame([{'Pair': 'S_total', 'Lists': 'up×up + down×down − up×down − down×up',
                                     'Score = −log10 (smallest P) (bar)': round(st, 4)}])], ignore_index=True)
    d_notes = [f'S_total of the example = {st:.4f} (positive: the two lists change in the same direction; ▼ in panel d). '
               'The bars of d have no y scale in the figure; the values are those above.']

    write(OUT / 'Source_Data_Supplementary_FigS1.xlsx', [
        'Source Data for Supplementary Fig. S1 (the Running Fisher test; a toy example).',
        'Toy gene lists: background U = 26 genes (GeneA-GeneZ); each list ranks all 26 genes, 13 increased (up) and 13 decreased (down).',
        'Every number is computed from the toy lists with the one-sided Fisher exact test (hypergeometric upper tail, U = 26), '
        'exactly as in the figure script.',
        'FigS1a_lists: the four ranked lists. FigS1b_one_cutoff: the 2 × 2 table of panel b. FigS1c_grid: −log10 P at every (G, K) for '
        'up × up. FigS1d_scores: the four scores and S_total.',
        'Source: the lists and cut-offs are read from figs/Supp_S1_RunningFisher/analysis/figS1_running_fisher.py.'],
        {'FigS1a_lists': (a, a_notes), 'FigS1b_one_cutoff': b, 'FigS1c_grid': (cg, c_notes), 'FigS1d_scores': (d, d_notes)})
    return {'b_P': P, 'b_M': M, 'best': (bg, bk, H.max()), 'scores': d['Score = −log10 (smallest P) (bar)'].tolist()}


# ------------------------------------------------------------------ Supplementary Fig. S2 (choice of the per-cell value)
LAB2 = {'simpleFC': 'simpleFC', 'log2diff': 'log2diff', 'zscore': 'z-score', 'P1': 'P1', 'P2': 'P2', 'P3_NB': 'P3 (NB)'}
KIND = {'UP': 'UMI ≥ 1, above the average (red)', 'DOWN-1': 'UMI ≥ 1, below the average (light blue)',
        'DOWN-2': 'UMI = 0, expected count < ln 2 (grey)', 'DOWN-3': 'UMI = 0, expected count ≥ ln 2 (blue)'}


def figS2():
    D = S2 / 'output/data'
    c = constants(S2 / 'analysis/figS2_panels.py', ['PLO', 'PHI', 'TORDER'])
    one = pd.read_csv(D / 'figS2b_one_cell.csv')
    genes = pd.DataFrame({'Gene': one.gene, 'UMI count (x)': one.x, 'Expected count (e)': one.e.round(5),
                          **{LAB2[m]: one[m].round(5) for m in c['TORDER']},
                          'Group (colour)': one.kind.map(KIND)})
    # what is drawn: numbers of genes per group, zeros of simpleFC, genes outside the axis of z-score and P1-P3
    rows = []
    for k in ['UP', 'DOWN-1', 'DOWN-2', 'DOWN-3']:
        s = one[one.kind == k]
        r = {'Group (colour)': KIND[k], 'Genes (n)': len(s), 'simpleFC = 0 (bar left of the log axis)': int((s.simpleFC == 0).sum())}
        for m in ['zscore', 'P1', 'P2', 'P3_NB']:
            r[f'{LAB2[m]} < {c["PLO"]:g} (left end bar)'] = int((s[m] < c['PLO']).sum())
            r[f'{LAB2[m]} > {c["PHI"]:g} (right end bar)'] = int((s[m] > c['PHI']).sum())
        rows.append(r)
    summ = pd.DataFrame(rows)
    tot = summ.drop(columns='Group (colour)').sum().to_dict(); tot['Group (colour)'] = 'All genes'
    summ = pd.concat([summ, pd.DataFrame([tot])], ignore_index=True)
    stat = pd.DataFrame([{'Candidate': LAB2[m], 'Min': one[m].min(), 'Median': one[m].median(), 'Max': one[m].max(),
                          'Median, UMI = 0 genes': one.loc[one.x == 0, m].median(),
                          'Median, UMI = 1 genes': one.loc[one.x == 1, m].median()} for m in c['TORDER']]).round(4)

    bench = pd.read_csv(D / 'figS2_delta_bench.csv')
    perc = pd.read_csv(D / 'figS2_delta_bench_percluster.csv', dtype={'cluster': str})
    order = bench.sort_values('composite', ascending=False).model.tolist()
    mean = perc.groupby('model')[['marker_top', 'negmarker_bot', 'composite']].mean()
    e1 = pd.DataFrame({'Candidate (order of the figure)': [LAB2[m] for m in order],
                       'Share of markers in UP, mean of 13 clusters (red bar)': bench.set_index('model').loc[order, 'marker_top'].values,
                       'Share of negative markers in DOWN, mean of 13 clusters (blue bar)': bench.set_index('model').loc[order, 'negmarker_bot'].values,
                       'Sum of the two shares, mean of 13 clusters': bench.set_index('model').loc[order, 'composite'].values,
                       'Sum, unrounded': mean.loc[order, 'composite'].round(6).values,
                       'Clusters (n)': perc.groupby('model').cluster.nunique().loc[order].values,
                       'Cells (n, all clusters)': perc.groupby('model').n.sum().loc[order].values})
    perc['rank'] = perc.groupby('cluster').composite.rank(ascending=False, method='min').astype(int)
    CL = ['1', '2a', '2b', '2c', '3', '4', '5', '6', '7a', '7b', '7c', '8', '9']
    perc['co'] = perc.cluster.map(CL.index); perc['mo'] = perc.model.map(order.index)
    perc = perc.sort_values(['mo', 'co'])
    e2 = pd.DataFrame({'Candidate': perc.model.map(LAB2), 'Cluster': perc.cluster, 'Cells (n)': perc.n,
                       'Share of markers in UP': perc.marker_top.round(4), 'Share of negative markers in DOWN': perc.negmarker_bot.round(4),
                       'Sum of the two shares (point)': perc.composite.round(4),
                       'Rank of the candidate in the cluster (1 = highest sum)': perc['rank']})
    q = perc.groupby('model').composite
    e3 = pd.DataFrame({'Candidate': [LAB2[m] for m in order], 'Median of the 13 clusters (bar)': q.median().loc[order].round(4).values,
                       '25th percentile (shade, left)': q.quantile(0.25).loc[order].round(4).values,
                       '75th percentile (shade, right)': q.quantile(0.75).loc[order].round(4).values,
                       'Clusters in which the candidate ranks first (n)': perc.groupby('model')['rank'].apply(lambda r: int((r == 1).sum())).loc[order].values})

    write(OUT / 'Source_Data_Supplementary_FigS2.xlsx', [
        'Source Data for Supplementary Fig. S2 (choice of the per-cell value; microglia of Hammond et al., 2019).',
        'Panels a and d are schematics and panel b a table of definitions (no data); they have no sheet.',
        'FigS2c_genes: every gene analysed (8,960) in the example cell (the first cell, column 1, of '
        'data/Hammond2019/output/counts_filtered.rds; its barcode can be read only with R): UMI count, expected count and the six '
        'candidate values; group = colour of panel c (upper row UMI ≥ 1, lower row UMI = 0). The histograms of c are binned from these values.',
        'FigS2c_summary: genes per group and the genes drawn as the separate bars (zeros of simpleFC; z-score and P1-P3 outside −4 to 8). '
        'FigS2c_candidates: range and medians of each candidate in this cell.',
        'FigS2e_mean: left plot (mean over the 13 clusters, monocytes/macrophages excluded; up to 1,000 cells per cluster). '
        'FigS2e_clusters: right plot (one point per cluster and candidate). FigS2e_summary: median (bar) and interquartile range (shade) '
        'of the right plot, and the number of clusters in which each candidate ranks first.',
        'Lists of 500 genes per direction (cell UP/DOWN, markers and negative markers of the cluster).',
        'Source tables: figs/Supp_S2_Delta/output/data/figS2b_one_cell.csv, figS2_delta_bench.csv and '
        'figS2_delta_bench_percluster.csv (figS2_delta_benchmark.R).'],
        {'FigS2c_genes': genes, 'FigS2c_summary': summ, 'FigS2c_candidates': stat, 'FigS2e_mean': e1, 'FigS2e_clusters': e2,
         'FigS2e_summary': e3})
    return {'n_genes': len(one), 'summ': summ, 'e1': e1, 'e3': e3}


# ------------------------------------------------------------------ Supplementary Fig. S3 (Perturb-seq knock-outs)
def load_pairs(ctx):
    d = np.load(S3 / f'output/data/rfind_allpairs_{ctx}.npz', allow_pickle=True)
    return d['S'].astype(np.float64), d['targets'].astype(str), d['kind'].astype(str), d['nuclei']


def shi_of(a, b):
    return [f for f, g in SHI if a in g and b in g]


def shi_ranks(ctx):
    """Every pair named by Shi et al. (both knock-outs in the data): score, rank among all knock-out pairs, top 1% or not
    (as shi_top1() of figS3e_celltypes.py)."""
    S, t, k, n = load_pairs(ctx)
    ko = np.where(k == 'KO')[0]; tk = t[ko]
    iu = np.triu_indices(len(ko), 1); Sk = S[np.ix_(ko, ko)][iu]
    order = np.argsort(-Sk, kind='stable'); rpos = np.empty_like(order); rpos[order] = np.arange(1, len(Sk) + 1)
    ix = {g: i for i, g in enumerate(tk)}; N = len(tk)
    seen, rows = set(), []
    for fig, genes in SHI:
        for a, b in combinations([g for g in genes if g in ix], 2):
            i, j = sorted((ix[a], ix[b]))
            if (i, j) in seen:
                continue
            seen.add((i, j)); q = i * N - i * (i + 1) // 2 + (j - i - 1)
            rows.append({'Knock-out 1': a, 'Knock-out 2': b, 'Statement of Shi et al. (first)': fig, 'RFind score': round(Sk[q], 2),
                         'Rank among all knock-out pairs': int(rpos[q]), 'Top % of all knock-out pairs': round(100 * rpos[q] / len(Sk), 4),
                         'In the top 1%': bool(100 * rpos[q] / len(Sk) <= 1)})
    return pd.DataFrame(rows), len(Sk), int(np.median(n[ko])), len(ko)


def figS3():
    pan = constants(S3 / 'analysis/figS3_rfind_panels.py', ['VLB', 'MODULES', 'COMPLEX'])
    S, tgt, kind, nuc = load_pairs('g151')
    ko, sf = kind == 'KO', kind == 'safe'
    tri = lambda m: S[np.ix_(m, m)][np.triu_indices(int(m.sum()), 1)]
    cls = {'Knock-out × knock-out': tri(ko), 'Knock-out × safe-targeting': S[np.ix_(ko, sf)].ravel(),
           'Safe-targeting × safe-targeting': tri(sf)}
    kk, ss = cls['Knock-out × knock-out'], cls['Safe-targeting × safe-targeting']
    smax = ss.max(); n_above, n_exp = int((kk > smax).sum()), len(kk) / (len(ss) + 1)

    # b: the 464 knock-outs with a partner above the highest safe x safe pair, in the clustered order of the figure
    Skk = S[np.ix_(ko, ko)].copy(); np.fill_diagonal(Skk, -np.inf)
    strong = np.flatnonzero((Skk > smax).any(1)); tks = tgt[ko][strong]
    Xs = S[np.ix_(np.flatnonzero(ko)[strong], np.flatnonzero(ko)[strong])].copy(); np.fill_diagonal(Xs, 0.0)
    Ds = np.clip(1 - np.corrcoef(Xs), 0, None); np.fill_diagonal(Ds, 0.0)
    ordr = hierarchy.leaves_list(hierarchy.linkage(squareform(Ds, checks=False), method='average'))
    Mb = Xs[np.ix_(ordr, ordr)].copy(); np.fill_diagonal(Mb, np.nan)
    names = tks[ordr]
    bmat = pd.DataFrame(np.round(Mb, 2), columns=names); bmat.insert(0, 'Knock-out (row = column order)', names)
    pos = {g: i for i, g in enumerate(names)}
    mods = '; '.join(f'{nm}: ' + ', '.join(f'{g} (row {pos[g] + 1})' for g in gs if g in pos)
                     for nm, gs in pan['MODULES'] + [('proteasome (ticks on the left)', ['Psmb4', 'Pomp', 'Psmc5', 'Psmc1'])])
    b_notes = [f'{len(strong)} of {int(ko.sum()):,} knock-outs (thalamic glutamatergic neurons) with an RFind score above the highest '
               f'safe × safe pair ({smax:.2f}) with at least one other knock-out; rows and columns in the order of the heat map '
               '(average linkage on 1 − Pearson correlation of the rows, diagonal set to 0 for the clustering; blank diagonal).',
               f'Colour scale ±{pan["VLB"]:g}: values beyond are drawn in the end colours (the values here are not clipped).',
               'Complexes labelled in b (rows counted from 1 at the top): ' + mods + '.']

    # c: summaries by type of pair, the line plots (binned densities) and the safe x safe ticks
    imax = np.argmax(ss); sfi = np.flatnonzero(sf); iu = np.triu_indices(len(sfi), 1)
    top_ss = (tgt[sfi[iu[0][imax]]], tgt[sfi[iu[1][imax]]])
    c1 = pd.DataFrame([{'Type of pair': k_, 'Pairs (n)': len(v), 'Mean': round(v.mean(), 3), 'Median': round(np.median(v), 3),
                        'SD': round(v.std(), 3), 'Share > 0': round((v > 0).mean(), 4), 'Min': round(v.min(), 2),
                        'Max': round(v.max(), 2), f'Pairs above {smax:.2f} (n)': int((v > smax).sum()),
                        'Pairs >= 30 (n)': int((v >= 30).sum())} for k_, v in cls.items()])
    c_notes = [f'Highest safe × safe pair (dashed line) = {smax:.4f} ({top_ss[0]} × {top_ss[1]}).',
               f'Knock-out × knock-out pairs above it: {n_above:,} of {len(kk):,} ({100 * n_above / len(kk):.2f}%); expected if '
               f'distributed like the control pairs = n KO pairs / (n control pairs + 1) = {len(kk):,} / {len(ss) + 1:,} = {n_exp:.1f}; '
               f'ratio {n_above / n_exp:.1f}.',
               f'Knock-outs {int(ko.sum()):,}, safe-targeting groups {int(sf.sum())}; median nuclei per knock-out {int(np.median(nuc[ko]))}, '
               f'per safe-targeting group {int(np.median(nuc[sf]))}.',
               'The full data (about 2.1 million pairs) are not listed; the plotted densities are in FigS3c_density_main and '
               'FigS3c_density_right_end.']

    def dens(v, edges):
        h, _ = np.histogram(v, bins=edges); return h, h / (len(v) * np.diff(edges))
    blocks = {}
    for nm, edges, types in [('main', np.arange(-40, 47, 1.0), list(cls)),
                             ('right', np.arange(30, 152.5, 2.5), ['Knock-out × knock-out', 'Knock-out × safe-targeting'])]:
        df = pd.DataFrame({'Bin start': edges[:-1], 'Bin end': edges[1:], 'Bin centre (x of the line)': (edges[:-1] + edges[1:]) / 2})
        for k_ in types:
            h, dn = dens(cls[k_], edges)
            df[f'{k_}: pairs (n)'] = h; df[f'{k_}: density (y)'] = dn.round(8)
        out = {k_: (int((cls[k_] < edges[0]).sum()), int((cls[k_] >= edges[-1]).sum())) for k_ in types}
        blocks[nm] = (df, [f'Bins [start, end) as numpy.histogram (last bin closed); density = pairs in the bin / (all pairs of the type × bin width). '
                           'Pairs outside the bins (not drawn): ' + '; '.join(f'{k_}: {a:,} below, {b:,} above' for k_, (a, b) in out.items()) + '.'])
    vv = ss[ss >= 30]; o = np.argsort(-vv)
    ticks = pd.DataFrame({'Safe-targeting group 1': [tgt[sfi[iu[0][i]]] for i in np.flatnonzero(ss >= 30)[o]],
                          'Safe-targeting group 2': [tgt[sfi[iu[1][i]]] for i in np.flatnonzero(ss >= 30)[o]],
                          'RFind score (orange tick)': np.round(vv[o], 4)})

    # d: the 20 most similar knock-out pairs
    iuk = np.triu_indices(int(ko.sum()), 1); Sk = S[np.ix_(ko, ko)][iuk]; tk = tgt[ko]
    top = np.argsort(-Sk, kind='stable')[:20]
    d = pd.DataFrame([{'Rank': r, 'Knock-out 1': tk[iuk[0][k_]], 'Knock-out 2': tk[iuk[1][k_]], 'RFind score (bar)': round(Sk[k_], 2),
                       'Reported by Shi et al. (red)': 'yes' if shi_of(tk[iuk[0][k_]], tk[iuk[1][k_]]) else 'no',
                       'Statement of Shi et al.': ', '.join(shi_of(tk[iuk[0][k_]], tk[iuk[1][k_]])),
                       'Complex or pathway (right of the bar)': pan['COMPLEX'].get(tuple(sorted((tk[iuk[0][k_]], tk[iuk[1][k_]]))), '')}
                      for r, k_ in enumerate(top, 1)])
    n_shi = int((d['Reported by Shi et al. (red)'] == 'yes').sum())
    d_notes = [f'{n_shi} of the 20 pairs are pairs that Shi et al. describe as similar in the main text; '
               f'among {len(Sk):,} knock-out pairs of the thalamic glutamatergic neurons.']

    # e: shares in the top 1% by cell type (the three drawn; no single hippocampal cell types, which the paper does not quote),
    # and every pair behind them
    CT = [('g151', 'Thalamus (TH Prkcd Grin2c Glut)', 'bar'), ('HIP', 'Hippocampus (Delta of DG, CA1, CA3 summed)', 'bar'),
          ('GABA', 'GABAergic neurons', 'bar')]
    erows, epairs = [], []
    for key, name, shown in CT:
        pr, npairs, _, nko = shi_ranks(key)
        a, n = int(pr['In the top 1%'].sum()), len(pr)
        erows.append({'Cell type': name, 'Drawn': shown, 'Pairs of Shi et al. in the data (n; under the bar)': n,
                      'Of them in the top 1% (n)': a, 'Share in the top 1% (%, bar)': round(100 * a / n, 1),
                      'Expected by chance (%)': 1, 'Knock-outs (n)': nko, 'Knock-out pairs (n)': npairs})
        epairs.append(pr.assign(**{'Cell type': name}))
    e = pd.DataFrame(erows)
    ep = pd.concat(epairs, ignore_index=True)
    ep = ep[['Cell type'] + [c_ for c_ in ep.columns if c_ != 'Cell type']]

    write(OUT / 'Source_Data_Supplementary_FigS3.xlsx', [
        'Source Data for Supplementary Fig. S3 (in vivo Perturb-seq of Shi et al., doi:10.64898/2026.03.16.711480).',
        'RFind score = S_total (4-term Running Fisher) between the UP/DOWN lists (500 genes each, by Delta) of two pseudo-bulk columns.',
        'Panel a is a schematic (no sheet). b-d: thalamic glutamatergic neurons (TH Prkcd Grin2c Glut; 1,947 knock-outs, '
        '100 safe-targeting groups, 12,286 genes); e: the thalamus and the cell types of Fig. 4.',
        'FigS3b_matrix: RFind scores among the 464 knock-outs of the heat map, in its order. FigS3c_summary: the three types of pair '
        '(with the dashed line, the counts above it and the number expected). FigS3c_density_main / _right_end: the plotted densities '
        '(binned; the 2.1 million single pairs are not listed). FigS3c_safe_ticks: the safe × safe pairs above +30 (orange ticks).',
        'FigS3d_top20: the 20 most similar knock-out pairs. FigS3e_celltypes: the bars. FigS3e_pairs: every pair of Shi et al. with its '
        'score and rank in each cell type.',
        'Source: figs/Supp_S3_Perturbseq/output/data/rfind_allpairs_{g151,HIP,GABA}.npz (figS3_rfind_allpairs.py); '
        'pairs of Shi et al. from figS3_shi_pairs.py; complexes of b and d read from figS3_rfind_panels.py. '
        'Summaries recomputed with the code of figS3_rfind_panels.py and figS3e_celltypes.py.'],
        {'FigS3b_matrix': (bmat, b_notes), 'FigS3c_summary': (c1, c_notes), 'FigS3c_density_main': blocks['main'],
         'FigS3c_density_right_end': blocks['right'], 'FigS3c_safe_ticks': ticks, 'FigS3d_top20': (d, d_notes),
         'FigS3e_celltypes': e, 'FigS3e_pairs': ep})
    return {'smax': smax, 'n_above': n_above, 'n_exp': n_exp, 'n_strong': len(strong), 'd': d, 'e': e, 'c1': c1}


# ------------------------------------------------------------------ Supplementary Fig. S5 (human PBMC of Kang et al., IFN-beta)
S5 = Path('figs/Supp_S5_Kang/output')
IDX5 = {'IFNb': 'IFNb-index', 'PGE2': 'PGE2-index', 'Interferon': 'Interferon-index', 'Chemokine': 'Chemokine-index'}
CT5 = {'CD14+ Monocytes': 'cM', 'FCGR3A+ Monocytes': 'ncM', 'Dendritic cells': 'DC', 'Megakaryocytes': 'Mkc',
       'CD4 T cells': 'TH', 'B cells': 'B', 'CD8 T cells': 'Tc', 'NK cells': 'NK'}       # abbreviations of the figure (TH = T_H)
COND5 = {'ctrl': 'Unstimulated', 'stim': 'IFN-β'}


def figS5():
    st = pd.read_csv(S5 / 'figS5_violin_stats.csv')
    nc = pd.read_csv(S5 / 'FigS5_celltype_ncells.csv').set_index('cell_type')
    vio = pd.DataFrame({
        'Panel': np.where(st.axis.isin(['IFNb', 'PGE2']), 'S5c', 'S5d'), 'Index': st.axis.map(IDX5), 'Cell type': st.cell_type.map(CT5),
        'Unstimulated cells (n)': st.cell_type.map(nc.n_Unstimulated), 'IFN-β cells (n)': st.cell_type.map(nc.n_IFNb),
        'Unstimulated patients (n)': st.n_donors_Unstimulated, 'IFN-β patients (n)': st.n_donors_IFNb,
        "Cohen's d, IFN-β - unstimulated (per-patient means; number above the violins and point of the summary panels of e, f)": st.d.round(3),
        'P, two-sided Wilcoxon rank-sum on per-patient means (nominal)': st.p.round(4)})

    s = pd.read_csv('data/Kang2018/output/scores_kang.csv')
    pm = (s.melt(id_vars=['cell_type', 'condition', 'donor'], value_vars=list(IDX5), var_name='axis', value_name='v')
            .groupby(['axis', 'cell_type', 'condition', 'donor']).v.agg(['mean', 'size']).reset_index())
    pm = pm[pm['size'] >= 10]
    pm['ax'] = pd.Categorical(pm.axis, list(IDX5)); pm['ct'] = pd.Categorical(pm.cell_type, list(CT5))
    pm = pm.sort_values(['ax', 'ct', 'condition', 'donor'])
    pts = pd.DataFrame({'Index': pm.axis.map(IDX5), 'Cell type': pm.cell_type.map(CT5), 'Condition': pm.condition.map(COND5),
                        'Patient': pm.donor, 'Cells (n)': pm['size'],
                        'Mean index of the patient (point of c, d; with the other index of the pair, point of the small panels of e, f)':
                            pm['mean'].round(4)})

    rg = pd.read_csv(S5 / 'figS5_scatter_axis_ranges.csv'); ou = pd.read_csv(S5 / 'figS5_scatter_outside.csv')
    w = ou.merge(rg, on=['pair', 'ct'], how='left')
    win = pd.DataFrame({'Panel': np.where(w.pair.str.startswith('IFNb'), 'S5e', 'S5f'), 'Pair (x x y)': w.pair,
                        'Cell type': w.ct.str.replace('T[H]', 'TH', regex=False),
                        'x from': w.xlo, 'x to': w.xhi, 'y from': w.ylo, 'y to': w.yhi, 'Cells (n)': w.n_cells,
                        'Cells outside the window (n; not drawn)': w.n_outside, 'Cells outside the window (%)': w.pct_outside})

    write(OUT / 'Source_Data_Supplementary_FigS5.xlsx', [
        'Source Data for Supplementary Fig. S5 (human PBMC of Kang et al., 2018; eight patients, unstimulated vs IFN-β 6 h; 8 cell types).',
        "FigS5c-d_stats: Cohen's d (pooled SD) of IFN-β minus unstimulated over per-patient means (patients with >= 10 cells), as written above "
        'the violins; the summary panels of e and f plot the same d. P = nominal two-sided Wilcoxon rank-sum on the per-patient means '
        '(for reference).',
        'FigS5c-f_points: per-patient means drawn as points in c and d; the small panels of e and f plot, for each patient, the means of the '
        'two indices of the pair.',
        'FigS5e-f_windows: the window of each small panel of e and f (all panels of a pair have the same width) and the cells outside it '
        '(not drawn; the shares in the legend). Row "all": all cells.',
        'Panels a and b (the t-SNE of the authors and the table of markers) have no sheet.',
        'Source tables: figs/Supp_S5_Kang/output/figS5_violin_stats.csv, FigS5_celltype_ncells.csv, figS5_scatter_axis_ranges.csv, '
        'figS5_scatter_outside.csv; per-cell scores: data/Kang2018/output/scores_kang.csv.'],
        {'FigS5c-d_stats': vio, 'FigS5c-f_points': pts, 'FigS5e-f_windows': win})
    return {'vio': vio, 'win': win, 'n_points': len(pts)}


if __name__ == '__main__':
    OUT.mkdir(exist_ok=True)
    r1 = figS1(); print('S1:', r1)
    r2 = figS2(); print('S2:', r2['e1'].to_string(), r2['e3'].to_string(), sep='\n')
    r3 = figS3(); print('S3:', {k: r3[k] for k in ('smax', 'n_above', 'n_exp', 'n_strong')}); print(r3['e'].to_string()); print(r3['c1'].to_string())
    r5 = figS5(); print('S5:', r5['n_points'], 'donor points'); print(r5['vio'].to_string()); print(r5['win'].to_string())
