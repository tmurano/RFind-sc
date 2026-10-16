# -*- coding: utf-8 -*-
"""Fig. 4 check quoted in Methods ("Analysis of Perturb-seq", References and stress response): how much the knock-out rankings
change with a lower nuclei cut-off of 5 instead of 10.

Spearman correlation of the RFind scores of the knock-outs present in both rankings.
Input : analysis/perturbseq/results/axis_<axis>_vs_<context>.csv (the rankings of Fig. 4), and the same files in
        analysis/perturbseq/results/base_check/all_minc5/ (03_axis_vs_ko.py with PERTURB_MINC=5 and
        PERTURB_RES=analysis/perturbseq/results/base_check/all_minc5)
Output: figs/Fig4_Perturbseq/output/fig4_ranking_checks.txt (also on stdout)
Run   : python figs/Fig4_Perturbseq/analysis/fig4_ranking_checks.py   (from the project root; pandas, scipy)
"""
import os
import pandas as pd
from scipy.stats import spearmanr

R = 'analysis/perturbseq/results'
FILES = {'DG': 'axis_DG_vs_g017_DG.csv', 'CA1': 'axis_CA1_vs_CA1.csv', 'CA3': 'axis_CA3_vs_g017_CA3.csv', 'GABA': 'axis_GABA_vs_GABA3.csv'}
out = []


def scores(path):
    d = pd.read_csv(path)
    return d[d['kind'] == 'KO'].set_index('target')['rfind']


def compare(label, a, b):
    common = a.index.intersection(b.index)
    rho = spearmanr(a[common], b[common]).statistic
    out.append(f'{label}: {len(common):,} knock-outs in both; Spearman {rho:.4f}')


out.append('Fig. 4 ranking check (figs/Fig4_Perturbseq/analysis/fig4_ranking_checks.py)')
out.append('Lower limit of 5 nuclei per knock-out instead of 10 (base_check/all_minc5)')
for ct, f in FILES.items():
    compare(f'   {ct}', scores(os.path.join(R, f)), scores(os.path.join(R, 'base_check', 'all_minc5', f)))
txt = '\n'.join(out) + '\n'
print(txt, end='')
open(os.path.join('figs', 'Fig4_Perturbseq', 'output', 'fig4_ranking_checks.txt'), 'w').write(txt)
