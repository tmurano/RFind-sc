"""RFind-sc demo with the Python package rfindsc: scoring your own cells against your own reference.
The same 1,000 microglia and DAM reference as demo/run_demo.py, scored with the package's public function score_panel_pairs().
The package takes the expected counts from the cells it is given (here these 1,000 cells; the paper took them from all 75,572 cells
of the atlas) and keeps the genes detected in at least 1% of these cells, so the values differ slightly from the paper's DAM-index;
the script reports how closely they agree.
Needs rfindsc 1.0.1 (README, "Companion packages") with numpy, scipy and pandas. Run from the project root:
  python demo/run_demo_package.py
Output: demo/output_package_py.csv (cell_id, condition, DAM_index)
"""
import os
import time

t_start = time.perf_counter()
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.io import mmread  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402
import rfindsc  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, 'demo', 'data')

# ---- your data: raw UMI counts, genes x cells, with gene symbols and cell names ----
X = mmread(os.path.join(D, 'counts.mtx.gz')).tocsc()
genes = pd.read_csv(os.path.join(D, 'genes.tsv.gz'), sep='\t', dtype={'gene': str, 'baseline_share': float})
cells = pd.read_csv(os.path.join(D, 'labels.tsv'), sep='\t', dtype={'cell_id': str, 'condition': str, 'total_umi': float})
# The packages match gene names without regard to case and warn about rows whose names then coincide; this matrix has two
# (PISD and Pisd), which are summed into one row here
key = genes.gene.str.lower()
if key.duplicated().any():
    codes, _ = pd.factorize(key)
    X = (csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes))))) @ X).tocsc()
    genes = genes.loc[~key.duplicated()].reset_index(drop=True)

# ---- your reference: up- and down-regulated genes, most changed first ----
ref = pd.read_csv(os.path.join(D, 'dam_reference.tsv'), sep='\t', dtype={'direction': str, 'rank': int, 'gene': str})
panels = {'DAM': {'up': ref.loc[ref.direction == 'up', 'gene'].tolist(),
                  'down': ref.loc[ref.direction == 'down', 'gene'].tolist()}}

# ---- one call; libsize = each cell's total UMI over all genes, because this matrix holds only a subset of the genes ----
out = rfindsc.score_panel_pairs(X, panels, config=rfindsc.rfindsc_config(n_cores=1), genes=genes.gene.tolist(),
                                cells=cells.cell_id.tolist(), libsize=cells.total_umi.to_numpy())
scores = out['scores']
score = (scores.loc[cells.cell_id, 'DAM'] if isinstance(scores, pd.DataFrame) else np.asarray(scores)[:, 0])
score = np.asarray(score, dtype=float)
usize = out['universe_size']
usize = usize['DAM'] if isinstance(usize, (dict, pd.Series)) else np.asarray(usize).ravel()[0]
print(f'rfindsc {rfindsc.__version__}: {len(score)} cells scored against the DAM reference (scoring universe {int(usize)} genes)')


def auc(s, pos):
    r = rankdata(s)
    n1, n0 = pos.sum(), (~pos).sum()
    return (r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


lpc = (cells.condition == 'LPC').to_numpy()
pd.DataFrame({'cell_id': cells.cell_id, 'condition': cells.condition,
              'DAM_index': [f'{v:.17g}' for v in score]}).to_csv(os.path.join(ROOT, 'demo', 'output_package_py.csv'), index=False)
t_run = time.perf_counter() - t_start

print('\nfirst cells:')
print(pd.DataFrame({'cell_id': cells.cell_id, 'condition': cells.condition, 'DAM_index': np.round(score, 4)}).head(5).to_string(index=False))
print(f'\nmedian DAM-index: LPC {np.median(score[lpc]):.3f}, saline {np.median(score[~lpc]):.3f}')
print(f'AUC (LPC vs saline): {auc(score, lpc):.4f}')
ex = pd.read_csv(os.path.join(ROOT, 'demo', 'expected_output.csv'), dtype={'cell_id': str, 'condition': str, 'DAM_index': float})
assert ex.cell_id.tolist() == cells.cell_id.tolist()
print(f"agreement with the paper's DAM-index of these cells (baseline from all cells of the atlas): "
      f"Spearman {spearmanr(score, ex.DAM_index).statistic:.4f}, AUC there {auc(ex.DAM_index.to_numpy(), lpc):.4f}")
print(f'written demo/output_package_py.csv; run time {t_run:.1f} s (Python {".".join(map(str, __import__("sys").version_info[:3]))}, '
      f'numpy {np.__version__})')
