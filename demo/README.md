# RFind-sc demo: the DAM-index of 1,000 microglia

A small, self-contained run of RFind-sc that reproduces numbers of the paper exactly. It scores 1,000 microglia from the
Hammond et al. (2019) atlas against the disease-associated microglia (DAM) reference of Keren-Shaul et al. (2017) and
checks that the per-cell DAM-index equals the one the paper's pipeline stored for these cells (Figs. 2 and 5).

- **Cells**: the DAM task of Fig. 5 (task H1): microglia from white matter 7 days after injection of lysolecithin (LPC) or
  saline (P100). From the 2,459 LPC and 2,432 saline microglia that passed quality control (non-microglial Mono/Mac cells
  removed), 500 of each were drawn at random (`set.seed(1)`).
- **Genes**: the 5,369 genes of the background for matching of this cell-reference pair (the 8,960 genes analysed in the data
  of Hammond et al. that were also reported by Keren-Shaul et al.; U = 5,368 distinct names, because `PISD` and `Pisd` differ
  only in case and are counted once, as in the paper). The per-cell lists of the paper are restricted to these genes before
  they are cut, so these genes give the same lists as all the genes analysed.
- **Reference**: the 500 most up-regulated and the 500 most down-regulated genes of DAM vs homeostatic microglia (Keren-Shaul
  et al. 2017, Table S3, ranked by fold change), taken before the genes of the data; the 408 up- and 427 down-regulated genes
  among the 5,369 genes are used, as in the paper.

## Files

| file | content |
|---|---|
| `data/counts.mtx.gz` | raw UMI counts, genes x cells, Matrix Market (integer), 926,973 non-zero entries |
| `data/genes.tsv.gz` | gene symbol and its baseline share (the gene's mean share of a cell's UMIs over all 75,572 QC cells; `delta_P2_basis.rds` of the paper, 17 significant digits) |
| `data/barcodes.tsv.gz` | cell barcodes, in the column order of the matrix |
| `data/labels.tsv` | cell_id, condition (`LPC` / `saline`), total UMI of the cell over all genes (used for the expected count) |
| `data/dam_reference.tsv` | the DAM reference: direction (`up` / `down`), rank, gene (lower case) |
| `expected_output.csv` | cell_id, condition and the DAM-index stored by the paper's pipeline (`data/Hammond2019/output/scores_hammond_9axes.csv`) |
| `run_demo.R`, `run_demo.py` | the demo in R and in Python (the kernel of this repository; reproduces the paper's values) |
| `run_demo_package.R`, `run_demo_package.py` | the same cells scored with the R package RFindsc and the Python package rfindsc, as you would score your own data |
| `make_demo_data.R` | how the files in `data/` were made from the full data (not needed to run the demo) |

The data come from Hammond et al., Immunity 2019 (GEO GSE121654) and Keren-Shaul et al., Cell 2017; cite them when you use
these files.

## What the demo computes

For every cell (as `analysis/hammond/05_delta_P2.R` and `07_score_axes.R`):

1. **Delta** = (x − e) / √(e + 1), with x = UMI count and e = baseline share of the gene × total UMI of the cell;
2. **UP** = the 500 genes with the largest positive Delta, **DOWN** = the 500 most negative;
3. **DAM-index** = RF(UP, up) − RF(UP, down) − RF(DOWN, up) + RF(DOWN, down), where RF is the Running Fisher score
   (−log10 of the smallest hypergeometric *P* over every pair of cut-offs of the two lists; background U = 5,368).

Then the AUC of the DAM-index for LPC vs saline cells, and the largest difference from `expected_output.csv`.
The R version uses the kernel `analysis/running_fisher.R`; the Python version uses its copy `analysis/running_fisher.py`
(with `analysis/_hyper.py`; the kernel of the Python package rfindsc). Neither needs the packages RFindsc or rfindsc.

## How to run

From the project root (the folder that contains `analysis/` and `demo/`):

```bash
Rscript demo/run_demo.R          # R >= 4.1 with the Matrix package (recommended package, installed with R)
python demo/run_demo.py          # Python >= 3.9 with numpy, scipy and pandas
```

Outputs: `demo/output_R.csv` and `demo/output_py.csv` (cell_id, condition, DAM_index).

## Expected output

Both versions print (the Python table has the same numbers):

```
demo: 5369 genes (P = 5368) x 1000 cells; reference 408 up / 427 down genes

first cells:
                          cell_id condition DAM_index
 P100_M_LPC_A4_AACCATGTCGCCGTGA-1       LPC    7.0129
 P100_M_LPC_A4_AACGTTGGTTAAGAAC-1       LPC    5.2932
 P100_M_LPC_A4_AACTCAGGTATCTGCA-1       LPC   -6.2090
 P100_M_LPC_A4_AACTCCCAGTTAACGA-1       LPC   12.7692
 P100_M_LPC_A4_AAGCCGCAGCCAGTAG-1       LPC   23.3063

median DAM-index: LPC 10.353, saline -9.025
AUC (LPC vs saline): 0.9045
max |difference| from the paper's stored DAM-index: 5.68e-14 -> OK
```

The largest difference from the stored values was 5.7 × 10⁻¹⁴ (R) and 7.1 × 10⁻¹⁴ (Python); the R and Python scores
differ by at most 2.8 × 10⁻¹⁴. The AUC of these 1,000 cells (0.9045) is close to that of all 4,891 cells of the task in the
paper (0.906, 95% CI 0.898–0.914; Fig. 5 and Supplementary Table 3).


## Scoring your own data with the packages

`run_demo_package.R` and `run_demo_package.py` score the same 1,000 cells with the public function `score_panel_pairs()` of the
companion packages (RFindsc and rfindsc 1.0.1; installation in the main README, "Companion packages"). This is how
you would score your own cells: a gene x cell matrix of raw UMI counts and a reference of up- and down-regulated genes, most
changed first, in one call.

```bash
Rscript demo/run_demo_package.R      # RFindsc installed; set RFINDSC_LIB=<library> if it is in a separate library
python demo/run_demo_package.py      # rfindsc installed (e.g. the virtual environment of the main README)
```

Both scripts first sum the two rows of the matrix whose names differ only in case (PISD and Pisd), because the packages match gene
names without regard to case. The packages take the expected counts from the cells they are given (here these 1,000 cells, half of
them from the injured white matter; the paper took them from all 75,572 cells of the atlas) and keep the genes detected in at least
1% of these cells (5,345 of the 5,368 genes), so the values are shifted from the paper's, but the cells keep almost the same order.
Both versions print:

```
RFindsc 1.0.1: 1000 cells scored against the DAM reference (scoring universe 5345 genes)

first cells:
                          cell_id condition DAM_index
 P100_M_LPC_A4_AACCATGTCGCCGTGA-1       LPC    1.5301
 P100_M_LPC_A4_AACGTTGGTTAAGAAC-1       LPC   -5.0596
 P100_M_LPC_A4_AACTCAGGTATCTGCA-1       LPC  -21.5596
 P100_M_LPC_A4_AACTCCCAGTTAACGA-1       LPC    8.0342
 P100_M_LPC_A4_AAGCCGCAGCCAGTAG-1       LPC   14.1598

median DAM-index: LPC 2.976, saline -20.467
AUC (LPC vs saline): 0.9054
agreement with the paper's DAM-index of these cells (baseline from all cells of the atlas): Spearman 0.9877, AUC there 0.9045
```

The R and the Python package give the same scores (largest difference 1.4 × 10⁻¹⁴). Outputs: `demo/output_package_R.csv` and
`demo/output_package_py.csv`. To reproduce the paper's values exactly, use `run_demo.R` / `run_demo.py` above, which take the
baseline of the whole atlas from `data/genes.tsv.gz`.

## Run time

Measured on 2026-10-20 on an Apple M3 Mac with 16 GB RAM: the whole script including loading the data (the time each
script prints) and the peak memory (maximum resident set size):

| version | software | time | peak memory |
|---|---|---|---|
| `run_demo.R` | R 4.5.2, Matrix | 2.5 s | 0.25 GB |
| `run_demo.py` | Python 3.12.12, numpy 2.4.2, scipy 1.17.1, pandas 3.0.1 | 5.8 s | 0.22 GB |
| `run_demo_package.R` | R 4.5.2, Matrix, RFindsc 1.0.1 | 2.1 s | 0.36 GB |
| `run_demo_package.py` | Python 3.12.12, numpy 2.4.2, scipy 1.17.1, pandas 3.0.1, rfindsc 1.0.1 | 1.6 s | 0.22 GB |

No installation beyond R or Python with the packages above is needed.
