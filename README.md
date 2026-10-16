# RFind-sc — analysis code for the paper

Scripts that reproduce the analyses and figures of the RFind-sc manuscript from public data: the Hammond et al. 2019
microglia atlas, the De Donno et al. 2026 hippocampus stress atlas and the Shi et al. 2026 in vivo Perturb-seq atlas.

## For reviewers

**Quick check (under 1 minute, no download).** `demo/` holds 1,000 microglia of the paper's DAM task (white matter after
LPC vs saline injection, 500 each; raw UMI counts) and the DAM reference. From the project root:

```bash
Rscript demo/run_demo.R      # R >= 4.1 with Matrix; ~4 s on one core
python demo/run_demo.py      # Python >= 3.9 with numpy, scipy, pandas; ~4-5 s on one core
```

Each prints the first per-cell DAM-indices, the AUC of LPC vs saline cells (**0.9045**) and the largest difference from the
DAM-index stored by the paper's pipeline (`demo/expected_output.csv`; **< 1e-13**, reported as `OK`); details in
[`demo/README.md`](demo/README.md).

**System requirements.** Tested on macOS 27.0.1 (Apple M3, 16 GB RAM) with the versions in [Requirements](#requirements).
The demo needs < 0.5 GB of memory; the pipelines need 16 GB RAM (the Fig. 5 benchmark uses about 10 GB). No non-standard
hardware is needed.

## Companion packages

To score your own data, use the R package **RFindsc** (Seurat objects; <https://github.com/tmurano/RFindsc>) or the Python
package **rfindsc** (AnnData; <https://github.com/tmurano/rfindsc-py>), version 0.3.8 (tag `v0.3.8` in both). The scripts
here do not need them. Installation takes < 1 s for RFindsc when Matrix is installed and about 6 s for rfindsc, plus the
download:

```bash
mkdir -p /tmp/rfindsc_lib            # R, into a temporary library (needs the R package remotes)
Rscript -e 'remotes::install_github("tmurano/RFindsc@v0.3.8", lib = "/tmp/rfindsc_lib", upgrade = "never")'
python -m venv --system-site-packages /tmp/rfindsc_venv   # Python, into a temporary environment (needs git)
/tmp/rfindsc_venv/bin/pip install git+https://github.com/tmurano/rfindsc-py@v0.3.8
```

`analysis/check_rfindsc.R` compares RFindsc with the code here and runs its tests (after the Hammond and De Donno
pipelines): `git clone --branch v0.3.8 https://github.com/tmurano/RFindsc tools/RFindsc`, then
`RFINDSC_LIB=/tmp/rfindsc_lib Rscript analysis/check_rfindsc.R`.

## What RFind-sc computes

For each cell, genes are ranked by Δ<sub>ij</sub> = (x<sub>ij</sub> − e<sub>ij</sub>) / √(e<sub>ij</sub> + 1) (x, UMI count;
e, count expected from a baseline profile at the cell's total UMI count); the most positive and most negative genes (up to
500 each) are the cell's **UP** and **DOWN** lists. A **Running Fisher** test over every pair of cut-offs matches them with a
reference (a differential expression list or a gene set) and gives one signed **index** per cell and reference:

```
index = RF(UP, up) + RF(DOWN, down) − RF(UP, down) − RF(DOWN, up)        (gene set: RF(UP, set) − RF(DOWN, set))
```

A positive index means that the cell deviates from the baseline in the direction of the reference (a similarity, not a
test). Kernel: [`analysis/running_fisher.R`](analysis/running_fisher.R) (Python copy: `analysis/running_fisher.py`);
details in the paper's Methods and Supplementary Notes 1 and 2.

## Repository layout

```
demo/                a one-minute demo with 1,000 cells and its expected output
analysis/            the kernel (running_fisher.R; Python copy running_fisher.py with _hyper.py), check_rfindsc.R and,
                     one folder per dataset, the steps from public data to scores:
  hammond/           Hammond et al. 2019: QC -> doublets -> genes analysed -> Delta -> references -> scores
                     (data_acquisition/ = the assembly of the input h5ad)
  dedonno/           De Donno et al. 2026 (wild-type mice): QC/UMAP -> per-cell-type Delta -> references -> scores
  perturbseq/        Shi et al. 2026: pseudo-bulk per knock-out -> stress responses -> knock-out matching (results/)
figs/<item>/         analysis/ (scripts) and output/ of each figure and supplementary table; figs/run_all_figs.sh runs them
                     all; figs/make_source_data*.py write the Source Data workbooks (to submission/, not tracked)
data/*/provided/     fixed inputs; also tracked: data/Hammond2019/output/references/ and data/Perturbseq/urls.txt
logs/runs/           run logs (not tracked)
```

## Requirements

Versions used for the paper:

**R 4.5.2 (Bioconductor 3.22)** — data.table 1.18.2.1, Matrix 1.7.5, hdf5r 1.3.12, irlba 2.3.7, uwot 0.2.4, RANN 2.6.2,
ggplot2 4.0.3, ggrepel 0.9.8, patchwork 1.3.2, scales 1.4.0, gridExtra 2.3, gtable 0.3.6, ragg 1.5.2, readxl 1.4.5,
systemfonts, png, msigdbr 26.1.0, AnnotationDbi 1.72.0, org.Mm.eg.db 3.22.0, GO.db 3.22.0 and RFindsc 0.3.7; recount3 1.20.0
with SummarizedExperiment 1.40.0 (only to recompute the Seizure reference); for the Fig. 5 benchmark, UCell 2.14.0, AUCell
1.32.0, decoupleR 2.16.0 (with withr), singscore 1.30.0, GSEABase 1.72.0, VISION 3.0.2, Beyondcell 2.2.2 (with Seurat 5.5.0
and SeuratObject 5.4.0) and GSVA 2.4.9 (with BiocParallel); the others as installed with Bioconductor 3.22. VISION and
Beyondcell come from their GitHub tags (`YosefLab/VISION@v3.0.2`, `cnio-bu/beyondcell@v2.2.2`; Beyondcell needs DMwR 0.4.1
from the CRAN archive); `figs/Fig5_benchmark/analysis/install_r452_bench_work.sh` installs them, with GSVA, without upgrading
the packages already installed.

**Python 3.12.12** — numpy 2.4.2, scipy 1.17.1, pandas 3.0.1, matplotlib 3.10.8, h5py 3.15.1 and openpyxl 3.1.5; duckdb
(Perturb-seq steps 01 and 01b only); scrublet with scikit-learn (only to re-run the Hammond doublet step); anndata and
scanpy (only `analysis/hammond/data_acquisition/run_pipeline.py`).

The demo and the companion packages were also run with R 4.6.0 and Python 3.13.9 on an Apple-silicon Mac with 8 GB RAM.

## Input data

Files under `data/` are not tracked and must be obtained from the sources below, except `data/*/provided/`,
`data/Hammond2019/output/references/` and `data/Perturbseq/urls.txt`.

| file | content | source |
|---|---|---|
| `data/Hammond2019/raw/Hammond2019_microglia_AllSamples_annotated.h5ad` | microglia atlas with the authors' annotation (raw counts in `layers/counts`) | GEO GSE121654 and the annotated object of S.E. Marsh (figshare file 37624073), assembled by `$PYTHON analysis/hammond/data_acquisition/run_pipeline.py data/Hammond2019/raw` |
| `data/Hammond2019/provided/*.csv` | metadata and mitochondrial fraction of the 76,149 annotated cells, and the doublet calls of the paper | this study, from the sources above |
| `data/DeDonno2026/raw/dedonno2026.h5ad` | hippocampus atlas, Stress vs No Stress (only the wild-type mice are analysed) | CZ CELLxGENE Discover, dataset version <https://datasets.cellxgene.cziscience.com/9940f1a7-5fc8-4e52-aba4-6dbba8c90303.h5ad> (schema 7.0.0) |
| `data/DeDonno2026/external/GSE184924_Raw_gene_count_matrix.tsv.gz` | ventral hippocampus, corticosterone vs vehicle (Glucocorticoid reference; read only when the fixed table below is deleted) | GEO GSE184924 |
| `data/DeDonno2026/provided/umap_261015.csv`, `DeDonno_GR_glucocorticoid_261016.csv`, `DeDonno_KA6h_261016.csv`, `fig1a_sim_umap_261015.csv` | fixed UMAP of the wild-type cells, fold changes of the Glucocorticoid reference and of the Seizure reference (hippocampal excitatory neurons after kainate-induced seizures; GEO GSE125068 via recount3, SRP179090) and UMAP of the simulated cells of Fig. 1a | this study ([Fixed inputs](#fixed-inputs)) |
| `data/Matcovitch2016/raw/GSE79812_Timeline.xlsx` | microglial development (Maturity reference) | GEO GSE79812 |
| `data/KerenShaul2017/provided/Mouse_DAM_KerenShaul2017_vs_homeostatic.csv` | DAM vs homeostatic microglia (DAM reference) | Keren-Shaul et al. 2017, Table S3 |
| `data/GSE169505/raw/GSE169505_Stressome_GeneCounts.txt.gz`, `GSE169505-GPL21103_series_matrix.txt.gz`, `GSE169505-GPL24247_series_matrix.txt.gz` | hippocampus after a 6-min forced swim (Forced swim reference: experiment TimeSeriesFST, male ventral hippocampus 3 h after the swim, 7 mice, vs home cage, 8 mice) | GEO GSE169505 |
| `data/GSE75431/raw/E-GEOD-75431-raw-counts.tsv` | cortical microglia, PS2APP vs non-transgenic, 13 months (PS2APP reference) | EMBL-EBI Expression Atlas E-GEOD-75431 ("raw counts") |
| downloaded by msigdbr | Hallmark gene sets (Interferon, G2M, NF-κB, UPR) | MSigDB via msigdbr 26.1.0 |
| `data/Perturbseq/urls.txt`, `gene_metadata.parquet` | URLs of the 317 parquet files of the atlas (tracked) and its gene table | Hugging Face `perturbai/wholebrain_crispr_atlas` (CC BY 4.0): `curl -L -o data/Perturbseq/gene_metadata.parquet https://huggingface.co/datasets/perturbai/wholebrain_crispr_atlas/resolve/a7c65dc0a64da4bd47cf6ef5f4dec6c7ef745e87/metadata/gene_metadata.parquet` |

The Perturb-seq analyses (Fig. 4, Supplementary Fig. S3) read five pseudo-bulk files in `data/Perturbseq/`: `pb_g017_DG.npz` and
`pb_g017_CA3.npz` (made by `01_pseudobulk_subclass.py` from the group-017 parquet shards, extracted from the atlas outside this
repository; put the shards in one folder and set `PERTURB_DIR` to it), `pb_CA1.npz` and `pb_GABA3.npz` (`01b_scan_subclasses.py`
and `01c_combine_gaba.py`) and `pb_g151.npz` (thalamus; from the same extraction). All five are available at Zenodo
(https://doi.org/10.5281/zenodo.23409227, CC BY 4.0); put them in `data/Perturbseq/` to skip the `01` steps below.

## How to run

Everything runs from the project root with relative paths (set `R_HOME` and `PYTHON`;
[Environment variables](#environment-variables)). In total, about 10 min for the Hammond and De Donno pipelines, 2.5-4.5 h
for Perturb-seq and about 5 h for the figures and tables.

```bash
bash analysis/hammond/run_all.sh                        # ~4 min
bash analysis/dedonno/run_all.sh                        # ~5 min

# Perturb-seq (Fig. 4, Supplementary Fig. S3)
PERTURB_DIR=/path/to/g017_shards $PYTHON analysis/perturbseq/01_pseudobulk_subclass.py   # PERTURB_DIR is required
$PYTHON analysis/perturbseq/01b_scan_subclasses.py      # 2-4 h, reads the atlas over HTTPS (81 GB; restartable);
                                                        # writes ~10 GB of part files (SCAN_PARTS)
$PYTHON analysis/perturbseq/01c_combine_gaba.py
$PYTHON analysis/perturbseq/02_stress_axis.py           # seconds
$PYTHON analysis/perturbseq/03_axis_vs_ko.py            # under a minute

bash figs/run_all_figs.sh                               # ~5 h (Fig. 5 benchmark ~3 h, always recomputed;
                                                        # rest of Fig. 5 ~1 h; Supplementary Fig. S3 ~1 h)
```

`figs/run_all_figs.sh` continues after a failed step and writes a summary to `logs/runs/run_all_figs_summary.tsv`; the Fig. 5
benchmark (`run_fig5_redo.sh`) skips the reference lengths it has finished when it is started again (`FIG5_LOGD` resumes a run
of another day). The figures of the paper were assembled from the panel PNGs.

## Which script makes what

Scripts are in `figs/<item>/analysis/` and write to the `output/` folder next to it; `figs/run_all_figs.sh` runs them in the
order of their inputs.

| item | script(s) |
|---|---|
| Fig. 1 | `fig1ab_simulated.R` and `fig1b_3d.py` (a, b; simulated cells), `fig1c_formula.py` (c), `fig1d_matching.R` (d) |
| Fig. 2 | `fig2_tsne_original.R` (a), `fig2_cluster_table.R` (b), `fig2_violins.R` (c, d; statistics in `Fig2_stats_table.csv`, where `type` B = differential expression list and C = gene set), `fig2_scatter_pairs.R` (e, f) |
| Fig. 3 | `fig3_umap.R` (a), `fig3_celltype_table.R` with `fig3_celltype_table_png.R` (b), `fig3_violins.R` (c, d), `fig3_scatter.R` (e, f); cell types from `fig3_celltypes.R` |
| Fig. 4 | `fig4a_schematic.py` (a), `fig4b_schematic.py` (b), `fig4c_rank_tables.py` (c; CA1 Glut), `fig4d_go_set.R` and `fig4d_heatmap.R` (d, e), `fig4f_gr_expression.py` (f); `fig4_ranking_checks.py` (the check of the rankings in Methods) |
| Fig. 5 | `run_fig5_redo.sh` (`fig5_bench_redo.R` at every reference length), `fig5_extra_refs.R` and `fig5_extra_score.R` (the four added type A tasks), then `fig5_ucell_zero.R`, `fig5_best_length.py`, `fig5_best_length_8tasks.py`, `fig5_balanced_accuracy.py`, `fig5_redo_figure.py` (b) and `fig5_rebuild_figure_draft.py` (other panels) |
| Supplementary Figs. S1-S4 | `figS1_running_fisher.py`; `figS2_delta_benchmark.R`, `figS2_panels.py`; `figS3_rfind_allpairs.py`, `figS3_rfind_panels.py`, `figS3e_celltypes.py`, `figS3_schematic.R`; `figS4_capability_map.py` (after the Fig. 5 steps) |
| Supplementary Tables 1-4 | `supp_table_ko_rankings.py` (1); `fig4_random_set_null.py`, `fig4_other_methods.py` (after the Supplementary Fig. S3 steps) and `supp_table_random_sets.py` (2); `fig5_redo_tables.py` (3; after the Fig. 5 steps); `supp_table_references.R` and `supp_table_references.py` (4; in `figs/Supp_Table4_References/`) |
| Source Data | `figs/make_source_data.py` (Figs. 2-5 and Supplementary Fig. S4) and `figs/make_source_data_supp.py` (Supplementary Figs. S1-S3), written to `submission/` |

Run alone: `figS3_shi_pairs.py` (the pairs named by Shi et al., imported by the Supplementary Fig. S3 scripts).

## Environment variables

- `R_HOME` — R with the packages above (default `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources`); install
  packages with `"$R_HOME/bin/exec/R"`, not `R CMD INSTALL`, whose wrapper may point to another R version.
- `PYTHON` — Python with the packages above (default: the first of `python3` and `python` on the `PATH` that imports numpy,
  scipy, pandas, matplotlib and openpyxl).
- `PERTURB_DIR` (required by `01_pseudobulk_subclass.py`: the folder with the group-017 parquet shards), `SCAN_PARTS` (part
  files of `01b_scan_subclasses.py`, about 10 GB, default `data/Perturbseq/scan_parts`; finished parts are skipped) and
  `FIG5_LOGD` (log folder of the Fig. 5 benchmark, default `logs/runs/fig5_redo_<yymmdd>`).
- For checks only: `FIG5_TEST_DIR` (output folder of the test mode of the Fig. 5 scripts) and `PERTURB_RES`, `PERTURB_MINC`,
  `PERTURB_DROP` and `PERTURB_RUNS` (output folder, minimum number of nuclei, genes left out of the matching, and the stress
  responses and cell types to run in `analysis/perturbseq/03_axis_vs_ko.py`).

## Reproducibility

### Processing values

| | Hammond et al. 2019 | De Donno et al. 2026 | Shi et al. 2026 (Perturb-seq) |
|---|---|---|---|
| fields read | raw counts in `layers/counts`; clusters from `cell_meta_hammond.csv` | raw counts in `raw/X` (genes, `raw/var/feature_name`); `obs/line`, `obs/condition`, `obs/donor_id`, `obs/cell_type_short` | `passes_qc`, `num_guides`, `predicted_subclass`, `gene_target` |
| mice and cells used | 41 mice, one sample each (the sample name is the prefix of the cell ID, e.g. `E14_F_B10`) | the 12 wild-type mice (`line` = `Ctrl`; `donor_id`) | nuclei with `passes_qc` and `num_guides` = 1 |
| mitochondrial threshold | 1.81% (the 5-MAD limit; median 0.44%, MAD 0.27%) | 5% (the cap; lower than the 5-MAD limit) | — |
| cells | 76,149 released; 327 removed by the mitochondrial and UMI rules and 250 more as doublets (253 Scrublet calls); 75,572 kept | 40,674 cells of the wild-type mice released; 26,246 after quality control; 26,235 after removing cell types with < 20 cells (Choroid Plexus, 11 cells) | — |
| doublets | Scrublet 0.2.3 on each sample: `expected_doublet_rate=0.06, random_state=0`; `scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)`; threshold 0.25 when the automatic threshold fails | removed by the authors before release | — |
| cell types | the authors' 13 microglial clusters (monocytes/macrophages kept in the baseline, not shown) | the authors' `cell_type_short` with numbered clusters merged (e.g. CA1 Glut 1–3 → CA1 Glut); 21 cell types | subclasses 037 DG Glut, 016 CA1-ProS Glut, 017 CA3 Glut; GABA = 052 Pvalb Gaba + 053 Sst Gaba + 049 Lamp5 Gaba |
| genes analysed | 8,960 | 11,647 | genes with CPM ≥ 1 in the non-targeting column |

**Gene names excluded** (regular expressions, as in the code): `Rik[0-9]?$`, `^Gm[0-9]`, `^mt-`, `^A[CFLI][0-9]` and `^[0-9]`
from the genes analysed and from the differential expression lists; the six lists also exclude `\.[0-9]+$`, genes without a
symbol and haemoglobin genes `^Hb[ab]-` (from blood cells in the samples). Ribosomal protein genes are kept.

**Differential expression lists** (one rule for all six): counts summed per group and converted to CPM over the nuclear genes
(mitochondrial genes left out of the library size; DAM and Maturity keep the form of the published values); the top 500 genes
per direction (Seizure: the top 200 up-regulated genes, without an expression floor) are stored and restricted to the genes
analysed when the cells are scored; the background of the matching is the genes analysed with a non-zero change in the list.
The numbers of genes before and after this restriction are in Supplementary Table 4.

### Fixed inputs

Steps that depend on software versions or threads use fixed copies of their outputs when present (delete a file to
recompute it):

- `data/Hammond2019/provided/scrublet_doublets_260914.csv` — Scrublet 0.2.3 with Python 3.10.19, numpy 1.26.4, scipy 1.15.2,
  scikit-learn 1.7.2 and annoy 1.17.3 reproduce every call.
- `data/DeDonno2026/provided/umap_261015.csv` (UMAP of the wild-type cells) — depends on uwot / irlba / BLAS.
- `data/DeDonno2026/provided/DeDonno_GR_glucocorticoid_261016.csv` (fold changes of the Glucocorticoid reference, from the GSE184924
  counts) — depends on the org.Mm.eg.db release (3.22.0) used to map Ensembl identifiers to gene names.
- `data/DeDonno2026/provided/fig1a_sim_umap_261015.csv` — UMAP of the simulated cells of Fig. 1a (uwot 0.2.4).
- `data/DeDonno2026/provided/DeDonno_KA6h_261016.csv` — fold changes of the Seizure reference from recount3 (genes with tied
  log2 fold changes can come in another order in a new download).
- `data/Hammond2019/provided/cell_meta_hammond.csv`, `mt_frac_hammond.csv` — fixed inputs that no script here rebuilds.

Some file names and data keys use short names of the references: CORT = Glucocorticoid, KA6h = Seizure, AcuteStress =
Forced swim, DAM_PS2APP = PS2APP.

### Random ties, gene names and data versions

- `AUCell_buildRankings` breaks tied ranks at random; the benchmark scripts set the seed just before each call.
- The annotated Hammond h5ad stores "-" as "." (e.g. `H2.D1` for `H2-D1`); `analysis/hammond/03_extract_qc.R` restores "-"
  where the dashed form is an MGI symbol and the dotted form is not (org.Mm.eg.db).
- Hugging Face, CELLxGENE and MSigDB are fixed by the repository commit, the dataset-version URL and the package version
  ([Input data](#input-data)); the CELLxGENE collection was revised on 2026-06-11 (schema 7.1.0), and the paper uses the
  schema 7.0.0 file. Downloads: CELLxGENE 2026-05-12; GEO GSE121654 and figshare 2026-05-27; Hugging Face 2026-08-14,
  2026-08-21 to 24 and 2026-09-15; GSE169505 and E-GEOD-75431 2026-10-05; recount3 2026-10-16.

## Citation

To be added (manuscript in preparation).

## License

MIT (see `LICENSE`). The data under `data/*/provided/`, `data/Hammond2019/output/references/`,
`figs/Fig5_benchmark/output/extra_refs/fig5_extra_refs.rds` and `demo/data/` are derived from the public datasets listed under
[Input data](#input-data) (for example the reference gene lists, including the MSigDB gene sets in `references_panel9.rds`,
the cell metadata, the recount3 table of the Seizure reference and the Hammond counts of the demo) and remain subject to the
terms of those sources; cite the original publications when you use them.
