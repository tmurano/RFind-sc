# RFind-sc — analysis code for the paper

> RFind-sc has three repositories: this analysis code of the paper and two packages that implement the method, the R package [RFindsc](https://github.com/tmurano/RFindsc) (for Seurat objects) and the Python package [rfindsc](https://github.com/tmurano/rfindsc-py) (for AnnData; the same scores). To use RFind-sc on your own data, install one of the packages.

Scripts that reproduce the analyses and figures of the RFind-sc manuscript from public data: the Hammond et al. 2019
microglia atlas, the De Donno et al. 2026 hippocampus stress atlas, the Shi et al. 2026 in vivo Perturb-seq atlas and the
Kang et al. 2018 human PBMC data (IFN-β). RFind-sc ranks the genes of each cell by their deviation from the expected counts
(Delta) and matches the ranking with external references by the Running Fisher test (kernel:
[`analysis/running_fisher.R`](analysis/running_fisher.R), Python copy `analysis/running_fisher.py`; details in the paper's
Methods and Supplementary Notes 1 and 2).

## For reviewers

**Quick check (under 1 minute, no download).** `demo/` holds 1,000 microglia of the paper's DAM task (white matter after
LPC vs saline injection, 500 each; raw UMI counts) and the DAM reference. From the project root:

```bash
Rscript demo/run_demo.R      # R >= 4.1 with Matrix; ~3 s on one core
python demo/run_demo.py      # Python >= 3.9 with numpy, scipy, pandas; ~6 s on one core
```

Each prints the AUC of LPC vs saline cells (**0.9045**) and the largest difference from the DAM-index stored by the paper's
pipeline (**< 1e-13**, reported as `OK`); details in [`demo/README.md`](demo/README.md).

**System requirements.** Tested on macOS 27.0.1 (Apple M3, 16 GB RAM) with the versions in [Requirements](#requirements):
the whole pipeline, without the Perturb-seq steps that read the atlas, ran from the public data in about 5 h (2026-10-19).
A run on a second computer (Apple M3, 8 GB RAM) with other versions of some R packages gave the same output tables, except
the run times and the simulated cells of Fig. 1b ([Fixed inputs](#fixed-inputs)). The demo needs < 0.5 GB of memory and the
pipelines 16 GB RAM; no non-standard hardware is needed.

## Companion packages

To score your own data, use the R package **RFindsc** (Seurat; <https://github.com/tmurano/RFindsc>) or the Python package
**rfindsc** (AnnData; <https://github.com/tmurano/rfindsc-py>), version 1.0.1; the scripts here do not need them. Each
installs in a few seconds plus the download:

```bash
Rscript -e 'remotes::install_github("tmurano/RFindsc@v1.0.1", upgrade = "never")'
pip install git+https://github.com/tmurano/rfindsc-py@v1.0.1
```

`analysis/check_rfindsc.R` compares RFindsc with the code here and runs its tests (after the Hammond and De Donno pipelines):
`git clone --branch v1.0.1 https://github.com/tmurano/RFindsc tools/RFindsc`, then
`RFINDSC_LIB=<library with RFindsc> Rscript analysis/check_rfindsc.R`.

## Repository layout

- `demo/`: the one-minute demo.
- `analysis/`: the kernel, `check_rfindsc.R` and one folder per dataset (`hammond/`, `dedonno/`, `kang/`, `perturbseq/`) with
  the steps from public data to scores.
- `figs/<item>/`: the scripts (`analysis/`) and outputs (`output/`) of each figure and supplementary table
  ([Which script makes what](#which-script-makes-what)); `figs/make_source_data*.py` write the Source Data workbooks.
- `data/*/provided/`: fixed inputs; `logs/runs/`: run logs (not tracked).

## Requirements

**R 4.5.2 (Bioconductor 3.22)**: data.table 1.18.2.1, Matrix 1.7.5, hdf5r 1.3.12, irlba 2.3.7, uwot 0.2.4, RANN 2.6.2,
ggplot2 4.0.3, ggrepel 0.9.8, patchwork 1.3.2, scales 1.4.0, gridExtra 2.3, gtable 0.3.6, ragg 1.5.2, readxl 1.4.5,
systemfonts, png, msigdbr 26.1.0, AnnotationDbi 1.72.0, org.Mm.eg.db 3.22.0, GO.db 3.22.0, RFindsc 1.0.1, and recount3 1.20.0
with SummarizedExperiment 1.40.0 (Seizure reference only). For the Fig. 5 benchmark: UCell 2.14.0, AUCell 1.32.0, decoupleR
2.16.0, singscore 1.30.0, GSEABase 1.72.0, VISION 3.0.2 and Beyondcell 2.2.2 (GitHub tags `YosefLab/VISION@v3.0.2` and
`cnio-bu/beyondcell@v2.2.2`; Beyondcell needs DMwR 0.4.1 from the CRAN archive), Seurat 5.5.0, SeuratObject 5.4.0 and GSVA
2.4.9; `figs/Fig5_benchmark/analysis/install_r452_bench_work.sh` installs them.

**Python 3.12.12**: numpy 2.4.2, scipy 1.17.1, pandas 3.0.1, matplotlib 3.10.8, h5py 3.15.1 and openpyxl 3.1.5; duckdb
(Perturb-seq steps 01 and 01b), scrublet with scikit-learn (only to re-run the Hammond doublet step) and anndata and scanpy
(only `analysis/hammond/data_acquisition/run_pipeline.py`).

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
| `data/Kang2018/raw/GSE96583_batch2.genes.tsv.gz`, `GSE96583_batch2.total.tsne.df.tsv.gz`, `GSM2560248_2.1.mtx.gz`, `GSM2560248_barcodes.tsv.gz`, `GSM2560249_2.2.mtx.gz`, `GSM2560249_barcodes.tsv.gz` | human PBMC of eight donors, unstimulated and IFN-β 6 h (batch 2): UMI counts, and the authors' t-SNE, demuxlet calls and cell types | GEO GSE96583 |
| `data/GSE191317/raw/GSE191317_Raw_gene_counts_matrix.txt.gz` | PBMC of seven donors, IFN-β vs mock, 18 h (IFNb reference of Supplementary Fig. S5) | GEO GSE191317 |
| `data/GSE272019/raw/GSE272019_rawcounts_PGE.txt.gz` | CD14+ monocytes of three donors, PGE2 vs unstimulated, 3 h (PGE2 reference of Supplementary Fig. S5) | GEO GSE272019 |
| `data/pbmc3k/raw/filtered_gene_bc_matrices/hg19/` | 2.7k PBMC of 10x Genomics, for the markers of Supplementary Fig. S5b (the Seurat tutorial that Kang et al. used to annotate the cells) | 10x Genomics dataset "3k PBMCs from a Healthy Donor", `pbmc3k_filtered_gene_bc_matrices.tar.gz` (SHA-256 847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5), unpacked in `data/pbmc3k/raw/` |
| downloaded by msigdbr | Hallmark gene sets (Interferon, G2M, NF-κB, UPR), and for Supplementary Fig. S5 Hallmark INTERFERON_ALPHA_RESPONSE and KEGG CHEMOKINE_SIGNALING_PATHWAY (C2, KEGG legacy) with human genes | MSigDB via msigdbr 26.1.0 |
| `data/Perturbseq/urls.txt`, `gene_metadata.parquet` | URLs of the 317 parquet files of the atlas (tracked) and its gene table | Hugging Face `perturbai/wholebrain_crispr_atlas` (CC BY 4.0): `curl -L -o data/Perturbseq/gene_metadata.parquet https://huggingface.co/datasets/perturbai/wholebrain_crispr_atlas/resolve/a7c65dc0a64da4bd47cf6ef5f4dec6c7ef745e87/metadata/gene_metadata.parquet` |

The Perturb-seq analyses (Fig. 4, Supplementary Fig. S3) read five pseudo-bulk files in `data/Perturbseq/`: `pb_g017_DG.npz` and
`pb_g017_CA3.npz` (made by `01_pseudobulk_subclass.py` from the group-017 parquet shards, extracted from the atlas outside this
repository; put the shards in one folder and set `PERTURB_DIR` to it), `pb_CA1.npz` and `pb_GABA3.npz` (`01b_scan_subclasses.py`
and `01c_combine_gaba.py`) and `pb_g151.npz` (thalamus; from the same extraction). All five are available at Zenodo
(https://doi.org/10.5281/zenodo.23409227, CC BY 4.0); put them in `data/Perturbseq/` to skip the `01` steps below.

## How to run

Everything runs from the project root with relative paths (set `R_HOME` and `PYTHON`;
[Environment variables](#environment-variables)). In total, about 7 min for the Hammond, De Donno and Kang pipelines,
2.5-4.5 h for Perturb-seq, about 1 h for the run times and about 4 h for the figures and tables (2026-10-19; Apple M3,
16 GB RAM).

```bash
bash analysis/hammond/run_all.sh                        # ~4 min
bash analysis/dedonno/run_all.sh                        # ~2 min
bash analysis/kang/run_all.sh                           # ~1 min (Supplementary Fig. S5)

# Perturb-seq (Fig. 4, Supplementary Fig. S3)
PERTURB_DIR=/path/to/g017_shards $PYTHON analysis/perturbseq/01_pseudobulk_subclass.py   # PERTURB_DIR is required
$PYTHON analysis/perturbseq/01b_scan_subclasses.py      # 2-4 h, reads the atlas over HTTPS (81 GB; restartable);
                                                        # writes ~10 GB of part files (SCAN_PARTS)
$PYTHON analysis/perturbseq/01c_combine_gaba.py
$PYTHON analysis/perturbseq/02_stress_axis.py           # seconds
$PYTHON analysis/perturbseq/03_axis_vs_ko.py            # under a minute

bash figs/Fig5_benchmark/analysis/fig5_runtime.sh       # ~1 h on an otherwise idle computer: the run times of
                                                        # Supplementary Table 3 (before run_all_figs.sh, which reads them)
bash figs/run_all_figs.sh                               # ~4 h (Fig. 5 benchmark ~2 h, always recomputed;
                                                        # rest of Fig. 5 ~1 h; Supplementary Fig. S3 ~1 h)
```

`figs/run_all_figs.sh` continues after a failed step and writes a summary to `logs/runs/run_all_figs_summary.tsv`; the Fig. 5
benchmark (`run_fig5_redo.sh`) skips the reference lengths it has finished when it is started again (`FIG5_LOGD` resumes a run
of another day). The figures of the paper were assembled from the panel PNGs.

## Which script makes what

Scripts are in `figs/<item>/analysis/`; `figs/run_all_figs.sh` runs them in the order of their inputs.

| item | script(s) |
|---|---|
| Fig. 1 | `fig1ab_simulated.R`, `fig1b_3d.py` (a, b), `fig1c_formula.py` (c), `fig1d_matching.R` (d) |
| Fig. 2 | `fig2_tsne_original.R` (a), `fig2_cluster_table.R` (b), `fig2_violins.R` (c, d), `fig2_scatter_pairs.R` (e, f) |
| Fig. 3 | `fig3_umap.R` (a), `fig3_celltype_table.R` (b), `fig3_violins.R` (c, d), `fig3_scatter.R` (e, f) |
| Fig. 4 | `fig4a_schematic.py` (a), `fig4b_schematic.py` (b), `fig4c_rank_tables.py` (c), `fig4d_go_set.R` (d), `fig4e_expression_vs_score.py` (e) |
| Fig. 5, Supplementary Fig. S4 | `run_fig5_redo.sh` (the benchmark), `fig5_extra_refs.R`, `fig5_extra_score.R` and the later steps of `figs/run_all_figs.sh` |
| Supplementary Figs. S1–S3, S5 | `figS1_running_fisher.py`; `figS2_*`; `figS3_*`; `figS5_*` (`figS5_compose.py` assembles S5) |
| Supplementary Tables 1–4 | `supp_table_ko_rankings.py` (1), `supp_table_random_sets.py` (2), `fig5_redo_tables.py` (3), `supp_table_references.R` and `.py` (4) |

## Environment variables

`R_HOME` (R with the packages above; default `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources`; install packages
with `"$R_HOME/bin/exec/R"`, whose wrapper cannot point to another R version), `PYTHON` (default: the first of `python3` and
`python` that imports numpy, scipy, pandas, matplotlib and openpyxl), `PERTURB_DIR` (required by `01_pseudobulk_subclass.py`),
`SCAN_PARTS` (part files of `01b_scan_subclasses.py`, about 10 GB) and `FIG5_LOGD` (log folder of the Fig. 5 benchmark). The
other variables, used only for checks, are described in the scripts that read them.

## Reproducibility

The processing (quality control, genes analysed and references) follows the paper's Methods; the numbers of genes of the
references are in Supplementary Table 4. Some file names and data keys use short names of the references: CORT =
Glucocorticoid, KA6h = Seizure, AcuteStress = Forced swim, DAM_PS2APP = PS2APP.

### Fixed inputs

Steps that depend on software versions or threads use fixed copies of their outputs when present (delete a file to
recompute it):

- `data/Hammond2019/provided/scrublet_doublets_260914.csv` — Scrublet 0.2.3 with Python 3.10.19, numpy 1.26.4, scipy 1.15.2,
  scikit-learn 1.7.2 and annoy 1.17.3 reproduce every call.
- `data/DeDonno2026/provided/umap_261015.csv` (UMAP of the wild-type cells) — depends on uwot / irlba / BLAS.
- `data/DeDonno2026/provided/DeDonno_GR_glucocorticoid_261016.csv` (fold changes of the Glucocorticoid reference, from the GSE184924
  counts) — depends on the org.Mm.eg.db release (3.22.0) used to map Ensembl identifiers to gene names.
- `data/DeDonno2026/provided/fig1a_sim_umap_261015.csv` — UMAP of the simulated cells of Fig. 1a (uwot 0.2.4). The simulated
  counts themselves depend on irlba and RANN (principal components and nearest neighbours of the real cells), so with other
  versions the cells of Fig. 1b, an illustration, take other indices (with irlba 2.4.1 and RANN 2.6.3, r = 0.79–0.88 with the
  indices in the paper).
- `data/DeDonno2026/provided/DeDonno_KA6h_261016.csv` — fold changes of the Seizure reference from recount3 (genes with tied
  log2 fold changes can come in another order in a new download).
- `data/Hammond2019/provided/cell_meta_hammond.csv`, `mt_frac_hammond.csv` — fixed inputs that no script here rebuilds.


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
