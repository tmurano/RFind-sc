# Hammond et al. 2019 microglia scRNA-seq: assembling the input h5ad

Procedure of 2026-05-27 (data retrieved on that date).

## 1. Purpose

Obtain the single-cell RNA-seq data of microglia of Hammond et al. (*Immunity* 2019) and put them in a form that the
analysis can use (an AnnData `.h5ad` file) with:

1. **gene expression**, both log-normalized and as raw UMI counts;
2. **the cell types (clusters) defined in the paper**.

The source of every product is recorded so that its provenance can be traced (PROVENANCE.txt, and `uns['provenance']`
in the h5ad).

Paper:
- Hammond TR, Dufort C, Dissing-Olesen L, ... Marsh SE, ... McCarroll SA, Stevens B.
  *Immunity* 50(1):253-271.e6 (2019).
- DOI: 10.1016/j.immuni.2018.11.004 / PMID: 30471926 / PMCID: PMC6655561
- 76,149 cells from development to old age and after brain injury (microglia and one cluster of
  monocytes/macrophages, Mono/Mac); the paper defines 9 microglial state clusters (14 annotated clusters with the
  subclusters and Mono/Mac).

## 2. Sources

| type | location |
|------|----------|
| raw counts | NCBI GEO **GSE121654** (`GSE121654_RAW.tar`, `.dge.txt.gz` files of 47 samples) |
| annotated object | released by the author S.E. Marsh, figshare file **37624073** (AnnData version) |
| repository of the annotated object | github.com/samuel-marsh/Hammond-et-al_2019_Microglia_scRNAseq |
| interactive portal | https://www.microgliasinglecell.com/ |

> Why two sources: the GEO raw data have **no cell-type annotation**, while the annotated object of S.E. Marsh holds
> **the paper's cluster annotation and the normalized expression** but no raw counts. The annotated object is therefore
> taken as the base, and the GEO raw counts are added to it.

## 3. Procedure (how to reproduce)

All steps are in `run_pipeline.py`. From the project root:

```bash
# environment (once)
python3 -m venv .venv
.venv/bin/pip install anndata scanpy pandas scipy matplotlib

# run (downloads from figshare and GEO and writes every product to the folder given as the argument)
.venv/bin/python analysis/hammond/data_acquisition/run_pipeline.py data/Hammond2019/raw
```

The pipelines of this repository read `data/Hammond2019/raw/Hammond2019_microglia_AllSamples_annotated.h5ad`, so give
`data/Hammond2019/raw` as the output folder (or move the h5ad there afterwards). The script needs `curl` and `tar`, and
downloads about 0.5 GB.

Steps inside the script:

1. **Download the annotated AnnData** (figshare 37624073, 366 MB).
2. **Clean up and reshape**:
   - remove the `X` prefix that R added to gene names (e.g. `X0610007P14Rik` -> `0610007P14Rik`); the original names
     are kept in `var['symbol_original']` (the key for matching with the GEO raw data);
   - store the paper's clusters as `obs['cell_type_paper']`;
   - embed the provenance in `uns['provenance']`;
   - write CSV files: cell metadata, mean expression per cluster, number of cells per cluster.
3. **UMAP and check figures**: HVG (2,000) -> scale -> PCA (50) -> neighbours -> UMAP. The paper's stored t-SNE and the
   new UMAP are drawn coloured by cluster (`QC_fig1`, `QC_fig2`); the UMAP coordinates are stored in the h5ad.
4. **Download and extract the raw counts** (GSE121654_RAW.tar, 142 MB -> 47 files).
5. **Add the raw counts**: each cell is matched to the annotated object by **GSM number + barcode** and its counts are
   stored in `layers['counts']`.

### Points to note in the implementation

- **Mismatched barcode prefixes**: the cell names of the annotated object come from `Sample_Barcode_Prefix`, but the GEO
  file names come from `Sample_GEO_Prefix`, and the numbers differ for the P30 and V2 samples (e.g. barcode `P30_Male_5`
  vs GEO file `P30_Male_1`). Matching through the GSM number (`GEO_ID_Number`) matches every cell correctly.
- **Gene names**: both the GEO raw data and the annotated object originally carry the `X` prefix; matching uses
  `symbol_original` (with the prefix).
- **The six percoll/nopercoll samples** (GSM3442047-3442052) are not in the paper's final object and are skipped
  automatically (their GSM numbers are absent from the annotated object).

### Checks

- After the merge, **the raw-count sum of every cell equals `obs['nCount_RNA']`** (asserted in the script; an exact
  match with the paper's metadata shows that the matching is correct).
- All 76,149 cells are matched.

## 4. Products

Written to the output folder given to `run_pipeline.py` (PROVENANCE.txt describes them):

| file | content |
|------|---------|
| **Hammond2019_microglia_AllSamples_annotated.h5ad** | main product: 76,149 cells x 20,579 genes; `X` = log-normalized, `layers['counts']` = raw counts, `obs['cell_type_paper']` = 14 clusters, `obsm` = t-SNE/UMAP, `uns['provenance']` = sources |
| Hammond2019_microglia_cell_metadata.csv | metadata of all cells (including the cell type) |
| Hammond2019_microglia_clusterMean_expression.csv | mean log-normalized expression, 14 clusters x genes |
| Hammond2019_microglia_cluster_counts.csv | number of cells per cluster |
| QC_fig1_paper_tSNE_by_cluster.png / QC_fig2_fresh_UMAP.png | check figures |

The script also leaves its downloads in the output folder (the figshare object `Hammond_AllSamples_anndata.h5ad`,
`GSE121654_RAW.tar` and the extracted files in `geo_raw/`).

This folder keeps `run_pipeline.py`, PROVENANCE.txt, this file and the two check figures of the run of 2026-05-27.
The h5ad and the CSV files are not tracked in git (README, Input data).

## 5. Properties of the data and quality control

- The expression in `X` is log-normalized (Seurat LogNormalize); `layers['counts']` holds the raw UMI counts.
- The cells were filtered by the authors before release (cells with fewer than 650 detected genes and contaminating
  cell types were removed; Hammond et al., 2019). The mitochondrial genes are not in the annotated object, so the
  mitochondrial fraction cannot be computed from it; the quality control of this study uses the fraction computed from
  the GEO DGE files (`data/Hammond2019/provided/mt_frac_hammond.csv`, copied to `data/Hammond2019/raw/` by
  `analysis/hammond/run_all.sh`; README, Input data) and is applied by the Hammond pipeline
  (`analysis/hammond/01_qc_mt_umi.R`, `02_scrublet.py`, `03_extract_qc.R`; Methods, Quality control), not by this script.
- The cells are microglia plus one cluster of monocytes/macrophages (Mono/Mac), so the cell-type annotation is the set
  of state clusters (`Paper_Cluster`: 1, 2a-2c, 3-6, 7a-7c, 8, 9, Mono/Mac).

## 6. Reproducibility check

- `run_pipeline.py` was run again from scratch in a separate folder, and its products were compared with those of the
  first run (shape, obs, var, X (log-normalized), raw counts, t-SNE, cell types and the CSV files); they matched.
- The intermediate files (the raw downloads from figshare and GEO, the extracted raw data and the copy made for the
  check) are not kept; they are regenerated by the script.
