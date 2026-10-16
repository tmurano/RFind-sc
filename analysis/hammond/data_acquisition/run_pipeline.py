#!/usr/bin/env python3
"""
Hammond et al. 2019 microglia scRNA-seq — reproducible pipeline (end-to-end)

Usage:
    .venv/bin/python run_pipeline.py <OUTDIR>

Regenerates every product from scratch in <OUTDIR>:
  1. download the authors' (S.E. Marsh) annotated anndata from figshare
  2. clean-up -> annotated h5ad + CSV files
  3. UMAP from the log-normalized data -> figures + coordinates saved in the h5ad
  4. download / extract the raw counts of GEO GSE121654
  5. merge the raw counts into layers['counts'] (checks sum == nCount_RNA)

Source: Hammond TR et al., Immunity 2019, DOI 10.1016/j.immuni.2018.11.004
  raw data: NCBI GEO GSE121654 / annotated object: figshare 37624073 (S.E. Marsh)
"""
import os, re, sys, glob, subprocess, datetime
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGSHARE_URL = "https://ndownloader.figshare.com/files/37624073"
GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE121654&format=file"

OUT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "rerun")
os.makedirs(OUT, exist_ok=True)
P = lambda *a: os.path.join(OUT, *a)

def sh(cmd):
    print("  $", cmd)
    subprocess.run(cmd, shell=True, check=True)

def step(msg):
    print(f"\n=== {msg} ===", flush=True)

# ---------------------------------------------------------------- 1. download
step("1. Download annotated anndata from figshare")
raw_h5ad = P("Hammond_AllSamples_anndata.h5ad")
if not (os.path.exists(raw_h5ad) and os.path.getsize(raw_h5ad) == 366187535):
    sh(f'curl -L --retry 5 --retry-delay 3 -o "{raw_h5ad}" "{FIGSHARE_URL}"')
print("  size:", os.path.getsize(raw_h5ad))

# ---------------------------------------------------------------- 2. clean/build
step("2. Clean & build annotated h5ad + CSVs")
a = ad.read_h5ad(raw_h5ad)
orig = a.var_names.to_list()
cleaned = [re.sub(r"^X(?=\d)", "", g) for g in orig]
a.var["symbol_original"] = orig
a.var["symbol"] = cleaned
a.var_names = pd.Index(cleaned)
a.var_names_make_unique()
a.obs["cell_type_paper"] = a.obs["Paper_Cluster"].astype(str)
a.uns["provenance"] = {
    "dataset_title": "Single-Cell RNA Sequencing of Microglia throughout the Mouse Lifespan and in the Injured Brain Reveals Complex Cell-State Changes",
    "journal": "Immunity 50(1):253-271.e6 (2019)",
    "doi": "10.1016/j.immuni.2018.11.004", "pmid": "30471926", "pmcid": "PMC6655561",
    "raw_data_repository": "NCBI GEO GSE121654",
    "processed_object_source": "S.E. Marsh, github.com/samuel-marsh/Hammond-et-al_2019_Microglia_scRNAseq, figshare file 37624073",
    "expression_in_X": "log-normalized (Seurat LogNormalize)",
    "cell_type_column": "cell_type_paper / Paper_Cluster (14 clusters from paper)",
    "retrieved_date": datetime.date.today().isoformat(),
    "assembled_by": "Claude Code for T. Murano (Miyakawa Lab)",
}
a.uns["paper_cluster_levels"] = sorted(a.obs["Paper_Cluster"].astype(str).unique().tolist())
OUT_H5AD = P("Hammond2019_microglia_AllSamples_annotated.h5ad")
a.write_h5ad(OUT_H5AD, compression="gzip")
print("  wrote", os.path.basename(OUT_H5AD), a.shape)

meta = a.obs.copy(); meta.index.name = "cell_barcode"
meta.to_csv(P("Hammond2019_microglia_cell_metadata.csv"))
clusters = a.obs["Paper_Cluster"].astype(str)
levels = sorted(clusters.unique().tolist())
Xcsr = a.X.tocsr()
rows = {lv: np.asarray(Xcsr[np.where(clusters.values == lv)[0]].mean(axis=0)).ravel() for lv in levels}
mean_df = pd.DataFrame(rows, index=a.var_names).T; mean_df.index.name = "Paper_Cluster"
mean_df.to_csv(P("Hammond2019_microglia_clusterMean_expression.csv"))
comp = clusters.value_counts().rename_axis("Paper_Cluster").reset_index(name="n_cells")
comp.to_csv(P("Hammond2019_microglia_cluster_counts.csv"), index=False)
print("  wrote 3 CSVs")

# ---------------------------------------------------------------- 3. UMAP/plots
step("3. Compute UMAP + verification figures")
a.obs["Paper_Cluster"] = a.obs["Paper_Cluster"].astype("category")
sc.pl.embedding(a, basis="tsne", color="Paper_Cluster",
                title="Hammond 2019 - paper tSNE (Paper_Cluster)",
                legend_loc="on data", legend_fontsize=7, size=3, show=False)
plt.savefig(P("QC_fig1_paper_tSNE_by_cluster.png"), dpi=150, bbox_inches="tight"); plt.close()
a.raw = a
sc.pp.highly_variable_genes(a, n_top_genes=2000)
ah = a[:, a.var.highly_variable].copy()
sc.pp.scale(ah, max_value=10)
sc.tl.pca(ah, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(ah, n_neighbors=15, n_pcs=50)
sc.tl.umap(ah)
a.obsm["X_umap"] = ah.obsm["X_umap"]
fig, ax = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(a, color="Paper_Cluster", legend_loc="on data", legend_fontsize=7, size=3, ax=ax[0], show=False, title="Fresh UMAP - Paper_Cluster")
sc.pl.umap(a, color="Age", size=3, ax=ax[1], show=False, title="Fresh UMAP - Age")
plt.savefig(P("QC_fig2_fresh_UMAP.png"), dpi=150, bbox_inches="tight"); plt.close()
a.write_h5ad(OUT_H5AD, compression="gzip")
print("  wrote figures + UMAP coords")

# ---------------------------------------------------------------- 4. GEO raw
step("4. Download & extract GEO GSE121654 raw counts")
tar = P("GSE121654_RAW.tar"); rawdir = P("geo_raw"); os.makedirs(rawdir, exist_ok=True)
if not (os.path.exists(tar) and os.path.getsize(tar) > 1e8):
    sh(f'curl -L --retry 5 --retry-delay 3 -o "{tar}" "{GEO_URL}"')
sh(f'tar -xf "{tar}" -C "{rawdir}"')
print("  raw files:", len(glob.glob(P("geo_raw", "*.dge.txt.gz"))))

# ---------------------------------------------------------------- 5. merge counts
step("5. Merge raw counts into layers['counts']")
a = ad.read_h5ad(OUT_H5AD)
obs_names = a.obs_names.to_numpy(); obs_set = set(obs_names)
genes = a.var["symbol_original"].to_numpy()
gene_pos = {g: i for i, g in enumerate(genes)}; n_genes = len(genes)
gsm_to_prefix = (a.obs.groupby("GEO_ID_Number", observed=True)["Sample_Barcode_Prefix"]
                 .first().astype(str).to_dict())
blocks = []
for fn in sorted(glob.glob(P("geo_raw", "*.dge.txt.gz"))):
    gsm = re.match(r"^(GSM\d+)_", os.path.basename(fn)).group(1)
    if gsm not in gsm_to_prefix:
        continue
    prefix = gsm_to_prefix[gsm]
    df = pd.read_csv(fn, sep="\t", index_col=0)
    cols = prefix + "_" + df.columns.astype(str)
    keep = np.array([c in obs_set for c in cols])
    if keep.sum() == 0:
        continue
    df = df.loc[:, keep]; kept_bc = cols[keep].to_numpy()
    gidx = np.array([gene_pos.get(g, -1) for g in df.index.to_numpy()]); valid = gidx >= 0
    sub = sp.csr_matrix(df.values[valid].T.astype(np.float32)); coo = sub.tocoo()
    full = sp.coo_matrix((coo.data, (coo.row, gidx[valid][coo.col])),
                         shape=(sub.shape[0], n_genes)).tocsr()
    blocks.append((kept_bc, full))
all_bc = np.concatenate([b for b, _ in blocks])
all_mat = sp.vstack([m for _, m in blocks]).tocsr()
order = {bc: i for i, bc in enumerate(all_bc)}
counts = all_mat[np.array([order[bc] for bc in obs_names])]
assert counts.shape == a.shape
a.layers["counts"] = counts
cs = np.asarray(counts.sum(axis=1)).ravel()
ok = np.allclose(cs, a.obs["nCount_RNA"].values)
print("  matched cells:", counts.shape[0], "| sum==nCount_RNA:", ok)
assert ok, "raw count sum mismatch vs nCount_RNA"
a.uns["provenance"]["raw_counts_layer"] = "layers['counts'] = raw UMI counts from GEO GSE121654"
a.write_h5ad(OUT_H5AD, compression="gzip")
print("\nPIPELINE DONE ->", OUT_H5AD)
print("layers:", list(a.layers.keys()), "| obsm:", list(a.obsm.keys()))
