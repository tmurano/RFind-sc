"""Hammond step 2 — doublet removal (Scrublet, per sample)
Input: layers/counts of the annotated h5ad (raw counts, CSR, cells x genes).
   Doublets arise within a sample, so Scrublet is run separately for each Sample_GEO_Prefix.
Output: data/Hammond2019/output/scrublet_doublets.csv (cell_id, sample, doublet_score, predicted_doublet)
"""
import h5py, numpy as np
from scipy.sparse import csr_matrix
import scrublet as scr

H5 = "data/Hammond2019/raw/Hammond2019_microglia_AllSamples_annotated.h5ad"
OUT = "data/Hammond2019/output/scrublet_doublets.csv"

def read_cat(f, path):
    obj = f[path]
    if isinstance(obj, h5py.Group):  # categorical: categories + codes
        cats = obj["categories"][:].astype(str); codes = obj["codes"][:]
        return np.array([cats[c] if c >= 0 else "NA" for c in codes])
    v = obj[:]
    return v.astype(str) if v.dtype == object or v.dtype.kind in "SU" else v

f = h5py.File(H5, "r")
shape = f["layers/counts"].attrs["shape"]
X = csr_matrix((f["layers/counts/data"][:], f["layers/counts/indices"][:],
                f["layers/counts/indptr"][:]), shape=tuple(shape))  # cells x genes
cell_id = f["obs/_index"][:].astype(str)
sample = read_cat(f, "obs/Sample_GEO_Prefix")
f.close()
print(f"loaded {X.shape[0]} cells x {X.shape[1]} genes; {len(np.unique(sample))} samples")

scores = np.full(X.shape[0], np.nan)
calls  = np.zeros(X.shape[0], dtype=bool)
for s in np.unique(sample):
    idx = np.where(sample == s)[0]
    sub = X[idx]
    try:
        sc = scr.Scrublet(sub, expected_doublet_rate=0.06, random_state=0)
        ds, pred = sc.scrub_doublets(min_counts=2, min_cells=3,
                                     min_gene_variability_pctl=85, n_prin_comps=30,
                                     verbose=False)
        if pred is None:  # auto threshold failed -> fallback 0.25
            pred = ds > 0.25
        scores[idx] = ds; calls[idx] = pred
    except Exception as e:
        print(f"  [{s}] n={len(idx)} FAILED: {e}")
        calls[idx] = False
    print(f"  [{s}] n={len(idx):5d}  doublets={int(np.nansum(calls[idx])):4d} "
          f"({100*np.nanmean(calls[idx]):.2f}%)")

import csv
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["cell_id", "sample", "doublet_score", "predicted_doublet"])
    for i in range(len(cell_id)):
        w.writerow([cell_id[i], sample[i], f"{scores[i]:.4f}", int(calls[i])])

print(f"\n=== TOTAL doublets: {int(calls.sum())} / {len(calls)} "
      f"({100*calls.mean():.2f}%) ===")
print("wrote", OUT)
