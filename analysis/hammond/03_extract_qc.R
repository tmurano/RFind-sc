## Hammond step 3 — raw-count extraction + QC subset
## (reads the raw counts from the h5ad, restores dashed gene names, drops all-zero genes, applies the flags of 01 + 02)
## Input:
##   annotated h5ad (layers/counts = CSR, cells x genes)
##   qc_flags_mt_umi.csv from 01 (mt/UMI flags) + scrublet_doublets.csv from 02 (doublets)
## Output:
##   counts_hammond.rds : raw genes x cells (after removing all-zero genes = 19,785 x 76,149)
##   counts_qc.rds      : raw genes x cells after removing the mt/UMI/doublet cells (input of 04, 05 and 08)
## Note: a CSR matrix (cells x genes) has the layout of a CSC matrix (genes x cells), so the dgCMatrix is built without a transpose.
suppressPackageStartupMessages({ library(hdf5r); library(Matrix); library(data.table) })
H5  <- "data/Hammond2019/raw/Hammond2019_microglia_AllSamples_annotated.h5ad"
OUTP<- "data/Hammond2019/output"

f <- H5File$new(H5, "r")
shape <- f[["layers/counts"]]$attr_open("shape")$read()   # [76149 cells, 20579 genes]
nC <- shape[1]; nG <- shape[2]
dat <- f[["layers/counts/data"]]$read()
idx <- f[["layers/counts/indices"]]$read()                # 0-based gene index
ptr <- f[["layers/counts/indptr"]]$read()                 # length nC+1 (per cell)
genes <- f[["var/_index"]]$read(); cells <- f[["obs/_index"]]$read()
f$close_all()
cat(sprintf("h5ad: %d cells x %d genes, nnz=%d\n", nC, nG, length(dat)))
## Gene names: the annotated h5ad stores "-" as "." (e.g. H2.D1 for H2-D1), so these genes would not match the references.
## Restore "-" where the dashed form is an MGI symbol and the dotted form is not (org.Mm.eg.db).
sym  <- AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype = "SYMBOL")
dash <- gsub(".", "-", genes, fixed = TRUE)
fix  <- grepl(".", genes, fixed = TRUE) & !(genes %in% sym) & (dash %in% sym)
cat(sprintf("gene names: %d of %d dotted names restored to '-' (e.g. %s)\n", sum(fix), sum(grepl(".", genes, fixed = TRUE)),
            paste(head(dash[fix], 8), collapse = ", ")))
genes[fix] <- dash[fix]
stopifnot(!anyDuplicated(genes))

## CSR(cell×gene) == CSC(gene×cell): i=gene idx, p=cell ptr → gene×cell dgCMatrix
M <- new("dgCMatrix", i = as.integer(idx), p = as.integer(ptr),
         x = as.numeric(dat), Dim = c(as.integer(nG), as.integer(nC)))
rownames(M) <- genes; colnames(M) <- cells
rm(dat, idx, ptr); gc()

## remove all-zero genes
nz <- Matrix::rowSums(M) > 0
M <- M[nz, ]
cat(sprintf("removed %d all-zero genes → %d genes × %d cells\n", sum(!nz), nrow(M), ncol(M)))
## (colSums equal nCount_RNA of the metadata; asserted when the counts layer was built, analysis/hammond/data_acquisition/run_pipeline.py)
saveRDS(M, file.path(OUTP, "counts_hammond.rds"))

## --- QC subset (apply the flags of 01 + 02) ---
qc  <- fread(file.path(OUTP, "qc_flags_mt_umi.csv"))[, .(cell_id, flag_mt, flag_umi_lo, flag_umi_hi)]
dbl <- fread(file.path(OUTP, "scrublet_doublets.csv"))[, .(cell_id, predicted_doublet)]
m <- merge(qc, dbl, by = "cell_id", all.x = TRUE)
setkey(m, cell_id); m <- m[colnames(M)]
keep <- !(m$flag_mt | m$flag_umi_lo | m$flag_umi_hi | as.logical(m$predicted_doublet))
keep[is.na(keep)] <- TRUE
cat(sprintf("QC: mt/UMI removed %d + doublets %d → keep %d / %d (%.2f%%)\n",
            sum(m$flag_mt|m$flag_umi_lo|m$flag_umi_hi, na.rm=TRUE),
            sum(as.logical(m$predicted_doublet), na.rm=TRUE),
            sum(keep), length(keep), 100*mean(keep)))
counts_qc <- M[, keep]
saveRDS(counts_qc, file.path(OUTP, "counts_qc.rds"))
saveRDS(colnames(counts_qc), file.path(OUTP, "cells_kept.rds"))
cat(sprintf("\nwrote counts_hammond.rds (%d×%d) + counts_qc.rds (%d×%d) + cells_kept.rds\n",
            nrow(M), ncol(M), nrow(counts_qc), ncol(counts_qc)))
