## Hammond step 4 — gene universe (hygiene filters)
## The same filters define one gene universe for both the per-cell Delta and the references (design decision).
##   (1) det >= 1%        : detected in >= 1% of the cells after QC, microglia and monocytes/macrophages (removes ambient / rare noise)
##   (2) relative SE < 5% : (sd/sqrt(n))/mean; stability of the FC/Delta estimate (~ total UMI >= 400)
##   (3) protein-coding   : Rik/Gm/mt-/clone-name genes removed (interpretability)
##   Sex genes are NOT removed (design decision): the reference contrasts are restricted to males, so sex is not a confounder
##     there, and in the per-cell Delta sex genes are kept as real expression.
##   (1) and (2) overlap by ~92%, but det also catches ambient genes with enough total UMI concentrated in <1% of cells (e.g. Hba-a2).
## Ambient RNA is not removed from individual cells (DecontX is not used; Methods, "Quality control and genes analysed").
## Input : counts_qc.rds (genes x cells, raw counts after mt/UMI/doublet removal)
## Output: genes_universe.rds (gene universe shared by the per-cell Delta and the references) + gene_universe.csv + counts_filtered.rds
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
OUTP <- "data/Hammond2019/output"
DET_MIN <- 0.01; RELSE_MAX <- 0.05

counts <- readRDS(file.path(OUTP, "counts_qc.rds"))   # gene×cell, raw counts after QC
n <- ncol(counts); G <- nrow(counts); gn <- rownames(counts)
sx  <- Matrix::rowSums(counts); sx2 <- Matrix::rowSums(counts*counts)
mu  <- sx / n; v <- (sx2 - sx^2/n) / (n-1)
relSE <- sqrt(pmax(v,0)) / sqrt(n) / mu
det <- Matrix::rowSums(counts > 0) / n
bad_pc <- grepl("Rik[0-9]?$",gn)|grepl("^Gm[0-9]",gn)|grepl("^mt-",gn)|grepl("^A[CFLI][0-9]",gn)|grepl("^[0-9]",gn)|grepl("Rik$",gn)
f_det <- det>=DET_MIN; f_se <- is.finite(relSE)&relSE<RELSE_MAX; f_pc <- !bad_pc

d <- data.table(gene=gn, total_UMI=sx, det=round(100*det,3), relSE=round(relSE,4),
                pass_det=f_det, pass_relSE=f_se, pass_pc=f_pc)
d[, keep := f_det & f_se & f_pc]
fwrite(d[order(-keep, relSE)], file.path(OUTP, "gene_universe.csv"))
genes_universe <- d[keep==TRUE]$gene
saveRDS(genes_universe, file.path(OUTP, "genes_universe.rds"))
saveRDS(counts[genes_universe, ], file.path(OUTP, "counts_filtered.rds"))

cat(sprintf("all %d genes → gene universe %d (%.1f%%) [sex genes kept]\n", G, length(genes_universe), 100*mean(d$keep)))
cat(sprintf("  passing each filter alone: det≥1%%=%d / relSE<5%%=%d / protein-coding=%d\n", sum(f_det), sum(f_se), sum(f_pc)))
cat(sprintf("  sequential: det(%d) +relSE(%d) +pc(%d)\n", sum(f_det), sum(f_det&f_se), sum(f_det&f_se&f_pc)))
cat(sprintf("  universe total UMI: min %.0f / median %.0f / min detection rate %.2f%%\n",
            min(sx[genes_universe]), median(sx[genes_universe]), 100*min(det[genes_universe])))
cat("\nwrote genes_universe.rds + gene_universe.csv + counts_filtered.rds\n")
