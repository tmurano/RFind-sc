## Hammond step 5 — per-cell Delta lists (P2 = Pearson residual with +1)
## Delta(cell_i, gene_j) = (obs - expc)/√(expc+1)
##   obs  = raw UMI counts (gene×cell, gene universe; no ambient-RNA correction)
##   expc = share[i] * cs[j]   (expected count)
##   share = rowMeans(cf %*% diag(1/cs)),  cs = total UMI over all genes (raw)
## Output:
##   delta_P2_basis.rds          : list(genes, share, cs) = basis to reconstruct the Delta of any cell
##   percell_signature_top500.rds: per-cell UP/DOWN gene indices, up to 500 each (the b1 input of RFind)
##   example_delta_cell{1,2,3}.csv: full Delta rank tables of three cells, for checking
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
OUTP <- "data/Hammond2019/output"; NCELL <- 500L; NDEEP <- 1500L

## DecontX not used (see 04) -> the raw QC'd counts_qc are used as obs.
all_g <- readRDS(file.path(OUTP, "counts_qc.rds"))         # all genes, raw counts after QC (gene×cell)
cf    <- readRDS(file.path(OUTP, "counts_filtered.rds"))   # gene universe (after the relSE 5% filter) (gene×cell)
cs <- Matrix::colSums(all_g)                               # per-cell total UMI (all genes)
stopifnot(identical(colnames(cf), colnames(all_g)))
genes <- rownames(cf); G <- nrow(cf); nC <- ncol(cf)

## share[i] = mean over cells of (count / total UMI of that cell) = the gene's average share of a cell's UMIs.
## Not scaled to CPM: a factor 1e6 would cancel exactly and only make the step look like a normalisation.
## The expected count is share[i] x N[j] — see the Delta definition in the header.
share <- Matrix::rowMeans(cf %*% Diagonal(x = 1 / cs))
## cf stays sparse; columns are densified in chunks (a full dense matrix would be 10259×75572×8 B ≈ 6 GB)

cat(sprintf("computing Delta P2 for %d genes × %d cells (dense chunks)...\n", G, nC))
## per-cell signature: UP = the top min(NCELL, n) genes with residual > 0 / DOWN = the bottom min(NCELL, n) genes with residual < 0.
## If fewer than NCELL genes qualify, all of them are used (genes of the opposite sign are never added); empty slots are NA.
## Running Fisher is insensitive to low ranks, so the exact number is not critical: capped at 500, all genes when fewer.
upI <- matrix(NA_integer_, NCELL, nC); dnI <- matrix(NA_integer_, NCELL, nC)
## Deep per-cell ranking: the lists are cut to 500 only after they are restricted to the universe of each
## target-reference pair (analysis/running_fisher.R: pair_universe / cut_in_universe), so the ranking is kept to depth NDEEP. n_pos / n_neg
## record how many genes of each sign a cell has, so a truncated deep list can be detected downstream.
upD <- matrix(NA_integer_, NDEEP, nC); dnD <- matrix(NA_integer_, NDEEP, nC); npos <- integer(nC); nneg <- integer(nC)
expc_g <- share                                             # × cs[j] gives the expected count
CH <- 4000L
for (st in seq(1L, nC, by=CH)) {
  en <- min(st+CH-1L, nC); idx <- st:en
  obs  <- as.matrix(cf[, idx])                             # G×chunk dense
  expc <- outer(as.numeric(expc_g), cs[idx])               # G×chunk expected counts
  P2   <- (obs - expc) / sqrt(expc + 1)                    # P2 residuals
  for (k in seq_along(idx)) {
    p <- P2[, k]
    pos <- which(p > 0); pos <- pos[order(p[pos], decreasing = TRUE)]   # positive residuals = candidate positive markers
    neg <- which(p < 0); neg <- neg[order(p[neg], decreasing = FALSE)]  # most negative first = negative markers
    nu <- min(NCELL, length(pos)); nd <- min(NCELL, length(neg))
    if (nu) upI[seq_len(nu), idx[k]] <- pos[seq_len(nu)]
    if (nd) dnI[seq_len(nd), idx[k]] <- neg[seq_len(nd)]
    npos[idx[k]] <- length(pos); nneg[idx[k]] <- length(neg)
    nu2 <- min(NDEEP, length(pos)); nd2 <- min(NDEEP, length(neg))
    if (nu2) upD[seq_len(nu2), idx[k]] <- pos[seq_len(nu2)]
    if (nd2) dnD[seq_len(nd2), idx[k]] <- neg[seq_len(nd2)]
  }
}
cat("done.\n")

## `share` = the baseline of the Delta (CPM = share * 1e6); nothing else in the repository reads this file.
saveRDS(list(genes=genes, share=share, cs=cs), file.path(OUTP, "delta_P2_basis.rds"))
saveRDS(list(genes=genes, up=upI, down=dnI, n_cell=NCELL, cells=colnames(cf)),
        file.path(OUTP, "percell_signature_top500.rds"))
saveRDS(list(genes=genes, up=upD, down=dnD, n_deep=NDEEP, n_pos=npos, n_neg=nneg, cells=colnames(cf)),
        file.path(OUTP, "percell_signature_deep.rds"))

## check: full Delta rank tables of three cells + top genes of cell 1
for (j in 1:3) {
  cj <- as.numeric(cf[, j]); expc <- expc_g * cs[j]; delta <- (cj - expc) / sqrt(expc + 1)
  dt <- data.table(gene=genes, obs=cj, expc=round(expc,3), delta=round(delta,4))[order(-delta)]
  fwrite(dt, file.path(OUTP, sprintf("example_delta_cell%d.csv", j)))
}
e1 <- expc_g*cs[1]; d1 <- (as.numeric(cf[,1])-e1)/sqrt(e1+1)
cat(sprintf("\n=== cell1 (cluster=%s) Delta P2 sanity check ===\n", colnames(cf)[1]))
cat("  range:", sprintf("[%.2f, %.2f]", min(d1), max(d1)), "\n")
cat("  top10 UP:", paste(genes[order(d1,decreasing=TRUE)[1:10]], collapse=", "), "\n")
cat("  top10 DOWN:", paste(genes[order(d1)[1:10]], collapse=", "), "\n")
cat("\nwrote delta_P2_basis.rds + percell_signature_top500.rds + example_delta_cell{1,2,3}.csv\n")
