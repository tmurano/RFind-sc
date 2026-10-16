## De Donno step 2 — per-cell Delta(P2) signature with a PER-CELL-TYPE baseline.
## Rule:
##   each cell is compared with the mean CPM of ALL cells of its own cell type (both conditions; wild-type mice only, as every
##   cell of cell_meta_qc.csv, 01_qc_embed.R),
##   Delta = (obs - exp)/sqrt(exp + 1), exp = baseline CPM/1e6 x cell total UMI (all raw genes).
##   UP = top-500 positive residuals, DOWN = top-500 negative residuals.
## Why per-cell-type: De Donno has 21 author cell types in the wild-type cells (numbered subtypes merged), so a global baseline makes
##   cell identity dominate the signature; the baseline includes both conditions so the Stress-vs-No-Stress readout is not biased.
##
## Input : data/DeDonno2026/output/{counts_universe.rds, cell_meta_qc.csv}
## Output: data/DeDonno2026/output/percell_signature_celltype.rds

suppressPackageStartupMessages({ library(Matrix); library(data.table) })
OUT <- "data/DeDonno2026/output"; NCELL <- 500L; NDEEP <- 1500L; CH <- 4000L

X <- readRDS(file.path(OUT, "counts_universe.rds")); cf <- X$counts; cs <- X$cs; cells <- X$cells; rm(X); gc()
meta <- fread(file.path(OUT, "cell_meta_qc.csv"))
stopifnot(identical(meta$cell_id, cells))
uni <- rownames(cf); G <- nrow(cf); nC <- ncol(cf)
cat(sprintf("universe %d genes (%d unique), %d cells\n", G, length(unique(uni)), nC))
## shares, not CPM: the 1e6 of CPM would cancel against the /1e6 in the expected count (as in
## analysis/hammond/05_delta_P2.R)
shr <- cf %*% Diagonal(x = 1 / cs)

ct <- meta$cell_type; cts <- sort(unique(ct)); K <- length(cts); ctidx <- match(ct, cts)
baseline_ct <- matrix(0, G, K, dimnames = list(NULL, cts))
for (ti in seq_len(K)) { cols <- which(ctidx == ti); stopifnot(length(cols) >= 20)
  baseline_ct[, ti] <- Matrix::rowMeans(shr[, cols, drop = FALSE])
  cat(sprintf("  baseline %-18s %5d cells\n", cts[ti], length(cols))) }

upI <- matrix(NA_integer_, NCELL, nC); dnI <- matrix(NA_integer_, NCELL, nC)
## Deep per-cell ranking: the lists are cut to 500 only after they are restricted to the universe of each
## target-reference pair (analysis/running_fisher.R: pair_universe / cut_in_universe), so the ranking is kept to depth NDEEP. n_pos / n_neg
## record how many genes of each sign a cell has, so a truncated deep list can be detected downstream.
upD <- matrix(NA_integer_, NDEEP, nC); dnD <- matrix(NA_integer_, NDEEP, nC); npos <- integer(nC); nneg <- integer(nC)
t0 <- Sys.time()
for (st in seq(1L, nC, by = CH)) {
  en <- min(st + CH - 1L, nC); idx <- st:en
  obs  <- as.matrix(cf[, idx])
  expc <- sweep(baseline_ct[, ctidx[idx], drop = FALSE], 2, cs[idx], "*")
  P2   <- (obs - expc) / sqrt(expc + 1)
  for (k in seq_along(idx)) { p <- P2[, k]
    pos <- which(p > 0); pos <- pos[order(p[pos], decreasing = TRUE)]
    neg <- which(p < 0); neg <- neg[order(p[neg], decreasing = FALSE)]
    nu <- min(NCELL, length(pos)); nd <- min(NCELL, length(neg))
    if (nu) upI[seq_len(nu), idx[k]] <- pos[seq_len(nu)]
    if (nd) dnI[seq_len(nd), idx[k]] <- neg[seq_len(nd)]
    npos[idx[k]] <- length(pos); nneg[idx[k]] <- length(neg)
    nu2 <- min(NDEEP, length(pos)); nd2 <- min(NDEEP, length(neg))
    if (nu2) upD[seq_len(nu2), idx[k]] <- pos[seq_len(nu2)]
    if (nd2) dnD[seq_len(nd2), idx[k]] <- neg[seq_len(nd2)] }
}
cat(sprintf("built in %.0fs; NA in UP %d / DOWN %d\n", as.numeric(Sys.time() - t0, units = "secs"),
            sum(is.na(upI)), sum(is.na(dnI))))
saveRDS(list(genes = uni, up = upI, down = dnI, n_cell = NCELL, cells = cells,
             baseline = "per_celltype_all_cells", baseline_celltypes = cts),
        file.path(OUT, "percell_signature_celltype.rds"))
saveRDS(list(genes = uni, up = upD, down = dnD, n_deep = NDEEP, n_pos = npos, n_neg = nneg, cells = cells,
             baseline = "per_celltype_all_cells", baseline_celltypes = cts),
        file.path(OUT, "percell_signature_celltype_deep.rds"))
cat("wrote percell_signature_celltype.rds\n")
