## Supplementary Fig. S2 — choice of the Delta statistic: the data of panels b, c and e, in one pass over the cells.
## Six candidate per-cell statistics for gene i in cell j (obs = UMI count, cpmj = CPM, expc = baseline CPM x total UMI / 1e6):
##   simpleFC = cpmj / baseline ; log2diff = log2(cpmj + 1) - log2(baseline + 1) ; zscore (see below) ;
##   P1 = (obs - expc)/sqrt(expc) ; P2 = (obs - expc)/sqrt(expc + 1) (the Delta used by RFind-sc) ;
##   P3_NB = (obs - expc)/sqrt(expc + expc^2/theta) (theta per gene by the method of moments on raw counts).
##   baseline = mean CPM over all cells.
## Cells: up to PER = 1,000 per author cluster (set.seed(1), clusters in sorted order). The comparison with the clusters uses the
## 13 microglial clusters, the table all 14 (Mono/Mac sorts last, so the cells drawn for the other clusters are the same).
## 1. Comparison with the clusters of Hammond et al. (panel e; no Running Fisher involved): for each of the 13 microglial author
##    clusters, positive markers = top N genes by one-vs-rest log2((cluster mean CPM + 1)/(rest + 1)) among genes detected in
##    >= 20% of the cluster; negative markers = bottom N among genes detected in >= 20% of the rest (Mono/Mac cells stay in the
##    baseline and in the "rest" of each cluster). Per cell, marker_top = share of the cell's top-N genes that are positive
##    markers; negmarker_bot = share of its bottom-N genes that are negative markers; composite = sum.
## 2. Check of the table of the candidates (panel b: what each fixes, what remains), per statistic over the cells of all 14
##    clusters: min / max of the statistic (divergence, asymmetry); among the top-N genes of a cell: share with a single molecule
##    (x = 1), mean x, mean expected count e and median baseline CPM; among the bottom-N genes: share undetected (x = 0) and the
##    number of genes tied with the N-th-lowest value (a broken DOWN order); rank agreement of simpleFC and log2diff within a
##    cell (Spearman over all genes; overlap of their top-N and bottom-N sets).
## 3. The six statistics of one cell over all genes of the universe (panel c, histograms; the first cell of the atlas), with four
##    kinds of genes: UP = detected, share at or above the mean share; DOWN-1 = detected, below the mean; DOWN-2 = not detected
##    and expected count < ln 2 (a zero is more likely than not, so missing it says little); DOWN-3 = not detected although
##    expected count >= ln 2 (truly low in this cell).
## N = 500 (the list length of the pipeline). The cells are read in blocks of 500 as dense columns.
## Input : data/Hammond2019/output/{counts_filtered.rds, counts_qc.rds}, data/Hammond2019/raw/cell_meta_hammond.csv
## Output: figs/Supp_S2_Delta/output/data/{figS2_delta_bench.csv, figS2_delta_bench_percluster.csv, figS2_table_check.csv,
##         figS2b_one_cell.csv} (drawn by figS2_panels.py)
## Run   : Rscript figs/Supp_S2_Delta/analysis/figS2_delta_benchmark.R   (from project root; a few minutes)

suppressPackageStartupMessages({ library(Matrix); library(data.table) })
SRC <- "data/Hammond2019/output"
OUT <- "figs/Supp_S2_Delta/output/data"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
N <- 500L; PER <- as.integer(Sys.getenv("S1_PER", "1000")); BLK <- 500L; set.seed(1)
t0 <- Sys.time()

cf    <- readRDS(file.path(SRC, "counts_filtered.rds"))       # 8,960 genes x 75,572 cells (gene universe)
all_g <- readRDS(file.path(SRC, "counts_qc.rds"))             # all genes, for the total UMI per cell
meta  <- fread("data/Hammond2019/raw/cell_meta_hammond.csv"); setkey(meta, cell_id)
cl    <- meta[colnames(cf)]$cell_type_paper
gn <- rownames(cf); G <- nrow(cf); cs <- Matrix::colSums(all_g); rm(all_g)
## CPM is used HERE: simpleFC and log2diff are defined in CPM, so they need it. For the Delta itself the
## factor cancels (baseline/1e6 = the gene's mean share), so the pipeline computes the Delta from the shares directly.
cpm <- cf %*% Diagonal(x = 1e6 / cs); baseline <- Matrix::rowMeans(cpm); lbase <- log2(baseline + 1)
.m <- Matrix::rowMeans(cf); .v <- (Matrix::rowMeans(cf * cf) - .m^2) * ncol(cf) / (ncol(cf) - 1)
theta <- ifelse(.v > .m, .m^2 / (.v - .m), Inf)
## z-score: the log-normalised value ln(1 + CP10k), centred and
## scaled per gene over all cells, values above 10 set to 10 = Seurat NormalizeData + ScaleData with their defaults (ScaleData clips
## only above; the per-cell value of the scFOCAL-type method of the benchmark).
lnm <- cpm; lnm@x <- log1p(lnm@x / 100)
zmu <- Matrix::rowMeans(lnm); zsd <- sqrt((Matrix::rowMeans(lnm * lnm) - zmu^2) * ncol(lnm) / (ncol(lnm) - 1)); rm(lnm)
stopifnot(all(zsd > 0))

FC <- list(
  simpleFC = function(obs, cpmj, expc) cpmj / baseline,
  log2diff = function(obs, cpmj, expc) log2(cpmj + 1) - lbase,
  zscore   = function(obs, cpmj, expc) pmin((log1p(cpmj / 100) - zmu) / zsd, 10),
  P1       = function(obs, cpmj, expc) (obs - expc) / sqrt(expc),
  P2       = function(obs, cpmj, expc) (obs - expc) / sqrt(expc + 1),
  P3_NB    = function(obs, cpmj, expc) (obs - expc) / sqrt(expc + expc^2 / theta))
models <- names(FC); M <- length(models)

## ---- cells: up to PER per author cluster ----
clusters_all <- sort(unique(cl))
draw <- setNames(lapply(clusters_all, function(k) { ik <- which(cl == k); if (length(ik) > PER) sample(ik, PER) else ik }),
                 clusters_all)
clusters <- setdiff(clusters_all, "Mono/Mac")                  # the 13 microglial clusters, as in Fig. 2c, d
cells <- unlist(draw, use.names = FALSE); kc <- rep(names(draw), lengths(draw)); nc <- length(cells)
cat(sprintf("%d cells drawn (%d in the 13 microglial clusters)\n", nc, sum(kc %in% clusters)))

## ---- 1. positive and negative markers of each microglial cluster ----
MK <- list()
for (k in clusters) {
  ik <- which(cl == k); mk <- Matrix::rowMeans(cpm[, ik]); mr <- Matrix::rowMeans(cpm[, -ik])
  dk <- Matrix::rowMeans(cpm[, ik] > 0); dr <- Matrix::rowMeans(cpm[, -ik] > 0)
  spec <- log2((mk + 1) / (mr + 1))
  sp <- spec; sp[dk < 0.20] <- -Inf; POS <- gn %in% gn[order(sp, decreasing = TRUE)][1:N]
  sn <- spec; sn[dr < 0.20] <-  Inf; NEG <- gn %in% gn[order(sn, decreasing = FALSE)][1:N]
  MK[[k]] <- list(POS = POS, NEG = NEG)
}
cat(sprintf("markers of %d clusters (%.0f s)\n", length(MK), as.numeric(Sys.time() - t0, units = "secs")))

## ---- one pass over the drawn cells: the six statistics, the top-N and bottom-N genes of each ----
SH_POS <- SH_NEG <- matrix(NA_real_, nc, M, dimnames = list(NULL, models))   # marker shares (part 1)
nr <- nc * M; r_model <- character(nr); r_cell <- integer(nr); r_min <- r_max <- r_upx1 <- r_upx <- r_upe <- r_upb <- r_dnx0 <-
  numeric(nr); r_tied <- integer(nr)                                          # per cell and statistic (part 2)
A <- data.table(cell = cells, spearman = NA_real_, top_overlap = NA_real_, bottom_overlap = NA_real_)
r <- 0L
for (b0 in seq(1L, nc, by = BLK)) {
  bi <- b0:min(nc, b0 + BLK - 1L); O <- as.matrix(cf[, cells[bi]]); C <- as.matrix(cpm[, cells[bi]])
  for (q in seq_along(bi)) {
    ci <- bi[q]; j <- cells[ci]; o <- O[, q]; cj <- C[, q]; e <- baseline / 1e6 * cs[j]; mk <- MK[[kc[ci]]]   # NULL for Mono/Mac
    UP <- DN <- S <- setNames(vector("list", M), models)
    for (m in models) {
      s <- FC[[m]](o, cj, e); up <- order(s, decreasing = TRUE)[1:N]; dn <- order(s, decreasing = FALSE)[1:N]
      UP[[m]] <- up; DN[[m]] <- dn; S[[m]] <- s; r <- r + 1L
      r_model[r] <- m; r_cell[r] <- j; r_min[r] <- min(s); r_max[r] <- max(s); r_upx1[r] <- mean(o[up] == 1)
      r_upx[r] <- mean(o[up]); r_upe[r] <- mean(e[up]); r_upb[r] <- median(baseline[up]); r_dnx0[r] <- mean(o[dn] == 0)
      r_tied[r] <- sum(s == s[dn[N]])
      if (!is.null(mk)) { SH_POS[ci, m] <- sum(mk$POS[up]) / N; SH_NEG[ci, m] <- sum(mk$NEG[dn]) / N }
    }
    set(A, ci, c("spearman", "top_overlap", "bottom_overlap"),
        list(cor(S$simpleFC, S$log2diff, method = "spearman"), length(intersect(UP$simpleFC, UP$log2diff)) / N,
             length(intersect(DN$simpleFC, DN$log2diff)) / N))
  }
  cat(sprintf("cells %d-%d done (%.0f s)\n", b0, max(bi), as.numeric(Sys.time() - t0, units = "secs")))
}

## ---- 1. comparison with the clusters (panel e) ----
rows <- list()
for (k in clusters) {
  ix <- which(kc == k)                                             # the cells of the cluster, in the order drawn
  for (m in models) {
    pos <- 0; neg <- 0; for (i in ix) { pos <- pos + SH_POS[i, m]; neg <- neg + SH_NEG[i, m] }
    rows[[length(rows) + 1]] <- data.table(cluster = k, model = m, n = length(ix), marker_top = pos / length(ix),
                                           negmarker_bot = neg / length(ix))
  }
}
pc <- rbindlist(rows); pc[, composite := marker_top + negmarker_bot]
fwrite(pc, file.path(OUT, "figS2_delta_bench_percluster.csv"))
tab <- pc[, .(marker_top = round(mean(marker_top), 3), negmarker_bot = round(mean(negmarker_bot), 3),
              composite = round(mean(composite), 3)), by = model][order(-composite)]
fwrite(tab, file.path(OUT, "figS2_delta_bench.csv"))
cat(sprintf("=== Delta benchmark (%d cells x %d genes, top/bottom %d, <= %d cells per cluster) ===\n", ncol(cf), G, N, PER)); print(tab)

## ---- 2. check of the table of the candidates (panel b) ----
R <- data.table(model = r_model, cell = r_cell, min = r_min, max = r_max, up_x1 = r_upx1, up_x = r_upx, up_e = r_upe,
                up_base = r_upb, dn_x0 = r_dnx0, dn_tied = r_tied)
tab <- R[, .(min = round(min(min), 2), max = round(max(max), 1), median_cell_max = round(median(max), 1),
             top_single_molecule = round(mean(up_x1), 3), top_mean_x = round(mean(up_x), 2), top_mean_e = round(mean(up_e), 3),
             top_median_baseline_cpm = round(median(up_base), 1), bottom_undetected = round(mean(dn_x0), 3), list_length = N,
             bottom_tied_median = median(dn_tied)), by = model]
fwrite(tab, file.path(OUT, "figS2_table_check.csv"))
print(tab)
cat(sprintf("\nsimpleFC vs log2diff within a cell: Spearman median %.3f [%.3f-%.3f]; top-%d overlap %.2f; bottom-%d overlap %.2f\n",
            median(A$spearman), min(A$spearman), max(A$spearman), N, mean(A$top_overlap), N, mean(A$bottom_overlap)))

## ---- 3. one cell, all genes (panel c, histograms) ----
j <- 1L; o <- as.numeric(cf[, j]); cj <- as.numeric(cpm[, j]); e <- baseline / 1e6 * cs[j]
C1 <- data.table(gene = rownames(cf), x = o, e = e, as.data.table(lapply(FC, function(f) f(o, cj, e))))
C1[, kind := fifelse(x > 0 & cj >= baseline, "UP", fifelse(x > 0, "DOWN-1", fifelse(e < log(2), "DOWN-2", "DOWN-3")))]
fwrite(C1, file.path(OUT, "figS2b_one_cell.csv"))
cat(sprintf("\none cell = %s (total UMI %d): %d genes, %d detected\n", colnames(cf)[j], as.integer(cs[j]), nrow(C1), sum(o > 0)))
print(C1[, .N, by = kind][order(kind)])
print(C1[, lapply(.SD, function(v) round(c(min = min(v), max = max(v)), 2)), .SDcols = names(FC)])
cat(sprintf("\nwrote %s / figS2_delta_bench{,_percluster}.csv, figS2_table_check.csv, figS2b_one_cell.csv (%.1f min)\n",
            OUT, as.numeric(Sys.time() - t0, units = "mins")))
