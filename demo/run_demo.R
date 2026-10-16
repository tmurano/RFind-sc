## RFind-sc demo (R): the DAM-index of 1,000 microglia (Hammond et al. 2019; white matter after LPC vs saline injection, 500 each)
## with the DAM reference of Keren-Shaul et al. 2017 (the top 500 genes per direction, 408 up and 427 down in the data), as in Figs. 2 and 5 of the paper.
## Uses only the kernel of this repository (analysis/running_fisher.R) and base R + Matrix (no RFindsc package needed).
## Steps, for every cell (as analysis/hammond/05_delta_P2.R and 07_score_axes.R):
##   1. Delta = (x - e) / sqrt(e + 1), x = UMI count, e = baseline share of the gene x total UMI of the cell (all genes);
##   2. UP = the 500 genes with the largest positive Delta, DOWN = the 500 most negative;
##   3. S = RF(UP, up) - RF(UP, down) - RF(DOWN, up) + RF(DOWN, down) (Running Fisher; universe P = the genes of the demo).
## Then the AUC of LPC vs saline (Mann-Whitney), and the check against the paper's stored DAM-index (demo/expected_output.csv).
## Run from the project root:  Rscript demo/run_demo.R      (or R --no-save --no-restore -f demo/run_demo.R)
## Output: demo/output_R.csv (cell_id, condition, DAM_index)
suppressPackageStartupMessages(library(Matrix))
t_start <- proc.time()[["elapsed"]]
source("analysis/running_fisher.R")                     # running_fisher(), rfind_score_cell()
D <- "demo/data"; NTOP <- 500L; CLIP <- 300

## ---- data ----
X <- as(readMM(gzfile(file.path(D, "counts.mtx.gz"))), "CsparseMatrix")          # genes x cells, raw UMI counts
genes <- read.delim(gzfile(file.path(D, "genes.tsv.gz")), colClasses = c("character", "numeric"))
cells <- read.delim(file.path(D, "labels.tsv"), colClasses = c("character", "character", "numeric"))
stopifnot(identical(readLines(gzfile(file.path(D, "barcodes.tsv.gz"))), cells$cell_id), nrow(X) == nrow(genes), ncol(X) == nrow(cells))
ref <- read.delim(file.path(D, "dam_reference.tsv"), colClasses = c("character", "integer", "character"))
b2 <- list(up = ref$gene[ref$direction == "up"], down = ref$gene[ref$direction == "down"])  # ranked by fold change, lower case
glow <- tolower(genes$gene)
P <- length(unique(glow))                                 # universe of the pair cell x reference (pair_universe)
cat(sprintf("demo: %d genes (P = %d) x %d cells; reference %d up / %d down genes\n", nrow(X), P, ncol(X), length(b2$up), length(b2$down)))

## ---- per-cell Delta, UP / DOWN lists, signed Running Fisher score ----
score <- numeric(ncol(X))
for (j in seq_len(ncol(X))) {
  x <- as.numeric(X[, j])
  e <- genes$baseline_share * cells$total_umi[j]          # expected count
  d <- (x - e) / sqrt(e + 1)                              # Delta
  pos <- which(d > 0); pos <- pos[order(d[pos], decreasing = TRUE)]
  neg <- which(d < 0); neg <- neg[order(d[neg], decreasing = FALSE)]
  up <- glow[head(pos, NTOP)]; dn <- glow[head(neg, NTOP)]
  score[j] <- rfind_score_cell(up, dn, b2, universe_size = P, clip = CLIP)
}
auc <- function(s, pos) { r <- rank(s); n1 <- sum(pos); n0 <- sum(!pos); (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
A <- auc(score, cells$condition == "LPC")
out <- data.frame(cell_id = cells$cell_id, condition = cells$condition, DAM_index = sprintf("%.17g", score))
write.csv(out, "demo/output_R.csv", row.names = FALSE, quote = FALSE)
t_run <- proc.time()[["elapsed"]] - t_start

cat("\nfirst cells:\n"); print(head(data.frame(cell_id = cells$cell_id, condition = cells$condition, DAM_index = round(score, 4)), 5), row.names = FALSE)
cat(sprintf("\nmedian DAM-index: LPC %.3f, saline %.3f\n", median(score[cells$condition == "LPC"]), median(score[cells$condition == "saline"])))
cat(sprintf("AUC (LPC vs saline): %.4f\n", A))
if (file.exists("demo/expected_output.csv")) {
  ex <- read.csv("demo/expected_output.csv", colClasses = c("character", "character", "numeric"))
  stopifnot(identical(ex$cell_id, cells$cell_id))
  dmax <- max(abs(ex$DAM_index - score))
  cat(sprintf("max |difference| from the paper's stored DAM-index: %.3g -> %s\n", dmax, if (dmax < 1e-8) "OK" else "DIFFERS"))
}
cat(sprintf("written demo/output_R.csv; run time %.1f s (%s)\n", t_run, R.version.string))
