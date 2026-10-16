## RFind-sc demo with the R package RFindsc: scoring your own cells against your own reference.
## The same 1,000 microglia and DAM reference as demo/run_demo.R, scored with the package's public function score_panel_pairs().
## The package takes the expected counts from the cells it is given (here these 1,000 cells; the paper took them from all 75,572
## cells of the atlas) and keeps the genes detected in at least 1% of these cells, so the values differ slightly from the paper's
## DAM-index; the script reports how closely they agree.
## Needs RFindsc 0.3.8 (README, "Installing the companion packages"). Run from the project root:
##   Rscript demo/run_demo_package.R          (set RFINDSC_LIB=<library> if RFindsc was installed into a separate library)
## Output: demo/output_package_R.csv (cell_id, condition, DAM_index)
suppressPackageStartupMessages(library(Matrix))
t_start <- proc.time()[["elapsed"]]
lib <- Sys.getenv("RFINDSC_LIB")
suppressPackageStartupMessages(if (nzchar(lib)) library(RFindsc, lib.loc = lib) else library(RFindsc))
D <- "demo/data"

## ---- your data: raw UMI counts, genes x cells, with gene symbols and cell names ----
X <- as(readMM(gzfile(file.path(D, "counts.mtx.gz"))), "CsparseMatrix")
genes <- read.delim(gzfile(file.path(D, "genes.tsv.gz")), colClasses = c("character", "numeric"))
cells <- read.delim(file.path(D, "labels.tsv"), colClasses = c("character", "character", "numeric"))
dimnames(X) <- list(genes$gene, cells$cell_id)

## ---- your reference: up- and down-regulated genes, most changed first ----
ref <- read.delim(file.path(D, "dam_reference.tsv"), colClasses = c("character", "integer", "character"))
panels <- list(DAM = list(up = ref$gene[ref$direction == "up"], down = ref$gene[ref$direction == "down"]))

## ---- one call; libsize = each cell's total UMI over all genes, because this matrix holds only a subset of the genes ----
out <- score_panel_pairs(X, panels, config = rfindsc_config(n_cores = 1), libsize = cells$total_umi)
score <- as.numeric(out$scores[cells$cell_id, "DAM"])
cat(sprintf("RFindsc %s: %d cells scored against the DAM reference (scoring universe %d genes)\n",
            as.character(packageVersion("RFindsc")), length(score), out$universe_size[["DAM"]]))

auc <- function(s, pos) { r <- rank(s); n1 <- sum(pos); n0 <- sum(!pos); (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
lpc <- cells$condition == "LPC"
write.csv(data.frame(cell_id = cells$cell_id, condition = cells$condition, DAM_index = sprintf("%.17g", score)),
          "demo/output_package_R.csv", row.names = FALSE, quote = FALSE)
t_run <- proc.time()[["elapsed"]] - t_start

cat("\nfirst cells:\n"); print(head(data.frame(cell_id = cells$cell_id, condition = cells$condition, DAM_index = round(score, 4)), 5),
                                row.names = FALSE)
cat(sprintf("\nmedian DAM-index: LPC %.3f, saline %.3f\n", median(score[lpc]), median(score[!lpc])))
cat(sprintf("AUC (LPC vs saline): %.4f\n", auc(score, lpc)))
ex <- read.csv("demo/expected_output.csv", colClasses = c("character", "character", "numeric"))
stopifnot(identical(ex$cell_id, cells$cell_id))
cat(sprintf("agreement with the paper's DAM-index of these cells (baseline from all cells of the atlas): Spearman %.4f, AUC there %.4f\n",
            cor(score, ex$DAM_index, method = "spearman"), auc(ex$DAM_index, lpc)))
cat(sprintf("written demo/output_package_R.csv; run time %.1f s (%s)\n", t_run, R.version.string))
