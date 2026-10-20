## Run time of scoring cells against one reference (Supplementary Table 3, sheet "Run time and memory").
## One run = one method, one number of cells, one reference length and one number of cores; fig5_runtime.sh runs the grid under
## /usr/bin/time -l (peak memory of the process) on an otherwise idle computer and collects the results.
## Data: the 75,572 cells of the data of Hammond et al. after quality control (counts_qc.rds; raw UMI counts of all genes); the first
##   n cells of a fixed random order (seed 1). Reference: DAM (Keren-Shaul et al.), the top L genes per direction of the published list
##   (L = 500, as in the paper, or all = the whole list), the genes present in the counts.
## Methods: RFind-sc (R package RFindsc, score_panel_pairs(); its own gene filter and baseline from the cells given), ULM of decoupleR
##   with the log2 fold changes as weights and UCell with the signed signature (up genes and down genes with "-"), with the input of
##   fig5_bench_redo.R (log1p of counts per 10,000 in blocks of 3,000 cells for ULM; raw counts for UCell).
## Output: one line appended to figs/Fig5_benchmark/output/results/fig5_runtime_raw.tsv (method, cells, L, cores, scoring seconds,
##   genes up / down given; the time covers scoring only, not reading the counts).
## Run (from the project root): R --no-save --no-restore -f figs/Fig5_benchmark/analysis/fig5_runtime_score.R --args method=rfindsc n=1000 L=500 cores=1
pdf(NULL)
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
kv <- function(k, d) { a <- grep(paste0("^", k, "="), args, value = TRUE); if (length(a)) sub("^[^=]*=", "", a[1]) else d }
METHOD <- kv("method", "rfindsc"); N <- kv("n", "1000"); LARG <- kv("L", "500"); CORES <- as.integer(kv("cores", "1"))
stopifnot(METHOD %in% c("rfindsc", "ulm", "ucell"))
OUT <- "figs/Fig5_benchmark/output/results/fig5_runtime_raw.tsv"

X <- readRDS("data/Hammond2019/output/counts_qc.rds")
set.seed(1); ord <- sample(ncol(X))
n <- if (N == "all") ncol(X) else as.integer(N)
X <- X[, sort(ord[seq_len(n)]), drop = FALSE]
libsize <- Matrix::colSums(X)
tab <- fread("data/Hammond2019/output/references/dam_KerenShaul.csv")
LV <- if (LARG == "all") .Machine$integer.max else as.integer(LARG)
up <- head(tab[`Fold Change` > 0][order(-`Fold Change`)]$Symbol, LV); dn <- head(tab[`Fold Change` < 0][order(`Fold Change`)]$Symbol, LV)
up <- up[up %in% rownames(X)]; dn <- dn[dn %in% rownames(X)]
cat(sprintf("%s: %d cells x %d genes; DAM L = %s: up %d, down %d; %d core(s)\n", METHOD, ncol(X), nrow(X), LARG, length(up), length(dn), CORES))

t0 <- proc.time()[["elapsed"]]
if (METHOD == "rfindsc") {
  suppressPackageStartupMessages(library(RFindsc))
  out <- score_panel_pairs(X, list(DAM = list(up = up, down = dn)), config = rfindsc_config(n_cores = CORES, top_n = if (LARG == "all") Inf else LV),
                             libsize = libsize)   # top_n = Inf: the package would otherwise cut the whole list to 500 per side
  s <- as.numeric(out$scores[, "DAM"])
} else if (METHOD == "ulm") {
  fc <- setNames(tab$`Fold Change`, tab$Symbol)[c(up, dn)]
  w <- sign(fc) * log2(abs(fc)); keep <- w != 0
  net <- data.frame(source = "DAM", target = c(up, dn)[keep], mor = w[keep])
  s <- numeric(ncol(X))
  for (b in split(seq_len(ncol(X)), ceiling(seq_len(ncol(X)) / 3000))) {
    m <- log1p(sweep(as.matrix(X[, b, drop = FALSE]), 2, 1e4 / libsize[b], "*")); cl <- colnames(m)
    u <- decoupleR::run_ulm(m, net, minsize = 5); s[b] <- u$score[match(cl, u$condition)]
  }
} else {
  f <- c(up, paste0(dn, "-"))
  s <- as.numeric(UCell::ScoreSignatures_UCell(X, features = list(S = f), ncores = CORES, name = "",
                                               maxRank = max(1500L, length(f)))[colnames(X), "S"])
}
secs <- proc.time()[["elapsed"]] - t0
stopifnot(length(s) == ncol(X), sum(is.finite(s)) > 0)
cat(sprintf("scoring time %.1f s\n", secs))
new <- !file.exists(OUT)
cat(if (new) "method\tcells\tL\tcores\tseconds\tup\tdown\n" else "",
    sprintf("%s\t%d\t%s\t%d\t%.2f\t%d\t%d\n", METHOD, ncol(X), LARG, CORES, secs, length(up), length(dn)), sep = "", file = OUT, append = TRUE)
