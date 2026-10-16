## De Donno step 4 — per-cell RFind scores on the 4 axes (CORT / KA6h / UPR / NFkB).
## 4-term Running Fisher (analysis/running_fisher.R: smallest p over every top-G x top-K cut-off pair), clip +/-300.
## Background: the universe of each target-reference pair = De Donno gene universe ∩ genes the reference
## could contain (analysis/running_fisher.R: pair_universe); the per-cell lists are restricted to it before the top 500 are taken (cut_in_universe).
## Input : data/DeDonno2026/output/{percell_signature_celltype_deep.rds, references_dedonno.rds, cell_meta_qc.csv}
## Output: data/DeDonno2026/output/scores_dedonno.csv  (cell meta + UMAP + 4 axes)

suppressPackageStartupMessages({ library(data.table); library(parallel) })
OUT <- "data/DeDonno2026/output"
source("analysis/running_fisher.R")
CLIP <- 300; NTOP <- 500L; NCORES <- max(1L, detectCores() - 2L)

sig <- readRDS(file.path(OUT, "percell_signature_celltype_deep.rds"))
stopifnot(identical(sig$baseline, "per_celltype_all_cells"))
genes <- sig$genes; glow <- tolower(genes)
n <- ncol(sig$up); up <- sig$up; dn <- sig$down; NDEEP <- sig$n_deep
refs <- readRDS(file.path(OUT, "references_dedonno.rds")); axes <- names(refs)
## pair universe per reference; gene sets (type B of the paper, stored with type = "C") have eligible = NULL and no order (ranked = FALSE)
UA  <- lapply(refs, function(e) pair_universe(genes, e$eligible))
INU <- lapply(UA, function(u) glow %in% u)
P   <- vapply(UA, length, integer(1))
b2 <- lapply(axes, function(a) { e <- refs[[a]]
  list(up = intersect(tolower(e$up), UA[[a]]), down = intersect(tolower(e$down), UA[[a]]), ranked = !identical(e$type, "C")) })
names(b2) <- axes
cat(sprintf("%d cells; target universe %d genes\n", n, length(unique(glow))))
for (a in axes) cat(sprintf("  [%-4s] pair universe P=%d  up=%d down=%d ranked=%s\n", a, P[[a]], length(b2[[a]]$up), length(b2[[a]]$down), b2[[a]]$ranked))

sc <- function(idx){ m <- matrix(0, length(idx), length(axes), dimnames = list(NULL, axes)); short <- 0L
  for (r in seq_along(idx)){ j <- idx[r]
    for (a in axes) {
      uj <- cut_in_universe(up[, j], INU[[a]], NTOP); dj <- cut_in_universe(dn[, j], INU[[a]], NTOP)
      ## a list shorter than NTOP is fine when the cell has no more genes of that sign; it is a truncation only when the
      ## deep ranking ran out first (n_pos / n_neg > NDEEP)
      short <- short + (length(uj) < NTOP && sig$n_pos[j] > NDEEP) + (length(dj) < NTOP && sig$n_neg[j] > NDEEP)
      m[r, a] <- rfind_score_cell(glow[uj], glow[dj], b2[[a]], universe_size = P[[a]], clip = CLIP)
    } }
  attr(m, "short") <- short; m }
t0 <- Sys.time()
RES <- mclapply(split(seq_len(n), ceiling(seq_len(n)/250)), sc, mc.cores = NCORES)
n_short <- sum(vapply(RES, function(x) attr(x, "short"), integer(1)))
S <- do.call(rbind, RES)
cat(sprintf("lists cut short by the deep ranking (NDEEP=%d): %d  (must be 0)\n", NDEEP, n_short)); stopifnot(n_short == 0L)
cat(sprintf("scored in %.0fs (%d cores)\n", as.numeric(Sys.time() - t0, units = "secs"), NCORES))

meta <- fread(file.path(OUT, "cell_meta_qc.csv"))
stopifnot(identical(meta$cell_id, sig$cells))
out <- cbind(meta, as.data.table(S))
fwrite(out, file.path(OUT, "scores_dedonno.csv"))
cat(sprintf("wrote scores_dedonno.csv (%d cells x %d cols)\n", nrow(out), ncol(out)))
print(out[line == "Ctrl", lapply(.SD, function(x) round(median(x), 2)), by = .(cell_type, condition), .SDcols = axes][order(cell_type, condition)])
