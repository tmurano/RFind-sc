## Fig. 5 — the share of cells to which UCell with the signed input (up genes + down genes marked "-") gives 0 (legend and
## Source Data of Fig. 5, and the legend of Supplementary Fig. S4 for its four added type A tasks). UCell v2 scores a signed signature as
## max(0, U+ - w_neg * U-) (w_neg = 1; Andreatta and Carmona 2026, section 2.3): U+ = UCell score of the up genes, U- = UCell score of
## the down genes, each from the ranks of the genes in the cell (maxRank 1500). This script scores the signed signature, with U+ and U-
## to check that max(0, U+ - U-) equals it (the script stops otherwise), on the cells of every signed task at L = 500 (the genes and
## the input of fig5_bench_redo.R), and reports, per task and group, the number of cells and the share of cells with a score of 0.
## Output: figs/Fig5_benchmark/output/results/fig5_ucell_zero.csv. Run: R 4.5.2 from the project root (a few minutes).
pdf(NULL); options(warn = 1)
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
NC <- max(1L, parallel::detectCores() - 2L)
uc <- function(cnt, up, dn) {
  f <- list(up = up, dn = dn, signed = c(up, paste0(dn, "-")))
  s <- UCell::ScoreSignatures_UCell(cnt, features = f, ncores = NC, name = "")[colnames(cnt), c("up", "dn", "signed")]
  data.table(Uplus = s[, "up"], Uminus = s[, "dn"], signed = s[, "signed"])
}
summ <- function(id, x, grp, nu, nd) {
  d <- max(abs(pmax(0, x$Uplus - x$Uminus) - x$signed))
  if (!isTRUE(d < 1e-12)) stop(sprintf("%s: max(0, U+ - U-) differs from the signed UCell score (max |diff| %.3g)", id, d), call. = FALSE)
  x[, group := grp][, .(id = id, n_cells = .N, up_genes = nu, down_genes = nd, share_score0 = mean(signed == 0)), by = group]
}
OUT <- list()
## Hammond (DAM: H1; Maturity: H4)
sig <- readRDS("data/Hammond2019/output/percell_signature_top500.rds"); cells <- sig$cells; univ <- sig$genes; rm(sig)
meta <- fread("data/Hammond2019/raw/cell_meta_hammond.csv")[match(cells, cell_id)]
panel <- readRDS("data/Hammond2019/output/references/references_panel9.rds")
cf <- readRDS("data/Hammond2019/output/counts_filtered.rds"); stopifnot(identical(colnames(cf), cells))
mg <- meta$paper_cluster != "Mono/Mac"
TH <- list(H1 = list(key = "DAM", ix = which(mg & meta$age %in% c("P100_Injured_WM", "P100_Control_WM"))),
           H4 = list(key = "Maturity", ix = which(mg & meta$age %in% c("P100", "E14", "P4/P5"))))
for (id in names(TH)) {
  k <- TH[[id]]$key; up <- intersect(panel[[k]]$up, univ); dn <- intersect(panel[[k]]$down, univ); ix <- TH[[id]]$ix
  OUT[[id]] <- summ(id, uc(cf[, ix], up, dn), meta$age[ix], length(up), length(dn)); cat(id, "done\n")
}
rm(cf); invisible(gc())
## De Donno (Glucocorticoid: D1-D3)
dm <- fread("data/DeDonno2026/output/cell_meta_qc.csv"); ref <- readRDS("data/DeDonno2026/output/references_dedonno.rds")$CORT
X <- readRDS("data/DeDonno2026/output/counts_universe.rds"); cu <- X$counts; colnames(cu) <- X$cells; rm(X); invisible(gc())
up <- intersect(ref$up, rownames(cu)); dn <- intersect(ref$down, rownames(cu))
DT <- list(D1 = "Oligodendrocytes", D2 = "Astro", D3 = "Endothelial")
for (id in names(DT)) {
  tix <- which(dm$cell_type == DT[[id]] & dm$line == "Ctrl")
  OUT[[id]] <- summ(id, uc(cu[, dm$cell_id[tix]], up, dn), dm$condition[tix], length(up), length(dn)); cat(id, "done\n")
}
## The four type A tasks added in Supplementary Fig. S4 (PS2APP: DA1; Forced swim: AS1-AS3): the share of zeros of the signed UCell
## score stored by fig5_extra_score.R at L = 500 (the score drawn in Supplementary Fig. S4b; U+ and U- are not stored there, so the
## check of summ() is not repeated). Input: figs/Fig5_benchmark/output/extra_refs/fig5_extra_{percell,auc}_L500.*
EXD <- "figs/Fig5_benchmark/output/extra_refs"
ex <- fread(cmd = sprintf("gzip -dc '%s/fig5_extra_percell_L500.csv.gz'", EXD), select = c("id", "group", "UCell"))
ng <- unique(fread(file.path(EXD, "fig5_extra_auc_L500.csv"))[, .(id, ref_up, ref_down)])
for (tk in c("DA1", "AS1", "AS2", "AS3")) {
  g <- ng[id == tk]; stopifnot(nrow(g) == 1L, ex[id == tk, .N] > 0L)
  OUT[[tk]] <- ex[id == tk, .(id = tk, n_cells = .N, up_genes = g$ref_up, down_genes = g$ref_down, share_score0 = mean(UCell == 0)),
                  by = group]
  cat(tk, "done\n")
}
res <- rbindlist(OUT); setcolorder(res, "id")
dir.create("figs/Fig5_benchmark/output/results", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "figs/Fig5_benchmark/output/results/fig5_ucell_zero.csv"); print(res, digits = 3)
