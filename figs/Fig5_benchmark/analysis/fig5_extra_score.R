## Supplementary Fig. S4 / Fig. 5e — the four added type A tasks of the benchmark, scored with RFind-sc and the methods of the Fig. 5
## benchmark (the settings and functions of fig5_bench_redo.R, read from that script; paper settings). The tasks were fixed before any
## result was seen; each uses the cells of a task of the main figure with one of the two references of fig5_extra_refs.R:
##   DA1      the cells of H1 (Hammond: LPC- vs saline-injected white matter, microglia) with the PS2APP reference
##   AS1-AS3  the cells of D1-D3 (De Donno: stress vs no stress, control line; oligodendrocytes, astrocytes, endothelial cells) with the
##            Forced swim reference
## RFind-sc: as analysis/hammond/07_score_axes.R and analysis/dedonno/04_score.R (the deep per-cell rankings, the pair universe, the
## top 500 cell genes per direction within it, clip 300). Every other method gets RFind-sc's reference genes (up / down within the pair
## universe; De Donno names matched case-insensitively, as fig5_bench_redo.R does) and the log2FC weights of the fold-change table.
## As in fig5_bench_redo.R: VISION on all QC cells of the data set (Hammond) or of the cell type (De Donno), AUCell rankings over all
## QC Hammond cells and over the task cells in De Donno, AddModuleScore / scFOCAL-type / Beyondcell with all cells of the cell type,
## Beyondcell cells without a score set to 0 (counted in n_na), UCell maxRank raised above 1,500 signature genes (note column),
## AUC with the paired percentile bootstrap over cells (boot_task: 1,000 resamples, set.seed(1) per task).
## Reference length L (Fig. 5e): the top L genes per direction of the fold-change table ('all' = every gene of that sign), by the rules of
## fig5_extra_refs.R (the top L of the whole list, cut to the analysed genes when used, as all references);
## RFind-sc's cell lists stay at 500 (as fig5_bench_redo.R). At L = 500 the lists must equal the stored references (the script stops
## otherwise).
## Inputs: figs/Fig5_benchmark/output/extra_refs/fig5_extra_refs.rds, data/GSE75431/output/dam_PS2APP_FC.csv,
##         data/GSE169505/output/acutestress_GSE169505_FC.csv (fig5_extra_refs.R); the Hammond and De Donno pipeline outputs read by
##         fig5_bench_redo.R (per-cell signatures, counts, cell metadata); figs/Fig5_benchmark/analysis/fig5_bench_redo.R (functions)
## Outputs (<outdir>, default figs/Fig5_benchmark/output/extra_refs/; names get _<parts> when only some parts are run):
##   fig5_extra_auc_L<L>.csv          AUC with 95% CI and the difference from RFind-sc with its CI, tasks DA1, AS1, AS2, AS3 (the columns
##                                    of results/fig5_redo_auc_L<L>.csv)
##   fig5_extra_percell_L<L>.csv.gz   per-cell scores (Beyondcell raw, with its NaN); not in the repository
##   logs/runs/fig5_extra/fig5_extra_L<L>.txt   summary of the run
## Run from the project root (R 4.5 with the benchmark packages; about 5 min per length and about 22 min at all genes on 14 cores,
## about 10 GB of memory):
##   R --no-save --no-restore -f figs/Fig5_benchmark/analysis/fig5_extra_score.R --args [L=500] [parts=DA,AS] [outdir=...]
##   (L = 100 | 250 | 500 | 1000 | 2000 | all; parts = DA (Hammond, DA1) and / or AS (De Donno, AS1-AS3))
pdf(NULL); options(warn = 1)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(parallel) })
ARGS <- commandArgs(trailingOnly = TRUE)
argv <- function(key, default) { a <- grep(paste0("^", key, "="), ARGS, value = TRUE); if (length(a)) sub(paste0("^", key, "="), "", a[1]) else default }
PARTS <- strsplit(argv("parts", "DA,AS"), ",")[[1]]; stopifnot(length(PARTS) > 0, all(PARTS %in% c("DA", "AS")))
PARTS <- intersect(c("DA", "AS"), PARTS)                          # Hammond first
LARG <- argv("L", "500"); stopifnot(LARG %in% c("100", "250", "500", "1000", "2000", "all"))   # reference length (Fig. 5e)
OUT <- argv("outdir", "figs/Fig5_benchmark/output/extra_refs"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- the settings and functions of fig5_bench_redo.R (top-level assignments of a function or a single constant only) ----
BENCH <- "figs/Fig5_benchmark/analysis/fig5_bench_redo.R"
is_def <- function(e) is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
  ((is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) || (is.atomic(e[[3]]) && length(e[[3]]) == 1L))
DEFS <- Filter(is_def, as.list(parse(BENCH, keep.source = FALSE)))
for (e in DEFS) eval(e, envir = globalenv())
FCKEYS <- c(ULMfc = "ULM (log2FC)", WMEANfc = "WMEAN (log2FC)")
UCELL_NOTE <- character(0); ORA_CHK <- character(0); SSG_CHK <- character(0); LOG <- character(0)
NC <- max(1L, parallel::detectCores() - 2L)                       # as fig5_bench_redo.R
options(Seurat.object.assay.version = "v3")
t0 <- Sys.time()
source("analysis/running_fisher.R")
SUF <- paste0("_L", LARG, if (setequal(PARTS, c("DA", "AS"))) "" else paste0("_", paste(PARTS, collapse = "-")))
LOGD <- "logs/runs/fig5_extra"; dir.create(LOGD, showWarnings = FALSE, recursive = TRUE)
TXT <- file.path(LOGD, sprintf("fig5_extra%s.txt", SUF))
lg(sprintf("Added type A tasks of the benchmark - %s; L = %s; parts %s; %d cores; VISION_ROUTE = %s; %s", format(Sys.time(), "%Y-%m-%d %H:%M"),
           LARG, paste(PARTS, collapse = ", "), NC, VISION_ROUTE, R.version.string))
PKLOG <- c("UCell", "AUCell", "decoupleR", "singscore", "VISION", "beyondcell", "Seurat", "GSVA")
lg("Packages: ", paste(sprintf("%s %s", PKLOG, vapply(PKLOG, pv, "")), collapse = "; "))
REF <- readRDS("figs/Fig5_benchmark/output/extra_refs/fig5_extra_refs.rds")
FCT <- list(AcuteStress = fread("data/GSE169505/output/acutestress_GSE169505_FC.csv"), DAM_PS2APP = fread("data/GSE75431/output/dam_PS2APP_FC.csv"))
## reference lists at length L; at L = 500 they must equal the stored references
LVAL <- if (LARG == "all") .Machine$integer.max else as.integer(LARG)
for (k in names(REF)) {
  ft <- FCT[[k]]
  up <- head(ft[`Fold Change` > 0][order(-`Fold Change`)]$Symbol, LVAL); dn <- head(ft[`Fold Change` < 0][order(`Fold Change`)]$Symbol, LVAL)
  if (LARG == "500" && !(identical(up, REF[[k]]$up) && identical(dn, REF[[k]]$down))) stop(k, ": the L = 500 rebuild differs from the stored reference", call. = FALSE)
  REF[[k]]$up <- up; REF[[k]]$down <- dn
}
lg(sprintf("Reference lists at L = %s: %s", LARG, paste(sprintf("%s up %d / down %d", names(REF), lengths(lapply(REF, `[[`, "up")),
                                                          lengths(lapply(REF, `[[`, "down"))), collapse = "; ")))
METH <- c("UCell (up genes)", "UCell", "AUCell", "AddModuleScore", "ULM (log2FC)", "WMEAN (log2FC)", "ORA", "singscore", "VISION",
          "Beyondcell", "ssGSEA", "scFOCAL-type")

## ---- RFind-sc per cell (07_score_axes.R / 04_score.R) ----
rf_score <- function(sig, ref, cols) {
  glow <- tolower(sig$genes); U <- pair_universe(sig$genes, ref$eligible); inu <- glow %in% U; P <- length(U)
  b2 <- list(up = intersect(tolower(ref$up), U), down = intersect(tolower(ref$down), U), ranked = TRUE)
  out <- numeric(length(cols)); short <- 0L
  for (r in seq_along(cols)) { j <- cols[r]
    uj <- cut_in_universe(sig$up[, j], inu, 500L); dj <- cut_in_universe(sig$down[, j], inu, 500L)
    short <- short + (length(uj) < 500L && sig$n_pos[j] > sig$n_deep) + (length(dj) < 500L && sig$n_neg[j] > sig$n_deep)
    out[r] <- rfind_score_cell(glow[uj], glow[dj], b2, universe_size = P, clip = 300)
  }
  if (short != 0L) stop("lists cut short by the deep ranking: ", short, call. = FALSE)
  list(score = out, P = P, up = b2$up, dn = b2$down)
}

## ---- evaluation ----
RES <- list(); PC <- list()
bc0 <- function(S) { if ("Beyondcell" %in% colnames(S)) { na <- is.na(S[, "Beyondcell"]); lg(sprintf("  Beyondcell NaN %d of %d cells -> 0", sum(na), length(na))); S[na, "Beyondcell"] <- 0 }; S }
mk <- function(rfv, sc) { sc <- lapply(sc[intersect(METH, names(sc))], function(v) as.numeric(unname(v))); do.call(cbind, c(list(`RFind-sc` = as.numeric(rfv)), sc)) }
unote <- function(sc) {                                           # UCell rows: the note of a raised maxRank (as fig5_bench_redo.R)
  x <- vapply(c("UCell (up genes)", "UCell"), function(m) { a <- attr(sc[[m]], "note"); if (is.null(a)) "" else a }, "")
  x[nzchar(x)]
}
auc_task <- function(tid, task, refname, cell_id, S, pos, grp, n_up, n_dn, m_note) {
  stopifnot(nrow(S) == length(pos), !anyNA(S[, setdiff(colnames(S), "Beyondcell")]), any(pos), any(!pos))
  PC[[tid]] <<- data.table(id = tid, cell_id = cell_id, positive = pos, group = grp, S)
  nav <- colSums(is.na(S))                                        # cells without a score (Beyondcell only; 0 for the others)
  r <- boot_task(tid, task, "condition", refname, bc0(S), pos)
  r[, c("n_na", "L", "ref_up", "ref_down", "note") := list(unname(nav[r$method]), LARG, n_up, n_dn, "")]
  for (m in intersect(names(m_note), r$method)) r[method == m, note := m_note[[m]]]
  RES[[tid]] <<- r
  lg(sprintf("%s (%d cells, %d positive; reference %d / %d): %s", tid, length(pos), sum(pos), n_up, n_dn,
             paste(sprintf("%s %.3f", r$method, r$AUC), collapse = ", ")))
  fwrite(rbindlist(RES, fill = TRUE), file.path(OUT, sprintf("fig5_extra_auc%s.csv", SUF))); writeLines(LOG, TXT)   # partial table after every task
}

## ================================ Hammond: PS2APP, task DA1 ================================
if ("DA" %in% PARTS) {
  sig <- readRDS("data/Hammond2019/output/percell_signature_deep.rds"); cells <- sig$cells; univ <- sig$genes
  meta <- fread("data/Hammond2019/raw/cell_meta_hammond.csv")[match(cells, cell_id)]
  mg <- meta$paper_cluster != "Mono/Mac"
  ix <- which(mg & meta$age %in% c("P100_Injured_WM", "P100_Control_WM"))   # the cells of H1
  ref <- REF$DAM_PS2APP
  rr <- rf_score(sig, ref, ix)
  s <- list(up = intersect(ref$up, univ), dn = intersect(ref$down, univ))
  if (!identical(tolower(s$up), rr$up) || !identical(tolower(s$dn), rr$dn)) stop("DA: the genes given to the methods differ from RFind-sc's", call. = FALSE)
  lg(sprintf("DA (DAM_PS2APP): RFind-sc scored %d cells; pair universe %d; reference up %d / down %d", length(ix), rr$P, length(s$up), length(s$dn)))
  rm(sig); invisible(gc())
  libsize <- Matrix::colSums(readRDS("data/Hammond2019/output/counts_qc.rds")); invisible(gc())
  cf <- readRDS("data/Hammond2019/output/counts_filtered.rds")
  stopifnot(identical(colnames(cf), cells), identical(names(libsize), cells), identical(rownames(cf), univ))
  ln <- cf %*% Diagonal(x = 1e4 / libsize); ln@x <- log1p(ln@x); dimnames(ln) <- dimnames(cf)
  say(sprintf("Hammond: %d genes x %d cells", nrow(cf), ncol(cf)))
  sigs <- list(DAM_PS2APP = s)
  V <- vision_run(cf, libsize, ln, sigs); vision_log("Hammond (all QC cells)", V, sigs); VS <- V$scores; rm(V); invisible(gc())
  g <- c(s$up, s$dn); w <- fcw(FCT$DAM_PS2APP, g); ct <- cf[, ix]; lt <- libsize[ix]
  sc <- list()
  set.seed(1); rk <- AUCell::AUCell_buildRankings(cf, plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE)   # all QC cells
  sc$AUCell <- aucell_scores(rk, s$up, cells[ix]); rm(rk); invisible(gc())
  sc[["UCell (up genes)"]] <- ucell_scores(ct, s$up, character(0)); sc$UCell <- ucell_scores(ct, s$up, s$dn)
  sc$AddModuleScore <- ams_scores(ln, s$up, ix)
  sc$singscore <- dc_ss(ct, lt, s$up, s$dn)[, "singscore"]
  fc <- dc_fc(ct, lt, g, w); sc[["ULM (log2FC)"]] <- fc[, "ULM (log2FC)"]; sc[["WMEAN (log2FC)"]] <- fc[, "WMEAN (log2FC)"]
  sc$ORA <- ora_scores(ct, lt, s$up, "DAM_PS2APP")
  sc[["scFOCAL-type"]] <- cos_focal(ln, g, w, ix, focal = TRUE)[, 3]
  sc$VISION <- VS[cells[ix], "DAM_PS2APP"]
  sc$Beyondcell <- bc_scores(ln, s$up, s$dn, ix)$score
  sc$ssGSEA <- ssgsea_scores(ct, lt, s$up, s$dn, "DAM_PS2APP")
  age <- meta$age[ix]
  auc_task("DA1", "Hammond: DAM (PS2APP), LPC vs saline white matter (microglia)", "DAM_PS2APP", cells[ix], mk(rr$score, sc),
           age == "P100_Injured_WM", age, length(s$up), length(s$dn), unote(sc))
  rm(ct, sc, cf, ln, VS); invisible(gc()); say("DA1 done")
}

## ================================ De Donno: Forced swim, tasks AS1-AS3 ================================
if ("AS" %in% PARTS) {
  dm <- fread("data/DeDonno2026/output/cell_meta_qc.csv")
  DT <- c(AS1 = "Oligodendrocytes", AS2 = "Astro", AS3 = "Endothelial")
  TIX <- lapply(DT, function(ct) which(dm$cell_type == ct & dm$line == "Ctrl")); BIX <- lapply(DT, function(ct) which(dm$cell_type == ct))
  sig <- readRDS("data/DeDonno2026/output/percell_signature_celltype_deep.rds")
  stopifnot(identical(sig$cells, dm$cell_id), identical(sig$baseline, "per_celltype_all_cells"))
  ref <- REF$AcuteStress; RF <- lapply(TIX, function(ix) rf_score(sig, ref, ix)); rm(sig); invisible(gc())
  X <- readRDS("data/DeDonno2026/output/counts_universe.rds"); stopifnot(identical(X$cells, dm$cell_id))
  keep <- sort(unique(unlist(BIX))); cu <- X$counts[, keep]; colnames(cu) <- X$cells[keep]; lib <- setNames(X$cs, X$cells); rm(X); invisible(gc())
  U <- pair_universe(rownames(cu), ref$eligible)
  ## the genes given to every method = RFind-sc's genes, mapped case-insensitively to the De Donno names (map_genes of fig5_bench_redo.R;
  ## at L >= 1,000 names differing from the reference only in case enter the lists)
  s <- list(up = map_genes(ref$up, rownames(cu), U), dn = map_genes(ref$down, rownames(cu), U))
  lg(sprintf("AcuteStress: genes mapped case-insensitively (not exact-name matches: up %d, down %d)", length(setdiff(s$up, ref$up)), length(setdiff(s$dn, ref$down))))
  if (!identical(tolower(s$up), RF[[1]]$up) || !identical(tolower(s$dn), RF[[1]]$dn)) stop("AcuteStress: genes differ from RFind-sc's lists", call. = FALSE)
  up <- s$up; dn <- s$dn; g <- c(up, dn); w <- fcw(FCT$AcuteStress, g)
  lg(sprintf("AcuteStress: up %d / down %d in the De Donno genes (pair universe %d)", length(up), length(dn), RF[[1]]$P))
  for (id in names(DT)) {
    tix <- TIX[[id]]; bix <- BIX[[id]]; kk <- match(tix, bix); tc <- dm$cell_id[tix]
    lnb <- cu[, dm$cell_id[bix]] %*% Diagonal(x = 1e4 / lib[dm$cell_id[bix]]); lnb@x <- log1p(lnb@x); dimnames(lnb) <- list(rownames(cu), dm$cell_id[bix])
    ct <- cu[, tc]; lt <- lib[tc]; sc <- list()
    sc[["UCell (up genes)"]] <- ucell_scores(ct, up, character(0)); sc$UCell <- ucell_scores(ct, up, dn)
    set.seed(1); rk <- AUCell::AUCell_buildRankings(ct, plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE); sc$AUCell <- aucell_scores(rk, up, tc); rm(rk)
    sc$AddModuleScore <- ams_scores(lnb, up, kk)
    sc$singscore <- dc_ss(ct, lt, up, dn)[, "singscore"]
    fc <- dc_fc(ct, lt, g, w); sc[["ULM (log2FC)"]] <- fc[, "ULM (log2FC)"]; sc[["WMEAN (log2FC)"]] <- fc[, "WMEAN (log2FC)"]
    sc$ORA <- ora_scores(ct, lt, up, id)
    sc[["scFOCAL-type"]] <- cos_focal(lnb, g, w, kk, focal = TRUE)[, 3]
    sg <- list(AcuteStress = list(up = up, dn = dn))
    V <- vision_run(cu[, dm$cell_id[bix]], lib[dm$cell_id[bix]], lnb, sg); vision_log(sprintf("%s (%s)", id, DT[[id]]), V, sg); sc$VISION <- V$scores[tc, "AcuteStress"]; rm(V)
    sc$Beyondcell <- bc_scores(lnb, up, dn, kk)$score
    sc$ssGSEA <- ssgsea_scores(ct, lt, up, dn, id)
    auc_task(id, sprintf("De Donno: Acute stress, Stress vs No Stress, %s (control line)", DT[[id]]), "AcuteStress", tc, mk(RF[[id]]$score, sc),
             dm$condition[tix] == "Stress", dm$condition[tix], length(up), length(dn), unote(sc))
    rm(lnb, ct); invisible(gc()); say(sprintf("%s done", id))
  }
}

## ================================ tables ================================
fwrite(rbindlist(PC, fill = TRUE), file.path(OUT, sprintf("fig5_extra_percell%s.csv.gz", SUF)))
lgs("", "Checks"); if (length(SSG_CHK)) lgs(SSG_CHK); if (length(ORA_CHK)) lgs(ORA_CHK)
lgs("", "Comparisons with RFind-sc (* = RFind-sc higher, x = lower: 95% bootstrap interval of the difference excludes 0 and difference > 0.01):")
fmt <- function(r, v) paste(sprintf("%s %.3f%s", r$method, r[[v]], ifelse(r$method == "RFind-sc", "", ifelse(r$d_lo > 0 & r$diff_vs_RFind > 0.01, "*",
                         ifelse(r$d_hi < 0 & r$diff_vs_RFind < -0.01, "x", "")))), collapse = ", ")
for (r in RES) lg(sprintf("%-4s AUC      %s", r$id[1], fmt(r, "AUC")))
lg(sprintf("Output: %s, %s", file.path(OUT, sprintf("fig5_extra_auc%s.csv", SUF)), file.path(OUT, sprintf("fig5_extra_percell%s.csv.gz", SUF))))
lg(sprintf("Run time %.1f min.", as.numeric(Sys.time() - t0, units = "mins"))); writeLines(LOG, TXT); say("done; summary in ", TXT)
