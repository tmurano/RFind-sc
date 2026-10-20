## Check that the RFindsc package (1.0.1; its tests are taken from the clone in tools/RFindsc) runs and reproduces the RFind-sc paper pipeline
##   1  the package's own tests (minimal stand-ins for the testthat functions used, so testthat is not needed)
##   2  kernel: running_fisher() vs analysis/running_fisher.R on random ranked lists
##   3  Hammond: per-cell scores computed with the package's functions (P2 signature from the paper's baseline, bundled
##      4 axes restricted to the universe) vs analysis/hammond/07 scores (data/Hammond2019/output/scores_hammond_9axes.csv)
##   4  De Donno: the same with the per-cell-type baseline and the bundled 4 axes vs analysis/dedonno/04 scores
##   5  score_cells() end to end on 2,000 Hammond cells (baseline from those cells only, so values differ from the paper;
##      checks that it runs and agrees in rank)
## Run: git clone https://github.com/tmurano/RFindsc tools/RFindsc; install RFindsc into <lib> (README, "For reviewers");
##      RFINDSC_LIB=<lib> Rscript analysis/check_rfindsc.R   (from project root)
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
LIB <- Sys.getenv("RFINDSC_LIB"); suppressPackageStartupMessages(library(RFindsc, lib.loc = LIB))
cat("RFindsc", as.character(packageVersion("RFindsc", lib.loc = LIB)), "\n")
ns <- asNamespace("RFindsc"); set.seed(1)
## a name that differs from another only in case (PISD and Pisd in the data of Hammond et al.; three names in those of De Donno
## et al.) can occur twice in a cell's lists; RFindsc >= 0.3.9 warns about it in running_fisher(), which steps 3 and 4 call
## cell by cell, so the expected warning is muted there (the scores are computed as before)
quiet_dup <- function(expr) withCallingHandlers(expr, rfindsc_duplicated_names = function(w) invokeRestart("muffleWarning"))

## ---- 1 package tests ----
chk <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
expect_equal <- function(a, b, tolerance = 1e-8, ...) chk(isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE)),
                                                            sprintf("not equal: %s vs %s", format(a)[1], format(b)[1]))
expect_identical <- function(a, b) chk(identical(a, b), "not identical")
expect_gt <- function(a, b) chk(a > b, "not >"); expect_gte <- function(a, b) chk(a >= b, "not >="); expect_lte <- function(a, b) chk(a <= b, "not <=")
expect_type <- function(x, t) chk(typeof(x) == t, "wrong type"); expect_true <- function(x) chk(isTRUE(x), "not TRUE")
expect_setequal <- function(a, b) chk(setequal(a, b), "sets differ")
expect_false <- function(x) chk(identical(x, FALSE), "not FALSE")
expect_error <- function(expr, ...) chk(inherits(tryCatch({ expr; NULL }, error = function(e) e), "error"), "no error raised")
expect_s3_class <- function(x, cl) chk(inherits(x, cl), "wrong class"); expect_length <- function(x, n) chk(length(x) == n, "wrong length")
skip_if_not <- function(condition, message = NULL) if (!isTRUE(condition)) stop("skipped: ", message)
test_path <- function(...) file.path("tools/RFindsc/tests/testthat", ...)
expect_lt <- function(a, b) chk(a < b, "not <")
expect_message <- function(expr, regexp = NULL, ...) {            # a message (matching regexp) is emitted
  msgs <- character(0)
  withCallingHandlers(expr, message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  chk(length(msgs) > 0 && (is.null(regexp) || any(grepl(regexp, msgs))), "no matching message")
}
expect_null <- function(x) chk(is.null(x), "not NULL")
capture_warnings <- function(expr) {                            # the messages of the warnings given while evaluating expr
  w <- character(0)
  withCallingHandlers(expr, warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  w
}
expect_warning <- function(expr, regexp = NULL, class = NULL, ...) {   # a warning matching regexp and of class is given
  hit <- FALSE
  withCallingHandlers(expr, warning = function(w) {
    if ((is.null(regexp) || grepl(regexp, conditionMessage(w))) && (is.null(class) || inherits(w, class))) hit <<- TRUE
    invokeRestart("muffleWarning") })
  chk(hit, "no matching warning")
}
expect_silent <- function(expr) chk(!length(capture_warnings(expr)), "a warning was given")
skip_if_not_installed <- function(pkg, minimum_version = NULL)
  if (!requireNamespace(pkg, quietly = TRUE) || (!is.null(minimum_version) && utils::packageVersion(pkg) < minimum_version))
    stop("skipped: ", pkg, " not installed")
skip_if <- function(condition, message = NULL) if (isTRUE(condition)) stop("skipped: ", message)
undo <- list()                                                   # functions that undo local_mocked_bindings() after a test
local_mocked_bindings <- function(..., .package) {               # replaces functions of a namespace until the test ends
  e <- asNamespace(.package); new <- list(...); old <- mget(names(new), envir = e)
  put <- function(v) for (n in names(v)) { unlockBinding(n, e); assign(n, v[[n]], envir = e); lockBinding(n, e) }
  put(new); undo <<- c(undo, list(function() put(old)))
}
n_fail <- 0L; n_skip <- 0L
test_that <- function(desc, code) {
  on.exit({ for (f in rev(undo)) f(); undo <<- list() })
  res <- tryCatch({ eval(substitute(code), envir = new.env(parent = ns)); "ok" },
                  error = function(e) if (startsWith(conditionMessage(e), "skipped: ")) "skip" else {
                    cat("      ", conditionMessage(e), "\n"); "FAIL" })
  if (res == "FAIL") n_fail <<- n_fail + 1L
  if (res == "skip") n_skip <<- n_skip + 1L
  cat(sprintf("  [%s] %s\n", res, desc))
}
cat("\n1. package tests\n")
if (!dir.exists("tools/RFindsc/tests/testthat"))
  stop("tools/RFindsc/tests/testthat not found: clone the package into tools/RFindsc first (README, \"Companion software\")", call. = FALSE)
for (f in list.files("tools/RFindsc/tests/testthat", pattern = "^test-.*\\.R$", full.names = TRUE)) sys.source(f, envir = environment())
cat(sprintf("  -> %d failed, %d skipped\n", n_fail, n_skip))

## ---- 2 kernel ----
E <- new.env(); sys.source("analysis/running_fisher.R", envir = E)
pool <- paste0("g", 1:3000)
d <- replicate(300, { g1 <- sample(pool, sample(20:500, 1)); g2 <- unique(c(sample(g1, min(length(g1), sample(0:30, 1))), sample(pool, sample(20:500, 1))))
                      RFindsc::running_fisher(g1, g2, 3000) - E$running_fisher(g1, g2, 3000) })
cat(sprintf("\n2. kernel: running_fisher() vs analysis/running_fisher.R, 300 random pairs: max |diff| = %g\n", max(abs(d))))

## ---- 3 Hammond ----
cf <- readRDS("data/Hammond2019/output/counts_filtered.rds"); cs <- Matrix::colSums(readRDS("data/Hammond2019/output/counts_qc.rds"))
gl <- tolower(rownames(cf)); U <- length(gl)
rate <- Matrix::rowMeans(cf %*% Diagonal(x = 1e6 / cs)) / 1e6            # the paper's baseline: mean CPM of all cells
bs <- default_biosets(); ax <- setdiff(names(bs), "meta")
## 0.3.5: each reference in its pair universe (prepare_panel), each cell's lists cut inside it (score_cell_ranking)
b2 <- lapply(bs[ax], function(p) ns$prepare_panel(ns$to_b2(p), gl, 500))
cat("\n3. pair universe per axis (package): ", paste(sprintf("%s %d", ax, vapply(b2, `[[`, 0L, "universe_size")), collapse = ", "), "\n")
pipe <- fread("data/Hammond2019/output/scores_hammond_9axes.csv"); stopifnot(identical(pipe$cell_id, colnames(cf)))
cells <- sample(ncol(cf), 300)
S <- t(vapply(cells, function(j) { r <- ns$rank_cell_residuals(as.numeric(cf[, j]), rate * cs[j])
  vapply(ax, function(a) quiet_dup(ns$score_cell_ranking(r, b2[[a]], gl, n_cell = 500, clip = 300)), 0) }, numeric(length(ax))))
dh <- abs(S - as.matrix(pipe[cells, ..ax]))
cat("\n3. Hammond, 300 cells, package vs analysis/hammond/07 (max |diff| per axis):\n"); print(signif(apply(dh, 2, max), 3))

## ---- 4 De Donno ----
X <- readRDS("data/DeDonno2026/output/counts_universe.rds"); meta <- fread("data/DeDonno2026/output/cell_meta_qc.csv")
pdd <- fread("data/DeDonno2026/output/scores_dedonno.csv"); stopifnot(identical(pdd$cell_id, X$cells))
glD <- tolower(rownames(X$counts)); UD <- length(glD)
## names = the package's axis names (RFindsc >= 0.3.1), values = the pipeline's keys in scores_dedonno.csv
dd <- load_dedonno_biosets(); map <- c(Glucocorticoid = "CORT", Seizure = "KA6h", UPR = "UPR", `NF-kB` = "NFkB")
b2d <- lapply(dd[names(map)], function(p) ns$prepare_panel(ns$to_b2(p), glD, 500))
cat("\n4. pair universe per axis (package): ", paste(sprintf("%s %d", names(map), vapply(b2d, `[[`, 0L, "universe_size")), collapse = ", "), "\n")
res4 <- rbindlist(lapply(c("DG Glut", "Oligodendrocytes"), function(ctp) {
  idx <- which(meta$cell_type == ctp)
  rate_ct <- Matrix::rowMeans(X$counts[, idx] %*% Diagonal(x = 1e6 / X$cs[idx])) / 1e6   # per-cell-type baseline (02_signature.R)
  pick <- sample(idx, 100)
  Sd <- t(vapply(pick, function(j) { r <- ns$rank_cell_residuals(as.numeric(X$counts[, j]), rate_ct * X$cs[j])
    vapply(names(map), function(a) quiet_dup(ns$score_cell_ranking(r, b2d[[a]], glD, n_cell = 500, clip = 300)), 0) },
    numeric(length(map))))
  data.table(cell_type = ctp, axis = names(map), max_abs_diff = apply(abs(Sd - as.matrix(pdd[pick, unname(map), with = FALSE])), 2, max))
}))
cat("\n4. De Donno, 100 cells each of DG and oligodendrocytes, package vs analysis/dedonno/04:\n"); print(res4[, .(cell_type, axis, max_abs_diff = signif(max_abs_diff, 3))])

## ---- 5 score_cells end to end ----
sub <- sample(ncol(cf), 2000)
t0 <- Sys.time()
out <- score_cells(cf[, sub], default_biosets(), rfindsc_config(n_cores = max(1, parallel::detectCores() - 2)), libsize = cs[sub])
cat(sprintf("\n5. score_cells() on 2,000 cells: %.0f s, universe %d genes, scores %d x %d\n", as.numeric(Sys.time() - t0, units = "secs"),
            length(out$universe), nrow(out$scores), ncol(out$scores)))
cat("   Spearman with the paper's scores (baseline from these 2,000 cells only):\n")
print(round(vapply(ax, function(a) cor(out$scores[, a], pipe[[a]][sub], method = "spearman"), 0), 3))

## ---- 6 a full-length DEG table gives the paper's reference: with the genes that the paper leaves out of every reference
## removed first (genes without informative names, mitochondrial and haemoglobin genes; analysis/hammond/06_build_references.R),
## the package keeps the top 500 per side of the full table and then the genes of the pair universe (prepare_panel(), 0.3.8),
## which must equal the bundled DAM after the same step. Compared as gene lists on the paper's 8,960-gene universe (scoring a
## cell subset with score_cells() would recompute the universe from those cells).
bad <- function(g) grepl("Rik[0-9]?$", g) | grepl("^Gm[0-9]", g) | grepl("^mt-", g) | grepl("^Hb[ab]-", g) |
  grepl("^A[CFLI][0-9]", g) | grepl("^[0-9]", g) | grepl("\\.[0-9]+$", g) | is.na(g)
k <- fread("data/KerenShaul2017/provided/Mouse_DAM_KerenShaul2017_vs_homeostatic.csv")[is.finite(`Fold Change`) & `Fold Change` != 0][!bad(Symbol)]
k <- k[order(-abs(`Fold Change`))][!duplicated(Symbol)]
full <- build_panel_from_deg(data.frame(gene = k$Symbol, log2FC = sign(k$`Fold Change`) * log2(abs(k$`Fold Change`))), label = "DAM_full")
prep <- function(p) ns$prepare_panel(ns$to_b2(p), gl, 500L)
pf <- prep(full); pb <- prep(bs$DAM)
cat(sprintf("\n6. DAM from the full Keren-Shaul table (%d up / %d down genes) vs the bundled DAM, both cut by the package (top 500 per\n   side, then the pair universe: %d up / %d down genes): up identical %s, down identical %s\n",
            nrow(full$up), nrow(full$down), length(pb$up), length(pb$down), identical(pf$up, pb$up), identical(pf$down, pb$down)))
