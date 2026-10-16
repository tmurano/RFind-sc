## Builds the demo data set in demo/data/ from the full Hammond et al. 2019 data (not needed to run the demo; kept to show how the
## demo files were made). Needs the raw h5ad and the outputs of analysis/hammond/run_all.sh.
##   cells : 1,000 microglia of the DAM task of Fig. 5 (task H1: white matter 7 days after LPC injection vs saline injection, P100;
##           QC cells, non-microglial Mono/Mac cells removed), 500 per group drawn at random (set.seed(1)).
##   genes : the pair universe of the DAM reference (Hammond gene universe of 8,960 genes intersected with the genes reported by
##           Keren-Shaul et al. 2017; analysis/running_fisher.R: pair_universe). The per-cell UP / DOWN lists of the paper are restricted
##           to this universe before they are cut to 500 genes, so these genes give the same lists as the whole universe.
##   counts: raw UMI counts from layers/counts of the h5ad (gene names with "-" restored as in analysis/hammond/03_extract_qc.R).
##   baseline: the gene shares (share) and the per-cell total UMI over all genes (cs) of data/Hammond2019/output/delta_P2_basis.rds,
##           written with 17 significant digits, so that Delta is computed exactly as in analysis/hammond/05_delta_P2.R.
##   reference: the DAM reference of the paper (references_panel9.rds: the top 500 genes per direction of Keren-Shaul et al. 2017,
##           Table S3, by fold change, taken before the genes of the dataset), restricted to the pair universe (408 + 427 genes).
##   expected : the DAM-index of these cells stored by the paper's pipeline (data/Hammond2019/output/scores_hammond_9axes.csv).
## Run from the project root: R --no-save --no-restore -f demo/make_demo_data.R   (R >= 4.1 with data.table, Matrix, hdf5r; < 1 GB RAM)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(hdf5r) })
OUT <- "demo/data"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
H5 <- "data/Hammond2019/raw/Hammond2019_microglia_AllSamples_annotated.h5ad"
N_PER_GROUP <- 500L

basis <- readRDS("data/Hammond2019/output/delta_P2_basis.rds")          # genes, share (per gene), cs (per QC cell, all genes)
panel <- readRDS("data/Hammond2019/output/references/references_panel9.rds"); e <- panel$DAM
meta  <- fread("data/Hammond2019/raw/cell_meta_hammond.csv")
stored <- fread("data/Hammond2019/output/scores_hammond_9axes.csv", select = c("cell_id", "DAM"))

## ---- cells (task H1 of figs/Fig5_benchmark/analysis/fig5_bench_redo.R) ----
qc <- names(basis$cs)
m <- meta[match(qc, cell_id)]
grp <- c(P100_Injured_WM = "LPC", P100_Control_WM = "saline")
ok <- m$paper_cluster != "Mono/Mac" & m$age %in% names(grp)
set.seed(1)
pick <- unlist(lapply(names(grp), function(a) sort(sample(which(ok & m$age == a), N_PER_GROUP))))
cells <- qc[pick]; lab <- unname(grp[m$age[pick]])
cat(sprintf("task H1 cells: %d LPC, %d saline (QC microglia); demo: %d + %d\n",
            sum(ok & m$age == "P100_Injured_WM"), sum(ok & m$age == "P100_Control_WM"), sum(lab == "LPC"), sum(lab == "saline")))

## ---- genes: pair universe of the DAM reference (as analysis/running_fisher.R: pair_universe) ----
glow <- tolower(basis$genes)
inU <- glow %in% unique(tolower(e$eligible))
genes <- basis$genes[inU]; share <- as.numeric(basis$share)[inU]
dupl <- genes[tolower(genes) %in% tolower(genes)[duplicated(tolower(genes))]]
if (length(dupl)) cat("genes that differ only in case (kept, as in the paper; P counts them once):", paste(dupl, collapse = ", "), "\n")
ref_up <- intersect(tolower(e$up), tolower(genes)); ref_dn <- intersect(tolower(e$down), tolower(genes))
cat(sprintf("pair universe: %d genes, P = %d distinct lower-case names (Hammond universe %d); DAM reference in it: %d up, %d down\n",
            length(genes), length(unique(tolower(genes))), length(glow), length(ref_up), length(ref_dn)))
stopifnot(length(ref_up) > 0L, length(ref_dn) > 0L)

## ---- raw counts of these cells and genes from the h5ad (CSR, cells x genes) ----
f <- H5File$new(H5, "r")
h5g <- f[["var/_index"]]$read(); h5c <- f[["obs/_index"]]$read(); ptr <- f[["layers/counts/indptr"]]$read()
dat <- f[["layers/counts/data"]]; ind <- f[["layers/counts/indices"]]
h5name <- ifelse(genes %in% h5g, genes, gsub("-", ".", genes, fixed = TRUE))   # names with "-" are stored with "." in the h5ad
gi <- match(h5name, h5g); stopifnot(!anyNA(gi), !anyDuplicated(gi))
ci <- match(cells, h5c); stopifnot(!anyNA(ci))
ii <- list(); jj <- list(); xx <- list(); tot <- numeric(length(ci))
for (k in seq_along(ci)) {
  a <- ptr[ci[k]] + 1; b <- ptr[ci[k] + 1]
  if (b < a) next
  x <- dat[a:b]; g <- ind[a:b] + 1L
  tot[k] <- sum(x)
  r <- match(g, gi); keep <- !is.na(r)
  ii[[k]] <- r[keep]; jj[[k]] <- rep(k, sum(keep)); xx[[k]] <- x[keep]
}
f$close_all()
X <- sparseMatrix(i = unlist(ii), j = unlist(jj), x = unlist(xx), dims = c(length(genes), length(cells)))
cs <- as.numeric(basis$cs[cells])
stopifnot(all(tot == cs))                                          # total UMI over all genes = cs of the paper
stopifnot(all(X@x == round(X@x)))
cat(sprintf("counts: %d genes x %d cells, %d non-zero entries\n", nrow(X), ncol(X), length(X@x)))

## ---- write ----
mtx <- file.path(OUT, "counts.mtx")
invisible(writeMM(as(X, "CsparseMatrix"), mtx))
invisible(system2("gzip", c("-f", "-9", "-n", shQuote(mtx))))
gz <- function(p) gzfile(p, "w")
con <- gz(file.path(OUT, "barcodes.tsv.gz")); writeLines(cells, con); close(con)
con <- gz(file.path(OUT, "genes.tsv.gz"))
writeLines(c("gene\tbaseline_share", sprintf("%s\t%.17g", genes, share)), con); close(con)
fwrite(data.table(cell_id = cells, condition = lab, total_umi = sprintf("%.17g", cs)), file.path(OUT, "labels.tsv"), sep = "\t")
fwrite(data.table(direction = rep(c("up", "down"), c(length(ref_up), length(ref_dn))), rank = c(seq_along(ref_up), seq_along(ref_dn)),
                  gene = c(ref_up, ref_dn)),
       file.path(OUT, "dam_reference.tsv"), sep = "\t")
dam <- stored$DAM[match(cells, stored$cell_id)]; stopifnot(!anyNA(dam))
fwrite(data.table(cell_id = cells, condition = lab, DAM_index = sprintf("%.17g", dam)), "demo/expected_output.csv")
cat("written:", paste(list.files(OUT), collapse = ", "), "and demo/expected_output.csv\n")
