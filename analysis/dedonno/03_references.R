## De Donno step 3 — the 4 reference axes (external DEG or named gene set).
##   CORT  signed  GSE184924 corticosterone 3 mg/kg 3 h vs vehicle (ventral hippocampus bulk)   [display name: Glucocorticoid-index]
##                 -> counts summed per group -> CPM over the nuclear genes -> log2(CPM + 1) difference; top 500 up / top 500 down
##   KA6h  up only GSE125068 hippocampal excitatory neurons, kainate SE 6 h vs saline (recount3)  [display name: Seizure-index]
##                 -> the same fold changes; the top 200 up-regulated genes, no expression floor (no DOWN list, so it is scored
##                    with two terms: the match with the cell's UP list minus the match with its DOWN list)
##   UPR   set     MSigDB Hallmark UNFOLDED_PROTEIN_RESPONSE (as-is)
##   NFkB  set     MSigDB Hallmark TNFA_SIGNALING_VIA_NFKB (as-is)
## One rule for the six ranked references of the paper: genes without informative names, mitochondrial genes and
## haemoglobin genes (from blood cells in the samples) are excluded, ribosomal protein genes are kept; the top N per direction are
## taken from the whole list and restricted to the genes analysed when the cells are scored (04_score.R); the eligible genes (the
## background of the matching) are the genes with a non-zero change, minus the excluded genes.
## CORT and KA6h are ranked lists (stored with type = "B"); UPR and NFkB are gene sets (type B of the paper; stored with type = "C").
## Fixed inputs: the two fold-change tables are computed once and kept in data/DeDonno2026/provided/ (CORT_FIX, KA6H_FIX), because
##   the Ensembl -> symbol mapping of CORT depends on the org.Mm.eg.db release (a newer release renamed a few genes and changed the
##   CORT index slightly, r = 0.997) and genes with tied fold changes can come in a different order in a new recount3 download.
##   A table is computed, and written there, only when it is absent (delete it to recompute): CORT from the GSE184924 count
##   matrix (data/DeDonno2026/external/), KA6h from recount3 (SRP179090; needs internet).
##
## Output: data/DeDonno2026/output/{DeDonno_GR_glucocorticoid.csv, DeDonno_KA6h_deferred.csv, references_dedonno.rds}

suppressPackageStartupMessages({ library(data.table); library(msigdbr) })
OUT <- "data/DeDonno2026/output"
CORT_TOPN <- 500L; KA6_TOPN <- 200L
bad_pc <- function(g) grepl("Rik[0-9]?$",g)|grepl("^Gm[0-9]",g)|grepl("^mt-",g)|grepl("^Hb[ab]-",g)|
  grepl("^A[CFLI][0-9]",g)|grepl("^[0-9]",g)|grepl("\\.[0-9]+$",g)|is.na(g)
cpm_nuc <- function(cnt, nuc) { s <- rowSums(cnt); s / sum(s[nuc]) * 1e6 }   # summed counts -> CPM over the nuclear genes

## ---------------- CORT (GSE184924) ----------------
CORT_FIX <- "data/DeDonno2026/provided/DeDonno_GR_glucocorticoid_261016.csv"   # computed 2026-10-16 with org.Mm.eg.db 3.22.0
if (file.exists(CORT_FIX)) {
  d <- fread(CORT_FIX); cat("CORT table: fixed input", CORT_FIX, sprintf("(%d genes)\n", nrow(d)))
} else {
  suppressPackageStartupMessages({ library(org.Mm.eg.db); library(AnnotationDbi) })
  M <- fread("data/DeDonno2026/external/GSE184924_Raw_gene_count_matrix.tsv.gz"); setnames(M, 1, "ens")
  veh  <- c("Sample_1","Sample_2","Sample_3","Sample_4")      # GSM5600427-30 vehicle
  cort <- c("Sample_5","Sample_6","Sample_7","Sample_8")      # GSM5600431-34 corticosterone
  M[, ens := sub("\\..*$", "", ens)]
  M[, Symbol := mapIds(org.Mm.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")]
  nuc <- !grepl("^mt-", M$Symbol)                             # genes without a symbol are nuclear and stay in the library size
  e <- cpm_nuc(as.matrix(M[, ..cort]), nuc); c0 <- cpm_nuc(as.matrix(M[, ..veh]), nuc)
  l2 <- log2(e + 1) - log2(c0 + 1)
  d <- data.table(Symbol = M$Symbol, `Fold Change` = round(sign(l2) * 2^abs(l2), 4), log2FC = round(l2, 4))[
    !is.na(Symbol) & is.finite(`Fold Change`) & `Fold Change` != 0]
  d <- d[order(-abs(`Fold Change`))][!duplicated(Symbol)]; setorder(d, -`Fold Change`)
  fwrite(d, CORT_FIX)
  cat(sprintf("CORT table: computed from GSE184924 with org.Mm.eg.db %s (%d genes; %d mitochondrial genes left out of the library size); wrote %s\n",
              as.character(packageVersion("org.Mm.eg.db")), nrow(d), sum(!nuc), CORT_FIX))
}
fwrite(d[, .(Symbol, `Fold Change`, log2FC)], file.path(OUT, "DeDonno_GR_glucocorticoid.csv"))
gr <- c("Fkbp5","Tsc22d3","Sgk1","Ddit4","Zbtb16","Klf9","Per1","Cebpd","Pdk4","Sult1a1")
cat("CORT reference: known GR targets (FC, rank):\n")
for (g in gr) { r <- match(g, d$Symbol); cat(sprintf("  %-8s %s\n", g, if (is.na(r)) "(absent)" else sprintf("%+.2f  rank %d/%d", d$`Fold Change`[r], r, nrow(d)))) }

## ---------------- KA6h (GSE125068 via recount3) ----------------
KA6H_FIX <- "data/DeDonno2026/provided/DeDonno_KA6h_261016.csv"   # recount3 table of SRP179090, downloaded 2026-10-16
ka_src <- "recount3 (SRP179090)"
if (file.exists(KA6H_FIX)) {
  kd <- fread(KA6H_FIX); cat("KA6h table: fixed input", KA6H_FIX, sprintf("(%d genes)\n", nrow(kd)))
} else {
  kd <- tryCatch({
    suppressMessages(library(recount3))
    ap  <- available_projects(organism = "mouse")
    row <- ap[ap$project == "SRP179090" & ap$project_type == "data_sources", ]
    rse <- create_rse(row[1, ])
    cd  <- as.data.frame(SummarizedExperiment::colData(rse))
    ids <- cd[[intersect(c("external_id","sra.run_acc","rail_id"), names(cd))[1]]]
    colFor <- function(v) which(vapply(ids, function(x) any(vapply(v, function(s) grepl(s, x), logical(1))), logical(1)))
    cSal <- colFor(c("SRR8441329","SRR8441330")); c6 <- colFor(c("SRR8441333","SRR8441334"))
    stopifnot(length(cSal) == 2, length(c6) == 2)
    cnt <- SummarizedExperiment::assay(rse, "raw_counts"); gn <- SummarizedExperiment::rowData(rse)$gene_name
    nuc <- !grepl("^mt-", gn)
    cat(sprintf("recount3 SRP179090: %d genes, %d mitochondrial (left out of the library size)\n", length(gn), sum(!nuc)))
    ee <- cpm_nuc(cnt[, c6, drop = FALSE], nuc); cc <- cpm_nuc(cnt[, cSal, drop = FALSE], nuc); ll <- log2(ee + 1) - log2(cc + 1)
    z <- data.table(Symbol = gn, `Fold Change` = round(sign(ll) * 2^abs(ll), 4), log2FC = round(ll, 6),
                    eCPM = round(ee, 3), cCPM = round(cc, 3))
    z[order(-abs(log2FC))][!duplicated(Symbol)]
  }, error = function(e) stop("recount3 download of SRP179090 failed (", conditionMessage(e), ") and the fixed input ", KA6H_FIX,
                              " is absent", call. = FALSE))
  fwrite(kd, KA6H_FIX)
  cat(sprintf("KA6h table: computed from recount3 %s (%d genes); wrote %s\n", as.character(packageVersion("recount3")), nrow(kd), KA6H_FIX))
}
fwrite(kd, file.path(OUT, "DeDonno_KA6h_deferred.csv"))
cat(sprintf("KA6h source: %s  (%d genes)\n", ka_src, nrow(kd)))

## ---------------- final 4 axes ----------------
r <- d[, .(Symbol, fc = `Fold Change`)][is.finite(fc) & fc != 0][!bad_pc(Symbol)][order(-abs(fc))][!duplicated(Symbol)]
k <- kd[is.finite(log2FC) & log2FC > 0][!bad_pc(Symbol)][order(-log2FC, Symbol)]
## eligible = genes each reference could contain (analysis/running_fisher.R: pair_universe): the genes with a finite, non-zero
## change, minus the excluded genes; sets = NULL.
k_elig <- kd[is.finite(log2FC) & log2FC != 0][!bad_pc(Symbol)]$Symbol
h <- as.data.table(msigdbr(species = "Mus musculus", collection = "H"))
hset <- function(nm) sort(unique(h[gs_name == nm]$gene_symbol))
refs <- list(
  CORT = list(up = head(r[fc > 0][order(-fc)]$Symbol, CORT_TOPN), down = head(r[fc < 0][order(fc)]$Symbol, CORT_TOPN), type = "B",
              eligible = unique(r$Symbol),
              desc = "GSE184924 corticosterone 3 mg/kg 3 h vs vehicle (hippocampus bulk); top-500 per direction by |FC|"),
  KA6h = list(up = head(k$Symbol, KA6_TOPN), down = character(0), type = "B", eligible = unique(k_elig),
              desc = sprintf("GSE125068 hippocampal excitatory neurons, kainate SE 6 h vs saline; up-only top-200, no expression floor [%s]", ka_src)),
  UPR  = list(up = hset("HALLMARK_UNFOLDED_PROTEIN_RESPONSE"), down = character(0), type = "C", eligible = NULL, desc = "MSigDB Hallmark UNFOLDED_PROTEIN_RESPONSE"),
  NFkB = list(up = hset("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),   down = character(0), type = "C", eligible = NULL, desc = "MSigDB Hallmark TNFA_SIGNALING_VIA_NFKB"))
saveRDS(refs, file.path(OUT, "references_dedonno.rds"))
for (a in names(refs)) cat(sprintf("  [%-4s] type %s  up=%d down=%d eligible=%s | %s\n", a, refs[[a]]$type, length(refs[[a]]$up),
                                   length(refs[[a]]$down), if (is.null(refs[[a]]$eligible)) "any gene" else length(refs[[a]]$eligible),
                                   paste(head(refs[[a]]$up, 6), collapse = ", ")))
cat(sprintf("msigdbr %s\nwrote references_dedonno.rds\n", as.character(packageVersion("msigdbr"))))
