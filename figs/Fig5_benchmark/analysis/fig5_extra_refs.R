## Supplementary Fig. S4 / Fig. 5e — the two additional type A references of the benchmark (data, groups and construction rules fixed
## before any benchmark result was seen). Each reference = ranked lists of up- and down-regulated genes by the signed linear fold change,
## built as the existing type A references of the paper:
##   Forced swim    von Ziegler et al., Nat. Commun. 13, 1824 (2022); GEO GSE169505 (bulk RNA-seq of the hippocampus). Experiment
##                  TimeSeriesFST, ventral hippocampus, male: 3 h after a 6-min forced swim (7 mice) vs home cage (8 mice). Built as the
##                  Glucocorticoid reference (analysis/dedonno/03_references.R): summed counts per group -> CPM over the nuclear
##                  genes -> log2(CPM + 1) difference; the genes excluded by bad_pc (below) left out; top 500 per direction among all
##                  genes (cut to the analysed genes when used). The genes allowed into the pair universe ('eligible') follow the rule
##                  of the Glucocorticoid reference; the rule is checked on the Glucocorticoid table (the script stops if it fails).
##                  Used for the tasks AS1-AS3 (the cells of the Glucocorticoid tasks D1-D3; fig5_extra_score.R).
##   PS2APP         Srinivasan et al., Nat. Commun. 7, 11295 (2016); EMBL-EBI Expression Atlas E-GEOD-75431 (= GEO GSE75431), raw
##                  counts per gene. Sorted cortical microglia at 13 months, PS2APP (5 mice) vs non-transgenic (5 mice) = the Atlas
##                  contrast "'PS2APP transgenic' vs 'non transgenic' in 'microglia cell' at '13 month'" (run accessions below).
##                  Built as the other references: summed per group -> CPM over the nuclear genes -> log2(CPM + 1) difference per
##                  gene row; one row per symbol (the largest |fold change|); the genes excluded by bad_pc left out; top 500 per
##                  direction of the whole list (cut to the Hammond universe when used, as all references). Mitochondrial rRNA
##                  (mt-Rnr2, mt-Rnr1) makes up most of the reads here, and more in PS2APP, so a CPM over all genes would shift every
##                  other gene down (median log2 fold change about -0.87 in the Hammond universe, about -0.19 over the nuclear genes);
##                  the all-gene version is only printed for comparison (median and Spearman correlation with the nuclear version).
##                  Used for the task DA1 (the cells of the DAM task H1; fig5_extra_score.R).
## In the reference lists, type = "B" is the pipeline's internal code for a ranked list with direction (the paper's type A).
## Keys in the outputs: AcuteStress = Forced swim, DAM_PS2APP = PS2APP.
## Inputs (data/<accession>/raw/, downloaded on 2026-10-05; sources in README.md, Input data):
##   data/GSE169505/raw/GSE169505_Stressome_GeneCounts.txt.gz, GSE169505-GPL21103_series_matrix.txt.gz, GSE169505-GPL24247_series_matrix.txt.gz
##   data/GSE75431/raw/E-GEOD-75431-raw-counts.tsv
##   data/DeDonno2026/provided/DeDonno_GR_glucocorticoid_261016.csv, data/DeDonno2026/output/references_dedonno.rds (rule check only)
##   data/Hammond2019/output/genes_universe.rds, data/Hammond2019/output/references/references_panel9.rds (overlap with DAM, printed)
## Outputs:
##   data/GSE169505/output/acutestress_GSE169505_FC.csv    fold-change table of the Forced swim reference (weights of ULM / WMEAN)
##   data/GSE75431/output/dam_PS2APP_FC.csv                fold-change table of the PS2APP reference
##   figs/Fig5_benchmark/output/extra_refs/fig5_extra_refs.rds   list(AcuteStress, DAM_PS2APP), each list(up, down, type, eligible, desc)
## Run from the project root (a few seconds):  R --no-save --no-restore -f figs/Fig5_benchmark/analysis/fig5_extra_refs.R
pdf(NULL)
suppressPackageStartupMessages(library(data.table))
TOPN <- 500L
## excluded from every reference: genes without informative names, mitochondrial and haemoglobin genes; ribosomal genes kept
bad_pc <- function(g) grepl("Rik[0-9]?$",g)|grepl("^Gm[0-9]",g)|grepl("^mt-",g)|grepl("^Hb[ab]-",g)|grepl("^A[CFLI][0-9]",g)|grepl("^[0-9]",g)|grepl("\\.[0-9]+$",g)|is.na(g)
to_lin <- function(l2) sign(l2) * 2^abs(l2)
OUTD <- "figs/Fig5_benchmark/output/extra_refs"
panel <- list()

## ================================ 1. Forced swim (GSE169505) ================================
## sample table from the two series-matrix files: one row per GEO sample with its characteristics; lib = the library name in the
## supplementary file name 'GSMxxxx_<lib>.h5.gz' (= the column names of the count table)
read_series <- function(f) {
  L <- readLines(gzfile(f))
  fld <- function(tag) lapply(grep(paste0("^!", tag, "\t"), L, value = TRUE), function(x) gsub('"', "", strsplit(x, "\t")[[1]][-1]))
  gsm <- fld("Sample_geo_accession")[[1]]; ttl <- fld("Sample_title")[[1]]; sup <- fld("Sample_supplementary_file_1")[[1]]
  d <- data.table(gsm = gsm, title = ttl, lib = sub("\\.h5\\.gz$", "", sub("^GSM[0-9]+_", "", basename(sup))))
  for (ch in fld("Sample_characteristics_ch1")) {                 # "key: value" per sample; a key may sit on different lines
    k <- sub(":.*$", "", ch); v <- sub("^[^:]*:\\s*", "", ch)
    for (key in unique(k[nzchar(k)])) { if (!key %in% names(d)) d[, (key) := NA_character_]; i <- which(k == key); set(d, i, key, v[i]) }
  }
  d
}
S <- rbindlist(lapply(c("data/GSE169505/raw/GSE169505-GPL21103_series_matrix.txt.gz", "data/GSE169505/raw/GSE169505-GPL24247_series_matrix.txt.gz"),
                      read_series), fill = TRUE)
S <- S[order(as.integer(sub("^GSM", "", gsm)))]
S <- S[experiment == "TimeSeriesFST" & region == "vHC" & Sex == "Male"]
sw <- S[condition == "Swim" & timepoint == "3h"]$lib; hc <- S[condition == "Homecage" & timepoint == "0min"]$lib
## the libraries of the reference, in the order in which their counts are summed (the order fixes the last bits of the CPM values)
SW <- paste0("20180806.A-", c(9, 19, 35, 39, 57, 67, 83))         # 3 h after the swim (7 mice)
HC <- paste0("20180806.A-", c(1, 23, 25, 43, 49, 71, 73, 87))     # home cage (8 mice)
stopifnot(identical(sw, SW), identical(hc, HC))
cat(sprintf("GSE169505: %d samples in the series matrices; TimeSeriesFST, vHC, male: swim 3 h %d, home cage %d\n",
            nrow(S), length(sw), length(hc)))
## the count table has no header for its first column (the gene symbols): fread warns and names it; it is renamed here
C <- fread("data/GSE169505/raw/GSE169505_Stressome_GeneCounts.txt.gz"); setnames(C, 1, "Symbol")
stopifnot(all(c(sw, hc) %in% names(C)))
nucC <- !grepl("^mt-", C$Symbol)                                    # CPM over the nuclear genes (mt- genes left out of the library size)
pb <- function(cols) { s <- rowSums(as.matrix(C[, ..cols])); s / sum(s[nucC]) * 1e6 }
l2 <- log2(pb(sw) + 1) - log2(pb(hc) + 1)
d <- data.table(Symbol = C$Symbol, `Fold Change` = round(to_lin(l2), 4), log2FC = round(l2, 4))[is.finite(`Fold Change`) & `Fold Change` != 0]
r <- d[!bad_pc(Symbol)][order(-abs(`Fold Change`))][!duplicated(Symbol)]
up <- head(r[`Fold Change` > 0][order(-`Fold Change`)]$Symbol, TOPN); dn <- head(r[`Fold Change` < 0][order(`Fold Change`)]$Symbol, TOPN)
dir.create("data/GSE169505/output", showWarnings = FALSE, recursive = TRUE); fwrite(r[order(-`Fold Change`)], "data/GSE169505/output/acutestress_GSE169505_FC.csv")
## eligible genes: the rule of the Glucocorticoid reference (all genes of the fold-change table without bad_pc); checked by applying it
## to the Glucocorticoid table, which must give the eligible genes stored with that reference
g0 <- fread("data/DeDonno2026/provided/DeDonno_GR_glucocorticoid_261016.csv")[is.finite(`Fold Change`) & `Fold Change` != 0]
e0 <- readRDS("data/DeDonno2026/output/references_dedonno.rds")$CORT$eligible
cat(sprintf("eligible rule check on the Glucocorticoid table: all genes %d, without bad_pc %d, stored %d\n", length(unique(g0$Symbol)),
            length(unique(g0[!bad_pc(Symbol)]$Symbol)), length(e0)))
if (!setequal(unique(g0[!bad_pc(Symbol)]$Symbol), e0)) stop("the eligible rule does not reproduce the Glucocorticoid reference", call. = FALSE)
elig <- unique(r$Symbol)
panel$AcuteStress <- list(up = up, down = dn, type = "B", eligible = elig,
                          desc = "Acute stress: ventral hippocampus 3 h after a 6-min swim (7) vs home cage (8), male (GSE169505, bulk); high = stress")
cat(sprintf("[AcuteStress] up %d / down %d (eligible %d; fold-change table %d rows); UP: %s | DOWN: %s\n", length(up), length(dn), length(elig),
            nrow(r), paste(head(up, 8), collapse = ", "), paste(head(dn, 8), collapse = ", ")))
gr <- c("Fkbp5", "Tsc22d3", "Sgk1", "Ddit4", "Zbtb16", "Klf9", "Per1", "Sult1a1", "Fos", "Arc", "Npas4", "Egr1")
cat("known stress / GR genes in the Acute stress reference (rank among up / down):",
    paste(sprintf("%s %s", gr, ifelse(gr %in% up, paste0("up", match(gr, up)), ifelse(gr %in% dn, paste0("down", match(gr, dn)), "-"))), collapse = "; "), "\n")

## ================================ 2. PS2APP (E-GEOD-75431 = GSE75431) ================================
univ <- readRDS("data/Hammond2019/output/genes_universe.rds")
TEST <- c("SRR2961958", "SRR2961957", "SRR2961970", "SRR2961953", "SRR2961939")   # PS2APP, microglia, 13 months
REFG <- c("SRR2961950", "SRR2961935", "SRR2961975", "SRR2961967", "SRR2961982")   # non-transgenic, microglia, 13 months
X <- fread("data/GSE75431/raw/E-GEOD-75431-raw-counts.tsv"); stopifnot(all(c(TEST, REFG) %in% names(X)))
mT <- rowSums(as.matrix(X[, ..TEST])); mR <- rowSums(as.matrix(X[, ..REFG]))   # counts summed per group
nuc <- !grepl("^mt-", X$`Gene Name`)
cat(sprintf("share of the reads in mt- genes: PS2APP %.2f, non-Tg %.2f\n", sum(mT[!nuc]) / sum(mT), sum(mR[!nuc]) / sum(mR)))
l2_all <- log2(mT / sum(mT) * 1e6 + 1) - log2(mR / sum(mR) * 1e6 + 1)        # all-gene CPM (printed for comparison only)
mT <- mT / sum(mT[nuc]) * 1e6; mR <- mR / sum(mR[nuc]) * 1e6                 # CPM over the nuclear genes
l2 <- log2(mT + 1) - log2(mR + 1)
cat(sprintf("median log2 fold change in the Hammond universe: nuclear CPM %.3f, all-gene CPM %.3f\n", median(l2[X$`Gene Name` %in% univ]), median(l2_all[X$`Gene Name` %in% univ])))
cat(sprintf("Spearman between the two versions over the universe genes: %.4f\n", cor(l2[X$`Gene Name` %in% univ], l2_all[X$`Gene Name` %in% univ], method = "spearman")))
d <- data.table(Symbol = X$`Gene Name`, `Fold Change` = round(to_lin(l2), 4), log2FC = round(l2, 4))[is.finite(`Fold Change`) & `Fold Change` != 0]
elig <- unique(d[!bad_pc(Symbol)]$Symbol)
r <- d[!bad_pc(Symbol)][order(-abs(`Fold Change`))][!duplicated(Symbol)]
up <- head(r[`Fold Change` > 0][order(-`Fold Change`)]$Symbol, TOPN); dn <- head(r[`Fold Change` < 0][order(`Fold Change`)]$Symbol, TOPN)
cat(sprintf("[DAM_PS2APP] top %d per direction: up %d / down %d in the Hammond universe\n", TOPN, sum(up %in% univ), sum(dn %in% univ)))
dir.create("data/GSE75431/output", showWarnings = FALSE, recursive = TRUE); fwrite(r[order(-`Fold Change`)], "data/GSE75431/output/dam_PS2APP_FC.csv")
panel$DAM_PS2APP <- list(up = up, down = dn, type = "B", eligible = elig,
                         desc = "DAM-type: PS2APP vs non-transgenic sorted microglia, 13 months (Srinivasan et al. 2016; E-GEOD-75431 / GSE75431, bulk, 5 + 5); high = disease")
k <- readRDS("data/Hammond2019/output/references/references_panel9.rds")$DAM
cat(sprintf("[DAM_PS2APP] up %d / down %d (eligible %d; fold-change table %d rows); overlap with the DAM reference (Keren-Shaul): up %d, down %d\n",
            length(up), length(dn), length(elig), nrow(r), length(intersect(up, k$up)), length(intersect(dn, k$down))))
cat("UP:", paste(head(up, 15), collapse = ", "), "\nDOWN:", paste(head(dn, 15), collapse = ", "), "\n")
core <- c("Apoe", "Lpl", "Spp1", "Cst7", "Itgax", "Clec7a", "Axl", "Trem2", "Tyrobp", "Ctsb", "Cd9", "P2ry12", "Tmem119", "Cx3cr1")
cat("DAM / homeostatic genes (rank among up or down):", paste(sprintf("%s %s", core, ifelse(core %in% up, paste0("up", match(core, up)),
    ifelse(core %in% dn, paste0("down", match(core, dn)), "-"))), collapse = "; "), "\n")

## ================================ output ================================
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
saveRDS(panel, file.path(OUTD, "fig5_extra_refs.rds"))
cat("wrote", file.path(OUTD, "fig5_extra_refs.rds"), ":", paste(names(panel), collapse = ", "), "\n")
