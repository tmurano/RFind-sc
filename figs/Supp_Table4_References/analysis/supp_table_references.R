## Supplementary Table 4, step 1 — the number of genes of each reference, before and after the restriction to the genes analysed.
## The lists are restricted exactly as in the scoring scripts (analysis/hammond/07_score_axes.R, analysis/dedonno/04_score.R,
## figs/Fig5_benchmark/analysis/fig5_extra_score.R): pair universe = genes analysed of the dataset ∩ genes the reference could contain
## (analysis/running_fisher.R: pair_universe; all genes analysed for gene sets), and the up / down lists of the reference are cut to
## it (gene names matched case-insensitively).
## In the reference files, the internal codes are: type "B" = ranked list with direction (type A of the paper), type "C" = gene set
## (type B of the paper); CORT = Glucocorticoid, KA6h = Seizure, NFkB = NF-kB, AcuteStress = Forced swim, DAM_PS2APP = PS2APP.
## Input : data/Hammond2019/output/genes_universe.rds, percell_signature_top500.rds (the genes analysed, checked to be the same),
##         data/Hammond2019/output/references/references_panel9.rds (DAM, Maturity, Interferon, G2M),
##         data/DeDonno2026/output/percell_signature_celltype.rds (the genes analysed), references_dedonno.rds (Glucocorticoid,
##         Seizure, UPR, NF-kB), figs/Fig5_benchmark/output/extra_refs/fig5_extra_refs.rds (Forced swim, PS2APP)
## Output: figs/Supp_Table4_References/output/SuppTable4_reference_sizes.csv (read by supp_table_references.py)
## Run from the project root (about 10 s):  R --no-save --no-restore -f figs/Supp_Table4_References/analysis/supp_table_references.R
pdf(NULL)
suppressPackageStartupMessages(library(data.table))
source("analysis/running_fisher.R")                                       # pair_universe
OUTD <- "figs/Supp_Table4_References/output"; dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

## genes analysed of each dataset (the genes of the per-cell rankings used for scoring)
gH <- readRDS("data/Hammond2019/output/genes_universe.rds")
stopifnot(identical(sort(gH), sort(readRDS("data/Hammond2019/output/percell_signature_top500.rds")$genes)))
gD <- readRDS("data/DeDonno2026/output/percell_signature_celltype.rds")$genes
GENES <- list(Hammond = gH, DeDonno = gD)

pH <- readRDS("data/Hammond2019/output/references/references_panel9.rds")
pD <- readRDS("data/DeDonno2026/output/references_dedonno.rds")
pX <- readRDS("figs/Fig5_benchmark/output/extra_refs/fig5_extra_refs.rds")
## reference name in the paper -> (reference object, dataset whose cells it scores)
REFS <- list(DAM = list(pH$DAM, "Hammond"), Maturity = list(pH$Maturity, "Hammond"), PS2APP = list(pX$DAM_PS2APP, "Hammond"),
             Interferon = list(pH$Interferon, "Hammond"), G2M = list(pH$G2M, "Hammond"),
             Glucocorticoid = list(pD$CORT, "DeDonno"), `Forced swim` = list(pX$AcuteStress, "DeDonno"), Seizure = list(pD$KA6h, "DeDonno"),
             `NF-kB` = list(pD$NFkB, "DeDonno"), UPR = list(pD$UPR, "DeDonno"))

rows <- lapply(names(REFS), function(nm) {
  e <- REFS[[nm]][[1]]; ds <- REFS[[nm]][[2]]; g <- GENES[[ds]]
  U <- pair_universe(g, e$eligible)
  up <- intersect(tolower(e$up), U); dn <- intersect(tolower(e$down), U)
  ## genes_analysed = rows of the gene table, as in Methods (De Donno: 11,647 rows, 3 symbols twice); background_U counts unique
  ## names, as the scoring does
  data.table(reference = nm, dataset = ds, type_code = e$type, genes_analysed = length(g),
             up_before = length(unique(e$up)), down_before = length(unique(e$down)),
             up_after = length(up), down_after = length(dn), background_U = length(U))
})
d <- rbindlist(rows)
print(d)
fwrite(d, file.path(OUTD, "SuppTable4_reference_sizes.csv"))
cat("wrote", file.path(OUTD, "SuppTable4_reference_sizes.csv"), "\n")
