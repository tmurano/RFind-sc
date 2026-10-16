## =====================================================================
## RFind-sc references of the Hammond data — builds the four references in one run
## Type A of the paper (ranked up / down lists of external fold changes: published for DAM, computed here from the published
##   timeline table for Maturity): DAM (Keren-Shaul et al., 2017) and Maturity (Matcovitch-Natan et al., 2016); stored with
##   type = "B" in the panel.
## Type B of the paper (gene sets): Interferon and G2M (MSigDB Hallmark); stored with type = "C" in the panel.
## Criteria: signed linear FC; p-values are not used; type A = the top 500 genes per direction of the published list (TOPN_B),
##   then those in the Hammond gene universe (8,960 genes); type B = membership.
## Output: one CSV per ranked list (DAM, Maturity) + data/Hammond2019/output/references/references_panel9.rds (all four)
## =====================================================================
suppressPackageStartupMessages({ library(Matrix); library(data.table); library(readxl); library(msigdbr) })
## Paths relative to the project root. Run: Rscript analysis/hammond/06_build_references.R (from the project root)
PIN <- "data/Hammond2019/output"; OUT <- "data/Hammond2019/output/references"; SRC <- "data/KerenShaul2017/provided"
TOPN_A <- 500L                      # not used
## Ranked lists (type A of the paper), one rule for the six references of the paper: the top 500 genes per
## direction are taken from the published list and stored as they are; those in the common universe are used when the cells are
## scored (07_score_axes.R), as for the references of De Donno (03_references.R; Seizure keeps the top 200 up-regulated genes).
## Excluded from every list: genes without informative names, mitochondrial genes and haemoglobin genes (from blood cells in the
## samples); ribosomal protein genes are kept.
TOPN_B <- 500L
bad_pc <- function(g) grepl("Rik[0-9]?$",g)|grepl("^Gm[0-9]",g)|grepl("^mt-",g)|grepl("^Hb[ab]-",g)|grepl("^A[CFLI][0-9]",g)|grepl("^[0-9]",g)|grepl("\\.[0-9]+$",g)|is.na(g)
to_lin <- function(l2fc) sign(l2fc)*2^abs(l2fc)
univ <- readRDS(file.path(PIN,"genes_universe.rds"))
cat(sprintf("common universe: %d genes\n\n", length(univ)))
panel <- list()

## =====================================================================
## Ranked lists, type A of the paper (type = "B" here) — external DEG (DAM: Keren-Shaul / Maturity: Matcovitch-Natan)
## =====================================================================
cat("\n===== Type-B (external DEG) =====\n")

## --- DAM Keren-Shaul 2017 Cell, Table S3 (5xFAD, DAM vs homeostatic, published FC) ---
k <- fread(file.path(SRC,"Mouse_DAM_KerenShaul2017_vs_homeostatic.csv"))[,.(Symbol,`Fold Change`,`P-Value`)]
k <- k[is.finite(`Fold Change`) & `Fold Change`!=0][!bad_pc(Symbol)]
k <- k[order(-abs(`Fold Change`))][!duplicated(Symbol)]
k_up <- head(k[`Fold Change`>0][order(-`Fold Change`)]$Symbol, TOPN_B)
k_dn <- head(k[`Fold Change`<0][order(`Fold Change`)]$Symbol, TOPN_B)
fwrite(k[order(-`Fold Change`), .(Symbol,`Fold Change`,`P-Value`)], file.path(OUT,"dam_KerenShaul.csv"))
## eligible = genes the reference study reported (finite, non-zero FC): with the target universe they form the pair
## universe, the Running Fisher background (analysis/running_fisher.R: pair_universe). Gene sets (type = "C") have eligible = NULL.
panel$DAM <- list(up=k_up, down=k_dn, type="B", eligible=unique(k$Symbol), desc="DAM vs homeostatic (Keren-Shaul 2017 Cell, 5xFAD, Table S3, real signed FC)")
cat(sprintf("[DAM] Keren-Shaul: up %d / down %d (in the universe %d / %d) | UP: %s\n", length(k_up), length(k_dn), sum(k_up %in% univ),
            sum(k_dn %in% univ), paste(head(k_up,5),collapse=", ")))


## --- Maturity Matcovitch-Natan 2016 Science (GSE79812, adult vs embryonic microglia; high = mature) ---
## Oriented so up = ADULT (mature) markers: FC = adult vs embryonic. Score is high in mature microglia
## and tracks real postnatal age (Spearman ~0.82 across the Hammond developmental timeline).
mt <- as.data.table(suppressMessages(read_excel("data/Matcovitch2016/raw/GSE79812_Timeline.xlsx", sheet="Timeline")))
setnames(mt,1,"Symbol"); cn <- colnames(mt)
emb <- cn[grepl("^mic_E1[346]",cn)]; adt <- cn[grepl("cortex|hipo",cn)]
numc <- function(cols) sapply(mt[,..cols], as.numeric)
E <- rowSums(numc(emb),na.rm=TRUE); A <- rowSums(numc(adt),na.rm=TRUE)   # samples summed per group, then scaled to one million
nuc <- !grepl("^mt-", mt$Symbol)                                   # per million over the nuclear genes (mitochondrial genes left out)
cat(sprintf("Maturity table: %d genes, %d mitochondrial (left out of the per-million scaling)\n", nrow(mt), sum(!nuc)))
E <- E/sum(E[nuc],na.rm=TRUE)*1e6; A <- A/sum(A[nuc],na.rm=TRUE)*1e6
l2m <- log2(A+1)-log2(E+1); linm <- to_lin(l2m)                    # adult vs embryonic -> high = mature
rm_dt <- data.table(Symbol=mt$Symbol, `Fold Change`=round(linm,4), log2FC=round(l2m,4))[is.finite(`Fold Change`) & `Fold Change`!=0][!bad_pc(Symbol)]
rm_dt <- rm_dt[order(-abs(`Fold Change`))][!duplicated(Symbol)]
m_up <- head(rm_dt[`Fold Change`>0][order(-`Fold Change`)]$Symbol, TOPN_B)
m_dn <- head(rm_dt[`Fold Change`<0][order(`Fold Change`)]$Symbol, TOPN_B)
## the whole table (not cut to the universe): fig5_bench_redo.R takes the top L per direction from it as here
fwrite(rm_dt[order(-`Fold Change`), .(Symbol,`Fold Change`,log2FC)], file.path(OUT,"maturity_Matcovitch_microglia.csv"))
m_elig <- data.table(Symbol=mt$Symbol, fc=linm)[is.finite(fc) & fc!=0][!bad_pc(Symbol)]$Symbol
panel$Maturity <- list(up=m_up, down=m_dn, type="B", eligible=unique(m_elig), desc="Maturity: adult vs embryonic (E13.5-16.5) microglia (Matcovitch-Natan 2016 Science, GSE79812, bulk); high = mature")
cat(sprintf("[Maturity] Matcovitch: up %d / down %d (in the universe %d / %d) | UP: %s\n", length(m_up), length(m_dn), sum(m_up %in% univ),
            sum(m_dn %in% univ), paste(head(m_up,5),collapse=", ")))

## =====================================================================
## Gene sets, type B of the paper (type = "C" here) — MSigDB Hallmark membership (Interferon / G2M)
## =====================================================================
cat("\n===== Type-C (MSigDB Hallmark membership) =====\n")
h <- as.data.table(msigdbr(species="Mus musculus", collection="H"))
hsets <- h[, .(genes=list(sort(unique(gene_symbol)))), by=gs_name]
get_set <- function(nm) hsets[gs_name==nm]$genes[[1]]
typeC_map <- c(Interferon="HALLMARK_INTERFERON_ALPHA_RESPONSE",
               G2M="HALLMARK_G2M_CHECKPOINT")
for(nm in names(typeC_map)){
  g <- get_set(typeC_map[[nm]]); stopifnot(length(g)>0)
  panel[[nm]] <- list(up=g, down=character(0), type="C", eligible=NULL, desc=sprintf("MSigDB Hallmark membership (%s)", typeC_map[[nm]]))
  cat(sprintf("[%s] %s: %d genes (%d in universe)\n", nm, typeC_map[[nm]], length(g), sum(g %in% univ)))
}

## =====================================================================
## Save the combined panel
## =====================================================================
## Four references: DAM, Maturity, Interferon, G2M.
## File name references_panel9.rds: the name that downstream readers use (RFindsc bundle).
ord <- c("DAM","Maturity","Interferon","G2M")
panel <- panel[ord]
saveRDS(panel, file.path(OUT,"references_panel9.rds"))
cat("\n===== combined 4-axis panel =====\n")
for(a in ord) cat(sprintf("  [%-11s] type%s  up=%5d  down=%5d  eligible=%s\n", a, panel[[a]]$type, length(panel[[a]]$up),
                          length(panel[[a]]$down), if (is.null(panel[[a]]$eligible)) "any gene" else length(panel[[a]]$eligible)))
cat(sprintf("\nwrote data/Hammond2019/output/references/references_panel9.rds (%d axes)\n", length(panel)))
