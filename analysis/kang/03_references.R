## Kang step 3 — the four references of Supplementary Fig. S5, by the rule of the paper (analysis/dedonno/03_references.R).
##   IFNb        type A  GSE191317 (Galbraith et al., 2022): PBMC of 7 donors, IFN-beta vs mock, 18 h; bulk 3' mRNA-seq, raw counts
##   PGE2        type A  GSE272019 (Sokhi et al., 2025): CD14+ monocytes of 3 donors, PGE2 vs unstimulated, 3 h; raw counts
##   Interferon  type B  MSigDB Hallmark INTERFERON_ALPHA_RESPONSE (human genes)
##   Chemokine   type B  KEGG CHEMOKINE_SIGNALING_PATHWAY (MSigDB C2, KEGG legacy; human genes)
## Type A: counts summed per group -> CPM over the nuclear genes -> log2(CPM + 1) difference; top 500 up / top 500 down taken
##   from the whole list (restricted to the genes analysed when the cells are scored, 04_score.R). Excluded: genes without
##   informative names, mitochondrial and haemoglobin genes (ribosomal protein genes kept). eligible = genes with a finite,
##   non-zero change, minus the excluded genes. Type B: the whole set (eligible = NULL, any gene).
## Gene names: the Kang matrix has older symbols (e.g. IRG1, FAM26F for ACOD1, CALHM6), so every reference gene is renamed to the
##   Kang symbol of its Ensembl ID where the ID is among the genes analysed (01_qc.R: gene_ensembl.csv) and keeps its own symbol
##   otherwise; the exclusions use the reference's own symbol. Genes that then share a name keep the largest |change|.
## Output: data/Kang2018/output/{IFNb_GSE191317.csv, PGE2_GSE272019.csv, references_kang.rds}
## Run from the project root: R -f analysis/kang/03_references.R
suppressPackageStartupMessages({ library(Matrix); library(data.table); library(msigdbr) })
OUT <- "data/Kang2018/output"; TOPN <- 500L
bad_hs <- function(g) is.na(g) | g == "" | grepl("^MT-", g) | grepl("^HB[ABDEGMQZ][0-9]?$", g) | grepl("\\.[0-9]+$", g) |
  grepl("^LOC[0-9]", g) | grepl("^[0-9]", g) | grepl("^ENSG", g)

ens <- fread(file.path(OUT, "gene_ensembl.csv"))                    # Ensembl ID -> Kang symbol (genes analysed)
kang_sym <- function(id, sym) { k <- ens$symbol[match(sub("\\..*$", "", id), ens$ens)]; ifelse(is.na(k), sym, k) }
## one type A list: id, sym = Ensembl IDs and the reference's own symbols; case, ctrl = count matrices of the two groups
type_a <- function(id, sym, case, ctrl) {
  nuc <- !grepl("^MT-", ifelse(is.na(sym), "", sym)); cpm <- function(x) { s <- rowSums(x); s / sum(s[nuc]) * 1e6 }
  e <- cpm(case); c0 <- cpm(ctrl); l2 <- log2(e + 1) - log2(c0 + 1)
  d <- data.table(Symbol = kang_sym(id, sym), Symbol_ref = sym, `Fold Change` = round(sign(l2) * 2^abs(l2), 4), log2FC = round(l2, 4),
                  caseCPM = round(e, 3), ctrlCPM = round(c0, 3))[is.finite(`Fold Change`) & `Fold Change` != 0]
  d <- d[!bad_hs(Symbol_ref)][order(-abs(log2FC))][!duplicated(Symbol)]; setorder(d, -log2FC)
  cat(sprintf("  %d genes; renamed to the Kang symbol of their Ensembl ID: %d (e.g. %s)\n", nrow(d), d[Symbol != Symbol_ref, .N],
              paste(head(d[Symbol != Symbol_ref][order(-abs(log2FC))][, paste0(Symbol_ref, "->", Symbol)], 6), collapse = ", ")))
  d
}
ref_a <- function(d, desc) list(up = head(d[log2FC > 0][order(-log2FC, Symbol)]$Symbol, TOPN),
                                down = head(d[log2FC < 0][order(log2FC, Symbol)]$Symbol, TOPN),
                                type = "B", eligible = unique(d$Symbol), desc = desc)   # type "B" in the code = type A of the paper

## IFNb: GSE191317, mock (D<n>m) vs IFN-beta (D<n>b) of the same 7 donors
M <- fread("data/GSE191317/raw/GSE191317_Raw_gene_counts_matrix.txt.gz")
mock <- grep("^D[0-9]+m$", names(M), value = TRUE); ifnb <- grep("^D[0-9]+b$", names(M), value = TRUE)
stopifnot(length(mock) == 7, length(ifnb) == 7, identical(sub("m$", "", mock), sub("b$", "", ifnb)))
cat(sprintf("IFNb (GSE191317): %s vs %s\n", paste(ifnb, collapse = ","), paste(mock, collapse = ",")))
dI <- type_a(M$Gene_ID, M$Symbol, as.matrix(M[, ..ifnb]), as.matrix(M[, ..mock]))
fwrite(dI, file.path(OUT, "IFNb_GSE191317.csv"))

## PGE2: GSE272019, unstimulated (CTRL_3h) vs PGE2 (PGE_3h) monocytes of the same 3 donors
E <- fread("data/GSE272019/raw/GSE272019_rawcounts_PGE.txt.gz")
ca <- grep("_PGE_3h_", names(E), value = TRUE); co <- grep("_CTRL_3h_", names(E), value = TRUE)
stopifnot(length(ca) == 3, length(co) == 3, identical(sub("_PGE_3h_", "_", ca), sub("_CTRL_3h_", "_", co)))
cat(sprintf("PGE2 (GSE272019): %s vs %s\n", paste(ca, collapse = ","), paste(co, collapse = ",")))
dP <- type_a(E$gene_id, E$gene_name, as.matrix(E[, ..ca]), as.matrix(E[, ..co]))
fwrite(dP, file.path(OUT, "PGE2_GSE272019.csv"))

## gene sets, renamed to the Kang symbols by Ensembl ID as above
h <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
k <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"))
gset <- function(tab, nm) { s <- tab[gs_name == nm]; stopifnot(nrow(s) > 0); sort(unique(kang_sym(s$ensembl_gene, s$gene_symbol))) }
cat(sprintf("msigdbr %s (%s)\n", as.character(packageVersion("msigdbr")), unique(h$db_version)))

refs <- list(
  IFNb       = ref_a(dI, "IFN-beta vs mock, PBMC of 7 donors, 18 h (GSE191317, bulk)"),
  PGE2       = ref_a(dP, "PGE2 vs unstimulated, CD14+ monocytes of 3 donors, 3 h (GSE272019, bulk)"),
  Interferon = list(up = gset(h, "HALLMARK_INTERFERON_ALPHA_RESPONSE"), down = character(0), type = "C", eligible = NULL,
                    desc = "MSigDB Hallmark INTERFERON_ALPHA_RESPONSE"),
  Chemokine  = list(up = gset(k, "KEGG_CHEMOKINE_SIGNALING_PATHWAY"), down = character(0), type = "C", eligible = NULL,
                    desc = "KEGG CHEMOKINE_SIGNALING_PATHWAY (MSigDB C2, KEGG legacy)"))
univ <- rownames(readRDS(file.path(OUT, "counts_universe.rds"))$counts)
for (a in names(refs)) cat(sprintf("  [%-10s] up %d / down %d (among the genes analysed %d / %d)\n", a, length(refs[[a]]$up),
                                   length(refs[[a]]$down), sum(refs[[a]]$up %in% univ), sum(refs[[a]]$down %in% univ)))
saveRDS(refs, file.path(OUT, "references_kang.rds"))
