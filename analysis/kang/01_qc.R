## Kang step 1 — cells, author cell-type annotation, gene universe (human PBMC, IFN-beta 6 h vs unstimulated; Supplementary Fig. S5).
## Data: Kang et al. 2018 Nat Biotechnol, GEO GSE96583, batch 2 (8 lupus donors pooled; GSM2560248 = ctrl, GSM2560249 = stim).
## Cells: all the authors' demuxlet singlets with an author cell type (24,673; the cells of the authors' released t-SNE,
##   Supplementary Fig. S5a); donor = "ind". No further cell QC: doublets are the authors' demuxlet calls, and the released
##   matrix has no mitochondrial UMIs (the 13 MT- genes are zero in every cell), so a mitochondrial fraction cannot be used.
## Genes (the rule of the paper, analysis/dedonno/01_qc_embed.R): symbols summed over duplicated Ensembl IDs; universe = detected
##   in >= 1% of the QC'd cells and relative SE of the mean UMI per cell < 5%, minus mitochondrial genes and genes without
##   informative names (clone names; ribosomal protein genes kept). The Ensembl IDs of each symbol are kept
##   (gene_ensembl.csv) so that references with newer symbols can be matched by Ensembl ID (03_references.R).
## Run from the project root: R -f analysis/kang/01_qc.R
## Output: data/Kang2018/output/{cell_meta_qc.csv, counts_universe.rds, gene_universe.csv, gene_ensembl.csv}
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
pdf(NULL)
RAW <- "data/Kang2018/raw"; OUT <- "data/Kang2018/output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

gn <- fread(file.path(RAW, "GSE96583_batch2.genes.tsv.gz"), header = FALSE, col.names = c("ens", "symbol"))
rd <- function(f) { m <- as(readMM(file.path(RAW, f)), "CsparseMatrix"); stopifnot(nrow(m) == nrow(gn)); m }
M1 <- rd("GSM2560248_2.1.mtx.gz"); M2 <- rd("GSM2560249_2.2.mtx.gz")
b1 <- fread(file.path(RAW, "GSM2560248_barcodes.tsv.gz"), header = FALSE)$V1
b2 <- fread(file.path(RAW, "GSM2560249_barcodes.tsv.gz"), header = FALSE)$V1
meta <- fread(file.path(RAW, "GSE96583_batch2.total.tsne.df.tsv.gz"), header = FALSE, skip = 1,
              col.names = c("cell_id", "tsne1", "tsne2", "donor", "condition", "cluster", "cell_type", "multiplets"))
## the authors renamed the stim barcodes that also occur in ctrl with the suffix "-11"; their matrix column has "-1"
meta[, bc := ifelse(condition == "stim", sub("-11$", "-1", cell_id), cell_id)]
meta[, col := ifelse(condition == "ctrl", match(bc, b1), length(b1) + match(bc, b2))]
stopifnot(!anyNA(meta$col), !anyDuplicated(meta$col))
M <- cbind(M1, M2)[, meta$col]; rm(M1, M2); invisible(gc())
cat(sprintf("batch 2: %d barcodes (%s)\n", nrow(meta), paste(names(table(meta$multiplets)), table(meta$multiplets), collapse = ", ")))

keep0 <- meta$multiplets == "singlet" & !is.na(meta$cell_type) & meta$cell_type != "NA"
meta <- meta[keep0]; M <- M[, keep0]
## duplicated symbols (several Ensembl IDs): counts summed
sym <- gn$symbol; f <- factor(sym, levels = unique(sym))
M <- as(fac2sparse(f) %*% M, "CsparseMatrix"); rownames(M) <- levels(f)
genes <- rownames(M); tot <- Matrix::colSums(M)

mtg <- grepl("^MT-", genes); mtpct <- as.numeric(Matrix::colSums(M[mtg, , drop = FALSE])) / tot * 100
stopifnot(sum(mtg) == 13, sum(Matrix::colSums(M[mtg, , drop = FALSE])) == 0)   # no mitochondrial UMI in the release
cs <- tot
cat(sprintf("cells: singlets with a cell type %d (all kept); total UMI per cell median %.0f [range %.0f-%.0f]; genes detected median %.0f\n",
            ncol(M), median(cs), min(cs), max(cs), median(Matrix::colSums(M > 0))))
meta[, `:=`(total_umi = cs, n_genes = as.integer(Matrix::colSums(M > 0)), mt_pct = mtpct)]

## cell types with >= 20 cells in each condition
nct <- meta[, .N, by = .(cell_type, condition)]
ok <- nct[, .(ok = all(N >= 20) && .N == 2), by = cell_type][ok == TRUE]$cell_type
cat("cell types dropped (< 20 cells in a condition):", paste(setdiff(unique(meta$cell_type), ok), collapse = ", "), "\n")
k2 <- meta$cell_type %in% ok; meta <- meta[k2]; M <- M[, k2]; cs <- cs[k2]

bad <- grepl("^MT-", genes) | grepl("\\.[0-9]+$", genes) | grepl("^LOC[0-9]", genes) | grepl("^[0-9]", genes) |
  grepl("^RP[0-9]+-", genes) | grepl("^CT[ABCD]-", genes) | grepl("^A[CFLP][0-9]{6}", genes)
n <- ncol(M); det <- Matrix::rowMeans(M > 0); sx <- Matrix::rowSums(M); sx2 <- Matrix::rowSums(M * M); mu <- sx / n
v <- (sx2 - sx^2 / n) / (n - 1); relSE <- ifelse(mu > 0, sqrt(pmax(v, 0)) / sqrt(n) / mu, Inf)
keepg <- det >= 0.01 & is.finite(relSE) & relSE < 0.05 & !bad
cat(sprintf("genes: %d detected in >= 1%% of cells; %d also relSE < 5%%; %d also informative names (universe)\n",
            sum(det >= 0.01), sum(det >= 0.01 & relSE < 0.05), sum(keepg)))
fwrite(data.table(gene = genes, total_UMI = sx, det = round(100 * det, 3), relSE = round(relSE, 4), bad_name = bad, keep = keepg)[order(-keep, relSE)],
       file.path(OUT, "gene_universe.csv"))
fwrite(gn[symbol %in% genes[keepg]], file.path(OUT, "gene_ensembl.csv"))
uni <- which(keepg)
cf <- M[uni, ]
saveRDS(list(counts = cf, cs = cs, cells = meta$cell_id), file.path(OUT, "counts_universe.rds"))
fwrite(meta[, .(cell_id, cell_type, condition, donor, cluster, tsne1, tsne2, total_umi, n_genes, mt_pct)],
       file.path(OUT, "cell_meta_qc.csv"))
cat(sprintf("QC'd %d cells x %d universe genes (HBB/HBA in universe: %s)\n", ncol(cf), nrow(cf),
            paste(intersect(c("HBB","HBA1","HBA2"), rownames(cf)), collapse = ",")))
print(dcast(meta[, .N, by = .(cell_type, condition)], cell_type ~ condition, value.var = "N"))
print(dcast(meta[, .N, by = .(donor, condition)], donor ~ condition, value.var = "N"))
print(meta[, .(median_umi = median(total_umi), median_genes = median(n_genes)), by = .(condition)])
