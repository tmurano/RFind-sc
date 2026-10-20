## Markers of the eight PBMC cell types of Kang et al. (2018) for Supplementary Fig. S5b.
## Kang et al. annotated their cells with the Seurat guided clustering tutorial on the 2.7k PBMC data of 10x Genomics (Zheng et
## al., 2017): FindAllMarkers gave the top 20 markers of each of the eight cell types, and each cell was assigned to the type whose
## mean raw counts of these markers had the highest Spearman correlation with the cell. The markers are not listed in the paper.
## This script repeats the tutorial (Seurat v2.4 version, compiled 2018-03-22; satijalab.org/seurat/archive/v2.4/pbmc3k_tutorial.html)
## with its parameters and computes FindAllMarkers as Seurat v2 did (natural-log avg_logFC = log(mean(expm1(x)) + 1) of each group;
## pct = fraction of cells with expression > 0; genes with max pct > 0.25 and avg_logFC > 0.25; Wilcoxon rank-sum test; P < 0.01).
## The result is checked against the top two markers of each cluster printed in the tutorial (16 genes with avg_logFC, pct.1 and
## pct.2), and the top three by avg_logFC are used in Supplementary Fig. S5b.
## Input : data/pbmc3k/raw/filtered_gene_bc_matrices/hg19/ (10x Genomics, pbmc3k_filtered_gene_bc_matrices.tar.gz)
## Output: figs/Supp_S5_Kang/output/{figS5_pbmc3k_markers.csv, figS5_pbmc3k_check.csv}
## Run from the project root: R -f figs/Supp_S5_Kang/analysis/figS5_pbmc3k_markers.R
pdf(NULL)
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table) })
OUT <- "figs/Supp_S5_Kang/output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

raw <- Read10X("data/pbmc3k/raw/filtered_gene_bc_matrices/hg19/")             # gene symbols made unique, as in v2
## CreateSeuratObject(min.cells = 3, min.genes = 200) of v2: genes detected in >= 3 cells, then cells with > 200 genes
raw <- raw[Matrix::rowSums(raw > 0) >= 3, ]; raw <- raw[, Matrix::colSums(raw > 0) > 200]
mito <- grep("^MT-", rownames(raw), value = TRUE)
pm <- Matrix::colSums(raw[mito, ]) / Matrix::colSums(raw); ng <- Matrix::colSums(raw > 0)
## FilterCells(nGene in (200, 2500), percent.mito in (-Inf, 0.05))
keep <- ng > 200 & ng < 2500 & pm > -Inf & pm < 0.05
pbmc <- CreateSeuratObject(counts = raw[, keep], project = "10X_PBMC")
pbmc$percent.mito <- pm[keep]
cat(sprintf("genes %d, cells %d (tutorial v2.4: 2,638 cells)\n", nrow(pbmc), ncol(pbmc)))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
## FindVariableGenes(ExpMean, LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5) = "mean.var.plot" in v3+
pbmc <- FindVariableFeatures(pbmc, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf),
                             verbose = FALSE)
cat(sprintf("variable genes %d (tutorial v2.4: 1,838)\n", length(VariableFeatures(pbmc))))
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = FALSE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 20, seed.use = 42, verbose = FALSE)
## shared nearest neighbours as in v2 (exact k-nearest neighbours with RANN, k = 30, prune 1/15); Louvain, resolution 0.6
pbmc <- FindNeighbors(pbmc, dims = 1:10, k.param = 30, prune.SNN = 1/15, nn.method = "rann", verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.6, algorithm = 1, n.start = 100, n.iter = 10, modularity.fxn = 1, random.seed = 0, verbose = FALSE)
cl <- as.character(Idents(pbmc)); print(table(cl))

## v2 FindAllMarkers(only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25), test "wilcox", return.thresh = 0.01
X <- GetAssayData(pbmc, layer = "data"); E <- X; E@x <- expm1(E@x)                 # expm1 of the log-normalized data
find_markers <- function(k) {
  i1 <- cl == k; i2 <- !i1
  pct1 <- round(Matrix::rowMeans(X[, i1] > 0), 3); pct2 <- round(Matrix::rowMeans(X[, i2] > 0), 3)
  lfc <- log(Matrix::rowMeans(E[, i1]) + 1) - log(Matrix::rowMeans(E[, i2]) + 1)
  g <- names(which(pmax(pct1, pct2) > 0.25 & lfc > 0.25))
  p <- vapply(g, function(z) { x <- X[z, ]; suppressWarnings(wilcox.test(x[i1], x[i2])$p.value) }, 0)
  data.table(cluster = k, gene = g, p_val = p, avg_logFC = lfc[g], pct.1 = pct1[g], pct.2 = pct2[g])[p_val < 0.01]
}
mk <- rbindlist(lapply(sort(unique(cl)), find_markers))
mk[, p_val_adj := pmin(1, p_val * nrow(pbmc))]

## clusters named as in the tutorial by its canonical markers (mean log-normalized expression of the marker set, highest cluster)
canon <- list(`CD4 T cells` = "IL7R", `CD14+ Monocytes` = c("CD14", "LYZ"), `B cells` = "MS4A1", `CD8 T cells` = "CD8A",
              `FCGR3A+ Monocytes` = c("FCGR3A", "MS4A7"), `NK cells` = c("GNLY", "NKG7"), `Dendritic cells` = c("FCER1A", "CST3"),
              Megakaryocytes = "PPBP")
avg <- sapply(sort(unique(cl)), function(k) Matrix::rowMeans(X[, cl == k, drop = FALSE]))
sc <- sapply(canon, function(gs) colMeans(avg[gs, , drop = FALSE]))     # clusters x types
## assign each type its best cluster, most specific types first (unique assignment)
ord <- names(sort(apply(sc, 2, function(v) diff(sort(v, decreasing = TRUE)[2:1])), decreasing = FALSE))
lab <- setNames(rep(NA_character_, nrow(sc)), rownames(sc))
for (t in ord) { free <- setdiff(rownames(sc), names(lab)[!is.na(lab)]); lab[free[which.max(sc[free, t])]] <- t }
print(lab)
mk[, cell_type := lab[cluster]]
setorder(mk, cell_type, -avg_logFC); mk[, rank := seq_len(.N), by = cell_type]
fwrite(mk, file.path(OUT, "figS5_pbmc3k_markers.csv"))

## check against the tutorial's printed top two per cluster (cluster ids of the tutorial: 0 CD4 T ... 7 Megakaryocytes)
tut <- data.table(cell_type = rep(c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells",
                                    "Dendritic cells", "Megakaryocytes"), each = 2),
  gene_tut = c("LDHB", "IL7R", "S100A9", "S100A8", "CD79A", "TCL1A", "CCL5", "GZMK", "FCGR3A", "LST1", "GNLY", "GZMB", "FCER1A", "HLA-",
               "PPBP", "PF4"),
  lfc_tut = c(1.15, 1.07, 3.83, 3.79, 2.98, 2.49, 2.16, 2.11, 2.28, 2.15, 3.76, 3.33, 2.73, 1.97, 5.94, 5.02),
  pct1_tut = c(0.924, 0.662, 0.996, 0.973, 0.936, 0.623, 0.974, 0.588, 0.962, 1.00, 0.961, 0.955, 0.844, 1.00, 1.00, 1.00),
  pct2_tut = c(0.482, 0.202, 0.216, 0.122, 0.042, 0.022, 0.230, 0.050, 0.137, 0.316, 0.131, 0.068, 0.011, 0.513, 0.024, 0.010))
tut[, rank := rep(1:2, 8)]
chk <- merge(tut, mk[rank <= 2, .(cell_type, rank, gene, avg_logFC = round(avg_logFC, 2), pct.1, pct.2)], by = c("cell_type", "rank"))
chk[, same_gene := gene == gene_tut | (gene_tut == "HLA-" & startsWith(gene, "HLA-"))]
fwrite(chk, file.path(OUT, "figS5_pbmc3k_check.csv")); print(chk)
cat(sprintf("top-2 genes reproduced: %d of 16; max |avg_logFC difference| %.3f\n", sum(chk$same_gene),
            max(abs(chk$avg_logFC - chk$lfc_tut))))
cat("\nTop three markers per cell type (avg_logFC):\n"); print(mk[rank <= 3, .(cell_type, rank, gene, avg_logFC = round(avg_logFC, 2), pct.1, pct.2)])
