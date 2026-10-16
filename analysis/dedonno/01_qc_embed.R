## De Donno step 1 — QC, author cell-type annotation, gene universe, UMAP.
## Cell metadata (condition / line / knockout / donor_id) is read straight from the h5ad obs.
## Only the wild-type mice (line == "Ctrl", 12 mice) are kept, before QC; the GR / MR knock-out lines enter no step.
##
## QC: mt% = mt-genes/total, cap = min(5-MAD upper, 5%); log10-UMI two-sided 5-MAD
##   (doublets pre-removed upstream). This mt cap is stricter than the authors' filtering, so fewer cells are
##   kept than in the published object (see log) and the author types separate less cleanly on our UMAP.
##   Cell types = author `cell_type_short` with numeric suffixes collapsed ("CA1 Glut 2" -> "CA1 Glut");
##   types with <20 cells dropped (in the wild-type cells: Choroid Plexus, 11 cells); the other 21 author types are kept as they
##   are (no further grouping).
## Gene universe: detection>=1%, relSE<0.05, protein-coding.
## Embedding: HVG 2000 (binned dispersion z) -> PCA 50 (scaled) -> UMAP (n_neighbors=30, min_dist=0.3).
##
## Input : data/DeDonno2026/raw/dedonno2026.h5ad
## Output: data/DeDonno2026/output/{cell_meta_qc.csv, counts_universe.rds}
## Run   : from project root (see analysis/dedonno/run_all.sh)

suppressPackageStartupMessages({ library(hdf5r); library(Matrix); library(data.table); library(irlba); library(uwot) })
set.seed(1)
H5 <- "data/DeDonno2026/raw/dedonno2026.h5ad"; OUT <- "data/DeDonno2026/output"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

f <- H5File$new(H5, "r")
rc <- function(p){ o<-f[[p]]; if(inherits(o,"H5Group")){c<-o[["categories"]]$read(); cd<-o[["codes"]]$read(); c[cd+1L]} else o$read() }
cat("reading raw/X...\n")
indptr <- f[["raw/X/indptr"]]$read(); xdata <- f[["raw/X/data"]]$read(); xidx <- f[["raw/X/indices"]]$read()
genes <- rc("raw/var/feature_name")
cell_id <- tryCatch(rc("obs/_index"), error=function(e) rc("obs/index"))
obs <- data.table(cell_id = cell_id, condition = rc("obs/condition"), line = rc("obs/line"),
                  knockout = rc("obs/knockout"), donor_id = rc("obs/donor_id"),
                  cell_type_short = rc("obs/cell_type_short"))
f$close_all()
ng <- length(genes); nc <- length(indptr)-1L
M <- new("dgCMatrix", i=as.integer(xidx), p=as.integer(indptr), x=as.numeric(xdata), Dim=c(ng,nc))
rownames(M) <- genes; colnames(M) <- cell_id; rm(xdata,xidx); gc()
cat(sprintf("h5ad: %d genes x %d cells, %d donors\n", ng, nc, uniqueN(obs$donor_id)))
## All analyses of these data use only the wild-type mice (line == "Ctrl"; GR/MR flox/flox without Cre), so that no
## knock-out cell enters the QC thresholds, the genes analysed, the embedding or the per-cell-type baseline (02_signature.R).
wt <- obs$line == "Ctrl"
M <- M[, wt]; obs <- obs[wt]; nc <- ncol(M); gc()
cat(sprintf("wild-type mice only: %d cells of %d mice kept\n", nc, uniqueN(obs$donor_id)))

## --- QC (deterministic) ---
mth <- function(x,k=5){ m<-median(x); s<-mad(x); c(m-k*s, m+k*s) }
mtg <- grepl("^mt-", genes, ignore.case=TRUE); tot <- Matrix::colSums(M)
mtpct <- as.numeric(Matrix::colSums(M[mtg,,drop=FALSE]))/tot*100
mt_hi <- min(mth(mtpct)[2], 5); lu <- log10(tot); ut <- mth(lu)
keep <- !(mtpct>mt_hi | lu<ut[1] | lu>ut[2])
cat(sprintf("QC: mt%% cap=%.2f, log10-UMI [%.3f, %.3f] -> keep %d/%d\n", mt_hi, ut[1], ut[2], sum(keep), nc))
M <- M[,keep]; cs <- tot[keep]; ob <- obs[keep]
ob[, celltype_fine := gsub(" [0-9]+$", "", cell_type_short)]    # collapse author number suffixes
tb <- table(ob$celltype_fine); keep2 <- ob$celltype_fine %in% names(tb)[tb>=20]
M <- M[,keep2]; cs <- cs[keep2]; ob <- ob[keep2]; n <- ncol(M)
## Cell type = the authors' type with its number suffix removed (21 types in the wild-type cells).
## The types are not grouped further, to stay close to the authors' definitions. The per-cell-type
## baseline in 02_signature.R is built on this column.
ob[, cell_type := celltype_fine]

## --- gene universe ---
det <- rowSums(M>0)/n; sx <- rowSums(M); sx2 <- rowSums(M*M); mu <- sx/n; v <- (sx2-sx^2/n)/(n-1)
relSE <- ifelse(mu>0, sqrt(pmax(v,0))/sqrt(n)/mu, Inf)
bad <- grepl("Rik[0-9]?$",genes)|grepl("^Gm[0-9]",genes)|grepl("^mt-",genes)|grepl("^A[CFLI][0-9]",genes)|grepl("^[0-9]",genes)|grepl("Rik$",genes)
gk <- det>=0.01 & is.finite(relSE) & relSE<0.05 & !bad; cf <- M[gk,]; G <- nrow(cf); rm(M); gc()
cat(sprintf("QC'd %d cells x %d universe genes, %d fine / %d grouped cell types\n",
            n, G, uniqueN(ob$celltype_fine), uniqueN(ob$cell_type)))

## --- normalize, HVG, PCA, UMAP ---
cpm <- cf %*% Diagonal(x=1e6/cs); lc <- cpm; lc@x <- log1p(lc@x)
mug <- Matrix::rowMeans(cpm); v2 <- (Matrix::rowMeans(cpm*cpm)-mug^2)*n/(n-1); disp <- ifelse(mug>0, v2/mug, 0)
lmn <- log1p(mug); ld <- log(disp+1e-8); bins <- cut(lmn, quantile(lmn, seq(0,1,0.05), na.rm=TRUE), include.lowest=TRUE)
zz <- ave(ld, bins, FUN=function(x){ s<-sd(x); if(is.na(s)||s==0) rep(0,length(x)) else (x-mean(x))/s })
hvg <- order(zz, decreasing=TRUE)[1:min(2000,G)]
## UMAP coordinates depend on the uwot/irlba/BLAS versions and threads (r = 0.97 between runs on two computers),
## so the paper uses the coordinates of 2026-10-15 (wild-type cells), kept as a fixed input. Delete that file to recompute the
## embedding.
UMAP_FIX <- "data/DeDonno2026/provided/umap_261015.csv"
if (file.exists(UMAP_FIX)) {
  u <- fread(UMAP_FIX)[match(ob$cell_id, cell_id)]; stopifnot(!anyNA(u$umap_x))
  um <- cbind(u$umap_x, u$umap_y); cat("UMAP: fixed coordinates from", UMAP_FIX, "\n")
} else {
  cat("PCA/UMAP...\n")
  pca <- prcomp_irlba(t(lc[hvg,]), n=50, center=TRUE, scale.=TRUE)
  um <- umap(pca$x, n_neighbors=30, min_dist=0.3, n_threads=4)
}
ob[, `:=`(umap_x = um[,1], umap_y = um[,2], total_umi = cs)]

setcolorder(ob, c("cell_id","cell_type","celltype_fine","cell_type_short","condition","line","knockout","donor_id","total_umi","umap_x","umap_y"))
fwrite(ob, file.path(OUT, "cell_meta_qc.csv"))
saveRDS(list(counts = cf, cs = cs, cells = ob$cell_id), file.path(OUT, "counts_universe.rds"))
cat(sprintf("wrote cell_meta_qc.csv (%d cells) + counts_universe.rds (%d x %d)\n", nrow(ob), G, n))
print(ob[, .N, by = cell_type][order(-N)])
print(ob[line == "Ctrl", .(n_donors = uniqueN(donor_id), n_cells = .N), by = condition])
