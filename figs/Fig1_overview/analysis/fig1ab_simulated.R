## Fig. 1a,b — illustration drawn from SIMULATED single-cell counts whose parameters are estimated from the De Donno et al.
## control-line cells, so that the virtual cells keep the structure of real cell types without copying any real cell.
## Model: for each real cell (wild-type mice, seven main cell types, 22,526 cells; 11,647 genes analysed), its expression share profile is smoothed
## over the cell and its 14 nearest neighbours (50 PCs); one virtual cell is drawn per real cell with a total UMI taken at
## random from real cells of the same type and negative binomial counts around total UMI x smoothed share. No real count
## is copied. The virtual counts are then processed exactly as real data:
##   a  log-normalised counts -> 2,000 variable genes -> 50 PCs -> UMAP (n_neighbors 30, min_dist 1.0; see the fixed input below)
##   b  Delta against the mean share of all virtual cells -> per-cell UP/DOWN (top 500 inside the background of each reference,
##      as in the paper) -> 4-term Running Fisher against
##      three references: Maturity (Matcovitch-Natan 2016), Interferon (MSigDB Hallmark, used whole) and
##      Glucocorticoid (Fig. 3); 8,000 cells. The axes of b carry generic names (Reference A-C-index; fig1b_axes.csv).
## Input : data/DeDonno2026/output/{counts_universe.rds, cell_meta_qc.csv, references_dedonno.rds},
##         data/Matcovitch2016/raw/GSE79812_Timeline.xlsx, data/Hammond2019/output/references/references_panel9.rds (Interferon set),
##         data/DeDonno2026/provided/fig1a_sim_umap_261015.csv (fixed UMAP coordinates of a, from the uwot 0.2.4 run; written by
##         this script when it is absent)
## Output: figs/Fig1_overview/output/{fig1_sim_cells.csv, fig1b_axes.csv, Fig1a_UMAP.png (pasted size)};
##         then run fig1b_3d.py for b
## Run   : <R> -f figs/Fig1_overview/analysis/fig1ab_simulated.R   (from project root)

pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(Matrix); library(irlba); library(uwot)
  library(RANN); library(parallel); library(readxl) })
source("analysis/running_fisher.R")
DD <- "data/DeDonno2026/output"; OUT <- "figs/Fig1_overview/output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
K_NN <- 15L; NB_SIZE <- 10; TOPN <- 500L; CLIP <- 300; NPICK <- 8000L; NCORES <- max(1L, detectCores() - 2L)
MAIN7 <- c("DG Glut", "CA1 Glut", "CA3 Glut", "Astro", "Microglia", "Oligodendrocytes", "Endothelial")

## ---- model estimated from the real cells ----
meta <- fread(file.path(DD, "cell_meta_qc.csv"))[line == "Ctrl" & cell_type %in% MAIN7]
L3 <- readRDS(file.path(DD, "counts_universe.rds"))
idx <- match(meta$cell_id, L3$cells); stopifnot(!anyNA(idx))
Y <- L3$counts[, idx, drop = FALSE]; cs <- L3$cs[idx]; genes <- rownames(Y); rm(L3); gc()
pca_of <- function(counts, depth) {                                   # log(1 + CP10k), 2,000 genes of highest variance/mean, 50 PCs
  lc <- log1p(Matrix::t(Matrix::t(counts) / (depth / 1e4)))
  m <- Matrix::rowMeans(lc); v <- Matrix::rowMeans(lc * lc) - m^2
  hvg <- order(v / (m + 1e-8), decreasing = TRUE)[1:2000]
  prcomp_irlba(t(as.matrix(lc[hvg, ])), n = 50, center = TRUE, scale. = TRUE)$x }
set.seed(11)
nn <- nn2(pca_of(Y, cs), k = K_NN)$nn.idx                            # the cell itself and its 14 nearest neighbours
Pshare <- Y %*% Diagonal(x = 1 / cs)                                 # genes x cells, share of each gene in each cell
set.seed(2026)
Nsim <- ave(cs, meta$cell_type, FUN = function(x) x[sample.int(length(x))])   # total UMI drawn from the same type
chunks <- split(seq_len(ncol(Y)), ceiling(seq_len(ncol(Y)) / 1000))
Ysim <- do.call(cbind, lapply(chunks, function(ii) {
  W <- sparseMatrix(i = rep(seq_along(ii), each = K_NN), j = as.vector(t(nn[ii, , drop = FALSE])), x = 1 / K_NN,
                    dims = c(length(ii), ncol(Y)))
  mu <- as.matrix(Pshare %*% Matrix::t(W)) %*% diag(Nsim[ii], length(ii))
  Matrix(matrix(rnbinom(length(mu), mu = as.vector(mu), size = NB_SIZE), nrow(mu)), sparse = TRUE) }))
rownames(Ysim) <- genes; rm(Pshare, Y); gc()
cs_sim <- Matrix::colSums(Ysim)
cat(sprintf("virtual cells: %d x %d genes; median total UMI %.0f (real %.0f)\n", ncol(Ysim), nrow(Ysim), median(cs_sim), median(cs)))
type <- factor(meta$cell_type)                                        # alphabetical; the order sets the colours

## ---- a: UMAP of the virtual cells ----
## The layout depends on the uwot version and threads (uwot 0.2.4 and 0.2.5 placed cells up to 2 units apart), so the
## coordinates of the uwot 0.2.4 run are a fixed input, as for the UMAP of the data of De Donno et al.
## Delete the file to recompute (it is then written again). The RNG is re-seeded below, so reading the file changes nothing else.
## The file holds the run on the wild-type data analysed here (11,647 genes).
UMAP_FIXED <- "data/DeDonno2026/provided/fig1a_sim_umap_261015.csv"
if (file.exists(UMAP_FIXED)) {
  fx <- fread(UMAP_FIXED)
  stopifnot(identical(fx$cell_id, sprintf("virtual%05d", seq_len(ncol(Ysim)))))
  um <- as.matrix(fx[, .(umap_x, umap_y)]); cat("UMAP coordinates read from", UMAP_FIXED, "\n")
} else {
  set.seed(11)
  um <- umap(pca_of(Ysim, cs_sim), n_neighbors = 30, min_dist = 1.0, spread = 1.0, repulsion_strength = 0.1, n_threads = 4)
  fwrite(data.table(cell_id = sprintf("virtual%05d", seq_len(ncol(Ysim))), umap_x = um[, 1], umap_y = um[, 2]), UMAP_FIXED)
  cat("UMAP coordinates written to", UMAP_FIXED, "\n")
}

## ---- b: indices of 8,000 virtual cells ----
ul <- tolower(genes)
baseline <- Matrix::rowMeans(Ysim %*% Diagonal(x = 1 / cs_sim))        # mean share over all virtual cells
set.seed(1); pick <- sort(sample(ncol(Ysim), NPICK))
## the rule of the paper's references: genes without informative names, mitochondrial and haemoglobin genes excluded;
## the top TOPN genes per direction of the whole list, restricted to the genes analysed below (R)
bad <- function(g) grepl("Rik[0-9]?$", g) | grepl("^Gm[0-9]", g) | grepl("^mt-", g) | grepl("^Hb[ab]-", g) | grepl("^A[CFLI][0-9]", g) |
  grepl("^[0-9]", g) | grepl("\\.[0-9]+$", g) | is.na(g)
signed_top <- function(sym, fc) {                                     # top TOPN up / down genes by FC of the whole list
  d <- data.table(g = sym, fc = fc)[is.finite(fc) & fc != 0 & !bad(g)][, g := tolower(g)]
  d <- d[order(-abs(fc))][!duplicated(g)]
  list(up = head(d[fc > 0][order(-fc)]$g, TOPN), down = head(d[fc < 0][order(fc)]$g, TOPN)) }
mt <- as.data.table(suppressMessages(read_excel("data/Matcovitch2016/raw/GSE79812_Timeline.xlsx", sheet = "Timeline")))
setnames(mt, 1, "Symbol"); cn <- colnames(mt); numc <- function(cols) sapply(mt[, ..cols], as.numeric)
E <- rowMeans(numc(cn[grepl("^mic_E1[346]", cn)]), na.rm = TRUE); A <- rowMeans(numc(cn[grepl("cortex|hipo", cn)]), na.rm = TRUE)
nuc <- !grepl("^mt-", mt$Symbol)                                     # per million over the nuclear genes, as the Maturity reference
E <- E / sum(E[nuc], na.rm = TRUE) * 1e6; A <- A / sum(A[nuc], na.rm = TRUE) * 1e6
dd <- readRDS(file.path(DD, "references_dedonno.rds"))
R <- list(Maturity = signed_top(mt$Symbol, log2(A + 1) - log2(E + 1)),
          Interferon = readRDS("data/Hammond2019/output/references/references_panel9.rds")$Interferon,   # MSigDB Hallmark IFN-alpha response (gene set, stored as type "C"), saved locally
          Glucocorticoid = dd$CORT)
## the matching as in the paper: the background of each reference = the genes analysed with a non-zero change in it
## (pair universe; all genes analysed for the gene set), and each cell's UP / DOWN lists = its top TOPN genes inside that background
mat_el <- unique(data.table(g = mt$Symbol, fc = log2(A + 1) - log2(E + 1))[is.finite(fc) & fc != 0 & !bad(g)]$g)
EL <- list(Maturity = mat_el, Interferon = NULL, Glucocorticoid = dd$CORT$eligible)
UA <- lapply(EL, function(e) pair_universe(genes, e)); INU <- lapply(UA, function(u) ul %in% u)
R <- lapply(setNames(nm = names(R)), function(a) { e <- R[[a]]
  list(up = intersect(tolower(e$up), UA[[a]]), down = intersect(tolower(e$down), UA[[a]]), ranked = !identical(e$type, "C")) })
score <- function(ii) {
  obs <- as.matrix(Ysim[, ii]); expc <- outer(baseline, cs_sim[ii]); d <- (obs - expc) / sqrt(expc + 1)
  t(vapply(seq_along(ii), function(k) { p <- d[, k]
    pos <- which(p > 0); pos <- pos[order(p[pos], decreasing = TRUE)]
    neg <- which(p < 0); neg <- neg[order(p[neg])]
    vapply(names(R), function(a) rfind_score_cell(ul[cut_in_universe(pos, INU[[a]], TOPN)], ul[cut_in_universe(neg, INU[[a]], TOPN)],
                                                 R[[a]], universe_size = length(UA[[a]]), clip = CLIP), numeric(1)) }, numeric(length(R)))) }
M <- do.call(rbind, mclapply(split(pick, ceiling(seq_along(pick) / 200)), score, mc.cores = NCORES))
colnames(M) <- names(R)

## ---- outputs ----
pal <- setNames(scales::hue_pal(l = 60, c = 110)(nlevels(type)), levels(type))
cells <- data.table(cell_id = sprintf("virtual%05d", seq_along(type)), type = as.character(type), umap_x = um[, 1], umap_y = um[, 2],
                    Maturity = NA_real_, Interferon = NA_real_, Glucocorticoid = NA_real_, in_b = FALSE)
cells[pick, `:=`(Maturity = M[, "Maturity"], Interferon = M[, "Interferon"], Glucocorticoid = M[, "Glucocorticoid"], in_b = TRUE)]
cells[, colour := pal[type]]
fwrite(cells, file.path(OUT, "fig1_sim_cells.csv"))
## b shows simulated values, so the axes carry generic names (Reference A-C-index), not the names of the references.
fwrite(data.table(column = names(R), label = paste0("Reference ", LETTERS[seq_along(R)], "-index")),
       file.path(OUT, "fig1b_axes.csv"))
print(cells[in_b == TRUE, lapply(.SD, function(x) round(median(x), 1)), by = type, .SDcols = names(R)])

## ---- a: the panel, drawn at the size it is pasted into Fig. 1 (43.9 x 33.8 mm, the size on the slide) ----
## Arial, axis titles 7 pt, axis lines 0.5 pt; point size 0.8 x 0.2416 x 0.894 (small dots for 22,526 cells on a 43.9-mm map).
## The draw order (sample(.N)) follows the RNG state left by set.seed(1) above.
PT_LW <- 0.5 / (ggplot2::.pt * 72.27 / 96)                            # ggplot linewidth that gives 0.5 pt
pa <- ggplot(cells[sample(.N)], aes(umap_x, umap_y, colour = type)) +
  geom_point(size = 0.8 * 0.2416 * 0.894, alpha = 0.85, stroke = 0) +
  scale_colour_manual(values = pal, guide = "none") +
  coord_fixed(ratio = 0.72) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 7, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 1)), axis.title.y = element_text(margin = margin(r = 1)),
        axis.line = element_line(colour = "grey25", linewidth = PT_LW, lineend = "square"),
        plot.margin = margin(1, 1, 1, 1, unit = "pt"))
OUT_PPTX <- OUT   # pasted-size panels go to output/
px <- function(mm) round(mm / 25.4 * 600)                            # ragg truncates mm sizes, so give whole pixels
ragg::agg_png(file.path(OUT_PPTX, "Fig1a_UMAP.png"), width = px(43.9), height = px(33.8), units = "px", res = 600, background = "white")
print(pa); invisible(dev.off())
cat("wrote Fig1a_UMAP.png, fig1_sim_cells.csv, fig1b_axes.csv\n")
