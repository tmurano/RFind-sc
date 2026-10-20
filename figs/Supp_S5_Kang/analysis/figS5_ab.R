## Supplementary Fig. S5a, b (Kang et al. 2018 PBMC, IFN-beta 6 h versus unstimulated).
## Panel a: the t-SNE of the most recent release of the authors (t-SNE coordinates, singlet calls and cell-type annotation of GEO
##   GSE96583, GSE96583_batch2.total.tsne.df.tsv.gz, released 2018-03-12). The cells are the singlets with a cell type (24,673;
##   the six singlets without a cell type are not drawn). Left, the eight cell types (the key, outside the map, gives the full names
##   with the abbreviations of the original Fig. 3a used in b-f); right, condition (key inside the map at the top left, as in the
##   original). Layout of Fig. 2a of this paper (no frame or ticks, "t-SNE 1/2"); colours as in the original figure (RColorBrewer
##   "Paired" blue and red for the conditions, "Set2" assigned to the cell types in the alphabetical order of the authors' labels);
##   points drawn in the row order of the released table (unstimulated first).
## Panel b: the eight cell types (abbreviations of the original figure) and their markers: the top three by avg_logFC of
##   FindAllMarkers in the Seurat guided clustering tutorial on the 2.7k PBMC data, which the original study used to annotate the
##   cells (figS5_pbmc3k_markers.R, which reproduces the top two printed in the tutorial for all eight types).
## Input : data/Kang2018/raw/GSE96583_batch2.total.tsne.df.tsv.gz, figs/Supp_S5_Kang/output/figS5_pbmc3k_markers.csv
## Output: figs/Supp_S5_Kang/output/{FigS5a_tSNE.png, FigS5b_celltype_table.png, figS5_celltype_table.csv}
## Run from the project root: R -f figs/Supp_S5_Kang/analysis/figS5_ab.R
pdf(NULL)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(grid); library(gridExtra); library(gtable) })
FIG <- "figs/Supp_S5_Kang/output"; OUT <- FIG; dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
FAM <- "Arial"; px <- function(mm) round(mm / 25.4 * 600)

## ---- cell types: authors' labels, abbreviations of the original Fig. 3a, colours, legend order ----
AUTH <- c("B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells")
stopifnot(identical(AUTH, sort(AUTH)))
PAL <- setNames(RColorBrewer::brewer.pal(8, "Set2"), AUTH)                    # Set2 in the alphabetical order of the labels
ABB <- c(`CD14+ Monocytes` = "cM", `FCGR3A+ Monocytes` = "ncM", `Dendritic cells` = "DC", Megakaryocytes = "Mkc",
         `CD4 T cells` = "T[H]", `B cells` = "B", `CD8 T cells` = "Tc", `NK cells` = "NK")   # order of the original key
COND_PAL <- c(ctrl = RColorBrewer::brewer.pal(12, "Paired")[2], stim = RColorBrewer::brewer.pal(12, "Paired")[6])
COND_LAB <- c(ctrl = "Unstimulated", stim = "IFN-β, 6 h")   # say the stimulus (the figure is read alone)

ts <- fread("data/Kang2018/raw/GSE96583_batch2.total.tsne.df.tsv.gz", header = FALSE, skip = 1,
            col.names = c("cell_id", "tsne1", "tsne2", "donor", "condition", "cluster", "cell_type", "multiplets"))
d <- ts[multiplets == "singlet" & cell_type %in% AUTH]                          # row order of the released table is kept
cat(sprintf("singlets %d; with a cell type %d (ctrl %d, stim %d)\n", ts[multiplets == "singlet", .N], nrow(d),
            d[condition == "ctrl", .N], d[condition == "stim", .N]))
d[, cond := factor(condition, levels = names(COND_LAB))]; d[, ct := factor(cell_type, levels = names(ABB))]

## ---- sizes: a and b share one 177-mm row of the composite; b is as narrow as its text, a takes the rest ----
ROW_W <- 177; GAP_AB <- 3; H_ROW <- 46; FS <- 6; PAD <- 1.0; lwd_pt <- function(pt) pt / (72.27 / 96)
mk <- fread(file.path(OUT, "figS5_pbmc3k_markers.csv"))[rank <= 3, .(Markers = paste(gene, collapse = ", ")), by = cell_type]
tab <- data.table(cell_type = unname(ABB), Markers = mk$Markers[match(names(ABB), mk$cell_type)],
  Description = c("Classical monocytes", "Non-classical monocytes", "Dendritic cells", "Megakaryocytes (platelets)",
                  "CD4+ helper T cells", "B cells", "CD8+ cytotoxic T cells", "Natural killer cells"),
  label_authors = names(ABB))
tab <- merge(tab, d[, .(n_unstim = sum(condition == "ctrl"), n_stim = sum(condition == "stim")), by = .(label_authors = cell_type)],
             by = "label_authors", sort = FALSE)
fwrite(tab, file.path(OUT, "figS5_celltype_table.csv"))   # the full names are given in the key of a, so the table shows cell type and markers
ragg::agg_png(tempfile(fileext = ".png"), width = 100, height = 100, units = "mm", res = 72)
txt_w <- function(s, face, parse = FALSE) max(vapply(s, function(z) convertWidth(grobWidth(textGrob(if (parse) parse(text = z) else z,
  gp = gpar(fontfamily = FAM, fontsize = FS, fontface = face))), "mm", valueOnly = TRUE), 0))
hdr <- c("Cell type", "Top markers")   # headers as Fig. 3b
need <- c(max(txt_w(paste0("bold(", tab$cell_type, ")"), "bold", TRUE), txt_w(hdr[1], "bold")), max(txt_w(tab$Markers, "italic"), txt_w(hdr[2], "bold"))) + 2 * PAD
invisible(dev.off())
W_B <- ceiling(sum(need)) + 1; W_A <- ROW_W - GAP_AB - W_B
cat(sprintf("row: a %.0f mm + gap %.0f + b %.0f mm (height %.0f mm)\n", W_A, GAP_AB, W_B, H_ROW))

## ---- panel a (layout of Fig. 2a of this paper: no frame or ticks, "t-SNE 1/2"; colours and abbreviations of the original;
## the cell-type key, outside, gives the full names with the abbreviations used in b-f; the condition key is inside at the top
## left, as the original, with names that say the stimulus) ----
PT <- 0.42; KEY <- 1.6
th <- theme_minimal(base_size = 7, base_family = FAM) +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), aspect.ratio = 1,
        axis.title = element_text(size = 7, colour = "black"), axis.title.x = element_text(margin = margin(t = 1)),
        axis.title.y = element_text(margin = margin(r = 1)),
        legend.position = "right", legend.justification = c(0, 0.5), legend.box.spacing = unit(0.8, "mm"),
        legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6, lineheight = 0.9, margin = margin(l = 0.5, unit = "mm")),
        legend.key.size = unit(2.4, "mm"), legend.key.spacing.y = unit(0, "mm"),
        legend.background = element_blank(), legend.key = element_blank(), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(0.5, 0.8, 0.5, 0.5, "mm"))
ys <- scale_y_continuous(expand = expansion(mult = c(0.03, 0.24)))   # room above the cells for the condition key and its title (both maps alike)
pa <- ggplot(d, aes(tsne1, tsne2, colour = cond)) + geom_point(size = PT, shape = 16, stroke = 0) + ys +
  scale_colour_manual(values = COND_PAL, labels = COND_LAB, name = "condition") +
  guides(colour = guide_legend(override.aes = list(size = KEY))) +
  labs(x = "t-SNE 1", y = "t-SNE 2") + th +
  theme(legend.position = "inside", legend.position.inside = c(0, 1), legend.justification = c(0, 1))
## the key gives the full names, with the abbreviations used in the other panels (first appearance; plotmath for T[H])
FULL <- c('"Classical monocytes (cM)"', '"Non-classical monocytes (ncM)"', '"Dendritic cells (DC)"', '"Megakaryocytes (Mkc)"',
          '"CD4+ helper T cells ("*T[H]*")"', '"B cells (B)"', '"CD8+ cytotoxic T cells (Tc)"', '"Natural killer cells (NK)"')
pb <- ggplot(d, aes(tsne1, tsne2, colour = ct)) + geom_point(size = PT, shape = 16, stroke = 0) + ys +
  scale_colour_manual(values = PAL, breaks = names(ABB), labels = parse(text = FULL), name = "cell type") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = KEY))) +
  labs(x = "t-SNE 1", y = "t-SNE 2") + th
ragg::agg_png(file.path(FIG, "FigS5a_tSNE.png"), width = px(W_A), height = px(H_ROW), units = "px", res = 600, background = "white")
print(pb + pa + plot_layout(ncol = 2)); invisible(dev.off()); cat("wrote FigS5a_tSNE.png\n")   # cell types left, condition right (as Fig. 3a)

## ---- panel b (table: cell type and markers, Arial 6 pt, as tall as a) ----
H_MM <- H_ROW - 2; NR <- nrow(tab)
ragg::agg_png(file.path(FIG, "FigS5b_celltype_table.png"), width = px(W_B), height = px(H_MM), units = "px", res = 600, background = "white")
wid <- need + (W_B - sum(need)) / 2; col_hj <- c(0.5, 0); col_x <- ifelse(col_hj == 0.5, 0.5, PAD / wid)
mk_theme <- function(j, parse) ttheme_minimal(padding = unit(c(0, 0), "mm"),
  core = list(fg_params = list(hjust = col_hj[j], x = col_x[j], fontsize = FS, fontfamily = FAM, parse = parse,
                               fontface = c("bold", "italic")[j], col = "grey12"),
              bg_params = list(fill = c("#F4F7FA", "white"), col = NA)),
  colhead = list(fg_params = list(fontface = "bold", hjust = col_hj[j], x = col_x[j], fontsize = FS, fontfamily = FAM, col = "white"),
                 bg_params = list(fill = "#34495E", col = "#34495E", lwd = lwd_pt(0.5))))
cols <- list(tab[, .(`Cell type` = paste0("bold(", cell_type, ")"))], tab[, .(`Top markers` = Markers)])   # plotmath: bold, T[H]
gs <- lapply(1:2, function(j) { g <- tableGrob(cols[[j]], rows = NULL, theme = mk_theme(j, j == 1)); g$widths <- unit(wid[j], "mm"); g })
g <- do.call(gtable_combine, c(gs, list(along = 1)))
g$heights <- unit(rep((H_MM - 0.3) / (NR + 1), NR + 1), "mm")
g <- gtable_add_grob(g, segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(col = "#34495E", lwd = lwd_pt(1.0))), t = 1, b = 1, l = 1, r = ncol(g))
g <- gtable_add_grob(g, segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(col = "grey70", lwd = lwd_pt(0.5))), t = nrow(g), b = nrow(g), l = 1, r = ncol(g))
grid.newpage(); pushViewport(viewport(y = 1, height = unit(H_MM - 0.3, "mm"), just = "top")); grid.draw(g); invisible(dev.off())
cat("wrote FigS5b_celltype_table.png\n"); print(tab)
cat(sprintf("compose: b at x = %.1f mm\n", W_A + GAP_AB))
