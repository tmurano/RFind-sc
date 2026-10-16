## Fig. 3a (De Donno) — UMAP of the control-line cells, coloured by cell type (left map; "a" in the code: mk_a, SZ$a) and by
##   condition (right map; "b" in the code), in the design of the t-SNE maps of Fig. 2a.
##   Only control-line mice (line == "Ctrl"; GR/MR flox/flox without Cre) are drawn. They are the only mice
##   kept by analysis/dedonno/01_qc_embed.R, so the embedding is computed on these cells (26,235) and the line filter below
##   keeps every cell.
##   cell-type map: the cell types analysed in Fig. 3c-f (the authors' seven major types and GABAergic neurons; = `shown` in
##          Fig3_celltype_ncells.csv, written by fig3_violins.R) are coloured and named on the map; the other author types
##          are light grey without names
##   condition map: the SAME UMAP, coloured by condition (No Stress vs Stress)
## Embedding computed by analysis/dedonno/01_qc_embed.R (own QC -> HVG/PCA/UMAP; not the authors' coordinates).
## Design as figs/Fig2_Hammond/analysis/fig2_tsne_original.R: the cell-type map without legend (in-map labels), the condition
##   map with a legend, identical fixed map size in both, no cell counts on the map (they go in the figure legend,
##   manuscript/Legends.md).
##
## Input : data/DeDonno2026/output/scores_dedonno.csv   (umap_x/y, cell_type, condition, line)
##         figs/Fig3_DeDonno/output/Fig3_celltype_ncells.csv   (run fig3_violins.R first)
## Output: figs/Fig3_DeDonno/output/fig3_umap_counts.csv
##         + {Fig3a_UMAP_celltype, Fig3a_UMAP_condition}.png (drawn at the pasted size)
## Run   : Rscript "figs/Fig3_DeDonno/analysis/fig3_umap.R"   (from project root)

pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
source("figs/Fig3_DeDonno/analysis/fig3_celltypes.R")
set.seed(1)
OUT <- "figs/Fig3_DeDonno/output"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

d <- fread("data/DeDonno2026/output/scores_dedonno.csv")
d <- d[is.finite(umap_x) & is.finite(umap_y)]
cat(sprintf("post-QC cells: %d\n", nrow(d)))
d <- d[line == "Ctrl"]
cat(sprintf("control-line cells drawn: %d / mice: %d / cell types: %d\n", nrow(d), uniqueN(d$donor_id), uniqueN(d$cell_type)))
stopifnot(all(d$cell_type %in% CT_LV))
shown <- fread(file.path(OUT, "Fig3_celltype_ncells.csv"))[shown == TRUE]$cell_type
stopifnot(length(shown) > 0, all(shown %in% CT_LV))
ct_lv <- CT_LV[CT_LV %in% unique(d$cell_type)]
d[, cell_type := factor(cell_type, levels = ct_lv)]
d[, is_shown := cell_type %in% shown]
cen <- d[, .(x=median(umap_x), y=median(umap_y), n=.N), by=.(cell_type, is_shown)]
fwrite(cen[order(cell_type)], file.path(OUT, "fig3_umap_counts.csv"))

## ---- drawn at the size pasted into Fig. 3 ----
## Each map is drawn at its pasted size (SZ below; cell-type map 42.7 x 39.1 mm, condition map 55.6 x 39.1 mm, 600 dpi): Arial,
## axis titles and legend 7 pt, cell-type names on the map 5 pt (LAB_MM below); dots (PTS, scaled with the map width) and
## keys scaled by SHRINK.
FAM <- "Arial"
SHRINK <- 0.19                                 # scale of the dot and key sizes below
PT  <- 0.9 * SHRINK                            # cell dot size
LAB_MM <- 5 / .pt                              # 5 pt, as the cluster names of Fig. 2a
base_theme <- theme_minimal(base_size = 7, base_family = FAM) +
  theme(panel.grid = element_blank(),          # no background grid on maps
        axis.title = element_text(size = 7, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 1)), axis.title.y = element_text(margin = margin(r = 1)),
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

## --- cell-type map (no legend; names of the analysed types on the map are the key) ---
## types not analysed are drawn first in light grey, so the analysed types sit on top
rx <- diff(range(d$umap_x)); ry <- diff(range(d$umap_y))
mk_a <- function(PT) ggplot(mapping = aes(umap_x, umap_y)) +
  geom_point(data = d[is_shown == FALSE], colour = "grey85", size = PT, alpha = 0.75, stroke = 0) +
  geom_point(data = d[is_shown == TRUE][sample(.N)], aes(colour = cell_type), size = PT, alpha = 0.75, stroke = 0) +
  scale_colour_manual(values = CT_PAL) +
  geom_label_repel(data = cen[is_shown == TRUE], aes(x, y, label = CT_LAB[as.character(cell_type)]), inherit.aes = FALSE,
                   ## box.padding 0.12 keeps the CA1 and DG boxes apart at the 39.1-mm size
                   size = LAB_MM, family = FAM, fontface = "bold", label.padding = 0.08, box.padding = 0.12,
                   label.size = 0.1, label.r = 0.08, fill = "white", alpha = 0.9,
                   colour = "black", segment.colour = "grey50", segment.size = 0.15, max.overlaps = Inf,
                   min.segment.length = 0, seed = 1) +
  ## key for the grey cells: "Others" = the author types not analysed, bottom-right corner (empty there)
  annotate("point", x = max(d$umap_x) - 0.27 * rx, y = min(d$umap_y) + 0.02 * ry, shape = 21, size = 5 * 1.5 * SHRINK,
           fill = "grey85", colour = "grey55", stroke = 0.15) +
  annotate("text", x = max(d$umap_x) - 0.23 * rx, y = min(d$umap_y) + 0.02 * ry, label = "Others", hjust = 0,
           size = LAB_MM, family = FAM, colour = "grey25") +
  ## clip off: "Others" may run into the blank strip right of the map
  coord_equal(clip = "off") + labs(x="UMAP 1", y="UMAP 2") + base_theme + theme(legend.position = "none")
pa <- mk_a(PT)

## --- condition map (same UMAP, legend on the right) ---
d[, condition := factor(condition, levels=c("No Stress","Stress"))]
cond_pal <- c("No Stress"="#4C78A8", "Stress"="#E45756")
## dots of the same size and alpha as in the cell-type map, so that the two maps are drawn alike
mk_b <- function(PT) ggplot(d[sample(.N)], aes(umap_x, umap_y, colour=condition)) +
  geom_point(size=PT, alpha=0.75, stroke=0) +
  scale_colour_manual(values=cond_pal, name="condition", drop=FALSE,
                      guide=guide_legend(override.aes=list(size=1.6, alpha=1))) +
  coord_equal() + labs(x="UMAP 1", y="UMAP 2") + base_theme +
  theme(legend.text = element_text(size = 7), legend.title = element_text(size = 7, face = "bold"),
        legend.key.height = unit(3.2, "mm"), legend.key.width = unit(2.2, "mm"),
        legend.box.spacing = unit(1, "mm"), legend.margin = margin(0, 0, 0, 0))
pb <- mk_b(PT)

## identical map size in the two maps (as in Fig. 2a; both are 39.1 mm high): the panel gets a fixed size, the largest one
## that fits both canvases with the data aspect kept (the condition map is limited by its legend, so the cell-type map keeps
## a blank strip on the right). SAME_MAP <- FALSE lets the cell-type map fill its own canvas instead (then the two maps
## differ in size).
SAME_MAP <- TRUE
OUT_PPTX <- OUT   # pasted-size panels go to output/
px <- function(mm) round(mm / 25.4 * 600)                                   # ragg truncates mm sizes, so give whole pixels
SZ <- list(a = c(42.7, 39.1), b = c(55.6, 39.1))                          # pasted size (mm)
nonpanel <- function(p) {   # width and height (mm) of everything around the panel (axis titles, legend)
  g <- ggplotGrob(p); pn <- g$layout[g$layout$name == "panel", ]
  w <- sum(grid::convertWidth(g$widths[-pn$l], "mm", valueOnly = TRUE))
  h <- sum(grid::convertHeight(g$heights[-pn$t], "mm", valueOnly = TRUE))
  c(w, h)
}
ragg::agg_png(tempfile(fileext = ".png"), width = 100, height = 100, units = "mm", res = 72)   # measure with the same device
npa <- nonpanel(pa); npb <- nonpanel(pb); invisible(dev.off())
fit <- function(sz, np) min(sz[1] - np[1], (sz[2] - np[2]) * rx / ry) - 0.3      # widest panel that fits (mm)
PWS <- if (SAME_MAP) rep(min(fit(SZ$a, npa), fit(SZ$b, npb)), 2) else c(fit(SZ$a, npa), fit(SZ$b, npb))
cat(sprintf("map panel a %.1f x %.1f mm, b %.1f x %.1f mm\n", PWS[1], PWS[1] * ry / rx, PWS[2], PWS[2] * ry / rx))
fix_panel <- function(p, PW) {   # panel fixed in mm; the spare width goes to the right (no left padding), spare height split
  g <- ggplotGrob(p); pn <- g$layout[g$layout$name == "panel", ]
  g$widths[pn$l] <- unit(PW, "mm"); g$heights[pn$t] <- unit(PW * ry / rx, "mm"); g$respect <- FALSE
  g <- gtable::gtable_add_cols(g, unit(1, "null"))
  gtable::gtable_add_rows(gtable::gtable_add_rows(g, unit(1, "null"), 0), unit(1, "null"))
}
save_pptx <- function(p, fname, sz) {
  ragg::agg_png(file.path(OUT_PPTX, fname), width = px(sz[1]), height = px(sz[2]), units = "px", res = 600, background = "white")
  grid::grid.newpage(); grid::grid.draw(p); invisible(dev.off())
  cat(sprintf("wrote %s (%.1f x %.1f mm)\n", fname, sz[1], sz[2]))
}
## dot size follows the map width (PT is the dot size for a 25-mm map), so that the maps look alike at any size
PTS <- PT * PWS / 25
set.seed(1); pa <- mk_a(PTS[1]); pb <- mk_b(PTS[2])
save_pptx(fix_panel(pa, PWS[1]), "Fig3a_UMAP_celltype.png", SZ$a)
save_pptx(fix_panel(pb, PWS[2]), "Fig3a_UMAP_condition.png", SZ$b)
print(d[, .N, by = condition])
print(cen[order(-n)])
