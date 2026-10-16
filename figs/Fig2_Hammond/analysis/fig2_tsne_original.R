## Fig. 2a (Hammond) — ORIGINAL t-SNE (authors' embedding + authors' clusters): coloured by cluster (left map) and by
## age / condition (right map)
##
## IMPORTANT: this panel uses Hammond et al. 2019's OWN t-SNE coordinates and OWN cluster annotations as published (not a
## re-embedding), drawn for the cells analysed in Fig. 2: the 74,189 microglia after our QC (monocytes/macrophages are not
## shown).
##
## Input : data/Hammond2019/raw/cell_meta_hammond.csv  (cols: tsne_x/tsne_y, paper_cluster)
## Output: figs/Fig2_Hammond/output/Fig2a_tSNE_clusters.png (33.3 x 35 mm) +
##         Fig2a_tSNE_conditions.png (52.2 x 35 mm), 600 dpi + output/fig2_tsne_cluster_counts.csv
## Run   : Rscript "figs/Fig2_Hammond/analysis/fig2_tsne_original.R"   (from project root)

pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
set.seed(1)
BASE <- "figs/Fig2_Hammond"; DATA <- file.path(BASE, "data"); OUT <- file.path(BASE, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

meta <- fread("data/Hammond2019/raw/cell_meta_hammond.csv")
meta <- meta[is.finite(tsne_x) & is.finite(tsne_y)]
qc <- fread("data/Hammond2019/output/scores_hammond_9axes.csv", select = c("cell_id", "cluster"))
meta <- meta[cell_id %in% qc[cluster != "Mono/Mac", cell_id]]   # the microglia after QC (Mono/Mac not shown)
stopifnot(nrow(meta) == 74189L)
cat(sprintf("cells: %d / author clusters: %d\n", nrow(meta), uniqueN(meta$paper_cluster)))

## authors' 14 clusters: order + 14 clearly-distinct colours
lev <- c("1","2a","2b","2c","3","4","5","6","7a","7b","7c","8","9","Mono/Mac")
lev <- lev[lev %in% unique(meta$paper_cluster)]   # Mono/Mac drops out with the QC filter above
## colours follow the violin signature colours of fig2_violins.R as far as 14 clusters allow:
## 2a/2b/2c = shades of its red, 3 = orange, 4 = blue, 7a/7b/7c = shades of its black, 9 = green;
## the clusters without a signature colour (1, 5, 6, 8, Mono/Mac) get colours outside those families.
pal <- c("1"="#F032E6","2a"="#E41A1C","2b"="#A50F15","2c"="#F7797B","3"="#E69F00",
         "4"="#2C7FB8","5"="#911EB4","6"="#9A6324","7a"="#1A1A1A","7b"="#5C5C5C",
         "7c"="#9A9A9A","8"="#808000","9"="#2CA02C","Mono/Mac"="#17BECF")
meta[, cluster := factor(paper_cluster, levels = lev)]
cen <- meta[, .(x = median(tsne_x), y = median(tsne_y), n = .N), by = cluster]
fwrite(cen[order(cluster)], file.path(OUT, "fig2_tsne_cluster_counts.csv"))

## ---- pasted-size rendering ----
## Each map is drawn at the exact size it is pasted into the figure (mm),
## 600 dpi with the dpi stored in the PNG, Arial throughout: axis titles / legend 7 pt, in-map cluster labels 5 pt.
FAM <- "Arial"; PT <- ggplot2::.pt
OUT_PPTX <- OUT   # pasted-size panels go to output/
W_CL <- 33.3; W_CO <- 52.2; H_TS <- 35.0                        # target sizes (mm): clusters map, conditions map, common height
PT_SIZE <- 0.12                                                  # cell dots

p <- ggplot(meta[sample(.N)], aes(tsne_x, tsne_y, colour = cluster)) +
  geom_point(size = PT_SIZE, alpha = 0.75, stroke = 0) +
  scale_colour_manual(values = pal, name = "Paper cluster", drop = FALSE,
                      guide = guide_legend(override.aes = list(size = 2.8, alpha = 1))) +
  ## cluster names on the map are the only key: 5 pt bold at the pasted size (14 labels on a 30-mm map)
  geom_label_repel(data = cen, aes(x, y, label = cluster), inherit.aes = FALSE, size = 5 / PT, family = FAM,
                   fontface = "bold", label.padding = unit(0.25, "mm"), box.padding = unit(0.15, "mm"),
                   label.r = unit(0.3, "mm"), fill = "white", alpha = 0.9, colour = "black", segment.colour = "grey50",
                   label.size = 0.08, segment.size = 0.15, max.overlaps = Inf, min.segment.length = 0, seed = 1) +
  ## cell / cluster counts are given in the figure legend, not on the map
  coord_equal() +
  labs(title = NULL, subtitle = NULL, x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal(base_size = 7, base_family = FAM) +
  ## no background grid
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 7),     # 7 pt at the pasted size
        axis.title.x = element_text(margin = margin(t = 1)), axis.title.y = element_text(margin = margin(r = 1)),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "none")

## identical map size in the two maps: the panel is given a fixed size, so the cluster map (no legend) and the
## age / condition map (legend) have the same map height and the same canvas height; only the second canvas is wider.
## The map size is chosen so that map + axis titles fill the 35 mm height exactly; the drawing is
## anchored top-left (no left padding) and any spare width falls at the right.
rx <- diff(range(meta$tsne_x)); ry <- diff(range(meta$tsne_y))
fix_panel <- function(g, pw, ph) {                               # ggplot -> gtable with the panel fixed at pw x ph mm
  gt <- ggplotGrob(g); gt$respect <- FALSE
  pl <- gt$layout[gt$layout$name == "panel", ]
  gt$widths[pl$l] <- unit(pw, "mm"); gt$heights[pl$t] <- unit(ph, "mm"); gt
}
save_topleft <- function(gt, fname, w_mm, h_mm) {
  ## size in pixels (mm / 25.4 * 600, rounded): units = "mm" truncates to one pixel short
  ragg::agg_png(file.path(OUT_PPTX, fname), width = round(w_mm / 25.4 * 600), height = round(h_mm / 25.4 * 600),
                units = "px", res = 600, background = "white")
  gw <- convertWidth(sum(gt$widths), "mm", valueOnly = TRUE); gh <- convertHeight(sum(gt$heights), "mm", valueOnly = TRUE)
  grid.newpage()
  pushViewport(viewport(x = 0, y = 1, width = unit(gw, "mm"), height = unit(gh, "mm"), just = c("left", "top")))
  grid.draw(gt); popViewport(); invisible(dev.off())
  cat(sprintf("wrote %s  (%.1f x %.1f mm canvas; drawing %.2f x %.2f mm)\n", fname, w_mm, h_mm, gw, gh))
}
library(grid)
## non-panel height (x title) measured on the device, then the map height is whatever is left of 35 mm
ragg::agg_png(tempfile(fileext = ".png"), width = W_CL, height = H_TS, units = "mm", res = 600)
g0 <- fix_panel(p, 10, 10)
PH <- H_TS - (convertHeight(sum(g0$heights), "mm", valueOnly = TRUE) - 10)
PW <- PH * rx / ry
PW <- min(PW, W_CL - (convertWidth(sum(g0$widths), "mm", valueOnly = TRUE) - 10))
PH <- PW * ry / rx; invisible(dev.off())
cat(sprintf("map panel: %.2f x %.2f mm\n", PW, PH))
save_topleft(fix_panel(p, PW, PH), "Fig2a_tSNE_clusters.png", W_CL, H_TS)
cat("wrote Fig2a_tSNE_clusters.png + counts csv\n")
print(cen[order(-n)])

## ============================================================================
## Second panel: SAME t-SNE coordinates, coloured by age / condition
## (developmental timeline E14->P540 + P100 LPC-demyelination vs Saline control)
## ============================================================================
meta[, age_grp := fifelse(age == "P100_Injured_WM",  "P100 LPC (injured)",
                   fifelse(age == "P100_Control_WM", "P100 Saline (control)", age))]
age_lv <- c("E14","P4/P5","P30","P100","P540","P100 Saline (control)","P100 LPC (injured)")
age_lv <- age_lv[age_lv %in% unique(meta$age_grp)]
meta[, age_grp := factor(age_grp, levels = age_lv)]
## E14 purple and P4/P5 teal, so that the two embryonic / early postnatal groups are distinct.
cond_pal <- c("E14"="#5E3C99","P4/P5"="#00A0B0","P30"="#33A02C","P100"="#FF7F00",
              "P540"="#B15928","P100 Saline (control)"="#1F78B4","P100 LPC (injured)"="#E7298A")
base_theme <- theme_minimal(base_size = 7, base_family = FAM) +
  theme(panel.grid = element_blank(),          # no background grid
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 7),          # 7 pt at the pasted size, as in the cluster map
        axis.title.x = element_text(margin = margin(t = 1)), axis.title.y = element_text(margin = margin(r = 1)),
        plot.margin = margin(0, 0, 0, 0),
        ## legend 7 pt at the pasted size
        legend.text = element_text(size = 7, lineheight = 0.9), legend.title = element_text(size = 7, face = "bold"),
        legend.key.height = unit(3.2, "mm"), legend.key.width = unit(2.2, "mm"), legend.key.spacing.x = unit(0.6, "mm"), legend.key.spacing.y = unit(0.5, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1.2, "mm"),
        legend.justification = c(0, 0.5))
## the two long condition names break before the bracket so the 7 pt legend fits the 52.2 mm canvas
## the embryonic age is written as in Hammond et al. 2019, "E14.5" (the data key stays "E14", the GEO sample name)
p2 <- ggplot(meta[sample(.N)], aes(tsne_x, tsne_y, colour = age_grp)) +
  geom_point(size = PT_SIZE, alpha = 0.75, stroke = 0) +
  scale_colour_manual(values = cond_pal, name = "age / condition", drop = FALSE,
                      labels = function(x) sub(" (", "\n(", sub("^E14$", "E14.5", x), fixed = TRUE),
                      guide = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  coord_equal() +
  labs(title = NULL, subtitle = NULL, x = "t-SNE 1", y = "t-SNE 2") + base_theme
save_topleft(fix_panel(p2, PW, PH), "Fig2a_tSNE_conditions.png", W_CO, H_TS)
cat("wrote Fig2a_tSNE_conditions.png\n")
print(meta[, .N, by = age_grp][order(-N)])
