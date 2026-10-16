## Fig4e (Perturb-seq) — heatmap of where the knock-outs of a gene set sit in the ranking against the stress response of each cell type.
##   A gene x cell-type comparison is shown as a heatmap. Rows = the knock-outs of the set (with their role); columns = DG, CA1,
##   CA3 and GABA. The paper uses the set GR-complex, the same five corticosteroid-receptor and GR-complex genes in the library
##   as Fig4d (Nr3c1, Nr3c2, Hsp90ab1, Ncoa1, Ncoa3; fig4d_go_set.R; Methods). Other sets (c-named, GR-lit) are kept as options.
##   Number and colour = the position from the opposed end as a percentage of the knock-outs ranked in that cell type (rank
##   counted from the most opposed knock-out, as in c, / number ranked x 100; percentages rather than ranks because the totals
##   differ between cell types). The RFind scores and the pipeline ranks (rank_among_KO, 1 = most aligned) are unchanged.
## Input : analysis/perturbseq/results/axis_<target>_vs_<ref>.csv   Output: figs/Fig4_Perturbseq/output/Fig4d_heatmap_<set>.png
##   (panel e of Fig. 4; the script and file names use "fig4d"). Drawn at the size the panel has in the final figure
##   (44.0 x 50.0 mm, 600 dpi, Arial; tick labels 6 pt, legend 7 pt).
## Run   : <R> -f figs/Fig4_Perturbseq/analysis/fig4d_heatmap.R --args GR-complex   (the panel of the paper; c-named or GR-lit
##         draw the other sets)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
pdf(NULL)
OUT <- "figs/Fig4_Perturbseq/output"; FONT <- "Arial"; dir.create(OUT, showWarnings = FALSE)
W_MM <- 44.0; H_MM <- 50.0; DPI <- 600   # size in the final figure (high enough for the two-line row names)
CT <- c(DG = "g017_DG", CA1 = "CA1", CA3 = "g017_CA3", GABA = "GABA3")
CT_NM <- c(DG = "DG\nGlut", CA1 = "CA1\nGlut", CA3 = "CA3\nGlut", GABA = "GABA\n")   # column labels = the names of De Donno et al., on two lines
## ("DG Glut" and the others do not fit the narrow columns on one line; GABA gets an empty second line so all labels align)
SETS <- list(
  ## the set of the paper: the corticosteroid receptors and the GR complex, every library member of this definition (5 rows)
  "GR-complex" = c(Nr3c1 = "GR", Nr3c2 = "MR", Hsp90ab1 = "chaperone", Ncoa1 = "coactivator", Ncoa3 = "coactivator"),
  ## other sets, kept as options: knock-outs of the DG lists of Fig4c, with Nr3c1 in GR-lit
  "c-named"    = c(Ncoa1 = "GR coactivator", Crebrf = "GR repressor", Gal = "galanin"),
  "GR-lit"     = c(Nr3c1 = "receptor", Ncoa1 = "coactivator", Kmt2d = "stabilises GR", Crebrf = "repressor"))
args <- commandArgs(TRUE); SET <- if (length(args)) args[1] else "GR-complex"
ROLE <- SETS[[SET]]; ROWS <- names(ROLE)
P <- rbindlist(lapply(names(CT), function(ct) {
  x <- fread(sprintf("analysis/perturbseq/results/axis_%s_vs_%s.csv", ct, CT[[ct]]))[kind == "KO"]
  setorder(x, -rfind); x[, rank := .I]; n <- nrow(x)
  stopifnot(x$rank == x$rank_among_KO)                       # the pipeline rank, 1 = most aligned
  x[, rank_opp := n - rank + 1L]                             # displayed rank, 1 = most opposed
  x[target %in% ROWS, .(ct, KO = target, rank, rank_opp, n_ko = n, pct = 100 * rank_opp / n)]
}))
## numbers = the position from the opposed end in % (rank counted from the most opposed knock-out, as in c, / number ranked;
## the scores do not change), one decimal; the colour shows the same position
P[, lab := sprintf("%.1f", pct)]
P[, ct_lab := factor(sprintf("%s\n%s\nKOs", CT_NM[ct], format(n_ko, big.mark = ",", trim = TRUE)),   # four lines: the columns are narrow
                     levels = unique(sprintf("%s\n%s\nKOs", CT_NM[ct], format(n_ko, big.mark = ",", trim = TRUE)))[order(match(unique(ct), names(CT)))])]
P[, ct := factor(ct, levels = names(CT))]
setorder(P, ct)
P[, ct_lab := factor(ct_lab, levels = unique(ct_lab))]
## the role goes on a second line, so that the four columns are wide enough for their 6-pt labels
P[, row := factor(sprintf("%s\n(%s)", KO, ROLE[KO]), levels = rev(sprintf("%s\n(%s)", ROWS, ROLE[ROWS])))]
print(dcast(P, KO ~ ct, value.var = "rank_opp")); print(dcast(P, KO ~ ct, value.var = "lab"))

## sizes: text in pt (geom_text takes mm, hence / .pt); tick-type labels 6 pt, legend 7 pt
p <- ggplot(P, aes(ct_lab, row, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab, colour = abs(pct - 50) > 35), size = 6 / .pt, family = FONT, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_fill_gradient2(low = "#1F5FA8", mid = "grey96", high = "#C0392B", midpoint = 50, limits = c(0, 100),
                       breaks = c(0, 100), labels = c("0\nopposed end", "100\naligned end"),   # the two ends only ("middle" crowds them)
                       name = "Position in that cell type\n(% from the opposed end)") +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL) +
  theme_minimal(7, base_family = FONT) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 6, face = "bold", colour = "grey10", lineheight = 1.0),
        axis.text.y = element_text(size = 6, face = "italic", colour = "grey10", lineheight = 0.95), legend.position = "bottom",
        legend.title = element_text(size = 7, margin = margin(b = 2.5)), legend.text = element_text(size = 7, margin = margin(t = 2), lineheight = 0.9),
        legend.key.width = unit(3.4, "mm"), legend.key.height = unit(1.8, "mm"),   # a short bar (17 mm), so that its end labels stay inside 44 mm
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1.5, "mm"), plot.margin = margin(1, 1, 1, 1)) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
## drawn directly at the final size; whole pixels = round(mm / 25.4 * 600)
ragg::agg_png(file.path(OUT, sprintf("Fig4d_heatmap_%s.png", SET)), width = round(W_MM / 25.4 * DPI), height = round(H_MM / 25.4 * DPI),
              units = "px", res = DPI, background = "white")
print(p); invisible(dev.off())
cat("wrote Fig4d_heatmap_", SET, "\n")
