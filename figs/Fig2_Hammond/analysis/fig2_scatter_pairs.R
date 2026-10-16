## Fig. 2e, f — 2D maps of two RFind indices per cell: Maturity x DAM (e) and DAM x Interferon (f).
## Each pair is drawn twice on identical axes:
##   dev : development series E14 -> P540 coloured by age, Saline/LPC grey;
##         r drawn for Early development (E14.5-P4/P5) and Maturation/ageing (P30-P540); r of the pooled series in the CSV only
##   inj : P100 white matter Saline / LPC coloured, development grey;
##         r drawn for Saline and LPC; r of Saline + LPC pooled in the CSV only
## Points = per-mouse means (mice with >= 20 cells in the highlighted groups).
## The age-group means are not drawn (no markers, no arrow); their positions anchor the group labels, which are the figure's
## colour key.
## What the marks mean is written in manuscript/Legends.md, not in the figure.
##
## Input : data/Hammond2019/output/scores_hammond_9axes.csv, data/Hammond2019/raw/cell_meta_hammond.csv
## Output: figs/Fig2_Hammond/output/Fig2{e,f}_scatter_{dev,inj}_{Maturity_vs_DAM,DAM_vs_Interferon}.png
##         (44 x 38 mm each, 600 dpi = the size pasted in the figure) + output/fig2_scatter_correlations.csv

pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })
OUT <- "figs/Fig2_Hammond/output"; set.seed(1)
d <- fread("data/Hammond2019/output/scores_hammond_9axes.csv")
m <- fread("data/Hammond2019/raw/cell_meta_hammond.csv")[, .(cell_id, age)]
d <- merge(d, m, by = "cell_id")
d <- d[cluster %in% c("1","2a","2b","2c","3","4","5","6","7a","7b","7c","8","9")]    # Mono/Mac excluded
d[, animal := sub("_[^_]+$", "", cell_id)]                                             # per-mouse sample ID
d[, grp := fifelse(age == "P100_Injured_WM", "LPC", fifelse(age == "P100_Control_WM", "Saline", age))]
devlv <- c("E14","P4/P5","P30","P100","P540"); injlv <- c("Saline","LPC"); grplv <- c(devlv, injlv)
d <- d[grp %in% grplv]; d[, grp := factor(grp, levels = grplv)]
pal <- c("E14"="#5E3C99","P4/P5"="#00A0B0","P30"="#33A02C","P100"="#FF7F00","P540"="#B15928",
         "Saline"="#1F78B4","LPC"="#E7298A")
BG <- "#8C8C8C"
DISP <- c(Maturity="Maturity-index", DAM="DAM-index", Interferon="Interferon-index")
MIN_CELLS_ANIMAL <- 20L
FAM <- "Arial"              # all text in these panels; geom_text/geom_label do not inherit it from the theme
## pasted-size rendering: each panel is drawn at the exact size it is pasted in the figure (44 x 38 mm,
## 600 dpi, dpi stored in the PNG). Tick labels 6 pt, axis titles 7 pt; r lines 5.5 pt (inj) and 5 pt (dev, whose longer
## group names must fit the panel width); group labels 5 pt (five labels in a 36-mm panel).
PT <- ggplot2::.pt; W_MM <- 44; H_MM <- 38
OUT_PPTX <- OUT   # pasted-size panels go to output/
R_SIZE <- 5.5 / PT; GRP_SIZE <- 5 / PT
th <- theme_minimal(7, base_family = FAM) + theme(axis.title = element_text(size = 7), axis.text = element_text(size = 6),
                                axis.title.x = element_text(margin = margin(t = 1)), axis.title.y = element_text(margin = margin(r = 1)),
                                panel.grid = element_blank(), legend.position = "none", plot.margin = margin(3.5, 4.5, 1.5, 0))
## plot.margin (pt): room above and to the right for half of the outermost tick labels (e.g. "50" at the top of e and the right of f),
## and below for the descenders of the x-axis title

## axis range = all-cells 0.5-99.5% UNION Saline/LPC 0.5-99.5%, identical for the dev and inj versions
lim <- function(ax){ q <- quantile(d[[ax]], c(.005, .995)); qi <- quantile(d[grp %in% injlv][[ax]], c(.005, .995))
  c(min(q[1], qi[1]), max(q[2], qi[2])) }

CORR <- list(); RTXT <- list()
mk <- function(xa, ya, mode){
  fg <- if (mode == "dev") devlv else injlv
  xl <- lim(xa); yl <- lim(ya); xr <- diff(xl); yr <- diff(yl)
  cr <- function(s) cor(s[[xa]], s[[ya]])
  if (mode == "dev") {
    ## the names of the two age groups as in Results and Methods
    labs_r <- c("development", "Early development (E14.5–P4/P5)", "Maturation/ageing (P30–P540)")
    vals <- c(cr(d[grp %in% devlv]), cr(d[grp %in% c("E14","P4/P5")]), cr(d[grp %in% c("P30","P100","P540")]))
    ## x1f = where the r value sits, as a fraction of the x range (in dev the name and the r value are written as one
    ## line instead, see below)
    cols <- c("grey30", "#5E3C99", "#B15928"); x1f <- 0.78
  } else {
    labs_r <- c("Saline/LPC", "Saline", "LPC")
    vals <- c(cr(d[grp %in% injlv]), cr(d[grp == "Saline"]), cr(d[grp == "LPC"]))
    ## 0.26: the r value starts right of the 5.5-pt "Saline" label
    cols <- c("grey30", unname(pal["Saline"]), unname(pal["LPC"])); x1f <- 0.26
  }
  CORR[[paste(mode, xa, ya)]] <<- data.table(mode = mode, pair = paste0(xa, "_vs_", ya), group = labs_r, r = round(vals, 3))
  ## pooled r (row 1) is computed and saved to the CSV but not drawn; only the two age groups (dev) or
  ## Saline and LPC (inj) are shown on the figure
  ## a value that rounds to 0 is written "0.00", not "-0.00" (e.g. LPC, Maturity x DAM: r = -0.002)
  rlab <- data.table(lab = labs_r, val = ifelse(round(vals, 2) == 0, "r = 0.00", sprintf("r = %+.2f", vals)), col = cols)[2:3]
  rlab[, `:=`(x0 = xl[1] + 0.03*xr, x1 = xl[1] + x1f*xr, y = yl[2] - 0.03*yr - 0.08*yr*(0:1))]
  ## dev: the group names ("Early development", "Maturation/ageing") are long, so the name and the r value are written as
  ## one line, right-aligned at the right edge of the panel, at 5 pt (the minimum of the Nature guide) so that it fits the
  ## 44-mm panel
  if (mode == "dev") rlab[, `:=`(lab = paste0(lab, " ", val), val = "", x0 = xl[2] - 0.005*xr)]
  rsz <- if (mode == "dev") 5 / PT else R_SIZE
  hj <- if (mode == "dev") 1 else 0
  ## dev: "Early development (E14.5–P4/P5) r = ..." is longer than the room left of its right end in panel e (its first
  ## letter would be cut by the left edge of the panel), so the two dev lines are not a layer: save_panel adds them to the
  ## laid-out plot outside the clipping of the panel, at the place, size and colour a geom_text layer would give them
  if (mode == "dev") RTXT[[paste(mode, xa, ya)]] <<- rlab[, .(lab, x0, y, col)]
  rtxt <- if (mode == "dev") NULL else list(
    geom_text(data = rlab, aes(x0, y, label = lab), colour = rlab$col, inherit.aes = FALSE, hjust = hj, vjust = 1, size = rsz, fontface = "bold", family = FAM),
    geom_text(data = rlab, aes(x1, y, label = val), colour = rlab$col, inherit.aes = FALSE, hjust = 0, vjust = 1, size = R_SIZE, fontface = "bold", family = FAM))
  cen <- d[grp %in% fg, .(x = mean(get(xa)), y = mean(get(ya))), by = grp][order(grp)]
  an  <- d[grp %in% fg, .(x = mean(get(xa)), y = mean(get(ya)), n = .N), by = .(grp, animal)][n >= MIN_CELLS_ANIMAL]
  ## point sizes and line widths for the 44-mm pasted width
  ptsz <- if (mode == "dev") 0.15 else 0.4; ptal <- if (mode == "dev") 0.5 else 0.85
  path <- cen[0, .(x, y, grp)]   # no path of the age-group means in any panel
  ggplot(mapping = aes(.data[[xa]], .data[[ya]])) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey75", linewidth = 0.15) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey75", linewidth = 0.15) +
    geom_point(data = d[!grp %in% fg][sample(.N)], colour = BG, size = 0.15, alpha = 0.35, stroke = 0) +
    geom_point(data = d[grp %in% fg][sample(.N)], aes(colour = grp), size = ptsz, alpha = ptal, stroke = 0) +
    rtxt +
    geom_point(data = an, aes(x, y, fill = grp), shape = 21, colour = "grey10", size = 1.1, stroke = 0.2, inherit.aes = FALSE) +
    geom_path(data = path, aes(x, y), inherit.aes = FALSE, colour = "black", linewidth = 0.4,
              arrow = arrow(length = unit(1.0, "mm"), type = "closed", ends = "last"), lineend = "round", linejoin = "round") +
    ## group labels only (no centroid marker, no leader line): they are the colour key
    ## the per-mouse points are passed to the repel layer with empty labels, so that the labels are pushed off the
    ## mouse points instead of covering them; empty labels are not drawn
    ## the embryonic age is written as in Hammond et al. 2019, "E14.5" (the data key stays "E14", the GEO sample name)
    geom_label_repel(data = rbind(cen[, .(x, y, grp, lab = sub("^E14$", "E14.5", as.character(grp)))], an[, .(x, y, grp, lab = "")]),
                     aes(x, y, label = lab, colour = grp), inherit.aes = FALSE, size = GRP_SIZE, fontface = "bold", family = FAM,
                     label.size = 0.08, label.padding = unit(0.12, "lines"), label.r = unit(0.3, "mm"), fill = "white", alpha = 0.92,
                     box.padding = unit(0.5, "mm"), point.padding = unit(0.3, "mm"), point.size = 1.1, max.iter = 1e5, force = 5, force_pull = 0.1,
                     max.overlaps = Inf, min.segment.length = Inf, seed = 1) +
    scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
    coord_cartesian(xlim = xl, ylim = yl) +
    labs(x = DISP[xa], y = DISP[ya]) + th   # no panel headers; left = development, right = injury (stated in the legend)
}

## drawn directly at the pasted size; pixels = mm / 25.4 * 600, rounded (units = "mm" truncates one pixel). rl = the r lines
## of dev (RTXT): drawn in the panel cell without clipping, data -> npc with the panel's ranges as geom_text does, 5 pt bold
save_panel <- function(p, file, rl = NULL) {
  ragg::agg_png(file, width = round(W_MM / 25.4 * 600), height = round(H_MM / 25.4 * 600), units = "px", res = 600,
                background = "white")
  gt <- ggplotGrob(p)
  if (!is.null(rl)) {
    pp <- ggplot_build(p)$layout$panel_params[[1]]; pc <- gt$layout[gt$layout$name == "panel", ]
    gt <- gtable::gtable_add_grob(gt, grid::textGrob(rl$lab, (rl$x0 - pp$x.range[1]) / diff(pp$x.range),
                                  (rl$y - pp$y.range[1]) / diff(pp$y.range), default.units = "native", hjust = 1, vjust = 1,
                                  gp = grid::gpar(col = scales::alpha(rl$col, NA), fontsize = 5 / PT * PT, fontfamily = FAM,
                                                  fontface = "bold", lineheight = 1.2)),
                                  t = pc$t, l = pc$l, clip = "off", name = "r_lines")
  }
  grid::grid.newpage(); grid::grid.draw(gt); invisible(dev.off())
}

PAIRS <- list(c("Maturity","DAM"), c("DAM","Interferon"))
plots <- list()
for (pr in PAIRS) for (mode in c("dev","inj")) {
  p <- mk(pr[1], pr[2], mode)
  ## panel e = Maturity x DAM, panel f = DAM x Interferon
  pl <- if (pr[1] == "Maturity") "e" else "f"
  fn <- file.path(OUT, sprintf("Fig2%s_scatter_%s_%s_vs_%s.png", pl, mode, pr[1], pr[2]))
  save_panel(p, file.path(OUT_PPTX, basename(fn)), RTXT[[paste(mode, pr[1], pr[2])]]); cat("wrote ", basename(fn), "\n", sep = "")
  plots[[length(plots) + 1]] <- p
}
## row order = Maturity x DAM (dev, inj), DAM x Interferon (dev, inj)
fwrite(rbindlist(CORR), file.path(OUT, "fig2_scatter_correlations.csv"))
cat("wrote fig2_scatter_correlations.csv\n"); print(rbindlist(CORR))
