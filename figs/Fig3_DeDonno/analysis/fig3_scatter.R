## Fig. 3e, f (De Donno) — per-cell 2D maps of two indices within each cell type (Ctrl line), in the style of Fig. 2e, f.
##   Panels: 3e = Glucocorticoid x Seizure, 3f = NF-kB x UPR (the data keys stay CORT and KA6h; the output file names use
##   the index names).
##   Small panels = the authors' seven major cell types (DG, CA1, CA3, Astro, Micro, Oligo, Endo) + GABA: four on top with
##   the colour key at the right end, 2 x 2 below next to the summary panel; cells coloured by condition.
##   Marks: cells + contours of each condition (30% and 60% of the maximum cell density) + outlined points = mean of each
##   mouse (>= 10 cells of that type); one No Stress / Stress colour key. No condition-mean points.
##   r: SHOW_R can write Pearson's r between the two indices within each condition (all cells of the cell type, the values
##   of fig3_axis_correlations.csv) above the cells of each small panel; it is off, as neither Results 3 nor the
##   Discussion cites the within-condition r.
##   Axes: on each axis every small panel has the same width (the median over the cell types of the width of their cells'
##   5th-95th percentile range), centred on the cell type's own median and widened where needed so that all per-mouse
##   means are inside. Cells outside the range are not drawn. Ranges written to fig3_scatter_axis_ranges.csv.
##   Summary panel: each cell type at its Cohen's d (Stress vs No Stress, per-mouse means) on the two indices.
##   r (within each condition, and pooled) is also written to fig3_axis_correlations.csv.
##   No titles or explanatory text on the figure; mark meanings go in manuscript/Legends.md.
##
## Input : data/DeDonno2026/output/scores_dedonno.csv
## Output: figs/Fig3_DeDonno/output/fig3_axis_correlations.csv + fig3_scatter_axis_ranges.csv
##         + {Fig3e_scatter_Glucocorticoid_vs_Seizure, Fig3f_scatter_NFkB_vs_UPR}.png (drawn at the pasted size)
## Run   : Rscript "figs/Fig3_DeDonno/analysis/fig3_scatter.R"   (from project root)

pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("figs/Fig3_DeDonno/analysis/fig3_celltypes.R")
set.seed(1)
OUT <- "figs/Fig3_DeDonno/output"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

d <- fread("data/DeDonno2026/output/scores_dedonno.csv")[line == "Ctrl"]
d[, condition := factor(condition, levels=c("No Stress","Stress"))]
cond_pal <- c("No Stress"="#4C78A8", "Stress"="#E45756")   # same colours as the condition map of Fig. 3a

DISP <- c(CORT="Glucocorticoid-index", KA6h="Seizure-index", NFkB="NF-κB-index", UPR="UPR-index")   # display names
PAIRS <- list(c("CORT","KA6h"), c("NFkB","UPR"))   # panels e and f
PANEL <- c(CORT_KA6h = "Fig3e_scatter_Glucocorticoid_vs_Seizure", NFkB_UPR = "Fig3f_scatter_NFkB_vs_UPR")
SUMMARY_STAT <- "d"   # summary panel: "d" = Cohen's d on per-mouse means, the numbers of the violins c, d (an effect size,
                      # so that all cell types share one pair of axes), "delta" = Stress minus No Stress mean
MIN_CELLS_ANIMAL <- 10L   # as the violins (MIN_CELLS of fig3_violins.R)
SHOW_R <- c(CORT_KA6h = FALSE, NFkB_UPR = FALSE)  # within-condition r on the small panels: off (the text does not cite it)
R_BAND <- 0.24            # height of the band above the cells that holds the r values, as a fraction of the y window
## The individual cells are drawn in the per-cell-type panels; the cell types are compared with each other on one pair of
## axes in the last (summary) panel.
## ZOOM, ZPAD, OFF, GRP and align() below are not used by the panels drawn: the ranges are set by win() in mk_facet.
ZOOM <- c("KA6h", "NFkB")   # these two axes are zoomed to the per-mouse means; the others follow the cells
ZPAD <- 0.4                 # padding on each side, as a fraction of the range of the per-mouse means
OFF  <- 0.5                 # a cell type is "far off" when its range centre is > OFF x median width from its group's median centre
GRP <- c(DG = "neuron", CA1 = "neuron", CA3 = "neuron", GABA = "neuron", Astro = "glia_vascular", Micro = "glia_vascular",
         Oligo = "glia_vascular", Endo = "glia_vascular")
## ---- drawn at the size pasted into Fig. 3 (W_MM x H_MM, 600 dpi) ----
## Arial, tick labels 6 pt, axis titles, panel titles, cell-type names and key 7 pt; line widths in pt (PT_LW). The sizes of
## cells and points are given for a 317.5-mm-wide drawing and scaled by SHRINK; the key slot is KEY_W wide, with a 2-pt gap
## between key and name.
W_MM <- 87.6; H_MM <- 60.1                              # pasted size (mm)
FAM <- "Arial"
SHRINK <- W_MM / 317.5                                  # scale of the point sizes below to the pasted width
## side of each point on which its name is written in the summary panels (r right, l left, a above, b below)
LAB_SIDE <- list(
  CORT_KA6h = c("DG Glut"="a", "CA1 Glut"="r", "CA3 Glut"="l", "GABA"="r", "Astro"="a", "Microglia"="a",
                "Oligodendrocytes"="a", "Endothelial"="b"),
  NFkB_UPR  = c("DG Glut"="a", "CA1 Glut"="l", "CA3 Glut"="b", "GABA"="r", "Astro"="b", "Microglia"="a",
                "Oligodendrocytes"="b", "Endothelial"="l"))
labs_sum <- function(d, side, rx, ry) {   # the points plus the origin (empty name); nx, ny = the start offset of each name (data units)
  d <- d[, .(x, y, ct = as.character(ct), lab = as.character(ct))]; s <- side[d$ct]
  d[, `:=`(nx = ifelse(s == "r", 0.12, ifelse(s == "l", -0.12, 0)) * diff(rx), ny = ifelse(s == "a", 0.09, ifelse(s == "b", -0.09, 0)) * diff(ry))]
  rbind(d, data.table(x = 0, y = 0, ct = NA_character_, lab = "", nx = 0, ny = 0))
}
Y_NAME_HJUST <- 1.2                                      # position of the y name along the lower block (0 bottom, 1 top)
KEY_W <- 0.8                                            # width of the key slot right of the top row (relative to a panel)
PT_LW <- function(pt) pt / (.pt * 72.27 / 96)           # ggplot linewidth that gives `pt` points
th <- theme_minimal(base_size = 7, base_family = FAM) +
  theme(axis.title=element_text(size=7, colour="black"), axis.text=element_text(size=6, colour="grey30"),
        axis.title.x=element_text(margin=margin(t=1)), axis.title.y=element_text(margin=margin(r=1)),
        panel.grid=element_blank(), plot.margin=margin(1, 1, 1, 1))

FOCUS <- c("DG Glut","CA1 Glut","CA3 Glut","GABA","Astro","Microglia","Oligodendrocytes","Endothelial")
dfc <- d[cell_type %in% FOCUS]
dfc[, ct := factor(CT_LAB[cell_type], levels=unname(CT_LAB[FOCUS]))]
CT_PAL_LAB <- setNames(unname(CT_PAL[FOCUS]), unname(CT_LAB[FOCUS]))   # cell-type colours for the summary panel

## per-cell-type ranges -> aligned ranges: common limits, except cell types far off, which keep their centre and the width
align <- function(lo, hi) {
  mid <- (lo + hi) / 2; off <- abs(mid - median(mid)) > OFF * median(hi - lo)
  clo <- min(lo[!off]); chi <- max(hi[!off]); w <- chi - clo
  list(lo = ifelse(off, pmin(lo, mid - w / 2), clo), hi = ifelse(off, pmax(hi, mid + w / 2), chi), off = off)
}

mk_facet <- function(xa, ya, show_r = FALSE){
  an <- dfc[, .(x=mean(get(xa)), y=mean(get(ya)), n=.N), by=.(ct, condition, donor_id)][n >= MIN_CELLS_ANIMAL]
  ## Every panel uses the SAME WIDTH on each axis, so that distances are comparable between panels; the window of each
  ## cell type is centred on its own median. The cell types are compared with each other in the summary panel.
  win <- function(v, lab) {
    ## width = median over the cell types of their 5th-95th percentile range; window POSITION = each cell type's median.
    ## A window common to all cell types would have to span the union of their ranges, much wider than the range of a
    ## single cell type; equal widths keep one scale for all panels, so displacements can be compared, while the
    ## absolute positions are not compared across cell types (each cell type has its own index offset).
    Q <- function(pr) vapply(split(v, lab), function(z) unname(quantile(z, pr)), 0)
    lo1 <- Q(.05); hi1 <- Q(.95); med <- Q(.5)
    w <- unname(median(hi1 - lo1))
    list(lo = med - w/2, hi = med + w/2, w = w)
  }
  wx <- win(dfc[[xa]], as.character(dfc$ct)); wy <- win(dfc[[ya]], as.character(dfc$ct))
  q <- data.table(ct = names(wx$lo), xlo = wx$lo, xhi = wx$hi, ylo = wy$lo[names(wx$lo)], yhi = wy$hi[names(wx$lo)])
  q[, ct := factor(ct, levels = levels(dfc$ct))]
  ## the per-mouse means must always be inside
  am <- an[, .(mxlo=min(x), mxhi=max(x), mylo=min(y), myhi=max(y)), by=ct]
  q <- merge(q, am, by="ct")[, `:=`(xlo=pmin(xlo,mxlo), xhi=pmax(xhi,mxhi), ylo=pmin(ylo,mylo), yhi=pmax(yhi,myhi))]
  pts <- merge(dfc[, .(ct, condition, x=get(xa), y=get(ya))], q, by="ct")[x>=xlo & x<=xhi & y>=ylo & y<=yhi]
  cat(sprintf("%s x %s: %.1f%% of cells drawn; common width x %.1f y %.1f (each panel centred on its own median)\n",
              xa, ya, 100*nrow(pts)/nrow(dfc), wx$w, wy$w))
  ranges[[paste(xa, ya)]] <<- q[, .(pair=paste(DISP[xa],"x",DISP[ya]), ct, xlo=round(xlo,2), xhi=round(xhi,2),
                                    ylo=round(ylo,2), yhi=round(yhi,2))]
  ## per-cell-type panels + a summary panel. The per-cell-type panels show the cells;
  ## the summary puts every cell type on ONE pair of axes so that they can be compared with each other.
  one <- function(k, i) {   # i = slot 1-8 (A-H of the design below)
    qq <- q[ct == k]
    ## with show_r, a band above the cells holds the within-condition r (all cells of the type, not only those drawn);
    ## the vertical zero line stops below the band so that it does not cross the numbers
    ytop <- if (show_r) qq$yhi + R_BAND * (qq$yhi - qq$ylo) else qq$yhi
    rr <- dfc[ct == k, .(r = cor(get(xa), get(ya))), by = condition]
    fmt_r <- function(r) sub("-", "−", sprintf("%+.2f", r))
    xmid <- (qq$xlo + qq$xhi) / 2; gap <- 0.03 * (qq$xhi - qq$xlo); ty <- qq$yhi + 0.5 * R_BAND * (qq$yhi - qq$ylo)
    ggplot(pts[ct == k][sample(.N)], aes(x, y)) +
      geom_hline(yintercept=0, linetype="dashed", colour="grey75", linewidth=PT_LW(0.35)) +
      annotate("segment", x=0, xend=0, y=-Inf, yend=if (show_r) qq$yhi else Inf, linetype="dashed", colour="grey75",
               linewidth=PT_LW(0.35)) +
      { if (show_r) list(
          annotate("text", x=xmid - gap, y=ty, label=fmt_r(rr[condition == "No Stress", r]), hjust=1, vjust=0.5,
                   colour=cond_pal[["No Stress"]], size=6 / .pt, family=FAM),
          annotate("text", x=xmid + gap, y=ty, label=fmt_r(rr[condition == "Stress", r]), hjust=0, vjust=0.5,
                   colour=cond_pal[["Stress"]], size=6 / .pt, family=FAM)) } +
      geom_point(aes(colour=condition), size=0.45 * SHRINK, alpha=0.32, stroke=0) +             # cells
      ## contours of each condition's cell cloud: they show the shift of the whole distribution of the cells, which is
      ## small compared with the spread of the cells
      stat_density_2d(aes(colour=condition), linewidth=PT_LW(0.45), bins=4, contour_var="ndensity",
                      breaks=c(0.3, 0.6), show.legend=FALSE) +
      geom_point(data=an[ct == k], aes(x, y, fill=condition), shape=21, colour="grey10", size=3.2 * SHRINK,
                 stroke=0.6 * SHRINK, inherit.aes=FALSE) +
      coord_cartesian(xlim=c(qq$xlo, qq$xhi), ylim=c(qq$ylo, ytop)) +
      scale_colour_manual(values=cond_pal, name=NULL) + scale_fill_manual(values=cond_pal, guide="none") +
      guides(colour=guide_legend(override.aes=list(size=2, alpha=1))) +
      ## about three ticks computed on the visible window (n.breaks on the scale would use the full data range and leave one tick)
      scale_x_continuous(breaks = scales::breaks_extended(3, Q = c(1, 2, 5))(c(qq$xlo, qq$xhi))) +   # steps of 1, 2, 5 (no 2.5)
      ## ticks only inside the cell window (with show_r, none in the band of the r values)
      scale_y_continuous(breaks = { b <- scales::breaks_extended(3, Q = c(1, 2, 5))(c(qq$ylo, qq$yhi))
                                    if (show_r) b[b >= qq$ylo & b <= qq$yhi] else b }) +
      ## axis names once for all small panels: y on the left column of the lower block (E, G; collected into one), x under
      ## its lowest row (G, H; collected into one); none in the top row
      labs(x = if (i %in% c(7, 8)) DISP[xa] else NULL, y = if (i %in% c(5, 7)) DISP[ya] else NULL,
           title=as.character(k)) + th +
      theme(plot.title=element_text(face="bold", size=7, hjust=0.5, margin=margin(b=1)),
            ## the one y name sits at the middle height of all small panels (top row + lower block), not of the lower block
            axis.title.y=element_text(hjust=Y_NAME_HJUST, margin=margin(r=1)),
            legend.position = if (i <= 4) "right" else "none")   # the key is collected from the top row only
  }
  ## summary panel: each cell type's stress-induced shift, centred on ITS OWN control mean. The centring is required —
  ## the absolute index value carries a per-cell-type offset, so DG/control and
  ## microglia/control cannot be compared directly; only the displacements can.
  gm  <- an[, .(x=mean(x), y=mean(y)), by=.(ct, condition)]
  sh  <- merge(gm, gm[condition=="No Stress", .(ct, x0=x, y0=y)], by="ct")[, .(ct, condition, x=x-x0, y=y-y0)]
  if (SUMMARY_STAT == "d") {   # Cohen's d on the per-mouse means (the numbers of the violins), not the raw difference
    cd <- function(a, b) (mean(b) - mean(a)) / sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
    sh <- an[, .(x = cd(x[condition == "No Stress"], x[condition == "Stress"]),
                 y = cd(y[condition == "No Stress"], y[condition == "Stress"]), condition = "Stress"), by = ct]
    print(sh[, .(ct, d_x = round(x, 2), d_y = round(y, 2))])
  }
  ld <- labs_sum(sh[condition=="Stress"], LAB_SIDE[[paste(xa, ya, sep = "_")]], range(c(0, sh$x)), range(c(0, sh$y)))
  p8 <- ggplot(sh[condition=="Stress"], aes(x, y)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey75", linewidth=PT_LW(0.35)) +
    geom_vline(xintercept=0, linetype="dashed", colour="grey75", linewidth=PT_LW(0.35)) +
    geom_point(aes(fill=ct), shape=21, colour="grey20", size=5.8 * SHRINK, stroke=0.8 * SHRINK, show.legend=FALSE) +
    geom_point(data=data.table(x=0, y=0), aes(x, y), shape=21, size=5.8 * SHRINK, stroke=1.0 * SHRINK,
               fill="white", colour="grey35", inherit.aes=FALSE) +
    ## each name starts on a fixed side of its point (LAB_SIDE; free placement would leave some names nearer to other points
    ## than to their own) and is then moved apart by the repel layer; a short grey line ties it to its point;
    ## the origin (white point) is passed with an empty name so that no name covers it
    ggrepel::geom_text_repel(data=ld, aes(label=lab, colour=ct), nudge_x=ld$nx, nudge_y=ld$ny, family=FAM, fontface="bold", size=7 / .pt,
                             seed=1, box.padding=0.3 * SHRINK, point.padding=0.3 * SHRINK, min.segment.length=0,
                             segment.colour="grey60", segment.size=PT_LW(0.35), max.overlaps=Inf, show.legend=FALSE,
                             max.time=5, max.iter=1e5) +
    scale_x_continuous(breaks = scales::breaks_extended(5)(range(c(0, sh$x))), expand = expansion(mult = c(0.22, 0.05))) +   # room left of the leftmost points for their names; ticks as without the room
    scale_fill_manual(values=CT_PAL_LAB) + scale_colour_manual(values=CT_PAL_LAB) +
    ## no title; the axis names say what is shown ("Cohen's d, <index>", or "Δ <index>" with SUMMARY_STAT = "delta");
    ## "Cohen's d," rather than "d," so that it is not taken for a panel letter
    labs(x = if (SUMMARY_STAT == "d") paste("Cohen's d,", DISP[xa]) else paste("\u0394", DISP[xa]),
         y = if (SUMMARY_STAT == "d") paste("Cohen's d,", DISP[ya]) else paste("\u0394", DISP[ya])) + th +
    theme(plot.margin=margin(t=1, r=4, b=1, l=1))   # room for the last x tick label
  ## layout: 4 small panels on top with the colour key in a narrow slot at the right end; below, 2 x 2 small panels on
  ## the left and the summary spanning two rows and the rest of the width on the right
  lv <- levels(droplevels(q$ct)); small <- lapply(seq_along(lv), function(i) one(lv[i], i)); names(small) <- LETTERS[seq_along(small)]
  ## nested so that a gap separates the summary from the small panels without opening gaps between the small panels of
  ## the top row; axis names collected within the lower block (one x name under Oligo-Endo, one y name left of it);
  ## no outer margin around the nested blocks or the whole panel (it is pasted at its drawn size)
  NOMARG <- plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
  ZM <- theme(plot.margin = margin(0, 0, 0, 0))   # the key slot and the spacer would otherwise add the default 5.5-pt margins
  top <- wrap_plots(c(small[1:4], list(guide_area() + ZM)), nrow = 1, widths = c(1, 1, 1, 1, KEY_W)) +
    plot_layout(guides = "collect", axis_titles = "collect") + NOMARG
  lower_left <- wrap_plots(small[5:8], ncol = 2) + plot_layout(axis_titles = "collect") + NOMARG
  bottom <- wrap_plots(list(lower_left, plot_spacer() + ZM, p8), nrow = 1, widths = c(2.25, 0.15, 2.15)) + NOMARG   # 2.25: room for the title 'Oligodendrocytes'
  (top / bottom + plot_layout(heights = c(1, 2)) + NOMARG) &
    theme(legend.direction = "vertical", legend.text = element_text(size = 7, margin = margin(l = 2)),   # 2-pt gap key-name
          legend.key.size = unit(2.6, "mm"), legend.key.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0))
}

OUT_PPTX <- OUT   # pasted-size panels go to output/
px <- function(mm) round(mm / 25.4 * 600)                            # ragg truncates mm sizes, so give whole pixels
corr_rows <- list(); plots <- list(); ranges <- list()
for (pr in PAIRS) {
  xa <- pr[1]; ya <- pr[2]; key <- paste(xa, ya, sep = "_")
  p <- mk_facet(xa, ya, show_r = SHOW_R[[key]]); plots[[key]] <- p
  ragg::agg_png(file.path(OUT_PPTX, paste0(PANEL[[key]], ".png")), width = px(W_MM), height = px(H_MM), units = "px",
                res = 600, background = "white")
  print(p); invisible(dev.off())
  cat("wrote ", PANEL[[key]], "\n", sep = "")
  for (ctv in FOCUS) { z <- dfc[cell_type == ctv]
    rc <- function(s) round(cor(s[[xa]], s[[ya]]), 3)
    corr_rows[[paste(key, ctv)]] <- data.table(pair=sprintf("%s_vs_%s", DISP[xa], DISP[ya]), cell_type=CT_LAB[[ctv]],
      n=nrow(z), r_NoStress=rc(z[condition == "No Stress"]), r_Stress=rc(z[condition == "Stress"]),
      r_pooled=rc(z), spearman_pooled=round(cor(z[[xa]], z[[ya]], method="spearman"), 3)) }
}
fwrite(rbindlist(corr_rows), file.path(OUT, "fig3_axis_correlations.csv"))
fwrite(rbindlist(ranges), file.path(OUT, "fig3_scatter_axis_ranges.csv"))
cat("\nwithin-cell-type correlations (not drawn; reference only):\n"); print(rbindlist(corr_rows))
