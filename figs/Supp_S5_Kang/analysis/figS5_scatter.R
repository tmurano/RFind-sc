## Supplementary Fig. S5e, f — e: IFNb-index x PGE2-index (type A), f: Interferon-index x Chemokine-index (type B), within each
## cell type, drawn with the code of figs/Fig3_DeDonno/analysis/fig3_scatter.R (same marks, windows, summary panel and sizes).
## Changes: the Kang data (the eight author cell types with the abbreviations, order and colours of Fig. 3a of the original
## study; figS5_celltypes.R), conditions unstimulated vs IFN-beta, means per donor (>= 10 cells) instead of per mouse; in the
## summary panel the shape of each point gives its lineage.
## Input : data/Kang2018/output/scores_kang.csv
## Output: figs/Supp_S5_Kang/output/{FigS5e_scatter_IFNb_vs_PGE2.png, FigS5f_scatter_Interferon_vs_Chemokine.png,
##         figS5_scatter_axis_ranges.csv (the window of each small panel), figS5_scatter_outside.csv (cells outside it)}
## Run from the project root: R -f figs/Supp_S5_Kang/analysis/figS5_scatter.R
## Ticks: as fig3_scatter.R (about three, steps of 1, 2, 5), but only those inside the window, with more tried when fewer than
## two fall inside.
tick_in <- function(lo, hi) { for (m in 3:6) { b <- scales::breaks_extended(m, Q = c(1, 2, 5))(c(lo, hi)); b <- b[b >= lo & b <= hi]
  if (length(b) >= 2) return(b) }; b }
pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("figs/Supp_S5_Kang/analysis/figS5_celltypes.R")
set.seed(1)
OUT <- FIG

## rows in the order of the cell IDs (the order in which the cells are drawn is a random permutation of it)
d <- fread("data/Kang2018/output/scores_kang.csv")[, .(cell_id, cell_type, condition, donor_id = donor, IFNb, PGE2, Interferon, Chemokine)][order(cell_id)]
d[, condition := factor(COND_LAB[condition], levels = COND_LAB)]
cond_pal <- COND_PAL   # same colours as the condition map (figS5_ab.R)

DISP <- c(IFNb="IFNb-index", PGE2="PGE2-index", Interferon="Interferon-index", Chemokine="Chemokine-index")   # display names
PAIRS <- list(c("IFNb", "PGE2"), c("Interferon", "Chemokine"))   # e = type A pair, f = type B pair
PANEL <- c(IFNb_PGE2 = "FigS5e_scatter_IFNb_vs_PGE2", Interferon_Chemokine = "FigS5f_scatter_Interferon_vs_Chemokine")
MIN_CELLS_DONOR <- 10L   # as the violins (MIN_CELLS of figS5_violins.R)
## The individual cells are drawn in the per-cell-type panels; the cell types are compared with each other on one pair of
## axes in the last (summary) panel, by Cohen's d on the per-donor means (the numbers of the violins c, d).
## ---- drawn at the size pasted into the figure (W_MM x H_MM, 600 dpi) ----
## Arial, tick labels 6 pt, axis titles, panel titles, cell-type names and key 7 pt; line widths in pt (PT_LW). The sizes of
## cells and points are given for a 317.5-mm-wide drawing and scaled by SHRINK; the condition key is a row KEY_H high above the whole
## panel, with a 2-pt gap between key and name.
W_MM <- 87.6; H_MM <- 60.1                              # pasted size (mm)
KEY_H <- 3.2                                            # height of the row of the condition key (mm)
FAM <- "Arial"
SHRINK <- W_MM / 317.5                                  # scale of the point sizes below to the pasted width
## side of each point on which its name is written in the summary panels (r right, l left, a above, b below)
LAB_SIDE <- list(IFNb_PGE2 = c("cM"="a", "ncM"="b", "DC"="a", "B"="r", "T[H]"="b", "Tc"="l", "NK"="a", "Mkc"="a"),
                 Interferon_Chemokine = c("cM"="a", "ncM"="r", "DC"="r", "B"="r", "T[H]"="b", "Tc"="l", "NK"="b", "Mkc"="a"))
labs_sum <- function(d, side, rx, ry) {   # the points plus the origin (empty name); nx, ny = the start offset of each name (data units)
  d <- d[, .(x, y, ct = as.character(ct), lab = as.character(ct))]; s <- side[d$ct]
  d[, `:=`(nx = ifelse(s == "r", 0.12, ifelse(s == "l", -0.12, 0)) * diff(rx), ny = ifelse(s == "a", 0.09, ifelse(s == "b", -0.09, 0)) * diff(ry))]
  rbind(d, data.table(x = 0, y = 0, ct = NA_character_, lab = "''", nx = 0, ny = 0))   # names are plotmath (T[H]); '' = empty
}
Y_NAME_HJUST <- 1.2                                      # position of the y name along the lower block (0 bottom, 1 top)
PT_LW <- function(pt) pt / (.pt * 72.27 / 96)           # ggplot linewidth that gives `pt` points
th <- theme_minimal(base_size = 7, base_family = FAM) +
  theme(axis.title=element_text(size=7, colour="black"), axis.text=element_text(size=6, colour="grey30"),
        axis.title.x=element_text(margin=margin(t=1)), axis.title.y=element_text(margin=margin(r=1)),
        panel.grid=element_blank(), plot.margin=margin(1, 1, 1, 1))

FOCUS <- CT_LV   # all eight author cell types (as Kang et al.)
dfc <- d[cell_type %in% FOCUS]
dfc[, ct := factor(CT_LAB[cell_type], levels=unname(CT_LAB[FOCUS]))]
CT_PAL_LAB <- setNames(unname(CT_PAL[FOCUS]), unname(CT_LAB[FOCUS]))   # cell-type colours for the summary panel
## summary panel: the shape of each point says its lineage (the groups split by the dotted lines of the violins), with a
## small key inside the panel where it is empty (e bottom left, f top left)
LIN <- c(cM = "Myeloid leukocytes", ncM = "Myeloid leukocytes", DC = "Myeloid leukocytes", Mkc = "Megakaryocytes",
         `T[H]` = "Lymphocytes", B = "Lymphocytes", Tc = "Lymphocytes", NK = "Lymphocytes")
LIN_SHAPE <- c(`Myeloid leukocytes` = 21, Megakaryocytes = 23, Lymphocytes = 24)   # circle, diamond, triangle
LIN_KEY <- list(e = list(pos = c(0.01, 0.01), just = c(0, 0)), f = list(pos = c(0.01, 0.99), just = c(0, 1)))

mk_facet <- function(xa, ya){
  an <- dfc[, .(x=mean(get(xa)), y=mean(get(ya)), n=.N), by=.(ct, condition, donor_id)][n >= MIN_CELLS_DONOR]
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
  ## the per-donor means must always be inside
  am <- an[, .(mxlo=min(x), mxhi=max(x), mylo=min(y), myhi=max(y)), by=ct]
  q <- merge(q, am, by="ct")[, `:=`(xlo=pmin(xlo,mxlo), xhi=pmax(xhi,mxhi), ylo=pmin(ylo,mylo), yhi=pmax(yhi,myhi))]
  pts <- merge(dfc[, .(ct, condition, x=get(xa), y=get(ya))], q, by="ct")[x>=xlo & x<=xhi & y>=ylo & y<=yhi]
  cat(sprintf("%s x %s: %.1f%% of cells drawn; common width x %.1f y %.1f (each panel centred on its own median)\n",
              xa, ya, 100*nrow(pts)/nrow(dfc), wx$w, wy$w))
  ranges[[paste(xa, ya)]] <<- q[, .(pair=paste(DISP[xa],"x",DISP[ya]), ct, xlo=round(xlo,2), xhi=round(xhi,2),
                                    ylo=round(ylo,2), yhi=round(yhi,2))]
  ## cells outside the window of their panel (not drawn; their share is given in the legend)
  no <- merge(dfc[, .(n_cells = .N), by = ct], pts[, .(n_drawn = .N), by = ct], by = "ct", all.x = TRUE)[is.na(n_drawn), n_drawn := 0L]
  no <- rbind(no[order(ct), .(ct = as.character(ct), n_cells, n_outside = n_cells - n_drawn)],
              data.table(ct = "all", n_cells = sum(no$n_cells), n_outside = sum(no$n_cells - no$n_drawn)))
  outside[[paste(xa, ya)]] <<- no[, .(pair = paste(DISP[xa], "x", DISP[ya]), ct, n_cells, n_outside,
                                      pct_outside = round(100 * n_outside / n_cells, 1))]
  ## per-cell-type panels + a summary panel. The per-cell-type panels show the cells;
  ## the summary puts every cell type on ONE pair of axes so that they can be compared with each other.
  one <- function(k, i) {   # i = slot 1-8 (A-H of the design below)
    qq <- q[ct == k]
    ggplot(pts[ct == k][sample(.N)], aes(x, y)) +
      geom_hline(yintercept=0, linetype="dashed", colour="grey75", linewidth=PT_LW(0.35)) +
      annotate("segment", x=0, xend=0, y=-Inf, yend=Inf, linetype="dashed", colour="grey75", linewidth=PT_LW(0.35)) +
      geom_point(aes(colour=condition), size=0.45 * SHRINK, alpha=0.32, stroke=0) +             # cells
      ## contours of each condition's cell cloud: they show the shift of the whole distribution of the cells, which is
      ## small compared with the spread of the cells
      stat_density_2d(aes(colour=condition), linewidth=PT_LW(0.45), bins=4, contour_var="ndensity",
                      breaks=c(0.3, 0.6), show.legend=FALSE) +
      geom_point(data=an[ct == k], aes(x, y, fill=condition), shape=21, colour="grey10", size=3.2 * SHRINK,
                 stroke=0.6 * SHRINK, inherit.aes=FALSE) +
      coord_cartesian(xlim=c(qq$xlo, qq$xhi), ylim=c(qq$ylo, qq$yhi)) +
      scale_colour_manual(values=cond_pal, labels=COND_SHOW, name=NULL) + scale_fill_manual(values=cond_pal, guide="none") +
      guides(colour=guide_legend(override.aes=list(size=2, alpha=1))) +
      ## about three ticks computed on the visible window (n.breaks on the scale would use the full data range and leave one tick)
      scale_x_continuous(breaks = tick_in(qq$xlo, qq$xhi)) +   # steps of 1, 2, 5 (no 2.5)
      scale_y_continuous(breaks = tick_in(qq$ylo, qq$yhi)) +
      ## axis names once for all small panels: y on the left column of the lower block (E, G; collected into one), x under
      ## its lowest row (G, H; collected into one); none in the top row
      labs(x = if (i %in% c(7, 8)) DISP[xa] else NULL, y = if (i %in% c(5, 7)) DISP[ya] else NULL,
           title=parse(text = paste0("bold(", as.character(k), ")"))[[1]]) + th +
      theme(plot.title=element_text(face="bold", size=7, hjust=0.5, margin=margin(b=1)),
            ## the one y name sits at the middle height of all small panels (top row + lower block), not of the lower block
            axis.title.y=element_text(hjust=Y_NAME_HJUST, margin=margin(r=1)),
            legend.position = "none")   # the condition key is drawn once, in its own row above the whole panel
  }
  ## summary panel: Cohen's d of each cell type on the per-donor means (the numbers of the violins), not the raw difference,
  ## so that all cell types share one pair of axes
  cd <- function(a, b) (mean(b) - mean(a)) / sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  sh <- an[, .(x = cd(x[condition == "Control"], x[condition == "IFN-β"]),
               y = cd(y[condition == "Control"], y[condition == "IFN-β"]), condition = "IFN-β"), by = ct]
  print(sh[, .(ct, d_x = round(x, 2), d_y = round(y, 2))])
  side <- LAB_SIDE[[paste(xa, ya, sep = "_")]]
  ld <- labs_sum(sh[condition=="IFN-β"], side, range(c(0, sh$x)), range(c(0, sh$y)))
  sh[, lin := factor(LIN[as.character(ct)], levels = names(LIN_SHAPE))]
  kp <- if (xa == "IFNb") LIN_KEY$e else LIN_KEY$f
  p8 <- ggplot(sh[condition=="IFN-β"], aes(x, y)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey75", linewidth=PT_LW(0.35)) +
    geom_vline(xintercept=0, linetype="dashed", colour="grey75", linewidth=PT_LW(0.35)) +
    geom_point(aes(fill=ct, shape=lin), colour="grey20", size=5.8 * SHRINK, stroke=0.8 * SHRINK) +
    geom_point(data=data.table(x=0, y=0), aes(x, y), shape=21, size=5.8 * SHRINK, stroke=1.0 * SHRINK,
               fill="white", colour="grey35", inherit.aes=FALSE) +
    ## each name starts on a fixed side of its point (LAB_SIDE; free placement would leave some names nearer to other points
    ## than to their own) and is then moved apart by the repel layer; a short grey line ties it to its point;
    ## the origin (white point) is passed with an empty name so that no name covers it
    ggrepel::geom_text_repel(data=ld, aes(label=paste0("bold(", lab, ")"), colour=ct), parse=TRUE, nudge_x=ld$nx, nudge_y=ld$ny, family=FAM, size=7 / .pt,
                             seed=1, box.padding=0.3 * SHRINK, point.padding=0.3 * SHRINK, min.segment.length=0,
                             segment.colour="grey60", segment.size=PT_LW(0.35), max.overlaps=Inf, show.legend=FALSE,
                             max.time=5, max.iter=1e5) +
    scale_x_continuous(breaks = scales::breaks_extended(5)(range(c(0, sh$x))), expand = expansion(mult = c(0.22, 0.05))) +   # room left of the leftmost points for their names; ticks as without the room
    scale_fill_manual(values=CT_PAL_LAB, guide="none") + scale_colour_manual(values=CT_PAL_LAB, guide="none") +
    scale_shape_manual(values=LIN_SHAPE, breaks=names(LIN_SHAPE), name=NULL, drop=FALSE) +
    guides(shape=guide_legend(override.aes=list(fill="grey80", colour="grey20", size=5.8 * SHRINK, stroke=0.8 * SHRINK))) +
    ## no title; the axis names say what is shown; "Cohen's d," rather than "d," so that it is not taken for a panel letter
    labs(x = paste("Cohen's d,", DISP[xa]), y = paste("Cohen's d,", DISP[ya])) + th +
    theme(plot.margin=margin(t=1, r=4, b=1, l=1),   # room for the last x tick label
          legend.position="inside", legend.position.inside=kp$pos, legend.justification=kp$just,
          legend.text=element_text(size=6, margin=margin(l=1)), legend.key.size=unit(2.3, "mm"),
          legend.key.spacing.y=unit(0, "mm"), legend.background=element_blank(), legend.key=element_blank(),
          legend.margin=margin(0, 0, 0, 0))
  ## layout: the condition key in a row at the top; 4 small panels below it; below them, 2 x 2 small panels on
  ## the left and the summary spanning two rows and the rest of the width on the right
  lv <- levels(droplevels(q$ct)); small <- lapply(seq_along(lv), function(i) one(lv[i], i)); names(small) <- LETTERS[seq_along(small)]
  ## nested so that a gap separates the summary from the small panels without opening gaps between the small panels of
  ## the top row; axis names collected within the lower block (one x name under its bottom row, one y name left of it);
  ## no outer margin around the nested blocks or the whole panel (it is pasted at its drawn size)
  NOMARG <- plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
  ZM <- theme(plot.margin = margin(0, 0, 0, 0))   # the key slot and the spacer would otherwise add the default 5.5-pt margins
  ## the condition key: one line at the top right of the whole panel, above the names of the cell types (the names
  ## "Unstimulated" and "IFN-β, 6 h" are too long for a slot right of the top row, which would narrow its panels)
  key <- ggplot(data.table(condition = factor(names(cond_pal), levels = names(cond_pal)), x = 0, y = 0), aes(x, y, colour = condition)) +
    geom_point(alpha = 0) + scale_colour_manual(values = cond_pal, labels = COND_SHOW, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) + theme_void(base_family = FAM) +
    theme(legend.position = "inside", legend.position.inside = c(1, 0.5), legend.justification = c(1, 0.5),
          legend.direction = "horizontal", legend.text = element_text(size = 7, margin = margin(l = 2)),   # 2-pt gap key-name
          legend.key.size = unit(2.6, "mm"), legend.key.spacing.x = unit(6, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(0, 0, 0, 0))
  top <- wrap_plots(small[1:4], nrow = 1) + plot_layout(axis_titles = "collect") + NOMARG
  lower_left <- wrap_plots(small[5:8], ncol = 2) + plot_layout(axis_titles = "collect") + NOMARG
  bottom <- wrap_plots(list(lower_left, plot_spacer() + ZM, p8), nrow = 1, widths = c(2.25, 0.15, 2.15)) + NOMARG
  key / top / bottom + plot_layout(heights = unit(c(KEY_H, 1, 2), c("mm", "null", "null"))) + NOMARG
}

px <- function(mm) round(mm / 25.4 * 600)                            # ragg truncates mm sizes, so give whole pixels
ranges <- list(); outside <- list()
for (pr in PAIRS) {
  xa <- pr[1]; ya <- pr[2]; key <- paste(xa, ya, sep = "_")
  p <- mk_facet(xa, ya)
  ragg::agg_png(file.path(FIG, paste0(PANEL[[key]], ".png")), width = px(W_MM), height = px(H_MM), units = "px",
                res = 600, background = "white")
  print(p); invisible(dev.off())
  cat("wrote ", PANEL[[key]], "\n", sep = "")
}
fwrite(rbindlist(ranges), file.path(OUT, "figS5_scatter_axis_ranges.csv"))
fwrite(rbindlist(outside), file.path(OUT, "figS5_scatter_outside.csv"))
print(rbindlist(outside))
