## Fig. 2c, d (Hammond) — RFind-sc index violins by Hammond cluster, with the effect size d0 above each violin.
##   QC'd data (our QC, scores_hammond_9axes.csv); the 13 microglial author clusters (Mono/Mac not shown).
##   Number above each violin = d0 = mean/SD of the cluster's scores on that axis (one-sample Cohen's d against 0, the
##     neutral point of a signed RFind score; see the stats block below). With ~74k cells a cell-level p would be pseudoreplicated and
##     uninformative, so the effect size is shown; the separation of cell groups is evaluated with bootstrap confidence
##     intervals in the Fig. 5 benchmark. Full table (with cluster-vs-rest d and AUC) -> output/Fig2_stats_table.csv.
##   Panels: 2c = Maturity-index [Matcovitch-Natan, high = mature] and DAM-index [Keren-Shaul] (ranked up/down lists, type A
##     of the paper; stored as type "B", which is also the type column of Fig2_stats_table.csv) |
##     2d = Interferon-index and G2M-index (MSigDB Hallmark gene sets, type B of the paper; stored as type "C").
##   Style: clean theme, thin violin + white median box, uniform annotation size.
##
## Input : data/Hammond2019/output/scores_hammond_9axes.csv
## Output: figs/Fig2_Hammond/output/Fig2c_violin_Maturity_DAM.png (panel 2c), Fig2d_violin_Interferon_G2M.png (panel 2d)
##         (83.8 x 78.6 mm each, 600 dpi = the size pasted in the figure)
##         + Fig2_stats_table.csv, Fig2_cluster_ncells.csv
## Run   : Rscript "figs/Fig2_Hammond/analysis/fig2_violins.R"   (from project root)

pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
BASE <- "figs/Fig2_Hammond"; DATA <- file.path(BASE,"data"); OUT <- file.path(BASE,"output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

d <- fread("data/Hammond2019/output/scores_hammond_9axes.csv")
## Mono/Mac (monocytes/macrophages, not microglia) is not shown.
ord <- c("1","2a","2b","2c","3","4","5","6","7a","7b","7c","8","9")   # 13 microglial author clusters
d <- d[cluster %in% ord]; d[, cluster := factor(cluster, levels = ord)]
## per-mouse sample ID (= Sample_GEO_Prefix; 41 mice). Animal-level points = per-mouse mean of the
## cells that mouse contributes to a cluster, drawn only when it contributes >= MIN_CELLS_ANIMAL cells.
d[, animal := sub("_[^_]+$", "", cell_id)]
MIN_CELLS_ANIMAL <- 20L
## record per-cluster cell counts -> Fig2_cluster_ncells.csv
NC <- d[, .(n_cells=.N), by=cluster][order(cluster)]
fwrite(NC, file.path(OUT, "Fig2_cluster_ncells.csv"))
cat(sprintf("%d cells across %d clusters (Mono/Mac excluded)\n", nrow(d), length(ord)))

## ---- GLOBAL per-cluster signature colours + fixed identity labels ----
## A violin is coloured in any panel iff BOTH: (1) the cluster has a signature colour, and
## (2) it passes the d0 gate of that panel (see mk below). Its colour is the cluster's fixed signature colour (same across panels).
## Clusters with NO signature colour (1, 5, 6, 8) stay grey everywhere, regardless of d0.
##   cl2-series (2a/2b/2c)=red · cl3=orange(metabolically active) · cl4=blue · cl7-series(7a/7b/7c)=black · cl9=green
CRED<-"#E41A1C"; CBLUE<-"#2C7FB8"; CBLACK<-"#000000"; CGREEN<-"#2CA02C"; CORANGE<-"#E69F00"
COL_OTHER <- "#D9D9D9"   # clusters without a signature colour = light grey
## pasted-size rendering: each sheet is drawn at the exact size it is pasted in the figure
## (83.8 x 78.6 mm; 600 dpi, dpi stored in the PNG), Arial; tick labels 6 pt,
## axis titles / title / legend 7 pt, d0 numbers 7 pt, reference subtitle and caption 6 pt.
## Title and subtitle start at the left edge of the plot (not of the panel) so the reference line fits the width.
FAM <- "Arial"; PT <- ggplot2::.pt
LW05 <- 0.5 / (PT * 72.27 / 96)   # linewidth for a true 0.5 pt line (grid lwd is in 1/96 inch, so 0.5/.pt alone gives ~0.375 pt)
W_MM <- 83.8; H_MM <- 78.6   # size of c and d in the figure (mm)
OUT_PPTX <- OUT   # pasted-size panels go to output/
LAB_SIZE <- 7 / PT   # the d0 numbers above the violins: 7 pt, the size of the numbers of Fig. 3c, d (in-panel text)
SIG <- list(
  list(cl=c("2a","2b","2c"), col=CRED,    lab="2a/2b/2c (proliferative)"),   # Hammond et al.'s word
  list(cl="3",               col=CORANGE, lab="cluster 3 (metabolically active)"),   # the names of Fig. 2b
  list(cl="4",               col=CBLUE,   lab="cluster 4 (ATM)"),   # cluster 4 = P4/P5 axon-tract-associated microglia
  list(cl=c("7a","7b","7c"), col=CBLACK,  lab="7a/7b/7c (juvenile/adult)"),
  list(cl="9",               col=CGREEN,  lab="cluster 9 (white-matter injury)"))
SIG_LAB    <- unlist(lapply(SIG, function(g) setNames(rep(g$lab, length(g$cl)), g$cl)))  # cluster -> identity label
LAB_LEVELS <- vapply(SIG, `[[`, "", "lab")                                              # legend order
LAB_COL    <- setNames(vapply(SIG, `[[`, "", "col"), LAB_LEVELS)                          # label -> colour
has_sig    <- function(cl) as.character(cl) %in% names(SIG_LAB)
## axis metadata (colour is the global rule above; no per-axis colour groups).
## nm = axis/y-axis display name. type = as stored in the panel: "B" = ranked up/down list (type A of the paper),
## "C" = gene set (type B of the paper). src = reference paper (ranked lists). hallmark = MSigDB set (gene sets).
## dir = case/control contrast for a subtitle (not set for these four axes).
AX <- list(
  Maturity   =list(type="B", nm="Maturity-index", src="Matcovitch-Natan et al., Science 2016 (adult vs embryonic; high = mature)"),
  DAM        =list(type="B", nm="DAM-index",  src="Keren-Shaul et al., Cell 2017"),
  Interferon =list(type="C", nm="Interferon-index", hallmark="INTERFERON_ALPHA_RESPONSE"),
  G2M        =list(type="C", nm="G2M-index", hallmark="G2M_CHECKPOINT",
                   lab="G[2]M-index"))   # G2M drawn with a subscript 2: "[2]" = the subscript (element_text_sub below)
## G2M drawn with a subscript 2, the subscript at SUB_PT = 5 pt, the minimum text size (Arial, the face of its
## label: bold in the title). plotmath draws a subscript at 0.7 x the label size (4.9 pt for 7 pt) and cannot change that
## scale, so the label is drawn in pieces by a text element of its own: "G" and "M-index" at the label size, and the "2" at
## 5 pt right after "G", its baseline lowered by plotmath's subscript drop (~1.3 pt for 7 pt). The label's box is that of the
## plain string "G2M-index", so the layout stays exactly as with the plain string; inside it the label is justified on its
## ink box, as plotmath does, so "G" and "M" sit where plotmath would put them (the top of "G" stays inside the image).
SUB_PT <- 5
element_text_sub <- S7::new_class("element_text_sub", parent = ggplot2::element_text,
                                  constructor = function(...) S7::new_object(ggplot2::element_text(...)))
registerS3method("element_grob", "element_text_sub", function(element, label = "", ...){
  m <- regmatches(label, regexec("^(.*)\\[(.*)\\](.*)$", label))[[1]]   # "G[2]M-index" -> "G", "2", "M-index"
  g <- getS3method("element_grob", "ggplot2::element_text")(element, label = if (length(m) == 4) paste(m[-1], collapse = "") else label, ...)
  if (length(m) != 4 || !inherits(g, "titleGrob")) return(g)
  tg <- g$children[[1]]   # the text of the plain string: its position, justification, angle and font
  grid::setChildren(g, grid::gList(grid::gTree(txt = m[-1], x = tg$x, y = tg$y, hjust = tg$hjust, vjust = tg$vjust,
                                               rot = tg$rot, gp = tg$gp, cl = "subtext")))
})
registerS3method("makeContent", "subtext", function(x){   # measured at drawing time, on the device of ggsave
  gp <- function(size){ g <- x$gp; g$fontsize <- size; g }
  sz <- c(x$gp$fontsize, SUB_PT, x$gp$fontsize)
  m  <- function(f, s, z){ u <- f(grid::textGrob(s, gp = gp(z)))                    # a metric in pt
    if (identical(f, grid::grobWidth)) grid::convertWidth(u, "pt", valueOnly = TRUE) else grid::convertHeight(u, "pt", valueOnly = TRUE) }
  w  <- mapply(m, x$txt, sz, MoreArgs = list(f = grid::grobWidth))                    # advance widths of "G", "2", "M-index"
  drop <- m(grid::grobDescent, expression(X["2"]), sz[1])                             # plotmath's subscript drop
  dy <- c(0, -drop, 0)                                                                # baselines of the pieces
  asc <- max(mapply(m, x$txt, sz, MoreArgs = list(f = grid::grobAscent)) + dy)        # ink box above / below the baseline
  dsc <- max(mapply(m, x$txt, sz, MoreArgs = list(f = grid::grobDescent)) - dy)
  u <- -x$hjust * sum(w) + c(0, cumsum(w)[1:2])                                       # along the text (pt)
  v <- dsc - x$vjust * (asc + dsc) + dy                                               # across the text (pt)
  r <- x$rot / 180 * pi
  grid::setChildren(x, do.call(grid::gList, lapply(1:3, function(i)
    grid::textGrob(x$txt[i], x = x$x + grid::unit(u[i]*cos(r) - v[i]*sin(r), "pt"), y = x$y + grid::unit(u[i]*sin(r) + v[i]*cos(r), "pt"),
                   hjust = 0, vjust = 0, rot = x$rot, gp = gp(sz[i])))))
}, envir = asNamespace("grid"))

## ---- stats: per axis x cluster ----
## Primary (shown on figure) = d0 = one-sample Cohen's d vs 0 = mean/SD: how far the cluster's
## score sits from 0 (the neutral point of a signed RFind score), in SD units. There is no natural
## "control" group here (unlike the Stress vs No Stress of the De Donno data), so the deviation from zero is the
## per-cluster readout. Also recorded (table only): cluster-vs-rest d and AUC (separation context).
cohend <- function(a,b){ n1<-length(a); n2<-length(b); (mean(a)-mean(b))/sqrt(((n1-1)*var(a)+(n2-1)*var(b))/(n1+n2-2)) }
auc_vs <- function(a,b){ r<-rank(c(a,b)); n1<-length(a); n2<-length(b); (sum(r[seq_len(n1)])-n1*(n1+1)/2)/(n1*n2) }
STAT <- rbindlist(lapply(names(AX), function(ax) rbindlist(lapply(ord, function(cl){
  a <- d[cluster==cl][[ax]]; b <- d[cluster!=cl][[ax]]
  na <- d[cluster==cl, .N, by=animal][N >= MIN_CELLS_ANIMAL, .N]
  data.table(axis=ax, type=AX[[ax]]$type, cluster=cl, is_target=has_sig(cl),
             n=length(a), n_animals=na, median=median(a), mean=mean(a), sd=sd(a),
             d0=mean(a)/sd(a), d_vs_rest=cohend(a,b), auc_vs_rest=auc_vs(a,b))
}))))
STAT[, cluster := factor(cluster, levels=ord)]
fwrite(STAT[order(axis,cluster)][, .(axis,type,cluster,is_target,n,n_animals,
          median=round(median,3), mean=round(mean,3), sd=round(sd,3),
          d0=round(d0,3), d_vs_rest=round(d_vs_rest,3), auc_vs_rest=round(auc_vs_rest,3))],
       file.path(OUT, "Fig2_stats_table.csv"))

## the note on the numbers is written once only, under the bottom-right violin plot
## (G2M, lower half of d); cap = 'note' | 'blank' (an empty line of the same height, under DAM, so the rows of c and d stay level)
## | 'none' (upper halves)
mk <- function(ax, cap = "none"){
  info <- AX[[ax]]; dd <- d[, .(cluster=as.character(cluster), val=get(ax))]
  v <- dd$val; qlo <- as.numeric(quantile(v,0.005)); qhi <- as.numeric(quantile(v,0.995))
  if (qhi-qlo < 1e-6){ qlo<-min(v); qhi<-max(v) }; span <- qhi-qlo
  A <- STAT[axis==ax]; A[, y := qhi + 0.06*span]; A[, lab := sprintf("%.1f", d0)]   # d0 = mean/SD; the row is labelled once at its left end
  ## GLOBAL colour rule: colour a cluster with its signature colour iff it HAS a
  ## signature colour AND passes the d0 gate; otherwise grey. (No-signature clusters always grey.)
  ##   ranked lists (type "B") = signed contrast -> gate = |d0| >= 1 (either extreme is interpretable);
  ##   gene sets (type "C") = one-sided membership -> gate = d0 >= 1 (positive enrichment only, so the
  ##              clusters with strongly negative values, e.g. on Interferon, are not coloured).
  ##   gene sets: if NO signature cluster reaches d0 >= 1, the single top-enriched (max d0) signature cluster is
  ##   coloured (Interferon: cluster 9, green). G2M has clusters with d0 >= 1.
  if (info$type == "C") {
    cand   <- A[has_sig(cluster)]
    colcl  <- as.character(cand[d0 >= 1]$cluster)
    if (!length(colcl) && nrow(cand)) colcl <- as.character(cand[which.max(d0)]$cluster)
  } else {
    colcl <- as.character(A[abs(d0) >= 1 & has_sig(cluster)]$cluster)
  }
  ## numbers in black bold = clusters that pass the d0 gate, with or without a signature colour (e.g. clusters 1, 5, 6
  ## on Maturity have |d0| >= 1 but no signature colour; they stay grey, only their numbers are black bold)
  boldcl <- if (info$type == "C") union(colcl, as.character(A[d0 >= 1]$cluster)) else as.character(A[abs(d0) >= 1]$cluster)
  ## the key names only the clusters coloured in this panel, e.g. "2b (proliferative)", "2a/2b (proliferative)"
  keylab <- character(0); keycol <- character(0)
  for (g in SIG) { m <- intersect(g$cl, colcl); if (!length(m)) next
    desc <- sub("^[^(]*", "", g$lab)                            # the part in brackets: "(proliferative)", "(metabolically active)", ...
    lab  <- paste0(paste(m, collapse = "/"), " ", desc)          # "3 (metabolically active)", "2a/2b (proliferative)": the x axis says Cluster (so the longer names of Fig. 2b fit)
    keylab[m] <- lab; keycol[lab] <- g$col }
  dd[, glab := ifelse(cluster %in% colcl, unname(keylab[cluster]), "other")]
  glevs <- c(unique(unname(keylab)), if ("other" %in% dd$glab) "other")
  gcols <- c(keycol, other=COL_OTHER)[glevs]
  dd[, cluster := factor(cluster, levels=ord)]; dd[, glab := factor(glab, levels=glevs)]
  ## subtitle = the reference only: the source paper (ranked lists) or the MSigDB set (gene sets)
  parts <- character(0)
  if (!is.null(info$src))      parts <- c(parts, paste0("Reference: ", info$src))
  ## gene sets: the same "Reference:" prefix as the ranked lists
  if (!is.null(info$hallmark)) parts <- c(parts, paste0("Reference: MSigDB Hallmark ", info$hallmark))
  if (!is.null(info$dir))      parts <- c(parts, info$dir)
  sub <- if (length(parts)) paste(parts, collapse="   —   ") else NULL
  ## G2M drawn with a subscript 2: title and y axis title drawn by element_text_sub (subscript 2 at 5 pt), with
  ## the font, face, size, margins and box of the theme's plain labels
  ylab <- if (is.null(info$lab)) info$nm else info$lab
  ttl  <- ylab
  fix <- if (!is.null(info$lab)) theme(plot.title=element_text_sub(), axis.title.y=element_text_sub())
  pa <- d[, .(val=mean(get(ax)), n=.N), by=.(cluster, animal)][n >= MIN_CELLS_ANIMAL]
  ggplot(dd, aes(cluster, val, fill=glab)) +
    ## line widths and point sizes for the pasted size
    geom_hline(yintercept=0, linetype="dashed", colour="grey70", linewidth=0.15) +
    geom_violin(scale="width", width=0.82, linewidth=0.1, colour="grey35") +
    geom_boxplot(aes(group=cluster), width=0.13, outlier.shape=NA, linewidth=0.1, colour="grey20", fill="white") +
    geom_point(data=pa, aes(cluster, val), inherit.aes=FALSE, shape=21, size=0.6, stroke=0.15,
               fill="white", colour="grey10", position=position_jitter(width=0.16, height=0, seed=1)) +
    ## black+bold numbers on clusters that pass the d0 gate (boldcl), grey otherwise
    geom_text(data=A[ as.character(cluster) %in% boldcl], aes(cluster, y, label=lab), inherit.aes=FALSE,
              size=LAB_SIZE, fontface="bold", colour="grey15", family=FAM) +
    geom_text(data=A[!as.character(cluster) %in% boldcl], aes(cluster, y, label=lab), inherit.aes=FALSE,
              size=LAB_SIZE, colour="grey60", family=FAM) +
    ## no row label for the numbers (they are named in the caption of d and in the figure legend); the same padding
    ## on both sides
    scale_x_discrete(expand=expansion(add=c(0.6, 0.6))) +
    scale_fill_manual(values=gcols, breaks=setdiff(glevs,"other"), name=NULL, drop=FALSE) +
    coord_cartesian(ylim=c(min(qlo-0.03*span,0), qhi+0.18*span)) +
    ## the caption (bottom right, under G2M only) names the numbers: effect size d0 = mean/SD
    labs(title=ttl, subtitle=sub,
         caption=switch(cap, note="Numbers above violins: effect size d0 (mean/SD)", blank=" ", none=NULL),
         x="Cluster", y=ylab) +          # what violins/points/numbers mean -> figure legend (manuscript/Legends.md), not in the figure
    ## text sizes in pt at the pasted size
    theme_classic(7, base_family=FAM) +
    theme(plot.title=element_text(face="bold", size=7, margin=margin(b=0.6)),
          plot.subtitle=element_text(size=6, colour="grey35", margin=margin(b=0.3)),
          plot.caption=element_text(size=6, colour="grey35", hjust=1, margin=margin(t=0.3)),
          axis.title=element_text(size=7),
          axis.title.x=element_text(margin=margin(t=0.8)), axis.title.y=element_text(margin=margin(r=1)),
          axis.text=element_text(size=6, colour="grey20"),
          axis.line=element_line(colour="grey55", linewidth=LW05),
          axis.ticks=element_line(colour="grey55", linewidth=LW05), axis.ticks.length=unit(0.8, "mm"),
          legend.position="top", legend.justification="right",
          plot.title.position="plot", plot.caption.position="plot",
          legend.text=element_text(size=7), legend.key.size=unit(2.4,"mm"),
          legend.key.spacing.x=unit(1.2, "mm"), legend.margin=margin(0, 0, 0, 0), legend.box.spacing=unit(0.6, "mm"),
          plot.margin=margin(3, 1, 0.5, 0)) +
    fix   # NULL except for a label with a subscript (G2M)
}

sheet <- function(axes, header, fname, h, caps){   # header and h are not used (header names the panel at each call below)
  P <- wrap_plots(Map(mk, axes, caps), ncol=1) +
    plot_annotation(theme=theme(plot.margin=margin(0, 0, 0, 0)))   # no outer patchwork margin
  ## drawn directly at the pasted size, given in pixels (mm / 25.4 * 600, rounded) so the PNG is exactly that size
  ## (units="mm" truncates)
  ggsave(file.path(OUT_PPTX, paste0(fname,".png")), P, width=round(W_MM/25.4*600), height=round(H_MM/25.4*600),
         units="px", dpi=600, bg="white", device=ragg::agg_png)
  cat(sprintf("wrote %s\n", fname))
}
## Panel 2c = the ranked lists (Maturity [Matcovitch-Natan] + DAM [Keren-Shaul]); panel 2d = the gene sets
## (Interferon / G2M).
sheet(c("Maturity","DAM"),   "Fig2c — Maturity-index / DAM-index",   "Fig2c_violin_Maturity_DAM",   9.4, c("none", "blank"))
sheet(c("Interferon","G2M"), "Fig2d — Interferon-index / G2M-index", "Fig2d_violin_Interferon_G2M", 9.4, c("none", "note"))
cat("\nwrote Fig2_stats_table.csv\n")
print(STAT[is_target==TRUE][order(axis,-d_vs_rest), .(axis, cluster, d=round(d_vs_rest,2), auc=round(auc_vs_rest,3))])
