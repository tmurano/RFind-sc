## Supplementary Fig. S5c, d — violins of the four indices per cell type, unstimulated vs IFN-beta, drawn with the code of
##   figs/Fig3_DeDonno/analysis/fig3_violins.R (same marks, sizes and colours of the conditions).
##   c: IFNb-index + PGE2-index (type A: IFN-beta-stimulated PBMC, GSE191317, and PGE2-stimulated CD14+ monocytes, GSE272019)
##   d: Interferon-index + Chemokine-index (type B: MSigDB Hallmark and KEGG gene sets)
## Changes from Fig. 3: the Kang data; the eight author cell types with the abbreviations, order and colours of Fig. 3a of the
##   original study (figS5_celltypes.R; Megakaryocytes included, as in the original paper); means per donor (>= 10 cells) instead
##   of per mouse; dotted lines = myeloid leukocytes | Megakaryocytes | lymphocytes; numbers with |d| >= 1.5 in black bold (the
##   PGE2-index falls), all others grey.
## Stats: per-donor means (>= MIN_CELLS cells), Cohen's d (pooled s.d.); the figure shows d only. The nominal two-sided Wilcoxon
##   rank-sum p is written to figS5_violin_stats.csv for reference only.
## Input : data/Kang2018/output/scores_kang.csv
## Output: figs/Supp_S5_Kang/output/{figS5_violin_stats.csv, FigS5_celltype_ncells.csv, FigS5c_violin_IFNb_PGE2.png,
##         FigS5d_violin_Interferon_Chemokine.png}
## Run from the project root: R -f figs/Supp_S5_Kang/analysis/figS5_violins.R
pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("figs/Supp_S5_Kang/analysis/figS5_celltypes.R")
OUT <- FIG
MIN_CELLS <- 10L; MIN_DONORS <- 3L   # >= 10 cells per donor (as Fig. 3), >= 3 donors per condition

## rows in the order of the cell IDs (the order in which the per-donor points are placed)
d <- fread("data/Kang2018/output/scores_kang.csv")[, .(cell_id, cell_type, condition, donor_id = donor, IFNb, PGE2, Interferon, Chemokine)][order(cell_id)]
d[, condition := factor(COND_LAB[condition], levels = COND_LAB)]

## ---- which cell types can be shown (>= MIN_DONORS donors per condition with >= MIN_CELLS cells) ----
pm_n <- d[, .N, by = .(cell_type, condition, donor_id)][N >= MIN_CELLS][, .(donors = .N), by = .(cell_type, condition)]
NC <- d[, .(n_Unstimulated = sum(condition == "Control"), n_IFNb = sum(condition == "IFN-β"), n_total = .N), by = cell_type]
NC <- merge(NC, dcast(pm_n, cell_type ~ condition, value.var = "donors", fill = 0), by = "cell_type", all.x = TRUE)
## donors_* = number of donors with >= MIN_CELLS cells of the type
setnames(NC, c("Control", "IFN-β"), c("donors_Unstimulated", "donors_IFNb"), skip_absent = TRUE)
for (k in c("donors_Unstimulated", "donors_IFNb")) { if (!k %in% names(NC)) NC[, (k) := 0L]; set(NC, which(is.na(NC[[k]])), k, 0L) }
SHOW <- CT_LV   # all eight author cell types (as Kang et al.)
NC[, shown := donors_Unstimulated >= MIN_DONORS & donors_IFNb >= MIN_DONORS & cell_type %in% SHOW]
stopifnot(sum(NC$shown) == length(SHOW))
NC <- NC[order(match(cell_type, CT_LV))]
fwrite(NC, file.path(OUT, "FigS5_celltype_ncells.csv"))
ord <- CT_LV[CT_LV %in% NC[shown == TRUE]$cell_type]
cat("shown in violins:", paste(CT_TXT[ord], collapse = ", "), "\n")

d <- d[cell_type %in% ord]; d[, cell_type := factor(cell_type, levels = ord)]
VLINES <- sort(unique(c(sum(ord %in% CT_MYELOID) + 0.5, if ("Megakaryocytes" %in% ord) match("Megakaryocytes", ord) + c(-0.5, 0.5))))
VLINES <- VLINES[VLINES > 0.5 & VLINES < length(ord) + 0.5]   # dotted lines = myeloid leukocytes | Megakaryocytes | lymphocytes
UNSTIM_GREY <- "#D9D9D9"
d[, fillcol := ifelse(condition == "Control", UNSTIM_GREY, CT_PAL[as.character(cell_type)])]

AX <- c("IFNb", "PGE2", "Interferon", "Chemokine")
META <- list(
  IFNb = list(nm="IFNb-index", sub="Reference: GSE191317 (IFN-β vs mock, PBMC, 18 h)"),
  PGE2 = list(nm="PGE2-index", sub="Reference: GSE272019 (PGE2 vs unstimulated, CD14+ monocytes, 3 h)"),
  Interferon = list(nm="Interferon-index", sub="Reference: MSigDB Hallmark INTERFERON_ALPHA_RESPONSE"),
  Chemokine  = list(nm="Chemokine-index", sub="Reference: KEGG CHEMOKINE_SIGNALING_PATHWAY"))
cohend <- function(a,b){ n1<-length(a); n2<-length(b); (mean(b)-mean(a))/sqrt(((n1-1)*var(a)+(n2-1)*var(b))/(n1+n2-2)) }
STAT <- rbindlist(lapply(AX, function(ax) rbindlist(lapply(ord, function(ct){
  pm <- d[cell_type==ct, .(m=mean(get(ax)), nn=.N), by=.(donor_id, condition)][nn>=MIN_CELLS]
  ns <- pm[condition=="Control"]$m; st <- pm[condition=="IFN-β"]$m
  data.table(axis=ax, cell_type=ct, n_donors_Unstimulated=length(ns), n_donors_IFNb=length(st), d=cohend(ns,st),
             p=tryCatch(suppressWarnings(wilcox.test(st,ns)$p.value), error=function(e) NA)) }))))
D_BLACK <- 1.5   # |d| >= 1.5 -> black bold (decreases too: the PGE2-index falls); otherwise grey
fwrite(STAT, file.path(OUT, "figS5_violin_stats.csv"))

## ---- drawn at the size pasted into the figure (W_MM x H_MM, 600 dpi) ----
## Arial, titles, axis titles and the d numbers 7 pt, tick labels, subtitles and caption 6 pt, axis lines and ticks 0.5 pt;
## the sizes of the per-donor points are given for a 254-mm-wide drawing and scaled by SHRINK.
W_MM <- 87.2; H_MM <- 73.0                              # pasted size (mm)
FAM <- "Arial"
SHRINK <- W_MM / 254                                    # scale of the point sizes below to the pasted width
PT_LW <- function(pt) pt / (.pt * 72.27 / 96)           # ggplot linewidth that gives `pt` points
LAB_SIZE <- 7 / .pt     # the d numbers: 7 pt (in-panel text, as the d0 numbers of Fig. 2c, d)
## the note is written once only, under the bottom-right panel (Chemokine, lower half of d);
## cap = 'note' | 'blank' (an empty line of the same height under the PGE2-index, so the rows of c and d stay level) | 'none' (upper halves)
mk <- function(ax, cap = "note"){
  m <- META[[ax]]; v <- d[[ax]]
  ## the y range covers all cells, so that no violin is cut (Fig. 3 uses the 0.5th-99.5th percentiles, which cut the long
  ## tails of the monocytes here); the numbers sit above the highest cell
  qlo <- min(v); qhi <- max(v)
  span <- qhi - qlo; ann_y <- qhi + 0.07*span   # numbers on white labels, so that violin tails do not cross them
  A <- STAT[axis==ax]; A[, cell_type := factor(cell_type, levels=ord)]; A[, y := ann_y]
  ## donor-level points = the per-donor means the stats are computed on (>= MIN_CELLS cells),
  ## placed over the dodged violin of their condition (unstimulated left, IFN-beta right)
  pm <- d[, .(m = mean(get(ax)), nn = .N), by = .(cell_type, condition, donor_id)][nn >= MIN_CELLS]
  set.seed(1)
  pm[, xpos := as.numeric(cell_type) + ifelse(condition == "Control", -0.195, 0.195) + runif(.N, -0.05, 0.05)]
  A[, lab := sprintf("%.1f", d)]
  A[lab == "-0.0", lab := "0.0"]   # a value that rounds to zero is written without a minus sign
  big <- A[!is.na(d) & abs(d) >= D_BLACK]; sml <- A[is.na(d) | abs(d) < D_BLACK]
  ggplot(d, aes(cell_type, get(ax), fill=fillcol, group=interaction(cell_type, condition))) +
    geom_hline(yintercept=0, linewidth=PT_LW(0.35), colour="grey70", linetype="dashed") +
    geom_vline(xintercept=VLINES, linewidth=PT_LW(0.45), linetype="dotted", colour="grey75") +
    geom_violin(scale="width", width=0.78, linewidth=PT_LW(0.3), colour="grey35", alpha=0.95, position=position_dodge(0.78)) +
    geom_boxplot(aes(group=interaction(cell_type, condition)), width=0.14, outlier.shape=NA, linewidth=PT_LW(0.3),
                 colour="grey20", fill="white", position=position_dodge(0.78), show.legend=FALSE) +
    geom_point(data=pm, aes(x=xpos, y=m), inherit.aes=FALSE, shape=21, size=1.9 * SHRINK * 1.2, stroke=0.35 * SHRINK,
               fill="white", colour="grey10") +
    geom_label(data=sml, aes(cell_type, y, label=lab), inherit.aes=FALSE, size=LAB_SIZE, family=FAM, colour="grey62", fill="white", border.colour=NA, label.padding=unit(0.2, "mm"), label.r=unit(0, "mm")) +
    geom_label(data=big, aes(cell_type, y, label=lab), inherit.aes=FALSE, size=LAB_SIZE, family=FAM, fontface="bold", colour="grey15", fill="white", border.colour=NA, label.padding=unit(0.2, "mm"), label.r=unit(0, "mm")) +
    ## no row label for the numbers (as in Fig. 2); the same padding on both sides
    scale_fill_identity() + scale_x_discrete(labels=parse(text=CT_LAB[ord]), expand=expansion(add=c(0.6, 0.6))) +
    coord_cartesian(ylim=c(min(qlo - 0.02*span, 0), qhi + 0.13*span)) +
    ## the meaning of the numbers is written once at the bottom right, as in Fig. 2
    labs(title=m$nm, subtitle=m$sub, x=NULL, y=m$nm,   # no x name (the cell-type names say it)
         caption=switch(cap, note="Numbers above violins: effect size d (IFN-β vs unstimulated)", blank=" ", none=NULL)) +   # other mark meanings -> manuscript/Legends.md
    theme_classic(base_size = 7, base_family = FAM) +
    theme(plot.title=element_text(face="bold", size=7, margin=margin(b=1)),
          plot.subtitle=element_text(size=6, colour="grey35", margin=margin(b=2)),
          plot.caption=element_text(size=6, colour="grey35", hjust=1, margin=margin(t=2)),
          axis.title=element_text(size=7, colour="black"),
          axis.title.x=element_text(margin=margin(t=1)), axis.title.y=element_text(margin=margin(r=1.5)),
          axis.text=element_text(size=6, colour="grey20"),
          ## the abbreviations under eight narrow columns are slanted by 30 degrees, as in Fig. 3
          axis.text.x=element_text(angle=30, hjust=1, vjust=1),
          axis.line=element_line(colour="grey55", linewidth=PT_LW(0.5)), axis.ticks=element_line(colour="grey55", linewidth=PT_LW(0.5)),
          axis.ticks.length=unit(0.7, "mm"),
          legend.position="none", plot.margin=margin(1, 1, 3, 0))
}
px <- function(mm) round(mm / 25.4 * 600)                            # ragg truncates mm sizes, so give whole pixels
sheet2 <- function(axes, fname, cap2){
  ## no outer margin (the sheet is pasted at its drawn size); 3 pt under the upper panel only
  P <- (mk(axes[1], cap = "none") / (mk(axes[2], cap = cap2) + theme(plot.margin = margin(1, 1, 0.5, 0)))) +
    plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
  ragg::agg_png(file.path(FIG, paste0(fname, ".png")), width = px(W_MM), height = px(H_MM), units = "px", res = 600,
                background = "white")
  print(P); invisible(dev.off())
  cat(sprintf("wrote %s\n", fname))
}
sheet2(c("IFNb", "PGE2"), "FigS5c_violin_IFNb_PGE2", "blank")                     # c: type A references
sheet2(c("Interferon", "Chemokine"), "FigS5d_violin_Interferon_Chemokine", "note")  # d: type B references
cat(sprintf("\n=== black numbers (|d| >= %.1f) ===\n", D_BLACK))
print(STAT[abs(d) >= D_BLACK][order(axis, -d), .(axis, cell_type, d = round(d, 2), p_nominal = signif(p, 2))])
