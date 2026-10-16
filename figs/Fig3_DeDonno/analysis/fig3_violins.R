## Fig. 3c, d (De Donno) — violins of the 4 indices per cell type, wild-type mice, No Stress vs Stress.
##   c: Glucocorticoid-index + Seizure-index (references = external DEG lists; data keys CORT and KA6h)
##   d: NF-kB-index + UPR-index (references = MSigDB Hallmark gene sets)
## Cell types = the 21 author types of the wild-type cells (fig3_celltypes.R; each keeps its own Delta baseline). Shown = the cell types of the
##   seven largest groups in the coarse annotation of the original study (Methods: Astro, DG, Oligo, CA1, Endothelial, Microglia, CA3)
##   plus GABAergic neurons
##   (the inhibitory neurons, also analysed in Fig. 4). Per-mouse rule: >= 3 mice per condition with
##   >= 10 cells (GABA has 6-26 cells per mouse; the seven major types have >= 20 in all 6 + 6 mice, so their numbers are the
##   same as with the >= 20 cells of Fig. 2). The other types are listed in the log and in Fig3_celltype_ncells.csv.
## Fill: No Stress = grey; Stress = the cell type's colour. Dotted lines = neurons | glia | vascular (endothelial) cells.
## Stats: per-mouse means (>= MIN_CELLS cells), Cohen's d (pooled s.d.); the figure shows d only (no stars, no FDR).
##   Numbers with d >= 1.5 (an increase under stress) are black bold, all others grey — a decrease such as -3.2 is grey
##   too. The nominal two-sided Wilcoxon rank-sum p is written to fig3_violin_stats.csv for reference only.
## Style follows figs/Fig2_Hammond/analysis/fig2_violins.R (title = index name only, subtitle = reference only,
##   same font sizes); mark meanings go in manuscript/Legends.md.
##
## Input : data/DeDonno2026/output/scores_dedonno.csv
## Output: figs/Fig3_DeDonno/output/fig3_violin_stats.csv + Fig3_celltype_ncells.csv
##         + {Fig3c_violin_Glucocorticoid_Seizure, Fig3d_violin_NFkB_UPR}.png (drawn at the pasted size)
## Run   : Rscript "figs/Fig3_DeDonno/analysis/fig3_violins.R"   (from project root)

pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("figs/Fig3_DeDonno/analysis/fig3_celltypes.R")
OUT <- "figs/Fig3_DeDonno/output"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
MIN_CELLS <- 10L; MIN_MICE <- 3L   # >= 10 cells per mouse, so that the GABAergic neurons (6-26 cells per mouse) can be shown

d <- fread("data/DeDonno2026/output/scores_dedonno.csv")[line == "Ctrl"]
d[, condition := factor(condition, levels = c("No Stress","Stress"))]

## ---- which cell types can be shown (>= MIN_MICE mice per condition with >= MIN_CELLS cells) ----
pm_n <- d[, .N, by = .(cell_type, condition, donor_id)][N >= MIN_CELLS][, .(mice = .N), by = .(cell_type, condition)]
NC <- d[, .(n_NoStress = sum(condition == "No Stress"), n_Stress = sum(condition == "Stress"), n_total = .N), by = cell_type]
NC <- merge(NC, dcast(pm_n, cell_type ~ condition, value.var = "mice", fill = 0), by = "cell_type", all.x = TRUE)
## mice20_* (columns of Fig3_celltype_ncells.csv) = number of mice with >= MIN_CELLS cells of the type
setnames(NC, c("No Stress", "Stress"), c("mice20_NoStress", "mice20_Stress"), skip_absent = TRUE)
for (k in c("mice20_NoStress", "mice20_Stress")) { if (!k %in% names(NC)) NC[, (k) := 0L]; set(NC, which(is.na(NC[[k]])), k, 0L) }
## GABA and the cell types of the seven largest groups in the coarse annotation of De Donno et al. (Fig. 1e, inset; DG = DG Glut, Oligo = Oligodendrocytes)
SHOW <- c("DG Glut", "CA1 Glut", "CA3 Glut", "GABA", "Astro", "Microglia", "Oligodendrocytes", "Endothelial")
NC[, shown := mice20_NoStress >= MIN_MICE & mice20_Stress >= MIN_MICE & cell_type %in% SHOW]
stopifnot(sum(NC$shown) == length(SHOW))
NC <- NC[order(match(cell_type, CT_LV))]
fwrite(NC, file.path(OUT, "Fig3_celltype_ncells.csv"))
ord <- CT_LV[CT_LV %in% NC[shown == TRUE]$cell_type]
cat("shown in violins:", paste(CT_LAB[ord], collapse = ", "), "\n")
cat("not shown      :", paste(CT_LAB[setdiff(CT_LV[CT_LV %in% NC$cell_type], ord)], collapse = ", "), "\n")

d <- d[cell_type %in% ord]; d[, cell_type := factor(cell_type, levels = ord)]
nNeuron <- sum(ord %in% CT_NEURON)
## a second dotted line before the vascular cells: endothelial cells are not glia
CT_VASCULAR <- c("Endothelial", "Vascular", "Pericytes")
VLINES <- c(nNeuron, if (any(ord %in% CT_VASCULAR)) min(which(ord %in% CT_VASCULAR)) - 1L) + 0.5
NOSTRESS_GREY <- "#D9D9D9"
d[, fillcol := ifelse(condition == "No Stress", NOSTRESS_GREY, CT_PAL[as.character(cell_type)])]

AX <- c("CORT","KA6h","NFkB","UPR")
META <- list(
  ## subtitle = reference only (gene-selection details go to Methods)
  ## display names say what each index measures: CORT -> Glucocorticoid-index, KA6h -> Seizure-index (data keys CORT / KA6h)
  CORT = list(nm="Glucocorticoid-index", sub="Reference: GSE184924 (corticosterone vs vehicle, 3 h)"),
  KA6h = list(nm="Seizure-index", sub="Reference: GSE125068 (kainate seizure vs saline, 6 h)"),
  ## the Hallmark-based axes start with "Reference:" as the DEG-based ones do (as in Fig. 2)
  NFkB = list(nm="NF-κB-index", sub="Reference: MSigDB Hallmark TNFA_SIGNALING_VIA_NFKB"),
  UPR  = list(nm="UPR-index",   sub="Reference: MSigDB Hallmark UNFOLDED_PROTEIN_RESPONSE"))
cohend <- function(a,b){ n1<-length(a); n2<-length(b); (mean(b)-mean(a))/sqrt(((n1-1)*var(a)+(n2-1)*var(b))/(n1+n2-2)) }
STAT <- rbindlist(lapply(AX, function(ax) rbindlist(lapply(ord, function(ct){
  pm <- d[cell_type==ct, .(m=mean(get(ax)), nn=.N), by=.(donor_id, condition)][nn>=MIN_CELLS]
  ns <- pm[condition=="No Stress"]$m; st <- pm[condition=="Stress"]$m
  data.table(axis=ax, cell_type=ct, n_mice_NoStress=length(ns), n_mice_Stress=length(st), d=cohend(ns,st),
             p=tryCatch(suppressWarnings(wilcox.test(st,ns)$p.value), error=function(e) NA)) }))))
D_BLACK <- 1.5   # d >= 1.5 -> black bold; otherwise grey (decreases included)
fwrite(STAT, file.path(OUT, "fig3_violin_stats.csv"))

## ---- drawn at the size pasted into Fig. 3 (W_MM x H_MM, 600 dpi) ----
## Arial, titles, axis titles and the d numbers 7 pt, tick labels, subtitles and caption 6 pt, axis lines and ticks 0.5 pt;
## the sizes of the per-mouse points are given for a 254-mm-wide drawing and scaled by SHRINK.
W_MM <- 87.2; H_MM <- 73.0                              # pasted size (mm)
FAM <- "Arial"
SHRINK <- W_MM / 254                                    # scale of the point sizes below to the pasted width
PT_LW <- function(pt) pt / (.pt * 72.27 / 96)           # ggplot linewidth that gives `pt` points
LAB_SIZE <- 7 / .pt     # the d numbers: 7 pt (in-panel text, as the d0 numbers of Fig. 2c, d)
## the note is written once only, under the bottom-right panel (UPR, lower half of d);
## cap = 'note' | 'blank' (an empty line of the same height under Seizure, so the rows of c and d stay level) | 'none' (upper halves)
mk <- function(ax, cap = "note"){
  m <- META[[ax]]; v <- d[[ax]]
  qlo <- as.numeric(quantile(v, 0.005)); qhi <- as.numeric(quantile(v, 0.995))
  if (qhi - qlo < 1e-6) { qlo <- min(v); qhi <- max(v) }
  span <- qhi - qlo; ann_y <- qhi + 0.05*span
  A <- STAT[axis==ax]; A[, cell_type := factor(cell_type, levels=ord)]; A[, y := ann_y]
  ## animal-level points = the per-mouse means the stats are computed on (>= MIN_CELLS cells),
  ## placed over the dodged violin of their condition (No Stress left, Stress right)
  pm <- d[, .(m = mean(get(ax)), nn = .N), by = .(cell_type, condition, donor_id)][nn >= MIN_CELLS]
  set.seed(1)
  pm[, xpos := as.numeric(cell_type) + ifelse(condition == "No Stress", -0.195, 0.195) + runif(.N, -0.05, 0.05)]
  A[, lab := sprintf("%.1f", d)]
  A[lab == "-0.0", lab := "0.0"]   # a value that rounds to zero is written without a minus sign
  big <- A[!is.na(d) & d >= D_BLACK]; sml <- A[is.na(d) | d < D_BLACK]
  ggplot(d, aes(cell_type, get(ax), fill=fillcol, group=interaction(cell_type, condition))) +
    geom_hline(yintercept=0, linewidth=PT_LW(0.35), colour="grey70", linetype="dashed") +
    geom_vline(xintercept=VLINES, linewidth=PT_LW(0.45), linetype="dotted", colour="grey75") +
    geom_violin(scale="width", width=0.78, linewidth=PT_LW(0.3), colour="grey35", alpha=0.95, position=position_dodge(0.78)) +
    geom_boxplot(aes(group=interaction(cell_type, condition)), width=0.14, outlier.shape=NA, linewidth=PT_LW(0.3),
                 colour="grey20", fill="white", position=position_dodge(0.78), show.legend=FALSE) +
    geom_point(data=pm, aes(x=xpos, y=m), inherit.aes=FALSE, shape=21, size=1.9 * SHRINK * 1.2, stroke=0.35 * SHRINK,
               fill="white", colour="grey10") +
    geom_text(data=sml, aes(cell_type, y, label=lab), inherit.aes=FALSE, size=LAB_SIZE, family=FAM, lineheight=0.8, colour="grey62") +
    geom_text(data=big, aes(cell_type, y, label=lab), inherit.aes=FALSE, size=LAB_SIZE, family=FAM, fontface="bold", lineheight=0.8, colour="grey15") +
    ## no row label for the numbers (as in Fig. 2); the same padding on both sides
    scale_fill_identity() + scale_x_discrete(labels=CT_LAB, expand=expansion(add=c(0.6, 0.6))) +
    coord_cartesian(ylim=c(min(qlo - 0.03*span, 0), qhi + 0.20*span)) +
    ## the meaning of the numbers is written once at the bottom right, as in Fig. 2
    labs(title=m$nm, subtitle=m$sub, x=NULL, y=m$nm,   # no x name (the cell-type names say it)
         caption=switch(cap, note="Numbers above violins: effect size d (Stress vs No Stress)", blank=" ", none=NULL)) +   # other mark meanings -> manuscript/Legends.md
    theme_classic(base_size = 7, base_family = FAM) +
    theme(plot.title=element_text(face="bold", size=7, margin=margin(b=1)),
          plot.subtitle=element_text(size=6, colour="grey35", margin=margin(b=2)),
          plot.caption=element_text(size=6, colour="grey35", hjust=1, margin=margin(t=2)),
          axis.title=element_text(size=7, colour="black"),
          axis.title.x=element_text(margin=margin(t=1)), axis.title.y=element_text(margin=margin(r=1.5)),
          axis.text=element_text(size=6, colour="grey20"),
          ## the names of the original study are too long for horizontal labels under eight narrow columns: slanted by 30 degrees
          axis.text.x=element_text(angle=30, hjust=1, vjust=1),
          axis.line=element_line(colour="grey55", linewidth=PT_LW(0.5)), axis.ticks=element_line(colour="grey55", linewidth=PT_LW(0.5)),
          axis.ticks.length=unit(0.7, "mm"),
          legend.position="none", plot.margin=margin(1, 1, 3, 0))
}
OUT_PPTX <- OUT   # pasted-size panels go to output/
px <- function(mm) round(mm / 25.4 * 600)                            # ragg truncates mm sizes, so give whole pixels
sheet2 <- function(axes, fname, cap2){
  ## no outer margin (the sheet is pasted at its drawn size); 3 pt under the upper panel only
  P <- (mk(axes[1], cap = "none") / (mk(axes[2], cap = cap2) + theme(plot.margin = margin(1, 1, 0.5, 0)))) +
    plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
  ragg::agg_png(file.path(OUT_PPTX, paste0(fname, ".png")), width = px(W_MM), height = px(H_MM), units = "px", res = 600,
                background = "white")
  print(P); invisible(dev.off())
  cat(sprintf("wrote %s\n", fname))
}
## file names use the index names; the data keys stay CORT / KA6h
sheet2(c("CORT","KA6h"), "Fig3c_violin_Glucocorticoid_Seizure", "blank")
sheet2(c("NFkB","UPR"),  "Fig3d_violin_NFkB_UPR", "note")

cat(sprintf("\n=== black numbers (d >= %.1f) ===\n", D_BLACK))
print(STAT[d >= D_BLACK][order(axis, -d), .(axis, cell_type, d = round(d, 2), p_nominal = signif(p, 2))])
