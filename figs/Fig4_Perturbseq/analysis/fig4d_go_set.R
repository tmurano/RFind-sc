## Fig4d (Perturb-seq) — where the knock-outs of a gene set sit in the ranking against the stress response.
##   Fig3c showed the Glucocorticoid-index rising with stress in DG and CA1 but not in CA3. If that rise is glucocorticoid
##   signalling, knocking out the genes of that pathway should move a cell in the direction opposite to the stress response
##   in DG and CA1, and not in CA3.
##   The set is either a public GO term or one of the named sets in PRESET below; every member of the set present in the
##   knock-out library is drawn whatever its position. The paper uses the preset GR-complex: the five corticosteroid-receptor
##   and GR-complex genes in the library (Nr3c1, Nr3c2, Hsp90ab1, Ncoa1, Ncoa3; Methods), a set defined in this study rather
##   than a public one. Supplementary Table 2 (fig4_random_set_null.py) compares its placement with random sets of five
##   knock-outs drawn from the same library.
##   One panel per cell type: all knock-outs ranked by RFind score (x = rank counted from the opposed end, 1 = most opposed,
##   as in c and the heatmap of e; below 0 = opposed), with the set's knock-outs marked and named. No significance calls.
##   With GR-complex, the names of Nr3c1 and Ncoa1 in CA3 swap sides so that their leader lines do not cross (SWAP below).
##
## Input : analysis/perturbseq/results/axis_<target>_vs_<ref>.csv; GO membership from org.Mm.eg.db (GOALL = the term and its children)
## Output: figs/Fig4_Perturbseq/output/Fig4d_<set>.png + figs/Fig4_Perturbseq/output/fig4d_<set>_positions.csv
##   (<set> = the argument with the colon removed, e.g. GR-complex or GO0042921). The figure is drawn at the size the panel
##   has in the final figure (W_MM x H_MM below, 600 dpi, Arial; tick labels 6 pt, titles / strips / legend 7 pt, the names
##   of the knock-outs 6 pt; axis lines and ticks 0.5 pt).
## Run   : Rscript "figs/Fig4_Perturbseq/analysis/fig4d_go_set.R" GR-complex   (from project root; the panel of the paper.
##         A GO id such as GO:0042921 draws that term instead; with no argument, GO:0042921 is drawn)

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel)
                                 library(AnnotationDbi); library(org.Mm.eg.db); library(GO.db) })
pdf(NULL)   # avoid Rplots.pdf in the project root
args <- commandArgs(TRUE)
GO   <- if (length(args) > 0) args[1] else "GO:0042921"
OUT  <- "figs/Fig4_Perturbseq/output"
OUT_PPTX <- OUT   # the panel at its final size goes straight into output/
FONT <- "Arial"
W_MM <- 89.5; H_MM <- 50.0; DPI <- 600   # size in the final figure (d, e and f share one row)
ENL  <- W_MM / 97.8                       # point sizes scale with the panel width (relative to 97.8 mm); text stays in pt
LW_AXIS <- 0.5 / (.pt * 72.27 / 96)       # 0.5 pt as a ggplot linewidth (lwd 1 = 1/96 inch)
CT   <- c(DG = "g017_DG", CA1 = "CA1", CA3 = "g017_CA3", GABA = "GABA3")   # target (De Donno) -> reference library (Perturb-seq)
## panel titles: the names of the cell types in De Donno et al. (Fig. 3b)
CT_LAB <- c(DG = "DG Glut", CA1 = "CA1 Glut", CA3 = "CA3 Glut", GABA = "GABA")

## The set is either a public GO term (GO:xxxxxxx, membership from org.Mm.eg.db) or one of the named sets below.
PRESET <- list(
  ## the set of the paper (Fig. 4d; Methods): the two corticosteroid receptors (Nr3c1, glucocorticoid receptor; Nr3c2,
  ## mineralocorticoid receptor) and the proteins that act on the glucocorticoid receptor directly, its chaperone and its
  ## coactivators; every library member of this definition (5 genes)
  "GR-complex" = list(name = "Knock-outs of corticosteroid receptors and their complexes",   # key: what the black dots are
                      genes = c("Nr3c1", "Nr3c2", "Hsp90ab1", "Ncoa1", "Ncoa3")),
  ## Other presets, not used in the paper: GR-lists (Nr3c1 with the GR regulators Ncoa1 and Crebrf) and c-named (Ncoa1,
  ## Crebrf and Gal). The knock-outs in marked_in_c are drawn as diamonds.
  "GR-lists" = list(name = "GR and the GR regulators Ncoa1 and Crebrf", genes = c("Nr3c1", "Ncoa1", "Crebrf"),
                    marked_in_c = c("Ncoa1", "Crebrf")),
  "c-named" = list(name = "Knock-outs Ncoa1, Crebrf and Gal", genes = c("Ncoa1", "Crebrf", "Gal"),
                   marked_in_c = c("Ncoa1", "Crebrf", "Gal")))
if (grepl("^GO:", GO)) {
  GO_NAME <- Term(GOTERM[[GO]])
  GENES <- sort(unique(suppressMessages(
    AnnotationDbi::select(org.Mm.eg.db, keys = GO, keytype = "GOALL", columns = "SYMBOL"))$SYMBOL))
  cat(sprintf("%s (%s): %d mouse genes; org.Mm.eg.db %s\n", GO, GO_NAME, length(GENES), packageVersion("org.Mm.eg.db")))
} else {
  stopifnot(GO %in% names(PRESET))
  GO_NAME <- PRESET[[GO]]$name; GENES <- PRESET[[GO]]$genes
  cat(sprintf("%s (%s): %d genes\n", GO, GO_NAME, length(GENES)))
}

R <- rbindlist(lapply(names(CT), function(ct) {
  x <- fread(sprintf("analysis/perturbseq/results/axis_%s_vs_%s.csv", ct, CT[[ct]]))[kind == "KO"]
  setorder(x, -rfind); x[, rank := .I]; x[, n_ko := .N]; x[, ct := ct]   # rank = the pipeline rank, 1 = most aligned
  ## only the displayed rank is counted from the opposed end (the scores are unchanged), as in c and the heatmap
  x[, rank_opp := n_ko - rank + 1L]                                        # displayed rank, 1 = most opposed
  x[, .(ct, target, rfind, z_vs_KO, n_cells, rank, rank_opp, n_ko)]
}))
R[, ct := factor(ct, levels = names(CT))]
H <- R[target %in% GENES]
MARKED <- if (!grepl("^GO:", GO) && !is.null(PRESET[[GO]]$marked_in_c)) PRESET[[GO]]$marked_in_c else character(0)
H[, grp := if (length(MARKED)) ifelse(target %in% MARKED, "marked", "receptor") else "set"]
cat(sprintf("in the library: %s\n", paste(sort(unique(H$target)), collapse = ", ")))
cat("rank shown on the x axis (1 = most opposed):\n")
print(dcast(H, target ~ ct, value.var = "rank_opp"))

## names are drawn only while they stay readable; for a large set only the receptor and the coactivator are named
KEY_GENES <- c("Nr3c1", "Ncoa1")
LAB <- if (uniqueN(H$target) <= 15) H else H[target %in% KEY_GENES]
cat(sprintf("labelled in the figure: %s\n", paste(sort(unique(LAB$target)), collapse = ", ")))

## which side of its point each name goes to (-1 = left, 1 = right): a point in the left half of the ranking has its name on the
## right, and the reverse. With the set GR-complex, in CA3 the names of Ncoa1 (rank 820, left half) and Nr3c1 (1181, right
## half) would face each other with crossing lines, so the two swap sides (Ncoa1 to the left, Nr3c1 to the right)
SIDE <- LAB[, ifelse(rank_opp > n_ko / 2, -1, 1)]
SWAP <- GO == "GR-complex" & LAB$ct == "CA3" & LAB$target %in% c("Ncoa1", "Nr3c1")
SIDE[SWAP] <- -SIDE[SWAP]
if (any(SWAP)) cat(sprintf("names moved to the other side of their point: %s\n", paste(LAB$ct[SWAP], LAB$target[SWAP], collapse = ", ")))

p <- ggplot(R, aes(rank_opp, rfind)) +
  geom_hline(yintercept = 0, linewidth = 0.15, colour = "grey70", linetype = "dashed") +
  geom_point(aes(colour = "Other knock-outs"), size = 0.12 * ENL) +
  geom_rug(data = H, aes(x = rank_opp), inherit.aes = FALSE, sides = "b", length = unit(0.05, "npc"),
           colour = "grey5", linewidth = 0.2) +
  geom_point(data = H, aes(fill = grp, shape = grp), size = 1.5 * ENL, stroke = 0.25, colour = "white") +
  ## labels: when the marked knock-outs pile up at one end (e.g. CA1 with GR-complex, all five within 6% of the opposed
  ## end), push the names to the other side and stack them vertically, so the points stay visible
  geom_text_repel(data = LAB, aes(label = target), size = 6 / .pt, fontface = "italic", family = FONT, colour = "grey15",
                  nudge_x = SIDE * 0.33 * LAB$n_ko, direction = "y",
                  box.padding = 0.3, point.padding = 0.08, min.segment.length = 0, segment.colour = "grey60", segment.size = 0.15,
                  max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = c("Other knock-outs" = "grey75"), labels = c("Other knock-outs" = "All other knock-outs (about 1,900)"), name = NULL) +
  scale_fill_manual(values = c(set = "grey5", receptor = "grey5", marked = "grey5"), name = NULL,
                    labels = c(set = if (grepl("^GO:", GO)) sprintf("%s (%s)", GO_NAME, GO) else GO_NAME,
                               receptor = "Nr3c1 (glucocorticoid receptor)", marked = if (GO == "c-named") "Ncoa1, Crebrf and Gal" else "GR regulators Ncoa1 and Crebrf")) +
  scale_shape_manual(values = c(set = 21, receptor = 21, marked = 23), name = NULL,
                     labels = c(set = if (grepl("^GO:", GO)) sprintf("%s (%s)", GO_NAME, GO) else GO_NAME,
                                receptor = "Nr3c1 (glucocorticoid receptor)", marked = if (GO == "c-named") "Ncoa1, Crebrf and Gal" else "GR regulators Ncoa1 and Crebrf")) +
  scale_x_continuous(breaks = c(1, 1000), expand = expansion(mult = 0.04)) +   # 1 and 1000 only: 500 and 1500 overlap
  facet_wrap(~ct, nrow = 1, labeller = as_labeller(CT_LAB)) +   # the y axis shared by the four panels; the titles name the neurons
  ## no panel title or subtitle: what the panel shows is given in the figure legend
  labs(x = "Knock-out rank (1 = most opposed)", y = "RFind score") +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 1.8)), shape = guide_legend(order = 1),
         colour = guide_legend(order = 2, override.aes = list(size = 0.9))) +
  theme_classic(7, base_family = FONT) +
  theme(plot.title = element_text(face = "bold", size = 7, margin = margin(b = 2)),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6, colour = "grey20"),
        axis.line = element_line(colour = "grey55", linewidth = LW_AXIS), axis.ticks = element_line(colour = "grey55", linewidth = LW_AXIS),
        strip.background = element_blank(), strip.text = element_text(size = 7, face = "bold", margin = margin(1, 0, 1, 0)),
        panel.spacing.x = unit(6, "pt"), legend.position = "bottom", legend.text = element_text(size = 7),
        ## the two legend entries one above the other, from the left edge of the figure: side by side, or centred under the
        ## panels, they are wider than the panel
        legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(0, "pt"),
        legend.location = "plot", legend.justification = "left",
        legend.key.size = unit(2.6, "mm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "pt"),
        plot.margin = margin(2, 3, 1, 2))

tag <- sub(":", "", GO)
## drawn directly at the final size; whole pixels = round(mm / 25.4 * 600)
ragg::agg_png(file.path(OUT_PPTX, sprintf("Fig4d_%s.png", tag)), width = round(W_MM / 25.4 * DPI), height = round(H_MM / 25.4 * DPI),
              units = "px", res = DPI, background = "white")
print(p); invisible(dev.off())

NUM <- H[, .(GO = GO, term = GO_NAME, cell_type = ct, KO = target, rfind, z_vs_KO, rank_from_aligned_end = rank,
             rank_from_opposed_end = rank_opp, n_ko, nuclei = n_cells)]
setorder(NUM, KO, cell_type)
fwrite(NUM, file.path(OUT, sprintf("fig4d_%s_positions.csv", tag)))
cat(sprintf("\n%% from the opposed end (100 * rank_from_opposed_end / n_ko; 0 = most opposed, 100 = most aligned)\n"))
print(dcast(NUM[, .(KO, cell_type, pct = round(100 * rank_from_opposed_end / n_ko, 1))], KO ~ cell_type, value.var = "pct"))
