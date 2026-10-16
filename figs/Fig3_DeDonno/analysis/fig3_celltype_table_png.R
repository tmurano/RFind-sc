## Fig. 3b table — rendering only, drawn at the size pasted into Fig. 3 (76.2 x 33.0 mm, 600 dpi).
##   Separate from fig3_celltype_table.R (which sources this file at its end) so that the panel can be redrawn from the saved
##   table without re-reading the h5ad and the count matrix. Layout as the Fig. 2b table: Arial 6 pt, the columns sized to
##   their longest entry and the rows filling the pasted height, no outer margin.
##
## Input : figs/Fig3_DeDonno/output/Fig3b_celltype_table.csv   (written by fig3_celltype_table.R)
## Output: figs/Fig3_DeDonno/output/Fig3b_celltype_table.png
## Run   : Rscript "figs/Fig3_DeDonno/analysis/fig3_celltype_table_png.R"   (from project root)

if (dev.cur() == 1L) pdf(NULL)   # no Rplots.pdf when run alone; when sourced by fig3_celltype_table.R its null device is kept
suppressPackageStartupMessages({ library(data.table); library(grid); library(gridExtra); library(gtable) })
OUT <- "figs/Fig3_DeDonno/output"
OUT_PPTX <- OUT   # pasted-size panels go to output/
tab <- fread(file.path(OUT, "Fig3b_celltype_table.csv"))

W_MM <- 76.2; H_MM <- 33.0                                            # pasted size
FAM <- "Arial"; FS <- 6                                               # font size (pt) of all cells
PAD <- 1.0                                                            # left/right padding inside a cell (mm)
px <- function(mm) round(mm / 25.4 * 600)                            # ragg truncates mm sizes, so give whole pixels
lwd_pt <- function(pt) pt / (72.27 / 96)                              # grid lwd that gives `pt` points

NR <- nrow(tab); NC <- ncol(tab)
col_hj   <- c(0.50, 0.00, 0.00)
col_face <- c("bold", "italic", "plain")                                 # gene symbols italic
ragg::agg_png(file.path(OUT_PPTX, "Fig3b_celltype_table.png"), width = px(W_MM), height = px(H_MM), units = "px",
              res = 600, background = "white")
## column widths: the longest entry (header included) + padding; the spare width is shared out evenly
txt_w <- function(s, face) max(vapply(s, function(z) convertWidth(grobWidth(textGrob(z, gp = gpar(fontfamily = FAM,
                                 fontsize = FS, fontface = face))), "mm", valueOnly = TRUE), 0))
need <- vapply(seq_len(NC), function(j) max(txt_w(tab[[j]], col_face[j]), txt_w(names(tab)[j], "bold")) + 2 * PAD, 0)
wid <- need + (W_MM - sum(need)) / NC
stopifnot(all(wid >= need))
col_x <- ifelse(col_hj == 0.5, 0.5, PAD / wid)                         # left-aligned text starts PAD from the cell edge
cat(sprintf("column widths (mm): %s (needed %s)\n", paste(round(wid, 1), collapse = ", "), paste(round(need, 1), collapse = ", ")))
th <- ttheme_minimal(
  padding = unit(c(0, 0), "mm"),
  core = list(
    fg_params = list(hjust = rep(col_hj, each = NR), x = rep(col_x, each = NR), fontsize = FS, fontfamily = FAM,
                     fontface = rep(col_face, each = NR), col = "grey12"),
    bg_params = list(fill = c("#F4F7FA", "white"), col = NA)),
  colhead = list(
    fg_params = list(fontface = "bold", hjust = col_hj, x = col_x, fontsize = FS, fontfamily = FAM, col = "white"),
    bg_params = list(fill = "#34495E", col = "#34495E", lwd = lwd_pt(0.5))))
g <- tableGrob(tab, rows = NULL, theme = th)
g$widths <- unit(wid, "mm")
g$heights <- unit(rep((H_MM - 0.3) / (NR + 1), NR + 1), "mm")         # rows fill the height (0.3 mm left for the rule)
g <- gtable_add_grob(g, segmentsGrob(x0=0, x1=1, y0=0, y1=0, gp=gpar(col="#34495E", lwd=lwd_pt(1.0))), t=1, b=1, l=1, r=NC)
g <- gtable_add_grob(g, segmentsGrob(x0=0, x1=1, y0=0, y1=0, gp=gpar(col="grey70", lwd=lwd_pt(0.5))), t=nrow(g), b=nrow(g), l=1, r=NC)
grid.newpage()
pushViewport(viewport(y = 1, height = unit(H_MM - 0.3, "mm"), just = "top"))
grid.draw(g)
invisible(dev.off())
cat("wrote Fig3b_celltype_table.png\n")
