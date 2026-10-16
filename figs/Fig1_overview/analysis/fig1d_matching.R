## Fig. 1d — Step 2 of RFind-sc: the per-cell ranked Delta-list is compared with an external reference (signed DEG list, or an
## unsigned gene set) by rank similarity (Running Fisher test); the four directional scores are combined into one signed index
## per cell and reference (schematic; the heading "Step 2: comparison with an external reference -> index" and the panel letter
## stay pptx text and are NOT drawn here).
## Drawn as in the figure layout (panel d of the figure slide; colours, words and arrangement as there) and at the size it
##   is pasted (180 mm wide, 600 dpi, Arial). Left part: labels 7 pt, small sub-labels (Gene-a ... Gene-z, "signed DEG list",
##   "unsigned gene set") 6 pt, subscripts 5 pt. Right part (formulas and signed axis): 7 pt, subscripts 5 pt, set as one block
##   centred on the ranked columns, the axis as wide as the longest formula. The + and − (left and right part) are 10 pt bold
##   and the vertical dots 8 pt: symbols, drawn larger than the 7-pt text so that they stay visible (as the 8-pt dots of
##   Supplementary Fig. S3a). Cell outlines, brackets and the axis 0.75 pt, the lines of the comparison 0.5 pt.
##   Same cell size as Fig. 1c.
## Content: the two ranked columns (#FF0000 #FF5555 #FFAAAA white white #AAAAFF #5555FF #0000FF: three steps of the pure
##   red and blue of Fig. 1c), the reference UP
##   (rows 1-3, red) and DOWN (rows 6-8, blue) brackets, the unsigned gene set (six cells #80B280, bracket and "set" #006400),
##   "or", the solid lines (UP-UP, DOWN-DOWN, "+") and the dotted crossing lines (UP-DOWN, DOWN-UP, "−"), the formulas
##   Score_(cell list, reference list) = −log10 P, index (S_total) = ..., unsigned set: S_total = ..., and the signed axis.
##   The crossing lines are drawn symmetric (row 2 <-> row 7 of the lists).
## Definitions follow analysis/running_fisher.R: S_total = Score(UP,UP) + Score(DOWN,DOWN) - Score(UP,DOWN) - Score(DOWN,UP);
##   unsigned set: Score(UP,set) - Score(DOWN,set). Schematic only: no data values.
## Written as two panels, the left part (the lists and the comparison) and the right part (the formulas and the signed
##   axis), cut in the empty gap between them (see the end of the script).
## Output: figs/Fig1_overview/output/Fig1d_matching_left.png, Fig1d_matching_right.png (PNG only)
## Run   : from the project root, R 4.5 with ragg:
##         $R_HOME/bin/exec/R --no-save --no-restore -f figs/Fig1_overview/analysis/fig1d_matching.R

pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(grid) })
OUT <- "figs/Fig1_overview/output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DPI <- 600; WMM <- 180; HMM <- 44.4                      # size in the pptx (mm)
WPX <- round(WMM / 25.4 * DPI); HPX <- round(HMM / 25.4 * DPI)   # whole pixels
FONT <- "Arial"; FS <- 7; FSS <- 6; FSUB <- 5            # pt: labels, small sub-labels, subscripts
DSUB <- 0.55                                             # drop of a subscript (mm)
FSR <- 7; FSUBR <- 5; DSUBR <- DSUB                     # right part (formulas, axis): 7 pt (Nature: text <= 7 pt), subscripts 5 pt
FSYM <- 10; FDOT <- 8                                    # pt: the + and − (bold) and the vertical dots (symbols)
LW  <- function(pt) pt * 96 / 72                         # line width in pt -> grid lwd (1 lwd = 1/96 in)
CAP <- function(pt) 0.716 * pt * 25.4 / 72               # cap height of Arial (mm)

## ---- colours of the pptx ----
RED <- "#FF0000"; BLUE <- "#0000FF"; SETC <- "#80B280"; SETD <- "#006400"
LIST <- c("#FF0000", "#FF5555", "#FFAAAA", "#FFFFFF", "#FFFFFF", "#AAAAFF", "#5555FF", "#0000FF")   # both ranked lists (three steps per colour)

## ---- layout (mm from the left / from the top) ----
cs <- 3.4; NR <- 8                                       # cell size (= Fig. 1c)
CT <- 9.5; CB <- CT + NR * cs                            # top / bottom of the ranked columns
rowc <- function(i) CT + (i - 0.5) * cs
xg  <- 9.8                                               # per-cell ranked Delta-list
xl0 <- xg + cs + 2.6; xl1 <- xl0 + 24                    # the lines of the comparison
xr  <- xl1 + 2.0                                         # reference: signed DEG list
xs  <- xr + cs + 18.4; ST <- CT + cs                     # unsigned gene set (six cells, centred on the lists)
xf  <- 96.0; YB <- 9.2                                   # formulas: left edge, centre of the first line (block centred on the columns)

draw <- function(vx0 = 0, vx1 = WMM) {                   # vx0-vx1: the part of the panel drawn on the device (mm)
  grid.newpage(); pushViewport(viewport(xscale = c(vx0, vx1), yscale = c(0, HMM)))
  U  <- function(v) unit(v, "native")
  Y  <- function(y) U(HMM - y)                            # y from the top
  Rect <- function(x, y, w, h, fill = NA, col = "black", lwd = 0.75)
    grid.rect(U(x), Y(y), U(w), U(h), just = c("left", "top"), gp = gpar(fill = fill, col = col, lwd = LW(lwd), linejoin = "mitre"))
  Ln <- function(xs, ys, col = "black", lwd = 0.5, lty = "solid", lineend = "butt")
    grid.lines(U(xs), Y(ys), gp = gpar(col = col, lwd = LW(lwd), lty = lty, lineend = lineend, linejoin = "mitre"))
  Tri <- function(xs, ys, col = "black") grid.polygon(U(xs), Y(ys), gp = gpar(fill = col, col = NA))
  ## rich text: pieces "text" or S(text, face, sub); y = centre of the capitals (vcentre) or the baseline
  S  <- function(t, face = "plain", sub = FALSE, size = NULL) list(t = t, face = face, sub = sub, size = size)
  pw <- function(p, base, ssub = FSUB) convertWidth(grobWidth(textGrob(p$t, gp = gpar(fontsize = if (p$sub) ssub else
                                         (p$size %||% base), fontfamily = FONT, fontface = p$face))), "mm", valueOnly = TRUE)
  RT <- function(pieces, x, y, hjust = 0, size = FS, col = "black", rot = 0, vcentre = TRUE, ssub = FSUB, dsub = DSUB) {
    pieces <- lapply(pieces, function(p) if (is.character(p)) S(p) else p)
    w <- vapply(pieces, pw, 0, base = size, ssub = ssub); tot <- sum(w)
    th <- rot * pi / 180; ux <- cos(th); uy <- -sin(th); dx <- sin(th); dy <- cos(th)
    if (vcentre) { x <- x + dx * CAP(size) / 2; y <- y + dy * CAP(size) / 2 }
    a <- -hjust * tot
    for (k in seq_along(pieces)) {
      p <- pieces[[k]]; d <- if (p$sub) dsub else 0
      grid.text(p$t, U(x + ux * a + dx * d), Y(y + uy * a + dy * d), just = c("left", "bottom"), rot = rot,
                gp = gpar(fontsize = if (p$sub) ssub else (p$size %||% size), fontfamily = FONT, fontface = p$face, col = col))
      a <- a + w[k]
    }
    invisible(tot)
  }
  sb <- function(t, face = "plain") S(t, face, sub = TRUE)
  column <- function(x, y0, cols) for (i in seq_along(cols)) Rect(x, y0 + (i - 1) * cs, cs, cs, fill = cols[i])
  bracket <- function(x, y0, y1, col, lab) {
    Ln(c(x, x), c(y0 + 0.35, y1 - 0.35), col = col, lwd = 0.75)
    RT(list(lab), x + 1.4, (y0 + y1) / 2, col = col)
  }

  ## ---- per-cell ranked Delta-list ----
  column(xg, CT, LIST)
  RT(list("Per-cell ranked"), xg + cs / 2, CT - 5.55, hjust = 0.5)
  RT(list("Delta-list"), xg + cs / 2, CT - 2.6, hjust = 0.5)
  for (g in list(c("a", 1), c("b", 2), c("c", 3), c("d", 4), c("z", 8)))
    RT(list("Gene-", S(g[1], "italic")), xg - 1.0, rowc(as.integer(g[2])), hjust = 1, size = FSS)
  xdots <- xg - 1.0 - pw(S("Gene-a"), FSS) / 2
  RT(list("\u2026\u2026"), xdots + 0.1, rowc(6), hjust = 0.5, size = FDOT, rot = 90, vcentre = FALSE)

  ## ---- rank similarity: UP-UP and DOWN-DOWN (solid, +), UP-DOWN and DOWN-UP (dotted, −) ----
  Ln(c(xl0, xl1), rep(rowc(1), 2)); Ln(c(xl0, xl1), rep(rowc(8), 2))
  Ln(c(xl0, xl1), c(rowc(2), rowc(7)), lty = "12", lineend = "round")
  Ln(c(xl0, xl1), c(rowc(7), rowc(2)), lty = "12", lineend = "round")
  xm <- (xl0 + xl1) / 2; L <- xl1 - xl0
  SYM <- function(t, x, y) RT(list(S(t, "bold")), x, y, hjust = 0.5, size = FSYM)
  SYM("+", xm, CT + 1.25 * cs); SYM("+", xm, CT + 6.75 * cs)
  SYM("\u2212", xl0 + 0.2 * L, CT + 4 * cs); SYM("\u2212", xl1 - 0.2 * L, CT + 4 * cs)
  RT(list("Rank similarity"), xm, CB + 3.0, hjust = 0.5)
  RT(list("(Running Fisher test)"), xm, CB + 5.95, hjust = 0.5)

  ## ---- reference gene list: signed DEG list (UP / DOWN) or unsigned gene set ----
  RT(list("Reference gene list"), (xr + xs + cs) / 2, 1.7, hjust = 0.5)
  column(xr, CT, LIST)
  RT(list("signed DEG list"), xr + cs / 2 + 2.2, CT - 2.1, hjust = 0.5, size = FSS)
  bracket(xr + cs + 1.6, CT, CT + 3 * cs, RED, "UP")
  bracket(xr + cs + 1.6, CB - 3 * cs, CB, BLUE, "DOWN")
  RT(list("or"), xr + cs + 14.4, CT + 4 * cs, hjust = 0.5)
  column(xs, ST, rep(SETC, 6))
  RT(list("unsigned"), xs + cs / 2, ST - 4.6, hjust = 0.5, size = FSS)
  RT(list("gene set"), xs + cs / 2, ST - 2.15, hjust = 0.5, size = FSS)
  bracket(xs + cs + 1.6, ST, ST + 6 * cs, SETD, "set")

  ## ---- right part: 7 pt, subscripts 5 pt; the formulas and the signed axis as one block centred on the ranked columns ----
  RR <- function(pieces, x, y, ...) RT(pieces, x, y, size = FSR, ssub = FSUBR, dsub = DSUBR, ...)
  ## the score of one comparison, the index (S_total, two lines) and the unsigned-set version
  y1 <- YB; y2 <- y1 + 4.6; y3 <- y2 + 3.35; y4 <- y3 + 4.6
  e1 <- RR(list("Score", sb("(cell list, reference list)"), " = −log", sb("10"), S("P", "italic")), xf, y1)
  w1 <- RR(list("index (", S("S", "italic"), sb("total", "italic"), ") = "), xf, y2)
  w2 <- RR(list("Score", sb("(UP, UP)"), " + Score", sb("(DOWN, DOWN)")), xf + w1, y2)
  wm <- pw(S("− "), FSR)
  w3 <- RR(list("− Score", sb("(UP, DOWN)"), " − Score", sb("(DOWN, UP)")), xf + w1 - wm, y3)
  w4 <- RR(list("unsigned set: ", S("S", "italic"), sb("total", "italic"), " = Score", sb("(UP, set)"), " − Score",
                sb("(DOWN, set)")), xf, y4)

  ## ---- signed axis under the formulas, as wide as the longest one: opposite to (−) / same as (+) the reference ----
  xa0 <- xf + 0.4; xa1 <- xf + max(e1, w1 + w2, w1 - wm + w3, w4); ya <- y4 + 9.0
  if (vx0 == 0) cat(sprintf("right ends of the formulas (mm): %.1f %.1f %.1f %.1f; axis %.1f-%.1f at %.1f\n",
                            xf + e1, xf + w1 + w2, xf + w1 - wm + w3, xf + w4, xa0, xa1, ya))
  ah <- 2.3; aw <- 0.8                                                          # arrow heads: length, half width (mm)
  Ln(c(xa0 + ah * 0.6, xa1 - ah * 0.6), c(ya, ya), lwd = 0.75)
  Tri(c(xa0, xa0 + ah, xa0 + ah), c(ya, ya - aw, ya + aw)); Tri(c(xa1, xa1 - ah, xa1 - ah), c(ya, ya - aw, ya + aw))
  xz <- (xa0 + xa1) / 2; Ln(c(xz, xz), c(ya - 1.1, ya + 1.1), lwd = 0.75)
  RR(list("0"), xz, ya + 3.0, hjust = 0.5)
  SYM("−", xa0 + 5.5, ya - 3.4); SYM("+", xa1 - 5.5, ya - 3.4)
  RR(list("opposite to"), xa0 + 6.8, ya + 3.0, hjust = 0.5)
  RR(list("the reference"), xa0 + 6.8, ya + 5.95, hjust = 0.5)
  RR(list("same as"), xa1 - 6.8, ya + 3.0, hjust = 0.5)
  RR(list("the reference"), xa1 - 6.8, ya + 5.95, hjust = 0.5)
  popViewport()
}

## ---- output: the left and the right part as separate panels ----
## The whole panel is drawn once to find the columns with drawing; it is cut at XCUT (an empty gap:
## the left part ends at about 74 mm, the right part starts at xf), and each part is drawn again on its own device, trimmed
## left and right to its drawing plus MARG and keeping the full height (so that the tops align in the pptx).
XCUT <- 85.0; MARG <- 1.0                                # mm
split_panel <- function(files) {
  mmpx <- WMM / WPX                                      # mm per pixel
  tmp <- tempfile(fileext = ".png")
  ragg::agg_png(tmp, width = WPX, height = HPX, units = "px", res = DPI, background = "white"); draw(); invisible(dev.off())
  im <- png::readPNG(tmp); unlink(tmp)
  ink <- which(colSums(im[, , 1] < 0.995 | im[, , 2] < 0.995 | im[, , 3] < 0.995) > 0)   # pixel columns with drawing
  cut <- round(XCUT / mmpx); m <- round(MARG / mmpx)
  if (any(abs(ink - cut) <= m)) stop("drawing within MARG of XCUT: move XCUT into the empty gap")
  for (k in 1:2) {
    cols <- if (k == 1) ink[ink < cut] else ink[ink > cut]
    p0 <- max(0, min(cols) - 1 - m); p1 <- min(WPX, max(cols) + m)                      # pixel edges of the part
    f <- file.path(OUT, files[k])
    ragg::agg_png(f, width = p1 - p0, height = HPX, units = "px", res = DPI, background = "white")
    draw(p0 * mmpx, p1 * mmpx); invisible(dev.off())
    cat(sprintf("wrote %s (%d x %d px = %.1f x %.1f mm at %d dpi; x = %.1f-%.1f mm of the one-piece panel)\n",
                f, p1 - p0, HPX, (p1 - p0) / DPI * 25.4, HPX / DPI * 25.4, DPI, p0 * mmpx, p1 * mmpx))
  }
}
split_panel(c("Fig1d_matching_left.png", "Fig1d_matching_right.png"))
