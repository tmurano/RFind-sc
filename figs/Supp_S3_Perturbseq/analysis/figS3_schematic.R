## Supplementary Fig. S3a — schematic of the procedure. The panel follows the code (figS3_rfind_allpairs.py = 03_axis_vs_ko.py):
##   nuclei grouped by guide (knock-outs, safe-targeting groups, non-targeting) -> pseudo-bulk = one column per group, one matrix
##   genes x groups -> Delta of every column against the summed counts of all columns -> genes sorted in each column (top = UP,
##   bottom = DOWN) -> Running Fisher score of every pair of columns (panel b; the all-pairs comparison is a plain Running
##   Fisher comparison of the columns, so the last step is labelled "Running Fisher").
##   Numbers = the thalamic glutamatergic neurons of panels b-d (TH Prkcd Grin2c Glut, pb_g151.npz).
## Red = higher / up, blue = lower / down (the paper's colours); knock-out identities by dot colour; safe-targeting orange (as b).
## Headings = the steps (what is done, e.g. "convert to Delta" rather than "Delta"); the method is in the grey line below each
##   heading.
## Drawn at its size in the whole figure (176 x 37.4 mm = 4157 x 883 px at 600 dpi; the whole
##   figure is 180 mm wide), Arial: headings 7 pt bold black, labels black, sub-captions dark grey 6 pt; line widths in pt (cell
##   outlines and nuclei 0.5 pt, UP/DOWN brackets 0.75 pt, step arrows 1.3 pt with a 3.0 x 2.8 mm head); other lengths (nuclei
##   1.9 mm) in drawing units, scaled with the drawing; the matrices have square cells and fill the height; the Delta and
##   sorted columns and the pair matrix on the colour map of panel b.
## Output: figs/Supp_S3_Perturbseq/output/FigS3a_schematic.png (at its size in the whole figure)
## Run   : Rscript figs/Supp_S3_Perturbseq/analysis/figS3_schematic.R   (from project root)

pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(grid) })
OUT <- "figs/Supp_S3_Perturbseq/output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)   # panel at its size in the whole figure
RED <- "#C0392B"; BLUE <- "#1F5FA8"                                   # red = UP (higher), blue = DOWN (lower)
div <- colorRampPalette(c(BLUE, "#AFC9E8", "white", "#F2B8B0", RED))(101)   # = the colour map of panel b
dcol <- function(v, m = 1) div[1 + round(50 * (pmax(pmin(v / m, 1), -1) + 1))]   # |v| >= m -> full red / blue
grey_ramp <- colorRampPalette(c("grey90", "grey15"))(101)              # pseudo-bulk greys
gcol <- function(v) grey_ramp[1 + round(100 * pmax(pmin(v, 1), 0))]
KO  <- c("#2A9D8F", "#7B4FB8", "#9A8C2B")     # three knock-outs (teal, purple, olive); not red/blue (= up/down) nor orange (= safe)
SAFE <- "#F0A030"; NON <- "grey55"; GFP <- "black"   # safe-targeting orange as panel b; nucleus outline black
SUBC <- "grey20"; LINE <- "grey30"            # sub-captions (dark grey); cell outlines

W <- 14.6; H <- 3.1                          # size of the drawing in drawing units (K below converts them to inches)
DPI <- 600; WMM <- 176.0; HMM <- 37.4        # size of the panel (mm)
WPX <- round(WMM / 25.4 * DPI); HPX <- round(HMM / 25.4 * DPI)   # 4157 x 883 px (whole pixels to ragg: units = "mm" truncates)
K  <- WPX / DPI / W                           # inches of the panel per drawing unit
LW <- function(pt) pt * 96 / 72.27            # line width in pt -> grid lwd (1 lwd = 1/96 in; checked: 0.5 pt = 4.15 px at 600 dpi)
MM <- function(mm) unit(mm / 25.4, "in")
FH <- 7; FT <- 6; FD <- 8                     # pt: headings, other text, the dots (⋮ …; drawn as marks)
draw <- function() {
  grid.newpage(); pushViewport(viewport(xscale = c(0, W), yscale = c(0, H)))
  U  <- function(v) unit(v, "native")
  R  <- function(x, yt, w, h, fill = NA, col = LINE, lwd = 0.5)            # lwd in pt
    grid.rect(U(x), U(H - yt), U(w), U(h), just = c("left", "top"), gp = gpar(fill = fill, col = col, lwd = LW(lwd), linejoin = "mitre"))
  Tx <- function(lab, x, yt, just = "centre", size = FT, face = "plain", col = "black", rot = 0)
    grid.text(lab, U(x), U(H - yt), just = just, rot = rot, gp = gpar(fontsize = size, fontfamily = "Arial", fontface = face, col = col, lineheight = 0.95))
  YH <- 0.19; YS <- 0.42                      # headings (bold) and the grey line below them
  Head <- function(lab, x, sub) { Tx(lab, x, YH, size = FH, face = "bold"); Tx(sub, x, YS, col = SUBC) }
  Ln <- function(xs, ys, col = "black", lwd = 0.75)                          # lwd in pt
    grid.lines(U(xs), U(H - ys), gp = gpar(col = col, lwd = LW(lwd), linejoin = "mitre", lineend = "butt"))
  ## the arrows between the steps: shaft 1.3 pt, head a filled triangle 3.0 mm long and 2.8 mm wide whose tip ends
  ## exactly at x1
  Step <- function(x0, x1, y = YA, hl = 3.0, hw = 1.4) {
    grid.lines(unit.c(U(x0), U(x1) - MM(hl * 0.6)), U(H - c(y, y)), gp = gpar(col = "black", lwd = LW(1.3), lineend = "butt"))
    grid.polygon(unit.c(U(x1), U(x1) - MM(hl), U(x1) - MM(hl)), unit.c(U(H - y), U(H - y) + MM(hw), U(H - y) - MM(hw)),
                 gp = gpar(col = NA, fill = "black"))
  }
  Ci <- function(x, yt, r, fill = NA, col = "black", lwd = 0.5) grid.circle(U(x), U(H - yt), unit(r * K, "in"), gp = gpar(fill = fill, col = col, lwd = LW(lwd)))
  nucleus <- function(x, yt, dot, r = 0.078) { Ci(x, yt, r, fill = "#EEF3F8", col = GFP); Ci(x, yt, r * 0.5, fill = dot, col = NA) }
  vdots <- function(x, yt) Tx("⋮", x, yt, size = FD, col = "grey15")
  ## a matrix of columns; NA in 'cols' = a gap drawn as "..." (the ~2,000 columns in between)
  mat <- function(x, yt, cols, heads, cs = 0.2, cw = 0.2, step = 0.25) {
    for (k in seq_along(cols)) {
      xk <- x + (k - 1) * step
      if (is.null(cols[[k]])) { Tx("…", xk + cw / 2, yt + 4 * cs, size = FD, col = "grey15"); next }
      R(xk, yt - 0.09, cw, 0.065, fill = heads[k], col = NA)
      for (i in seq_along(cols[[k]])) R(xk, yt + (i - 1) * cs, cw, cs, fill = cols[[k]][i])
    }
    vdots(x + (length(cols) - 1) * step / 2 + cw / 2, yt + 8 * cs + 0.17)
  }
  mini <- function(x, yt, m, cs = 0.2) for (i in 1:nrow(m)) for (j in 1:ncol(m))
    R(x + (j - 1) * cs, yt + (i - 1) * cs, cs, cs, fill = if (is.na(m[i, j])) "white" else dcol(m[i, j], 0.8))

  px <- 3.85; py <- 0.72; cs <- 0.2           # matrices: left of the first, top of the cells, cell height
  YA <- py + 4 * cs                           # the arrows at mid-height of the matrices

  ## ---- nuclei, grouped by the guide they carry (the safe-targeting groups drawn as several rows: there are 100 of them) ----
  XL <- 1.55; X0 <- 1.71; NP <- 0.215         # right end of the labels, first nucleus, nucleus pitch
  Tx("group nuclei by guide", 1.68, YH, size = FH, face = "bold")
  rows <- list(list("knock-out 1", KO[1], 0.50), list("knock-out 2", KO[2], 0.77), list("knock-out 1,947", KO[3], 1.30),
               list("safe-targeting 1", SAFE, 1.67), list("safe-targeting 100", SAFE, 2.20), list("non-targeting", NON, 2.58))
  for (r in rows) { Tx(r[[1]], XL, r[[3]], just = "right")
                    for (k in 0:2) nucleus(X0 + k * NP, r[[3]], r[[2]]) }
  for (k in 3:6) nucleus(X0 + k * NP, 2.58, NON)                         # non-targeting: more nuclei
  vdots(X0 + NP, 1.035); vdots(XL - 0.47, 1.035); vdots(X0 + NP, 1.935); vdots(XL - 0.47, 1.935)

  ## ---- pseudo-bulk: the counts of the nuclei of each group summed into one column -> one matrix ----
  set.seed(3)
  pbv <- function(v) gcol(pmin(v + runif(8, -0.05, 0.05), 1))
  prof <- c(0.55, 0.35, 0.8, 0.6, 0.25, 0.45, 0.1, 0.3)                # a common gene profile; groups differ a little
  cols_pb <- list(pbv(prof * c(1.1, 0.8, 1, 1, 1.2, 0.9, 1, 1)), pbv(prof * c(1.05, 0.8, 1, 1.1, 1.2, 0.9, 1, 0.9)), NULL,
                  pbv(prof * c(0.8, 1.2, 1, 0.9, 0.8, 1.1, 1, 1.2)), pbv(prof), NULL, pbv(prof), pbv(prof))
  mat(px, py, cols_pb, c(KO[1], KO[2], NA, KO[3], SAFE, NA, SAFE, NON), cs = cs)
  Head("sum per group", px + 0.975, "(pseudo-bulk; one column each)")
  Tx("genes", px - 0.14, py + 4 * cs, rot = 90)
  Tx("knock-outs, safe-targeting,\nnon-targeting", px + 0.975, py + 8 * cs + 0.48)
  Step(3.12, px - 0.3)

  ## ---- Delta of every column against the summed counts of all columns ----
  dx <- 6.5
  set.seed(5)
  base <- c(0.6, -0.4, 0.9, 0.3, -0.7, 0.2, -0.2, 0.5)
  ko12 <- base + rnorm(8, 0, 0.25)                                       # knock-outs 1 and 2 change the same genes
  dv <- list(ko12 + rnorm(8, 0, 0.15), ko12 + rnorm(8, 0, 0.15), NULL, c(-0.5, 0.7, -0.3, -0.6, 0.4, 0.2, 0.8, -0.9),
             rnorm(8, 0, 0.18), NULL, rnorm(8, 0, 0.18))
  hd <- c(KO[1], KO[2], NA, KO[3], SAFE, NA, SAFE)
  MD <- 0.7                                                              # |Delta| >= 0.7 drawn full red / blue
  mat(dx, py, lapply(dv, function(v) if (is.null(v)) NULL else dcol(v, MD)), hd, cs = cs)
  Head("convert to Delta", dx + 0.85, "(vs. all columns summed)")
  Step(px + 2.1, dx - 0.15)

  ## ---- sort the genes of each column: UP and DOWN (no list length drawn) ----
  sx <- 9.05
  sv <- lapply(dv, function(v) if (is.null(v)) NULL else dcol(sort(v, decreasing = TRUE), MD))
  mat(sx, py, sv, hd, cs = cs)
  Head("sort genes", sx + 0.85, "(in each column)")
  bx <- sx + 6 * 0.25 + 0.28
  Ln(c(bx, bx + 0.07, bx + 0.07, bx), c(py, py, py + 2 * cs, py + 2 * cs), col = RED)
  Ln(c(bx, bx + 0.07, bx + 0.07, bx), c(py + 6 * cs, py + 6 * cs, py + 8 * cs, py + 8 * cs), col = BLUE)
  Tx("UP", bx + 0.14, py + cs, just = "left", col = RED)
  Tx("DOWN", bx + 0.14, py + 7 * cs, just = "left", col = BLUE)
  Step(dx + 1.85, sx - 0.15)

  ## ---- compare every pair of columns (panel b; 8 x 8 columns, drawn as a matrix) ----
  mx <- 12.55; mc <- cs                                                 # same cell size and top as the columns
  g <- c(1, 1, 1, 2, 2, 3, 4, 4)                                      # knock-outs 1-3 alike, 4-5 alike and opposed to 1-3; 1,947; safe
  m <- matrix(0, 8, 8)
  set.seed(11)
  for (i in 1:8) for (j in 1:8) {
    m[i, j] <- if (g[i] == 1 && g[j] == 1) 0.9 else if (g[i] == 2 && g[j] == 2) 0.8 else
      if ((g[i] == 1 && g[j] == 2) || (g[i] == 2 && g[j] == 1)) -0.6 else 0 }
  m <- (m + t(m)) / 2 + matrix(rnorm(64, 0, 0.06), 8); m <- (m + t(m)) / 2; diag(m) <- NA
  mini(mx, py, m, cs = mc)
  hx <- c(KO[1], KO[2], "#8C8C8C", "#8C8C8C", "#8C8C8C", KO[3], SAFE, SAFE)
  for (k in 1:8) { R(mx + (k - 1) * mc, py - 0.09, mc, 0.065, fill = hx[k], col = NA)
                   R(mx - 0.09, py + (k - 1) * mc, 0.065, mc, fill = hx[k], col = NA) }
  Head("compare every pair", mx + 4 * mc, "(Running Fisher)")
  Tx("similar knock-outs\nform blocks", mx + 4 * mc, py + 8 * mc + 0.36)
  Step(bx + 0.62, mx - 0.25)
  popViewport()
}
ragg::agg_png(file.path(OUT, "FigS3a_schematic.png"), width = WPX, height = HPX, units = "px", res = DPI, background = "white"); draw(); invisible(dev.off())
## PNG only
cat(sprintf("wrote FigS3a_schematic.png (%d x %d px = %.1f x %.1f mm at %d dpi)\n", WPX, HPX, WPX / DPI * 25.4, HPX / DPI * 25.4, DPI))
