## Fig. 2b table — identity of the 14 author-defined Hammond clusters.
##
## The "Dominant age", "Top markers" and "Identity / state" text is a curated summary; this script CHECKS it against
## the data and writes the evidence next to the figure:
##   fig2_cluster_age_composition.csv  % of each cluster's cells per age/condition
##   fig2_cluster_markers_top30.csv    one-vs-rest top-30 genes per cluster
##                                     (mean log1p CPM in cluster - rest; genes detected in >=10% of the cluster)
##   fig2_cluster_marker_check.csv     rank of every listed marker in its cluster's one-vs-rest list
##
## Input : data/Hammond2019/output/{counts_qc.rds, scores_hammond_9axes.csv}, data/Hammond2019/raw/cell_meta_hammond.csv
## Output: figs/Fig2_Hammond/output/Fig2b_cluster_identity.png (72.6 x 45 mm, 600 dpi = the size pasted in the figure)
##         + output/Fig2b_cluster_identity.csv + the three check tables above

pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(grid); library(gridExtra); library(gtable) })
OUT <- "figs/Fig2_Hammond/output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lv <- c("1","2a","2b","2c","3","4","5","6","7a","7b","7c","8","9","Mono/Mac")

## ---------------- data: cluster sizes and age composition ----------------
sc   <- fread("data/Hammond2019/output/scores_hammond_9axes.csv")[, .(cell_id, cluster)]
meta <- fread("data/Hammond2019/raw/cell_meta_hammond.csv")[, .(cell_id, age)]
d <- merge(sc, meta, by = "cell_id")
## the embryonic age is written as in Hammond et al. 2019, "E14.5" (the data key stays "E14", the GEO sample name)
agelab <- c(E14 = "E14.5", `P4/P5` = "P4/P5", P30 = "P30", P100 = "P100", P540 = "P540",
            P100_Control_WM = "P100 saline WM", P100_Injured_WM = "P100 inj. WM")
d[, age := factor(agelab[age], levels = agelab)]
N <- d[, .N, by = cluster]; setkey(N, cluster)
comp <- dcast(d[, .N, by = .(cluster, age)][, pct := round(100 * N / sum(N), 1), by = cluster],
              cluster ~ age, value.var = "pct", fill = 0)
comp <- comp[match(lv, cluster)]
fwrite(comp, file.path(OUT, "fig2_cluster_age_composition.csv"))

## ---------------- data: one-vs-rest markers ----------------
X  <- readRDS("data/Hammond2019/output/counts_qc.rds")            # gene x cell, raw counts
cl <- d$cluster[match(colnames(X), d$cell_id)]; stopifnot(!anyNA(cl))
lc <- X %*% Diagonal(x = 1e6 / Matrix::colSums(X)); lc@x <- log1p(lc@x)
bin <- X; bin@x[] <- 1
tot_sum <- Matrix::rowSums(lc); n_all <- ncol(X)
MK <- rbindlist(lapply(lv, function(k) {
  j <- which(cl == k); nk <- length(j)
  m_in  <- Matrix::rowSums(lc[, j, drop = FALSE]) / nk
  m_out <- (tot_sum - m_in * nk) / (n_all - nk)
  det   <- Matrix::rowSums(bin[, j, drop = FALSE]) / nk
  z <- data.table(cluster = k, gene = rownames(X), diff = m_in - m_out, det_in = det)[det_in >= 0.10][order(-diff)]
  z[, rank := seq_len(.N)]
}))
fwrite(MK[rank <= 30, .(cluster, rank, gene, diff = round(diff, 3), det_in = round(det_in, 3))],
       file.path(OUT, "fig2_cluster_markers_top30.csv"))

## ---------------- curated summary (checked against the data below) ----------------
tab <- data.table(
  Cluster = lv,
  ## dominant age = the most frequent age or condition among the cells of the cluster (one or two; "mixed" where none
  ## predominates; fig2_cluster_age_composition.csv). The numbers of cells are not in the table (Fig2_cluster_ncells.csv).
  `Dominant age` = c("P4/P5","E14.5–P4/P5","E14.5–P4/P5","E14.5–P4/P5","E14.5","P4/P5","E14.5–P4/P5","E14.5–P4/P5",
                     "P100","P100/P30","P100","mixed","P100 LPC WM","mixed"),   # 7c = P100 (37.8% of its cells; others <= 17%)
  ## markers and identity follow Hammond et al. 2019's own wording:
  ## genes the authors name for each cluster (Results text and Figs 1-5) are listed first when they rank in the
  ## top ~16 of this dataset's one-vs-rest list; the rest are this dataset's top markers. 7a-7c: the authors say
  ## they "were not defined by any unique genes", so the shared top genes are listed.
  `Top markers` = c(
    "Arg1, Cd74, Ccl24",
    "Rrm2, Pcna, Mcm6",
    "Ube2c, Birc5, H2afx",
    "Cenpa, Hist1h2bc, Ccnb2",
    "Fabp5, Mif, Ldha, Tpi1",
    "Spp1, Igf1, Gpnmb",
    "Hmox1, Gclm, Prdx1",
    "Ms4a7, Ms4a6b, Ms4a6c, Ccr1",
    "Selplg, Tmem119, Mafb",
    "Selplg, Mafb, Tmem119",
    "Selplg, Tmem119, Mafb",
    "Ccl4, Ccl3, Cd83",
    "Ifi27l2a, Ifitm3, Rtp4",
    "F13a1, H2-Aa, Pf4, Ms4a7"),
  `Identity / state` = c(
    "Small Arg1+ population",
    "Proliferative microglia",
    "Proliferative microglia",
    "Proliferative microglia",
    "Metabolically active (Fabp5+ Mif+)",
    "Axon tract-associated microglia (ATM)",
    "Hmox1+ microglia",
    "Ms4a7+, resemble border macrophages",
    "Juvenile/adult; no unique genes",
    "Juvenile/adult; no unique genes",
    "Juvenile/adult; no unique genes",
    "Inflammatory (Ccl4+); expand in aging",
    "White-matter injury (Ifi27l2a+)",
    "Non-microglial monocytes/macrophages"))
stopifnot(sum(N$N) == nrow(d))
tab <- tab[Cluster != "Mono/Mac"]   # monocytes/macrophages are not shown in the figures; their markers stay in the CSVs

## ---------------- check listed markers against the one-vs-rest ranking ----------------
expand <- function(s) {                      # "Mcm2/5/6" -> Mcm2, Mcm5, Mcm6 ; "Ms4a6c/7" -> Ms4a6c, Ms4a7
  s <- gsub("\\)", "", gsub("\\s*\\(\\+mt,\\s*", ", ", s)); toks <- trimws(strsplit(s, ",")[[1]])
  unlist(lapply(toks, function(t) { p <- strsplit(t, "/")[[1]]
    if (length(p) == 1) return(p)
    stem <- sub("[0-9]+[a-z]*$", "", p[1]); c(p[1], paste0(stem, p[-1])) }))
}
CHK <- rbindlist(lapply(seq_len(nrow(tab)), function(i) {
  g <- setdiff(expand(tab$`Top markers`[i]), "histones")
  z <- MK[cluster == tab$Cluster[i]]
  data.table(cluster = tab$Cluster[i], marker = g, rank = z$rank[match(g, z$gene)])   # matrix gene names keep "-" (03_extract_qc.R)
}))
fwrite(CHK, file.path(OUT, "fig2_cluster_marker_check.csv"))
cat("\n=== listed markers: rank in their cluster's one-vs-rest list (NA = detected in <10% of cluster) ===\n")
print(CHK[, .(markers = paste(sprintf("%s:%s", marker, ifelse(is.na(rank), "NA", rank)), collapse = "  ")), by = cluster], right = FALSE)
cat("\n=== age composition (% of cluster) ===\n"); print(comp)

## ---------------- rendering at the pasted size ----------------
## The table is drawn directly at 72.6 x 45 mm (600 dpi, dpi stored in the PNG), Arial 5 pt (header 5 pt bold), and fills
## the canvas (no outer margin). Style: header #34495E with white bold text, zebra rows #F4F7FA / white,
## cluster names bold, gene symbols italic. At 5 pt the longest rows do not fit on one line in 72.6 mm, so a cell whose
## text is wider than its column is broken onto two lines at a space (words unchanged; the CSV keeps one line).
library(systemfonts)
FAM <- "Arial"; FS <- 5; W_MM <- 72.6; H_MM <- 45; PADX <- 0.4          # PADX = left/right text inset in a cell (mm)
NR <- nrow(tab); NC <- ncol(tab)
col_hj   <- c(0.50, 0.00, 0.00, 0.00)
col_face <- c("bold", "plain", "italic", "plain")               # gene symbols italic
txt_w <- function(x, face) string_width(x, family = FAM, size = FS, res = 600, weight = ifelse(face == "bold", "bold", "normal"),
                                         italic = face == "italic") / 600 * 25.4                       # mm
## column widths: cluster / dominant age / markers fit their longest entry (header included); markers ignore the few entries
## that are much longer than the rest (they wrap); the identity column takes the remaining width
hdr_w <- txt_w(names(tab), "bold")
cw <- numeric(NC)
for (j in 1:2) cw[j] <- max(hdr_w[j], txt_w(tab[[j]], col_face[j])) + 2 * PADX
w3 <- sort(txt_w(tab[[3]], col_face[3]), decreasing = TRUE)
cw[3] <- max(hdr_w[3], w3[4]) + 2 * PADX                         # the 3 longest marker lists wrap
cw[4] <- W_MM - sum(cw[1:3])
## greedy word wrap at spaces to the column width
wrap_cell <- function(x, face, width) {
  if (txt_w(x, face) <= width) return(x)
  w <- strsplit(x, " ")[[1]]; lines <- character(0); cur <- w[1]
  for (k in w[-1]) { cand <- paste(cur, k); if (txt_w(cand, face) <= width) cur <- cand else { lines <- c(lines, cur); cur <- k } }
  paste(c(lines, cur), collapse = "\n")
}
disp <- copy(tab)
for (j in seq_len(NC)) disp[[j]] <- vapply(tab[[j]], wrap_cell, "", face = col_face[j], width = cw[j] - 2 * PADX)
nl <- apply(as.matrix(disp), 1, function(r) max(lengths(strsplit(r, "\n"))))   # lines per row
cat(sprintf("column widths (mm): %s ; rows on two lines: %s\n", paste(round(cw, 2), collapse = " / "),
            paste(tab$Cluster[nl > 1], collapse = ", ")))
## row heights: one extra text line (5 pt x 1.05 line height) for a wrapped row, everything else shares the 45 mm evenly
LINE <- FS / 72.27 * 25.4 * 1.05
h1 <- (H_MM - sum(nl - 1) * LINE) / (NR + 1)
col_x <- ifelse(col_hj == 0.5, 0.5, PADX / cw)
th <- ttheme_minimal(
  core = list(
    fg_params = list(hjust = rep(col_hj, each = NR), x = rep(col_x, each = NR), fontsize = FS, fontfamily = FAM, lineheight = 1.05,
                     fontface = rep(col_face, each = NR), col = "grey12"),
    bg_params = list(fill = c("#F4F7FA", "white"), col = NA)),
  colhead = list(
    fg_params = list(fontface = "bold", hjust = col_hj, x = col_x, fontsize = FS, fontfamily = FAM, col = "white"),
    bg_params = list(fill = "#34495E", col = "#34495E", lwd = 0.5)))
g <- tableGrob(disp, rows = NULL, theme = th)
g$widths  <- unit(cw, "mm")
g$heights <- unit(c(h1, h1 + (nl - 1) * LINE), "mm")
g <- gtable_add_grob(g, segmentsGrob(x0=0, x1=1, y0=0, y1=0, gp=gpar(col="#34495E", lwd=0.8)), t=1, b=1, l=1, r=NC)
g <- gtable_add_grob(g, segmentsGrob(x0=0, x1=1, y0=0, y1=0, gp=gpar(col="grey70", lwd=0.5)), t=nrow(g), b=nrow(g), l=1, r=NC)

## No title on the figure: what the table shows is stated in the figure legend.
OUT_PPTX <- OUT   # pasted-size panels go to output/
ragg::agg_png(file.path(OUT_PPTX, "Fig2b_cluster_identity.png"), width = round(W_MM / 25.4 * 600), height = round(H_MM / 25.4 * 600),
              units = "px", res = 600, background = "white")
grid.newpage(); grid.draw(g); invisible(dev.off())
fwrite(tab, file.path(OUT, "Fig2b_cluster_identity.csv"))
cat("\nwrote Fig2b_cluster_identity.png + Fig2b_cluster_identity.csv + fig2_cluster_{age_composition,markers_top30,marker_check}.csv\n")
