## Fig. 3b table — the De Donno cell types analysed in Fig. 3c-f (same layout as the Fig. 2b table, fig2_cluster_table.R).
##
## Columns: cell type (label used in Fig. 3) | marker genes | name. Also computed, not shown in the table: the numbers of
##          cells (all cells after our QC, control line; checked against the QC'd metadata) and the author annotation that
##          each type merges (h5ad `author_cell_type`, number suffixes joined).
## Markers (the rule of the Fig. 2b table): the markers reported by the original study (ORIG below) are ranked in a
## one-vs-rest list computed here (mean log1p CPM in the type - rest; genes detected in >= 10% of the type), and the three
## ranked highest among those within the top 16 are shown:
##   fig3_celltype_markers_top30.csv   one-vs-rest top-30 genes per type
##   fig3_celltype_marker_check.csv    rank of every marker of the original study in its type's list
##
## Input : data/DeDonno2026/output/{counts_universe.rds, cell_meta_qc.csv}, data/DeDonno2026/raw/dedonno2026.h5ad
## Output: figs/Fig3_DeDonno/output/Fig3b_celltype_table.csv + the two check tables above
##         + Fig3b_celltype_table.png (via fig3_celltype_table_png.R; drawn at the pasted size)
## Run   : Rscript "figs/Fig3_DeDonno/analysis/fig3_celltype_table.R"   (from project root)

pdf(NULL)   # the default device without a file (no Rplots.pdf); the panels are drawn on their own devices below
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(hdf5r)
                                 library(grid); library(gridExtra); library(gtable) })
source("figs/Fig3_DeDonno/analysis/fig3_celltypes.R")
OUT <- "figs/Fig3_DeDonno/output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

meta <- fread("data/DeDonno2026/output/cell_meta_qc.csv")
stopifnot(all(meta$cell_type %in% CT_LV))
## rows = only the cell types shown in the violins c, d (the 8 representative types);
## the list is read from fig3_violins.R's output so the table and the violins always agree (run the violins first)
NCF <- file.path(OUT, "Fig3_celltype_ncells.csv")
stopifnot(file.exists(NCF))
shown <- fread(NCF)[shown == TRUE]$cell_type
lv <- CT_LV[CT_LV %in% shown]
N <- meta[, .(N_all = .N, N_ctrl = sum(line == "Ctrl")), by = cell_type]; setkey(N, cell_type)

## ---------------- author annotation, straight from the h5ad ----------------
f <- H5File$new("data/DeDonno2026/raw/dedonno2026.h5ad", "r")
rc <- function(p){ o <- f[[p]]; cats <- o[["categories"]]$read(); cd <- o[["codes"]]$read(); cats[cd + 1L] }
amap <- unique(data.table(cell_type_short = rc("obs/cell_type_short"), author = rc("obs/author_cell_type")))
f$close_all()
amap[, cell_type := gsub(" [0-9]+$", "", cell_type_short)]
AUTH <- amap[, {
  a <- sort(unique(author)); base <- unique(gsub(" [0-9]+$", "", a))
  num <- sub(".* ([0-9]+)$", "\\1", a[grepl(" [0-9]+$", a)])
  .(author = if (length(num)) sprintf("%s %s", base, paste(sort(as.integer(num)), collapse = ", ")) else base)
}, by = cell_type]
stopifnot(!anyDuplicated(AUTH$cell_type))
setkey(AUTH, cell_type)

## ---------------- one-vs-rest markers (same statistic as the Fig. 2b table) ----------------
X <- readRDS("data/DeDonno2026/output/counts_universe.rds"); cf <- X$counts; cs <- X$cs; rm(X); gc()
stopifnot(ncol(cf) == nrow(meta))
lc <- cf %*% Diagonal(x = 1e6 / cs); lc@x <- log1p(lc@x)
bin <- cf; bin@x[] <- 1
tot <- Matrix::rowSums(lc); n_all <- ncol(cf); genes <- rownames(cf)
MK <- rbindlist(lapply(lv, function(k) {
  j <- which(meta$cell_type == k); nk <- length(j)
  m_in <- Matrix::rowSums(lc[, j, drop = FALSE]) / nk; m_out <- (tot - m_in * nk) / (n_all - nk)
  det <- Matrix::rowSums(bin[, j, drop = FALSE]) / nk
  z <- data.table(cell_type = k, gene = genes, diff = m_in - m_out, det_in = det)[det_in >= 0.10][order(-diff)]
  z[, rank := seq_len(.N)]
}))
fwrite(MK[rank <= 30, .(cell_type, rank, gene, diff = round(diff, 3), det_in = round(det_in, 3))],
       file.path(OUT, "fig3_celltype_markers_top30.csv"))

## ---------------- markers reported by the original study, and names ----------------
## De Donno et al. (2026): the top 15 marker genes per cluster in its Supplementary Fig. 1 (Wilcoxon test; read from the figure;
## clusters with a number suffix joined into the cell types used here) and the markers named in its main text (Prox1 for the
## dentate gyrus and Neurod6 for CA; Fibcd1 and Spock1 to separate CA1 from CA3). Genes outside the genes analysed (e.g. mt-)
## have no rank.
ORIG <- list(
  `DG Glut` = c("mt-Nd1 Meg3 Gria2 mt-Cytb mt-Atp6 mt-Co3 mt-Nd4 mt-Co2 mt-Co1 Camk2b Dgkh Camk2a Celf2 Adcy1 mt-Nd2",      # DG Glut 1
                "Pcp4 Nrgn Ppp3ca Olfm1 Ncdn Chn1 Bex2 Snca Calm2 Map1b Sncb Snap25 Cplx2 Ywhah Caly", "Prox1"),           # DG Glut 2
  `CA1 Glut` = c("Meg3 Snhg11 Gria2 Atp2b1 Opcml Grin2b Cpne7 Celf2 Shank1 Camk2a R3hdm1 Miat mt-Nd1 Cttnbp2 Grm5",       # CA1 Glut 1
                 "Nrgn Cck Snhg11 Meg3 Atp6v1g2 Snrpn Calm2 Sncb Stmn3 Atp1b1 Fkbp1a Snap25 Snca Mllt11 Ywhah",            # CA1 Glut 2
                 "Atp1b1 Cck Rtn1 Atp6v0c Ly6h Ndfip1 Tspan13 Gpm6a Atp2b1 Caly Scg5 Meg3 Prnp Gria2 Snhg11",              # CA1 Glut 3
                 "Neurod6 Fibcd1 Spock1"),
  `CA3 Glut` = c("Meg3 Snhg11 Slc17a7 Gria2 Syt1 Nptxr mt-Nd1 Opcml Atp1b1 mt-Nd4 Atp2b1 Syp Dnm1 Rnf112 Camk2a",       # CA3 Glut 1
                 "Snap25 Atp6v1g2 Snrpn Snca Nrgn Selenow Gng3 Syt1 Snhg11 Meg3 Ywhah Chn1 Basp1 Hpca Mllt11",             # CA3 Glut 2
                 "Calm2 Snap25 Snca Ywhah Stmn3 Hpca Rtn1 Pcsk1n Selenow Snrpn Scn1b Calm1 Calm3 Atp1b1 Atp6v0c",          # CA3 Glut 3
                 "Neurod6 Fibcd1 Spock1"),
  GABA = "Snhg11 Gad2 Gad1 Syt1 Atp1b1 Meg3 Nap1l5 Nrxn3 Slc32a1 Ttc3 Snap25 Slc6a1 Ndrg4 Elavl2 Rab3c",
  Astro = c("Aldoc Clu Mt1 Apoe Mt3 Mt2 Slc1a2 Cst3 Slc1a3 Ndrg2 Dbi Gstm1 Ckb Sparcl1 Ntsr2",                            # Astro 1
            "Apoe Slc1a2 Atp1a2 Clu Slc1a3 Ntsr2 Gpr37l1 Cst3 Plpp3 Gm3764 Bcan Pla2g7 Htra1 Gja1 Dclk1"),                  # Astro 2
  Microglia = "Hexb C1qa C1qb C1qc Ctss Ctsd Lgmn Cst3 Csf1r P2ry12 Selplg Tyrobp Fcer1g Trem2 Laptm5",
  Oligodendrocytes = "Plp1 Mbp Cnp Cldn11 Trf Mag Aplp1 Sept4 Mal Mobp Cryab Enpp2 Car2 Mog Fth1",
  Endothelial = "Cldn5 Ly6c1 Bsg Itm2a Ly6a Igfbp7 Flt1 Slco1a4 Pltp Tsc22d1 Slc2a1 Ly6e Pglyrp1 Klf2 Id3")
ORIG <- lapply(ORIG, function(x) unique(unlist(strsplit(x, " "))))
TOP_RANK <- 16L; N_SHOW <- 3L
NAME <- c(`DG Glut` = "Dentate gyrus granule neurons", `CA1 Glut` = "CA1 pyramidal neurons",
          `CA2 Glut` = "CA2 pyramidal neurons", `CA3 Glut` = "CA3 pyramidal neurons",
          `Mossy Glut` = "Mossy glutamatergic neurons", GABA = "GABAergic neurons",
          `Cajal-Retzius` = "Cajal-Retzius cells", `Granule Neuroblasts` = "Granule-cell neuroblasts",
          NPCs = "Neural progenitor cells", Astro = "Astrocytes", Microglia = "Microglia",
          Macrophages = "Macrophages", Oligodendrocytes = "Oligodendrocytes",
          MFOLs = "Myelin-forming oligodendrocytes", OPCs = "Oligodendrocyte precursor cells",
          COPs = "Committed oligodendrocyte precursors", Endothelial = "Endothelial cells",
          Vascular = "Vascular cells (smooth-muscle markers)", Pericytes = "Pericytes",
          Meningeal = "Vascular leptomeningeal cells", Ependymal = "Ependymal cells",
          `Choroid Plexus` = "Choroid plexus epithelial cells")
stopifnot(all(lv %in% names(ORIG)), setequal(names(NAME), CT_LV))

CHK <- rbindlist(lapply(lv, function(k) {
  g <- ORIG[[k]]; z <- MK[cell_type == k]
  data.table(cell_type = k, marker = g, in_universe = g %in% genes, rank = z$rank[match(g, z$gene)])
}))
fwrite(CHK, file.path(OUT, "fig3_celltype_marker_check.csv"))
cat(sprintf("=== markers of the original study within the top %d of the type's one-vs-rest list (gene:rank) ===\n", TOP_RANK))
print(CHK[!is.na(rank) & rank <= TOP_RANK][order(rank), .(markers = paste(sprintf("%s:%d", marker, rank), collapse = "  ")), by = cell_type], right = FALSE)

## The table shows cell type, markers and name (as the Fig. 2b table); the cell numbers are in Fig3_celltype_ncells.csv
## (fig3_violins.R) and the authors' annotation in AUTH.
## the table shows the three markers of the original study ranked highest within the top 16 (three fit the pasted width)
TOP3 <- vapply(lv, function(k) { z <- CHK[cell_type == k & !is.na(rank) & rank <= TOP_RANK][order(rank)]
  if (nrow(z) < N_SHOW) stop(sprintf("%s: fewer than %d markers of the original study within the top %d", k, N_SHOW, TOP_RANK), call. = FALSE)
  paste(head(z$marker, N_SHOW), collapse = ", ") }, "")
tab <- data.table(`Cell type` = unname(CT_LAB[lv]),
                  `Top markers` = unname(TOP3[lv]), Description = unname(NAME[lv]))
stopifnot(sum(N$N_all) == nrow(meta))                    # N covers every QC'd cell; the table shows the analysed types

## ---------------- rendering (as the Fig. 2b table) ----------------
## the PNG is drawn at the size pasted into Fig. 3 by fig3_celltype_table_png.R, which reads the table
## written here (so the panel can be redrawn without re-reading the h5ad and the count matrix)
fwrite(tab, file.path(OUT, "Fig3b_celltype_table.csv"))
source("figs/Fig3_DeDonno/analysis/fig3_celltype_table_png.R")
cat("\nwrote Fig3b_celltype_table.png, Fig3b_celltype_table.csv + fig3_celltype_{markers_top30,marker_check}.csv\n")
