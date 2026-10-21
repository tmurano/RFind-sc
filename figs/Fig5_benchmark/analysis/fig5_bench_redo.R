## Fig. 5 / Supplementary Table 2 — the benchmark (adapted from fig5_bench_v2.R; run with figs/Fig5_benchmark/analysis/run_fig5_redo.sh).
## fig5_bench_v2.R, fig5_final.R, fig5_reflength_sweep.R, fig5_scfocal_type.R and fig5_ucell_up.R, named below, are the scripts of an
## earlier stored run of the benchmark; they are not in the repository. Of their outputs, results/fig5_final_auc.csv and
## results/fig5_reflength_sweep_auc.csv are kept in figs/Fig5_benchmark/output/ and read for the checks below (the tasks in CHANGED,
## H1, H4 and D1-D3, are not checked against them, nor are ULM / WMEAN of the signed tasks). Every output here is named fig5_redo_*.
## Reference genes (up / down): the top L of the whole list, then the genes in the pair universe; lengths in LEN_DAM, LEN_MAT and
##   LEN_CORT below; gene sets as they are.
## As in fig5_bench_v2.R / fig5_final.R: data loading; input = ln(1 + counts per 10,000) with the library size over ALL genes (raw
##   counts where a method wants counts); AddModuleScore bins and the centring groups (Hammond: all QC cells; De Donno: all cells
##   of the cell type); boot_task (AUC; RFind-sc = column 1).
## 95% CIs: paired percentile bootstrap over MICE (boot_task: all cells of a drawn mouse; mice drawn within the
##   positive and the negative group when every mouse is in one group, else from all mice, as in the tasks defined by a cluster;
##   set.seed(1) per task, 1,000 resamples; the draws do not depend on the methods). Mouse = the sample of a Hammond cell (the name
##   before its barcode) or the donor of a De Donno cell (donor_id).
## Differences from fig5_bench_v2.R:
##   - 'UCellUp' = UCell with the up genes only (the usual use of UCell; = fig5_ucell_up.R), at every L. 'UCell' = UCell with the up
##     genes and the down genes marked "-" (UCell's own signed input; 2nd string).
##   - 'ORA' (2nd string) = decoupleR run_ora with its defaults: per cell the top 5% of the genes by expression (ties at the boundary
##     broken at random, seed 42, as run_ora does), one-sided Fisher test of the overlap with the reference genes against a
##     background of 20,000 genes, score = -log10 P; a signed reference gives its up genes (as the other gene-set methods).
##     run_ora builds a long table of every gene x cell (far too slow for the atlas), so its computation is transcribed
##     here (ora_scores) and checked against decoupleR::run_ora on 200 cells of every reference (must agree exactly).
##   - H1 (DAM, LPC vs saline white matter): the 25 non-microglial cells (Mono/Mac) are removed, as in the other Hammond tasks.
##   - for signed references, 'ULM' / 'WMEAN' (decoupleR with weights +1 / -1) are not methods of the paper and are not written
##     (NOT_WRITTEN_PM1 below; 'ULM' / 'WMEAN' of the gene sets, weight 1, are the decoupleR of the type B tasks and are written);
##     they come out of dc_ss(), which gives ULM, WMEAN and singscore in one pass.
##   - every method scores each reference once, on the union of the cells of the tasks using it (H2 / H3 / H5 share
##     cells); every score is a per-cell value (ssGSEA's joint constant changes no AUC and no sign).
##   - Beyondcell cells without a score (NaN): the rule of the main table is BEYONDCELL_NAN (below); the number and fraction of
##     these cells and the AUC with CI with them set to 0 are written to results/fig5_redo_beyondcell_L500.csv (L = 500 only; at the other
##     L, to the summary in <logdir> only).
##   - checks against the stored results run for the tasks not in CHANGED only; outputs carry the L in their names.
##   - the De Donno data are restricted to the wild-type mice (analysis/dedonno/01_qc_embed.R: QC, genes analysed and per-cell-type
##     baseline from the wild-type cells only), so D1-D3 are in CHANGED and their reference lengths are checked against LEN_CORT below.
## Also:
##   - 'ULMfc' / 'WMEANfc' = decoupleR run_ulm / run_wmean with the log2 fold change of the reference as the weight (mor = the
##     weights of scFOCAL-type, fcw()), instead of +1 / -1: the form of decoupleR that uses the magnitude (as PROGENy's
##     continuous weights). Signed references only (a gene set has one weight). Same input and settings as 'ULM' / 'WMEAN'.
##   - the per-cell scores are written at every L (results/fig5_redo_percell_L<L>.csv.gz), read by fig5_best_length.py.
## Not computed (2nd string, recorded only): SingleR, scmap, CHETAH,
##   clustifyr, scDRS, the two modified versions.
## Methods:
##   VISION       v3.0.2 (GitHub YosefLab/VISION). Input counts x s / library size, NOT log (s by VISION_SCALE: 1e4 = CP10k, or
##                the median library size; VISION takes log2(1 + x) itself and z-scores each
##                cell over all genes). Signed reference = one signature (+1 up / -1 down), gene set = +1. Per-cell score of
##                calcSignatureScores() with the Vision() defaults (znorm_columns, sig_gene_threshold 0.001,
##                min_signature_genes 5), pool = FALSE; no analyze(), no smoothing. One object per data set (Hammond: all QC
##                cells; De Donno: all cells of the cell type). Checked against the closed form (max |diff| < 1e-8).
##   Beyondcell   v2.2.2 (GitHub cnio-bu/beyondcell, tag v2.2.2). GMT with REF_UP / REF_DOWN; perform.reversal = FALSE
##                (positive = the direction of the reference); bcScore(expr.thres = 0.1) on ln(1 + CP10k) (normalised,
##                unscaled) of the reference genes; score = bc@normalized (no bcRegressOut / bcSubset). NaN (cells that
##                express < 10% of the reference genes) handled by BEYONDCELL_NAN. Checked against a line-by-line
##                transcription of bcScore (bc_manual; must agree exactly).
##   ssGSEA       GSVA gsva(ssgseaParam()), alpha 0.25, minSize 1, maxSize Inf; ranks over all genes of the universe.
##                Signed reference and normalisation by SSGSEA_SIGNED / SSGSEA_NORMALIZE. Scored in blocks of cells with
##                normalize = FALSE; the GSVA constant is applied once over the task cells (= one unchunked call; checked).
##   scFOCAL-type as fig5_scfocal_type.R: Spearman rho between the reference log2FC and ln(1 + CP10k) z-scored per gene over the
##                centring group of cos_focal() (the atlas; De Donno: the cell type), clipped by SCFOCAL_CLIP (Seurat ScaleData clips only above +10;
##                fig5_scfocal_type.R clipped at +-10); signed references only (blank for gene sets).
## Checks (full runs), in the tasks not in CHANGED: at L = 500, AddModuleScore is always computed, selected or not, and its AUC must
##   reproduce results/fig5_final_auc.csv (the script stops otherwise); RFind-sc must reproduce AUC and CI of
##   results/fig5_final_auc.csv (L = 500) or results/fig5_reflength_sweep_auc.csv (other L). At L = 500 every method of
##   results/fig5_final_auc.csv is also compared with it, task by task (log only: 'reproduced' / 'DIFFERS').
## Arguments (after --args, any order):
##   methods=all (default) | a comma list of
##     UCellUp,UCell,AUCell,AddModuleScore,ULM,WMEAN,ULMfc,WMEANfc,ORA,singscore,VISION,Beyondcell,ssGSEA,scFOCAL
##     (RFind-sc = column 1)
##   L=500 (default = the paper) | 100 | 250 | 1000 | 2000 | all  — the signed reference lists are rebuilt by the rules of
##     fig5_reflength_sweep.R (lists_hammond / lists_cort, restricted to the pair universe) and given to EVERY method; RFind-sc is
##     recomputed at the same L as in fig5_reflength_sweep.R; at L != 500 only the signed tasks (H1, H4, D1-D3). L = 500 must
##     rebuild the paper's lists exactly (the script stops otherwise).
##   test — 600 cells per task (300 positive + 300 negative; random, seed 1), output to the folder FIG5_TEST_DIR
##     (environment variable; default logs/runs/fig5_test/, which git ignores); for timing and debugging.
##   logdir=<folder> — where the summary goes (default logs/runs/fig5_redo_<yymmdd>).
## Output (not test), figs/Fig5_benchmark/output/results/ (<suffix> = _L<L>, plus _<methods> for a subset run):
##   fig5_redo_auc<suffix>.csv             the columns of results/fig5_final_auc.csv + n_na, L, ref_up, ref_down, note
##   fig5_redo_percell<suffix>.csv.gz      per-cell scores of every task (Beyondcell raw, with its NaN); 'group' = age (Hammond) or
##                                         condition (De Donno)
##   fig5_redo_beyondcell<suffix>.csv      Beyondcell: number and fraction of NaN cells and the AUC with CI with them set to 0
##                                         (rule 'as_zero'); L = 500 only
##   <logdir>/fig5_redo<suffix>.txt        summary: settings, package versions, checks, AUC table
## Run   : R 4.5.2, from the project root; ~10 GB RAM like fig5_final.R (VISION on the Hammond atlas is
##   computed before the AUCell rankings so the two peaks do not add up). All six lengths:
##   bash figs/Fig5_benchmark/analysis/run_fig5_redo.sh (keeps the Mac awake; logs in logs/runs/fig5_redo_<yymmdd>/).

pdf(NULL); options(warn = 1)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(parallel) })

## ============= settings =============
BEYONDCELL_NAN   <- "as_zero"       # a NaN counts against Beyondcell (NaN cells are not excluded).
                                    #   Options: 'as_zero', 'as_min', 'worst'
##   Cells with NaN (fewer than 10% of the reference genes expressed, or none):
##   'as_zero'     NaN -> 0 (the sign-change point of a signed score; the bottom for a gene set, whose scores are >= 0)
##   'worst'       every NaN counted as an error: NaN positives ranked below every cell, NaN negatives above every cell
##   'as_min'      NaN cells ranked below every other cell (tied)
##   The NaN fraction and the Beyondcell AUC with CI with the NaN cells set to 0 are written to results/fig5_redo_beyondcell_L500.csv.
SSGSEA_SIGNED    <- "up_minus_down" # options: 'up_minus_down', 'up_only'
##   Signed references: 'up_minus_down' = ES(up) - ES(down) from one call (GenePattern ssGSEA 'combine' mode; ssGSEA2.0);
##   'up_only' = ES of the up genes (as AUCell / AddModuleScore use them here). Gene sets: the ES of the set in either case.
SSGSEA_NORMALIZE <- "joint"         # GSVA default (options: 'joint', 'separate', 'none')
##   'joint'    GSVA normalize = TRUE with up and down in ONE call: ES / (max - min over both sets x all task cells) = one
##              constant, so every AUC equals 'none' (the normalisation changes up - down only with 'separate')
##   'separate' up and down each divided by their own range (= two GSVA calls); changes up - down
##   'none'     raw ES
UCELL_LONG       <- "raise"         # the setting UCell recommends (options: 'raise', 'skip'; fairness mode L >= 1,000 only)
##   Needed only in the fairness mode with L >= 1000: UCell stops when a signature has more genes than maxRank (default 1500;
##   up and down genes count together). 'raise' = maxRank = the number of signature genes, for those signatures only (UCell's
##   own error message says to increase maxRank); 'skip' = UCell reported as not applicable at that L. At L = 500 the
##   default (1500) is kept and this setting is not used.
VISION_SCALE     <- "cp10k"         # options: 'cp10k', 'median_umi'
##   VISION input (VISION takes log2(1 + x) itself; its per-cell z-scores change non-linearly with the scale, so the AUC depends on it):
##   'cp10k'      counts x 1e4 / library size (all genes) = the normalisation of every other method
##   'median_umi' counts x median library size / library size, the median over the cells given to VISION (the vignette's example)
SCFOCAL_CLIP     <- "upper"         # Seurat ScaleData clips above +10 only (options: 'upper', 'both')
##   scFOCAL-type z-scores: 'upper' = z > 10 set to 10 only, as Seurat ScaleData (scale.max = 10; FastRowScale / FastSparseRowScale
##   clip only above); 'both' = also z < -10 set to -10 (not Seurat's default; the clip of fig5_scfocal_type.R).
## ------ further settings (from the documentation and defaults of each package):
##   duplicated gene names (De Donno: 3 names x 2 rows) made unique with make.unique for VISION and ssGSEA (reference genes = the
##   first row, as fig5_final.R); Beyondcell expr.thres 0.1, no bcRegressOut / bcSubset, score = bc@normalized; the routes below.
VISION_ROUTE     <- "package"   # 'package' = VISION v3.0.2 | 'closed_form' = the same formula without the package (report
                                #   it as a reimplementation). Decide the fallback BEFORE the run if VISION cannot be installed.
BEYONDCELL_ROUTE <- "package"   # 'package' = beyondcell v2.2.2 | 'source' = R/objects.R, Genesets.R, Score.R of the tag
                                #   v2.2.2 sourced from BEYONDCELL_SRC (needs Seurat, SeuratObject, scales, stringr) |
                                #   'transcription' = bc_manual only (a reimplementation). The package code runs in the first two.
BEYONDCELL_SRC   <- "tools/beyondcell_v2.2.2"   # for 'source': git clone --branch v2.2.2 https://github.com/cnio-bu/beyondcell
BC_EXPR_THRES    <- 0.1         # bcScore default
BC_BLOCK_ELEMS   <- 2e7         # Beyondcell: genes x cells per dense block (160 MB); blocks do not change any score
VISION_SIG_GENE_THRES <- 0.001; VISION_MIN_GENES <- 5L   # Vision() defaults (used by 'closed_form'; the package uses its own)
SSG_ALPHA        <- 0.25        # GSVA default (Barbie 2009)
SSG_BLOCK        <- 3000L       # cells per ssGSEA call; no effect on the scores (normalize = FALSE per block, constant once)
OTHER_NA         <- "stop"      # NA in any method other than Beyondcell is unexpected: 'stop' | 'as_min'

## ==================================================== arguments and checks ====================================================
ARGS <- commandArgs(trailingOnly = TRUE); TEST <- "test" %in% ARGS
argv <- function(key, default) { a <- grep(paste0("^", key, "="), ARGS, value = TRUE); if (length(a)) sub(paste0("^", key, "="), "", a[1]) else default }
bad_args <- setdiff(ARGS, c("test", grep("^(methods|L|logdir)=", ARGS, value = TRUE)))
if (length(bad_args)) stop("unknown argument(s): ", paste(bad_args, collapse = " "), call. = FALSE)
LAB_ALL <- c("100", "250", "500", "1000", "2000", "all")
LARG <- argv("L", "500"); if (!LARG %in% LAB_ALL) stop("L must be one of ", paste(LAB_ALL, collapse = ", "), call. = FALSE)
LVAL <- if (LARG == "all") .Machine$integer.max else as.integer(LARG); FAIR <- LARG != "500"
KEYS <- c("UCellUp", "UCell", "AUCell", "AddModuleScore", "ULM", "WMEAN", "ULMfc", "WMEANfc", "ORA", "singscore", "VISION", "Beyondcell",
          "ssGSEA", "scFOCAL")
MARG <- argv("methods", "all")
SEL <- if (MARG == "all") KEYS else strsplit(MARG, ",")[[1]]
if (!length(SEL) || !all(SEL %in% KEYS)) stop("methods= must be all or a comma list of: ", paste(KEYS, collapse = ","), call. = FALSE)
SEL <- KEYS[KEYS %in% SEL]
COLS <- list(UCellUp = "UCell (up genes)", UCell = "UCell", AUCell = "AUCell", AddModuleScore = "AddModuleScore", ULM = "ULM", WMEAN = "WMEAN",
             ULMfc = "ULM (log2FC)", WMEANfc = "WMEAN (log2FC)", ORA = "ORA", singscore = "singscore", VISION = "VISION", Beyondcell = "Beyondcell", ssGSEA = "ssGSEA",
             scFOCAL = "scFOCAL-type")
METH_ALL <- unlist(COLS, use.names = FALSE)                       # output order
SIGNED_ONLY <- c("scFOCAL", "ULMfc", "WMEANfc")         # keys that need the order / values of the reference (not gene sets)
FCKEYS <- c(ULMfc = "ULM (log2FC)", WMEANfc = "WMEAN (log2FC)")
NOT_WRITTEN_PM1 <- c("ULM", "WMEAN")                              # weights +1 / -1: written for gene sets only (header)
CHANGED <- c("H1", "H4", "D1", "D2", "D3")                        # tasks without a stored result to reproduce: cells or reference genes differ
                                                                  # from those of the stored results (H4: Maturity reference; D1-D3: wild-type De Donno cells)

need <- c(BEYONDCELL_NAN = "Beyondcell" %in% SEL, SSGSEA_SIGNED = "ssGSEA" %in% SEL, SSGSEA_NORMALIZE = "ssGSEA" %in% SEL,
          UCELL_LONG = any(c("UCell", "UCellUp") %in% SEL) && FAIR && LVAL >= 1000L, VISION_SCALE = "VISION" %in% SEL, SCFOCAL_CLIP = "scFOCAL" %in% SEL)
todo <- names(need)[need & vapply(names(need), function(k) grepl("^TODO", get(k)), logical(1))]
if (length(todo)) stop(sprintf(paste0("\n  Settings not decided yet (edit the top of this script before any run):\n  %s\n",
                                      "  The settings are described in Methods."),
                               paste(sprintf("%s <- \"%s\"", todo, vapply(todo, get, "")), collapse = "\n  ")), call. = FALSE)
chk_opt <- function(k, v) if (need[[k]] && !get(k) %in% v) stop(sprintf("%s = '%s' is not one of: %s", k, get(k), paste(v, collapse = ", ")), call. = FALSE)
chk_opt("BEYONDCELL_NAN", c("as_zero", "as_min", "worst")); chk_opt("SSGSEA_SIGNED", c("up_minus_down", "up_only"))
chk_opt("SSGSEA_NORMALIZE", c("joint", "separate", "none")); chk_opt("UCELL_LONG", c("raise", "skip"))
chk_opt("VISION_SCALE", c("cp10k", "median_umi")); chk_opt("SCFOCAL_CLIP", c("upper", "both"))
CHK_LOAD <- !TEST && !FAIR                                        # AddModuleScore must reproduce results/fig5_final_auc.csv
stopifnot(VISION_ROUTE %in% c("package", "closed_form"), BEYONDCELL_ROUTE %in% c("package", "source", "transcription"),
          OTHER_NA %in% c("stop", "as_min"))

PKG <- list(UCellUp = "UCell", UCell = "UCell", AUCell = "AUCell", ULM = "decoupleR", WMEAN = "decoupleR", ULMfc = "decoupleR", WMEANfc = "decoupleR",
            ORA = "decoupleR",                                    # decoupleR::run_ora: the check of the transcribed ORA
            singscore = c("singscore", "GSEABase"),
            VISION = if (VISION_ROUTE == "package") "VISION",
            Beyondcell = switch(BEYONDCELL_ROUTE, package = c("beyondcell", "Seurat", "SeuratObject"),
                                source = c("Seurat", "SeuratObject", "scales", "stringr"), transcription = character(0)),
            ssGSEA = "GSVA")
need_pkg <- unique(unlist(PKG[intersect(names(PKG), SEL)]))
miss <- need_pkg[!vapply(need_pkg, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop(sprintf("Package(s) not installed: %s. Install them first (figs/Fig5_benchmark/analysis/install_r452_bench_work.sh; README), or run the other methods with methods=...",
                               paste(miss, collapse = ", ")), call. = FALSE)
if ("VISION" %in% SEL && VISION_ROUTE == "package" && packageVersion("VISION") < "3.0.1")
  stop("VISION >= 3.0.1 is needed (3.0.0 returns scores without cell names); the benchmark uses v3.0.2", call. = FALSE)
if ("Beyondcell" %in% SEL && BEYONDCELL_ROUTE == "package" && packageVersion("beyondcell") != "2.2.2")
  stop("beyondcell ", packageVersion("beyondcell"), " is installed; the benchmark uses v2.2.2", call. = FALSE)
BCENV <- NULL
if ("Beyondcell" %in% SEL && BEYONDCELL_ROUTE == "source") {
  dsc <- read.dcf(file.path(BEYONDCELL_SRC, "DESCRIPTION"))
  if (dsc[1, "Package"] != "beyondcell" || dsc[1, "Version"] != "2.2.2") stop("BEYONDCELL_SRC is not beyondcell 2.2.2", call. = FALSE)
  BCENV <- new.env(parent = globalenv())
  ## toplevel.env = globalenv(): with the default (toplevel.env = BCENV) setClass() registers "geneset" / "beyondcell" under an invented
  ## package name and new("geneset") then fails in validObject (checked with the tag v2.2.2 files on toy data)
  for (f in c("objects.R", "Genesets.R", "Score.R")) sys.source(file.path(BEYONDCELL_SRC, "R", f), envir = BCENV, toplevel.env = globalenv())
}
if ("Beyondcell" %in% SEL && BEYONDCELL_ROUTE != "transcription") options(Seurat.object.assay.version = "v3")   # maintainers' advice (#158/#159)
bc_fun <- function(f) if (BEYONDCELL_ROUTE == "package") getExportedValue("beyondcell", f) else get(f, envir = BCENV)

TAG <- if (MARG == "all") "" else paste0("_", gsub("[^A-Za-z0-9]+", "-", MARG))
SUF <- paste0("_L", LARG, TAG)
OUT <- if (TEST) Sys.getenv("FIG5_TEST_DIR", "logs/runs/fig5_test") else "figs/Fig5_benchmark/output"   # test: a separate folder (same file names)
LOGDIR <- argv("logdir", sprintf("logs/runs/fig5_redo_%s", format(Sys.Date(), "%y%m%d")))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); if (!TEST) dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)
TXT <- if (TEST) file.path(OUT, sprintf("fig5_redo%s_test.txt", SUF)) else file.path(LOGDIR, sprintf("fig5_redo%s.txt", SUF))
f_out <- function(stem, ext = "csv") file.path(OUT, "results", sprintf("%s%s.%s", stem, SUF, ext))   # the tables in output/results/
dir.create(file.path(OUT, "results"), showWarnings = FALSE, recursive = TRUE)
NC <- max(1L, parallel::detectCores() - 2L)
t0 <- Sys.time(); say <- function(...) cat(sprintf("[%5.0fs] ", as.numeric(Sys.time() - t0, units = "secs")), ..., "\n", sep = "")
LOG <- character(0); lg <- function(...) { x <- paste0(...); LOG <<- c(LOG, x); cat(x, "\n", sep = "") }   # one line
lgs <- function(...) for (x in c(...)) lg(x)                  # one line per argument
pv <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
sha <- function(p) tryCatch({ x <- suppressWarnings(packageDescription(p))$RemoteSha; if (is.null(x)) "" else paste0(" (", substr(x, 1, 7), ")") },
                             error = function(e) "")                # GitHub commit of VISION / beyondcell installed with remotes
source("analysis/running_fisher.R")                                  # pair_universe, cut_in_universe, rfind_score_cell

lg(sprintf("Benchmark of the Fig. 5 redo - %s; %s; L = %s%s; methods: %s", format(Sys.time(), "%Y-%m-%d %H:%M"),
           if (TEST) "TEST (600 cells per task)" else "all task cells", LARG, if (FAIR) " (signed tasks H1, H4, D1-D3)" else " (the paper's length)",
           paste(sub("^(ULM|WMEAN)$", "\\1 (gene sets)", unlist(COLS[SEL])), collapse = ", ")))
lg("Script: figs/Fig5_benchmark/analysis/fig5_bench_redo.R; settings as described in Methods")
lg(sprintf("%s; %s", R.version.string, R.version$platform))
PKLOG <- c("Matrix", "data.table", "UCell", "AUCell", "decoupleR", "singscore", "GSEABase", "VISION", "beyondcell", "Seurat", "SeuratObject",
           "GSVA", "BiocVersion")                         # BiocVersion = the Bioconductor release
lg("Packages: ", paste(sprintf("%s %s%s", PKLOG, vapply(PKLOG, pv, ""), vapply(PKLOG, sha, "")), collapse = "; "))
lg(sprintf("Settings: BEYONDCELL_NAN = %s; SSGSEA_SIGNED = %s; SSGSEA_NORMALIZE = %s; UCELL_LONG = %s; VISION_SCALE = %s; SCFOCAL_CLIP = %s; ",
           BEYONDCELL_NAN, SSGSEA_SIGNED, SSGSEA_NORMALIZE, UCELL_LONG, VISION_SCALE, SCFOCAL_CLIP),
   sprintf("VISION_ROUTE = %s; BEYONDCELL_ROUTE = %s%s; ", VISION_ROUTE, BEYONDCELL_ROUTE, if (BEYONDCELL_ROUTE == "source") paste0(" (", BEYONDCELL_SRC, ")") else ""),
   sprintf("BC_EXPR_THRES = %g; SSG_ALPHA = %g; OTHER_NA = %s; UCell / AUCell / AddModuleScore / decoupleR / singscore as fig5_final.R",
           BC_EXPR_THRES, SSG_ALPHA, OTHER_NA))
lg("Fixed by the script: duplicated gene names -> make.unique (VISION, ssGSEA); Beyondcell without bcRegressOut / ",
   "bcSubset, score = bc@normalized")

## ============================================ scoring functions (fig5_final.R) ============================================
auc_v <- function(s, pos) { n1 <- sum(pos); n0 <- length(s) - n1; r <- rank(s); (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
DONOR_FILE <- "data/DeDonno2026/output/cell_meta_qc.csv"
mouse_of <- function(cell_id) {                                   # the mouse of each cell: Hammond = the sample prefix of the cell name
  dd <- !grepl("_", cell_id)                                      # (e.g. P100_M_SALINE_A3); De Donno (no '_' in the names) = donor_id
  m <- sub("_[^_]*$", "", cell_id)
  if (any(dd)) { d <- fread(DONOR_FILE, select = c("cell_id", "donor_id")); m[dd] <- d$donor_id[match(cell_id[dd], d$cell_id)] }
  stopifnot(!anyNA(m)); m
}
boot_task <- function(id, task, kind, ref, S, pos, mouse, B = 1000) {  # paired percentile bootstrap over MICE (all cells of a drawn mouse)
  est <- apply(S, 2, auc_v, pos = pos)                            # mice drawn within the positive and the negative group when every mouse
  stopifnot(length(mouse) == length(pos), !anyNA(mouse))          # is in one group, else from all mice (tasks defined by a cluster)
  cm <- split(seq_along(pos), mouse)
  nested <- all(vapply(cm, function(i) length(unique(pos[i])) == 1L, logical(1)))
  mp <- sort(unique(mouse[pos])); mn <- sort(unique(mouse[!pos])); ma <- sort(names(cm))
  set.seed(1)
  bs <- t(replicate(B, {
    repeat {
      d <- if (nested) c(mp[sample.int(length(mp), length(mp), replace = TRUE)], mn[sample.int(length(mn), length(mn), replace = TRUE)])
           else ma[sample.int(length(ma), length(ma), replace = TRUE)]
      i <- unlist(cm[d], use.names = FALSE)
      if (any(pos[i]) && any(!pos[i])) break                      # (a draw from all mice must contain both groups)
    }
    apply(S[i, , drop = FALSE], 2, auc_v, pos = pos[i]) }))
  d <- bs[, 1] - bs
  data.table(id = id, task = task, definition = kind, reference = ref, method = colnames(S), n_pos = sum(pos), n_neg = sum(!pos), AUC = est,
             lo = apply(bs, 2, quantile, 0.025), hi = apply(bs, 2, quantile, 0.975),
             diff_vs_RFind = est[1] - est, d_lo = apply(d, 2, quantile, 0.025), d_hi = apply(d, 2, quantile, 0.975))
}
ams_scores <- function(ln_all, genes, ix) {                       # Seurat AddModuleScore on log-normalised values (verbatim)
  avg <- Matrix::rowMeans(ln_all); u <- rownames(ln_all)
  bin <- ceiling(rank(avg, ties.method = "first") / length(avg) * 24)
  set.seed(1)
  ctrl <- unique(unlist(lapply(genes, function(g) { pool <- u[bin == bin[match(g, u)]]; sample(pool, min(100L, length(pool))) })))
  (Matrix::colMeans(ln_all[genes, ix, drop = FALSE]) - Matrix::colMeans(ln_all[ctrl, ix, drop = FALSE]))
}
dc_ss <- function(cnt, lib, up, dn, block = 3000) {                # ULM, WMEAN (mor +1/-1) and singscore (verbatim)
  signed <- length(dn) > 0
  net <- data.frame(source = "ref", target = c(up, dn), mor = c(rep(1, length(up)), rep(-1, length(dn))))
  n <- ncol(cnt); o <- matrix(NA_real_, n, 3, dimnames = list(NULL, c("ULM", "WMEAN", "singscore")))
  for (b in split(seq_len(n), ceiling(seq_len(n) / block))) {
    m <- log1p(sweep(as.matrix(cnt[, b, drop = FALSE]), 2, 1e4 / lib[b], "*")); cl <- colnames(m)
    u <- decoupleR::run_ulm(m, net, minsize = 5); o[b, 1] <- u$score[match(cl, u$condition)]
    w <- decoupleR::run_wmean(m, net, minsize = 5, times = 100, seed = 42); w <- w[w$statistic == "wmean", ]; o[b, 2] <- w$score[match(cl, w$condition)]
    rk <- singscore::rankGenes(m)
    o[b, 3] <- if (signed) singscore::simpleScore(rk, upSet = GSEABase::GeneSet(up), downSet = GSEABase::GeneSet(dn))$TotalScore
               else singscore::simpleScore(rk, upSet = GSEABase::GeneSet(up))$TotalScore
  }
  o
}
dc_fc <- function(cnt, lib, genes, w, block = 3000) {            # ULM and WMEAN with mor = the log2FC weights w
  stopifnot(length(genes) == length(w), !anyNA(w))
  keep <- w != 0                                                  # Fold Change exactly +-1 -> weight 0 (only at L = all: Maturity 1 gene,
  net <- data.frame(source = "ref", target = genes[keep], mor = w[keep])   # Glucocorticoid 5): left out, as a weight of 0 says nothing
  n <- ncol(cnt); o <- matrix(NA_real_, n, 2, dimnames = list(NULL, unname(FCKEYS)))
  for (b in split(seq_len(n), ceiling(seq_len(n) / block))) {     # the same blocks, input and settings as dc_ss
    m <- log1p(sweep(as.matrix(cnt[, b, drop = FALSE]), 2, 1e4 / lib[b], "*")); cl <- colnames(m)
    u <- decoupleR::run_ulm(m, net, minsize = 5); o[b, 1] <- u$score[match(cl, u$condition)]
    v <- decoupleR::run_wmean(m, net, minsize = 5, times = 100, seed = 42); v <- v[v$statistic == "wmean", ]; o[b, 2] <- v$score[match(cl, v$condition)]
  }
  o
}
UCELL_NOTE <- character(0)
ucell_scores <- function(cnt, up, dn) {                            # as fig5_final.R; maxRank only touched above 1500 genes (fairness)
  f <- if (length(dn)) c(up, paste0(dn, "-")) else up
  if (length(f) <= 1500L) return(as.numeric(UCell::ScoreSignatures_UCell(cnt, features = list(S = f), ncores = NC, name = "")[colnames(cnt), "S"]))
  if (UCELL_LONG == "skip") return(NULL)
  if (UCELL_LONG != "raise") stop("UCell: a signature has ", length(f), " genes > maxRank 1500 and UCELL_LONG is not set", call. = FALSE)
  nt <- sprintf("maxRank raised to %d (signature genes; UCELL_LONG = raise)", length(f)); UCELL_NOTE <<- c(UCELL_NOTE, nt)
  structure(as.numeric(UCell::ScoreSignatures_UCell(cnt, features = list(S = f), ncores = NC, name = "", maxRank = length(f))[colnames(cnt), "S"]),
            note = nt)                                            # -> the 'note' column of the UCell row of this task
}
aucell_scores <- function(rk, up, cells) {                        # verbatim
  as.numeric(AUCell::getAUC(AUCell::AUCell_calcAUC(list(S = up), rk, aucMaxRank = ceiling(0.05 * nrow(rk)), verbose = FALSE))["S", cells])
}
## ---- ORA = decoupleR::run_ora with its defaults, transcribed (decoupleR 2.17: R/statistic-ora.R) ----
## run_ora(mat, network): n_up = ceiling(0.05 * nrow(mat)), n_bottom = 0, n_background = 20000, seed = 42, minsize = 5. Per cell the
## n_up genes with the largest values are taken, ties broken by rnorm() drawn over the long table (gene-major: pivot_longer of the
## genes x cells matrix) under withr::with_seed(42) (.ora_slice_targets: arrange(condition, value, rand), tail(n_up)); then
## fisher.test(alternative = "greater") on TP = |top ∩ set|, FP = |set \ top|, FN = |top \ set|, TN = 20000 - TP - FP - FN
## (= phyper(TP - 1, |set|, 20000 - |set|, |top|, lower.tail = FALSE)); score = -log10(P). The random tie-break depends on the
## matrix given to run_ora, so each block of ORA_BLOCK cells is scored as one run_ora call would score it. The same transcription
## is checked on the first ORA_CHECK_N cells of every reference against decoupleR::run_ora itself (must agree exactly).
ORA_BLOCK <- 500L; ORA_BG <- 20000L; ORA_SEED <- 42L; ORA_CHECK_N <- 200L; ORA_CHK <- character(0)   # 500 cells per block: ~0.2 GB per worker (4 workers)
ora_block <- function(m, set) {                                   # m = dense genes x cells (ln(1 + CP10k)), row names = genes
  G <- nrow(m); B <- ncol(m); n_up <- ceiling(0.05 * G)
  if (n_up >= G) n_up <- G                                        # .ora_check_ns
  set.seed(ORA_SEED); R <- matrix(stats::rnorm(G * B), G, B, byrow = TRUE)
  rn <- rownames(m); S <- unique(set); k <- length(S)
  vapply(seq_len(B), function(j) {
    top <- order(m[, j], R[, j], decreasing = TRUE)[seq_len(n_up)]
    obs <- unique(rn[top]); tp <- sum(obs %in% S)                  # intersect / setdiff of the contingency table count unique genes
    -log10(stats::phyper(tp - 1, k, ORA_BG - k, length(obs), lower.tail = FALSE))
  }, numeric(1))
}
ora_scores <- function(cnt, lib, set, tid) {                       # cnt = raw counts genes x cells (all genes of the universe)
  stopifnot(sum(unique(set) %in% rownames(cnt)) >= 5L)              # run_ora minsize = 5
  lnm <- function(b) log1p(sweep(as.matrix(cnt[, b, drop = FALSE]), 2, 1e4 / lib[b], "*"))
  n <- ncol(cnt); blocks <- split(seq_len(n), ceiling(seq_len(n) / ORA_BLOCK))
  res <- mclapply(blocks, function(b) ora_block(lnm(b), set), mc.cores = min(4L, NC))
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) stop(sprintf("ORA %s: %d blocks failed: %s", tid, sum(bad), paste(unique(unlist(res[bad])), collapse = " | ")), call. = FALSE)
  o <- unlist(res, use.names = FALSE)
  k <- seq_len(min(ORA_CHECK_N, n)); m0 <- lnm(k)                  # check: run_ora itself on the first cells (one call on that matrix)
  d <- suppressMessages(decoupleR::run_ora(m0, data.frame(source = "ref", target = unique(set), mor = 1)))
  ref <- d$score[match(colnames(m0), d$condition)]; mine <- ora_block(m0, set)
  same <- (is.infinite(ref) & is.infinite(mine)) | abs(ref - mine) < 1e-12
  ORA_CHK <<- c(ORA_CHK, sprintf("  ORA %s: transcription vs decoupleR::run_ora %s on %d cells: %s (max |diff| %.2g; set %d genes; %d cells scored)",
                                 tid, pv("decoupleR"), length(k), if (isTRUE(all(same))) "identical" else "DIFFERS",
                                 max(c(0, abs(ref - mine)[is.finite(ref) & is.finite(mine)])), length(unique(set)), n))
  if (!isTRUE(all(same))) stop(sprintf("ORA %s: the transcription differs from decoupleR::run_ora", tid), call. = FALSE)
  o
}
fcw <- function(tab, genes) {                                     # fig5_final.R; case-insensitive fallback only for genes not matched exactly
  i <- match(genes, tab$Symbol); j <- is.na(i); i[j] <- match(tolower(genes[j]), tolower(tab$Symbol))
  fc <- tab$`Fold Change`[i]; stopifnot(!anyNA(fc)); sign(fc) * log2(abs(fc))
}
## cosine, cosine (gene-centred) [fig5_final.R cos_scores] and scFOCAL-type [fig5_scfocal_type.R score_block], in blocks of task cells
## (same values; the per-gene mean / sd over the centring group are taken block-wise so long reference lists fit in memory)
grp_stats <- function(lnm, genes, block = 5000L) {                # per-gene mean and sd (n - 1) over ALL columns of lnm
  n <- ncol(lnm); mu <- Matrix::rowMeans(lnm[genes, , drop = FALSE]); ss <- numeric(length(genes))
  for (b in split(seq_len(n), ceiling(seq_len(n) / block))) ss <- ss + rowSums((as.matrix(lnm[genes, b, drop = FALSE]) - mu)^2)
  list(mu = mu, sd = sqrt(ss / (n - 1)))
}
cos_focal <- function(lnm, genes, w, ix_task, focal = TRUE, block = 5000L) {   # lnm = ln(1 + CP10k) genes x cells, centring group
  st <- grp_stats(lnm, genes); rw <- rank(w); nw <- sqrt(sum(w^2))
  cs <- function(A) { v <- as.vector(crossprod(A, w)) / (sqrt(colSums(A^2)) * nw); v[!is.finite(v)] <- 0; v }   # 0/0 -> 0 (fig5_final.R)
  o <- matrix(NA_real_, length(ix_task), 3, dimnames = list(NULL, c("cosine", "cosine (gene-centred)", "scFOCAL-type")))
  for (b in split(seq_along(ix_task), ceiling(seq_along(ix_task) / block))) {
    Xt <- as.matrix(lnm[genes, ix_task[b], drop = FALSE]); Xc <- Xt - st$mu
    o[b, 1] <- cs(Xt); o[b, 2] <- cs(Xc)
    if (!focal) next                                              # scFOCAL-type not selected (cosine only, e.g. for the check)
    Z <- Xc / st$sd; Z[!is.finite(Z)] <- 0; Z[Z > 10] <- 10                           # Seurat ScaleData: z > scale.max (10) -> 10
    if (SCFOCAL_CLIP == "both") Z[Z < -10] <- -10                                     # +-10, as fig5_scfocal_type.R
    o[b, 3] <- as.vector(cor(apply(Z, 2, rank), rw))                                  # Spearman = Pearson of ranks
  }
  o
}

## ---- VISION v3.0.2 ----
vision_closed <- function(lnm, rn, kept, block = 10000L) {       # VISION's score from ln(1 + CP10k): per cell, z over all genes,
  W <- Matrix::sparseMatrix(i = unlist(lapply(kept, function(d) match(names(d), rn)), use.names = FALSE),   # (sum s_g z_g) / sum |s_g|
                            j = rep(seq_along(kept), lengths(kept)), x = as.numeric(unlist(kept, use.names = FALSE)),
                            dims = c(nrow(lnm), length(kept)))
  ng <- nrow(lnm); sw <- Matrix::colSums(W); sa <- Matrix::colSums(abs(W))
  o <- matrix(NA_real_, ncol(lnm), length(kept), dimnames = list(colnames(lnm), names(kept)))
  for (b in split(seq_len(ncol(lnm)), ceiling(seq_len(ncol(lnm)) / block))) {
    Lb <- lnm[, b, drop = FALSE]; mu <- Matrix::colMeans(Lb)
    sdv <- sqrt(pmax(Matrix::colSums(Lb^2) - ng * mu^2, 0) / (ng - 1)); sdv[sdv == 0] <- 1          # VISION: scale 1 if sd = 0
    num <- as.matrix(Matrix::crossprod(Lb, W))
    o[b, ] <- sweep((num - outer(mu, sw)) / sdv, 2, sa, "/")
  }
  o
}
vision_run <- function(cnt, lib, lnm, sigs) {                     # one object per data set / cell type; sigs = list(name = list(up, dn))
  stopifnot(identical(rownames(lnm), rownames(cnt)), identical(colnames(lnm), colnames(cnt)), !is.null(colnames(cnt)), inherits(lnm, "dgCMatrix"))
  rn <- make.unique(rownames(cnt))                                # De Donno has 3 duplicated names; unique names as Seurat would give
  for (nm in names(sigs)) stopifnot(!grepl("_(UP|DOWN)$", nm, ignore.case = TRUE), all(c(sigs[[nm]]$up, sigs[[nm]]$dn) %in% rn),
                                    !length(intersect(sigs[[nm]]$up, sigs[[nm]]$dn)))
  sig_vec <- function(s) c(setNames(rep(1, length(s$up)), s$up), setNames(rep(-1, length(s$dn)), s$dn))
  vs <- if (VISION_SCALE == "cp10k") 1e4 else median(lib)         # VISION input = counts x vs / library size (VISION_SCALE), NOT log
  lnv <- function() if (VISION_SCALE == "cp10k") lnm else {       # ln(1 + input) for the closed form (log base cancels in the z-score)
    x <- cnt %*% Diagonal(x = vs / lib); x@x <- log1p(x@x); dimnames(x) <- dimnames(cnt); x }
  if (VISION_ROUTE == "package") {
    cp <- cnt %*% Diagonal(x = vs / lib); dimnames(cp) <- list(rn, colnames(cnt))             # CP10k or median-scaled, NOT log
    S <- lapply(names(sigs), function(nm) VISION::createGeneSignature(name = nm, sigData = sig_vec(sigs[[nm]])))
    vis <- suppressMessages(VISION::Vision(cp, signatures = S, pool = FALSE)); rm(cp)         # everything else at defaults
    lost <- setdiff(names(sigs), names(vis@sigData))
    if (length(lost)) stop("VISION dropped signature(s) (< min_signature_genes after its filters): ", paste(lost, collapse = ", "), call. = FALSE)
    vis <- suppressMessages(VISION::calcSignatureScores(vis, sig_gene_importance = FALSE))   # no analyze(): no kNN, no pooling
    sc <- VISION::getSignatureScores(vis)                          # cells x (signatures, <signed>_UP, <signed>_DOWN)
    if (is.null(rownames(sc))) stop("VISION returned scores without cell names (VISION >= 3.0.1 needed)", call. = FALSE)
    stopifnot(identical(rownames(sc), colnames(cnt)), identical(rownames(vis@exprData), rn), all(names(sigs) %in% colnames(sc)))
    kept <- lapply(setNames(nm = names(sigs)), function(nm) vis@sigData[[nm]]@sigDict)          # genes VISION actually used
    sc <- sc[, names(sigs), drop = FALSE]; rm(vis); invisible(gc())
    dmax <- max(abs(vision_closed(lnv(), rn, kept) - sc))
    if (!is.finite(dmax) || dmax >= 1e-8) stop(sprintf("VISION: package scores differ from the closed form (max |diff| %.3g)", dmax), call. = FALSE)
  } else {                                                        # 'closed_form': VISION's gene filter and score, without the package
    det <- tabulate(lnm@i[lnm@x > 0] + 1L, nbins = nrow(lnm)); ok <- rn[det >= VISION_SIG_GENE_THRES * ncol(lnm)]
    kept <- lapply(sigs, function(s) { d <- sig_vec(s); d[names(d) %in% ok] })
    small <- names(kept)[lengths(kept) < VISION_MIN_GENES]
    if (length(small)) stop("VISION (closed form): fewer than 5 genes left in ", paste(small, collapse = ", "), call. = FALSE)
    sc <- vision_closed(lnv(), rn, kept); dmax <- NA_real_
  }
  list(scores = sc, kept = kept, dmax = dmax, scale = vs)
}
vision_log <- function(where, V, sigs) for (nm in names(sigs))
  lg(sprintf("  VISION %s %s: input counts x %.6g / library size; genes used up %d / down %d of %d / %d; package vs closed form max |diff| %s",
             where, nm, V$scale, sum(V$kept[[nm]] > 0), sum(V$kept[[nm]] < 0), length(sigs[[nm]]$up), length(sigs[[nm]]$dn),
             if (is.na(V$dmax)) "n/a (closed_form route)" else sprintf("%.2g", V$dmax)))

## ---- Beyondcell v2.2.2 ----
bc_gmt <- function(up, dn) {                                      # set names must end in _UP / _DOWN; 2nd field must not be a gene
  f <- tempfile(fileext = ".gmt")
  writeLines(c(paste(c("REF_UP", "S3-bench-v2", up), collapse = "\t"), if (length(dn)) paste(c("REF_DOWN", "S3-bench-v2", dn), collapse = "\t")), f)
  f
}
bc_quiet <- function(expr) withCallingHandlers(expr, warning = function(w)
  if (grepl("^Using count matrix as input", conditionMessage(w))) invokeRestart("muffleWarning"))
bc_manual <- function(X, up, dn, thres, digits = 2) {            # line-by-line transcription of bcScore (R/Score.R, v2.2.2)
  rn <- tolower(rownames(X)); up <- unique(tolower(up)); dn <- unique(tolower(dn)); al <- unique(c(up, dn))
  nexp <- colSums(X[rn %in% al, , drop = FALSE] > 0)
  below <- nexp < length(al) * thres | nexp == 0
  f <- function(s) { Y <- X[rn %in% s, , drop = FALSE]; m <- colMeans(Y); sdv <- apply(Y, 2, sd); m * ((colSums(Y) - sdv) / (m + sdv)) }
  if (length(dn)) {
    u <- f(up); d <- f(dn); nn <- is.na(u) & is.na(d); u[is.na(u)] <- 0; d[is.na(d)] <- 0; s <- u - d; s[nn | below] <- NaN
  } else {
    s <- f(up); s[below] <- NaN
  }
  if (is.null(digits)) s else round(s, digits)
}
bc_scores <- function(lnm, up, dn, ix) {                          # lnm = ln(1 + CP10k) genes x cells; ix = task cells (columns)
  g <- c(up, dn)
  stopifnot(length(up) >= 2L, length(dn) == 0L || length(dn) >= 2L, !anyDuplicated(tolower(g)), all(g %in% rownames(lnm)), !is.null(colnames(lnm)))
  gs <- NULL
  if (BEYONDCELL_ROUTE != "transcription") {
    gs <- suppressMessages(bc_fun("GenerateGenesets")(bc_gmt(up, dn), perform.reversal = FALSE))
    stopifnot(identical(gs@mode, if (length(dn)) c("up", "down") else "up"), identical(gs@inverse.score, FALSE), identical(names(gs@genelist), "REF"))
  }
  nb <- max(1L, ceiling(length(ix) * length(g) / BC_BLOCK_ELEMS))
  blocks <- split(seq_along(ix), sort(rep_len(seq_len(nb), length(ix))))       # contiguous blocks of near-equal size
  out <- rep(NA_real_, length(ix))
  for (b in blocks) {
    X <- as.matrix(lnm[g, ix[b], drop = FALSE])                  # dense base matrix of the reference genes only (all bcScore reads)
    m <- unname(bc_manual(X, up, dn, BC_EXPR_THRES))
    if (BEYONDCELL_ROUTE == "transcription") s <- m else {
      bc <- tryCatch(suppressMessages(bc_quiet(bc_fun("bcScore")(X, gs, expr.thres = BC_EXPR_THRES))),
                     error = function(e) if (all(is.na(m))) NULL else stop(e))      # bcScore stops if no cell passes expr.thres
      s <- if (is.null(bc)) rep(NaN, length(b)) else unname(bc@normalized["REF", colnames(X)])
      if (!identical(is.na(s), is.na(m)) || !isTRUE(all(abs(s - m) < 1e-9, na.rm = TRUE)))
        stop("Beyondcell: package scores differ from the transcription bc_manual()", call. = FALSE)
    }
    out[b] <- s
  }
  list(score = out)
}
bc_rule_table <- function(tid, raw, pos, rfv, mouse) {           # Beyondcell: NaN cells, and the AUC with CI with them set to 0 (paired with RFind-sc)
  na <- is.na(raw); s <- raw; s[na] <- 0
  o <- if (!(any(pos) && any(!pos))) data.table(rule = "as_zero", n_cells = length(pos), n_pos = sum(pos)) else {
    b <- boot_task(tid, "", "", "", cbind(`RFind-sc` = rfv, Beyondcell = s), pos, mouse)
    data.table(rule = "as_zero", n_cells = length(pos), n_pos = sum(pos), AUC = b$AUC[2], lo = b$lo[2], hi = b$hi[2], RFind_AUC = b$AUC[1],
               diff_vs_RFind = b$diff_vs_RFind[2], d_lo = b$d_lo[2], d_hi = b$d_hi[2])
  }
  cbind(data.table(id = tid, n_nan = sum(na), frac_nan = mean(na)), o)
}

## ---- ssGSEA (GSVA) ----
SSG_CHK <- character(0)
ssgsea_raw <- function(m, gs) {                                   # m = dense genes x cells (unique rownames); sets x cells, NOT normalised
  e <- suppressMessages(GSVA::gsva(GSVA::ssgseaParam(m, gs, alpha = SSG_ALPHA, normalize = FALSE, minSize = 1, maxSize = Inf), verbose = FALSE))
  stopifnot(all(names(gs) %in% rownames(e)), identical(colnames(e), colnames(m)))
  e[names(gs), , drop = FALSE]
}
ssgsea_scores <- function(cnt, lib, up, dn, tid) {                # cnt = raw counts genes x task cells (all genes of the universe)
  rn <- make.unique(rownames(cnt))                                # GSVA stops on duplicated row names (De Donno: 3)
  two <- length(dn) > 0 && SSGSEA_SIGNED == "up_minus_down"
  gs <- if (two) list(up = up, dn = dn) else list(up = up)
  stopifnot(all(unlist(gs) %in% rn), !anyDuplicated(up), !anyDuplicated(dn), !length(intersect(up, dn)))
  n <- ncol(cnt); es <- matrix(NA_real_, length(gs), n, dimnames = list(names(gs), colnames(cnt)))
  for (b in split(seq_len(n), ceiling(seq_len(n) / SSG_BLOCK))) {
    m <- log1p(sweep(as.matrix(cnt[, b, drop = FALSE]), 2, 1e4 / lib[b], "*")); rownames(m) <- rn
    es[, b] <- ssgsea_raw(m, gs)
    if (b[1] == 1L) {                                             # (a) GSVA normalize = TRUE = one constant; (b) counts = ln(1 + CP10k)
      k <- seq_len(min(500L, length(b))); e0 <- es[, b[k], drop = FALSE]
      e1 <- suppressMessages(GSVA::gsva(GSVA::ssgseaParam(m[, k, drop = FALSE], gs, alpha = SSG_ALPHA, normalize = TRUE, minSize = 1, maxSize = Inf),
                                        verbose = FALSE))[names(gs), , drop = FALSE]
      c1 <- max(abs(e1 - e0 / diff(range(e0))))
      cnt_k <- as.matrix(cnt[, b[k], drop = FALSE]); rownames(cnt_k) <- rn
      c2 <- max(abs(ssgsea_raw(cnt_k, gs) - e0))
      SSG_CHK <<- c(SSG_CHK, sprintf("  ssGSEA %s: GSVA normalize=TRUE vs one constant max |diff| %.2g; raw counts vs ln(1 + CP10k) max |diff| %.2g",
                                     tid, c1, c2))
      if (!(c1 < 1e-10 && c2 < 1e-10)) stop(sprintf("ssGSEA check failed in %s (%.3g, %.3g)", tid, c1, c2), call. = FALSE)
    }
  }
  if (SSGSEA_NORMALIZE == "joint") es <- es / diff(range(es))
  if (SSGSEA_NORMALIZE == "separate") for (r in rownames(es)) es[r, ] <- es[r, ] / diff(range(es[r, ]))
  if (two) es["up", ] - es["dn", ] else es["up", ]
}

## ---- reference lists at length L (fig5_reflength_sweep.R, copied) and the RFind-sc recomputation of the fairness mode ----
## the rule of the six references: the top L per direction of the whole list (genes without informative names,
## mitochondrial and haemoglobin genes excluded), restricted to the pair universe when used (make_ref_L, map_genes)
bad_pc <- function(g) grepl("Rik[0-9]?$",g)|grepl("^Gm[0-9]",g)|grepl("^mt-",g)|grepl("^Hb[ab]-",g)|grepl("^A[CFLI][0-9]",g)|grepl("^[0-9]",g)|grepl("\\.[0-9]+$",g)|is.na(g)
lists_hammond <- function(tab, L) {                              # 06_build_references.R: k_up / k_dn, m_up / m_dn (tables already without bad_pc)
  list(up = head(tab[`Fold Change` > 0][order(-`Fold Change`)]$Symbol, L),
       down = head(tab[`Fold Change` < 0][order(`Fold Change`)]$Symbol, L))
}
lists_cort <- function(L) {                                      # 03_references.R: CORT
  d <- fread("data/DeDonno2026/provided/DeDonno_GR_glucocorticoid_261016.csv")
  r <- d[, .(Symbol, fc = `Fold Change`)][is.finite(fc) & fc != 0][!bad_pc(Symbol)][order(-abs(fc))][!duplicated(Symbol)]
  list(up = head(r[fc > 0][order(-fc)]$Symbol, L), down = head(r[fc < 0][order(fc)]$Symbol, L))
}
make_ref_L <- function(e, lists, genes) {
  glow <- tolower(genes); U <- pair_universe(genes, e$eligible)
  list(inu = glow %in% U, P = length(U), b2 = list(up = intersect(tolower(lists$up), U), down = intersect(tolower(lists$down), U),
                                                  ranked = TRUE))
}
check500 <- function(nm, lists, e) {
  ok <- identical(lists$up, e$up) && identical(lists$down, e$down)
  lg(sprintf("  %s: rebuilt L = 500 lists identical to the stored reference: %s (up %d / down %d)", nm, ok, length(lists$up), length(lists$down)))
  if (!ok) stop(sprintf("%s: the rebuilt L = 500 lists differ from the stored reference", nm), call. = FALSE)
}
map_genes <- function(x, rn, U) { m <- rn[match(tolower(x), tolower(rn))]; unique(m[!is.na(m) & tolower(m) %in% U]) }   # = RFind-sc's genes
## up / down lengths of the lists given to every method (RFind-sc's lists at L in the pair universe; the top L of the whole list,
## then the genes analysed) in the signed tasks: DAM (H1) and Maturity (H4) in the 8,960 genes of Hammond et al., Glucocorticoid
## (D1-D3) in the 11,647 genes of the wild-type De Donno cells. The other tasks are gene sets (no length).
LEN_DAM <- list("100" = c(71L, 75L), "250" = c(192L, 205L), "500" = c(408L, 427L), "1000" = c(869L, 880L),
                "2000" = c(1701L, 1821L), "all" = c(1701L, 3667L))
LEN_MAT <- list("100" = c(71L, 70L), "250" = c(178L, 163L), "500" = c(353L, 287L), "1000" = c(722L, 596L),
                "2000" = c(1468L, 1212L), "all" = c(4040L, 4732L))
LEN_CORT <- list("100" = c(39L, 84L), "250" = c(113L, 205L), "500" = c(243L, 401L), "1000" = c(553L, 812L),
                 "2000" = c(1274L, 1592L), "all" = c(6027L, 5557L))
chk_len <- function(tid, s) {                                     # the genes given to the competitors = RFind-sc's lists at this L
  if (tid == "H1") e <- LEN_DAM[[LARG]] else if (tid == "H4") e <- LEN_MAT[[LARG]] else if (tid %in% c("D1", "D2", "D3")) e <- LEN_CORT[[LARG]] else {
    if (!LARG %in% SW$L) return(invisible(NULL))                  # L = 100 / 250 ... are all in the sweep; kept for safety
    r <- SW[SW$id == tid & SW$L == LARG, ]
    e <- if (nrow(r) == 1L) c(r$ref_up, r$ref_down) else c(NA, NA)
  }
  if (anyNA(e) || e[1] != length(s$up) || e[2] != length(s$dn))
    stop(sprintf("%s: reference lengths %d / %d differ from the expected lengths of RFind-sc's lists at L = %s", tid, length(s$up),
                 length(s$dn), LARG), call. = FALSE)
}
NCELL <- 500L; CLIP <- 300; BLOCK <- 250L; NCORES <- 4L          # as fig5_reflength_sweep.R
sweep_cells <- function(tix, cnt, expected, deep, glow, RL) {     # copied from fig5_reflength_sweep.R
  inu <- RL[[1]]$inu; P <- RL[[1]]$P
  blocks <- split(tix, ceiling(seq_along(tix) / BLOCK))
  res <- mclapply(blocks, function(idx) {
    obs <- as.matrix(cnt[, idx, drop = FALSE]); expc <- expected(idx)
    D <- (obs - expc) / sqrt(expc + 1); rm(obs, expc)                            # P2 residual, as the pipeline
    nb <- length(idx); sc <- matrix(NA_real_, nb, length(RL)); lu <- ld <- integer(nb); mis <- 0L
    for (k in seq_len(nb)) {
      p <- unname(D[, k])
      pos <- which(p > 0); pos <- pos[order(p[pos], decreasing = TRUE)]           # full ranking, most extreme first
      neg <- which(p < 0); neg <- neg[order(p[neg], decreasing = FALSE)]
      dc <- deep$col[idx[k]]; du <- deep$up[, dc]; dd <- deep$down[, dc]
      mis <- mis + !identical(head(pos, deep$n_deep), du[!is.na(du)]) + !identical(head(neg, deep$n_deep), dd[!is.na(dd)])
      uj <- cut_in_universe(pos, inu, NCELL); dj <- cut_in_universe(neg, inu, NCELL); lu[k] <- length(uj); ld[k] <- length(dj)
      for (i in seq_along(RL)) sc[k, i] <- rfind_score_cell(glow[uj], glow[dj], RL[[i]]$b2, universe_size = P, clip = CLIP)
    }
    list(sc = sc, lu = lu, ld = ld, mis = mis)
  }, mc.cores = NCORES)
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) stop(sprintf("%d blocks failed: %s", sum(bad), paste(unique(unlist(res[bad])), collapse = " | ")))
  list(sc = do.call(rbind, lapply(res, `[[`, "sc")), lu = unlist(lapply(res, `[[`, "lu")), ld = unlist(lapply(res, `[[`, "ld")),
       mis = sum(vapply(res, `[[`, integer(1), "mis")))
}
pick <- function(ix, posf) {                                      # test: 300 positive + 300 negative cells
  if (!TEST) return(ix)
  p <- posf(ix); a <- ix[p]; b <- ix[!p]; set.seed(1)
  sort(c(a[sample.int(length(a), min(300L, length(a)))], b[sample.int(length(b), min(300L, length(b)))]))
}

## ---- one task: NA rules, bootstrap, checks ----
RES <- list(); PC <- list(); BCN <- list(); CHK <- character(0)
BM <- fread("figs/Fig5_benchmark/output/results/fig5_final_auc.csv")
SW <- fread("figs/Fig5_benchmark/output/results/fig5_reflength_sweep_auc.csv", colClasses = list(character = "L"))
finish_task <- function(tid, task, kind, refname, cell_id, rfv, pos, sc, na_note, n_up, n_dn, ck = list(), m_note = character(0),
                        grp = NA_character_) {
  nw <- if (n_dn > 0) NOT_WRITTEN_PM1 else character(0)          # computed, not written (header): ULM / WMEAN of signed references
  sc <- lapply(sc[intersect(setdiff(METH_ALL, nw), names(sc))], function(v) as.numeric(unname(v)))
  na_note <- na_note[!names(na_note) %in% nw]; m_note <- m_note[!names(m_note) %in% nw]
  stopifnot(length(rfv) == length(pos), all(lengths(sc) == length(pos)), all(lengths(ck) == length(pos)), !anyNA(pos))
  if (tid %in% CHANGED) ck <- list()                              # no stored result for a task in CHANGED
  for (m in names(ck)) {                                          # data loading of the competitor side (normalisation, centring group,
    b0 <- BM[BM$id == tid & BM$method == m, ]; a <- auc_v(as.numeric(ck[[m]]), pos)   # reference genes, FC weights), all task cells
    tol <- max(1e-6, 3 / (sum(pos) * sum(!pos)))                  # 1e-6, or 3 cell pairs (near-ties that rounding may reorder)
    d <- if (nrow(b0) == 1L) abs(a - b0$AUC) else Inf
    CHK <<- c(CHK, sprintf("  %s %s AUC %.6f vs results/fig5_final_auc.csv: |diff| %.2g -> %s", tid, m, a, d, if (d < tol) "reproduced" else "NOT reproduced"))
    if (!(d < tol)) stop(sprintf("%s %s: AUC differs from results/fig5_final_auc.csv (|diff| %.3g); check the data loading", tid, m, d), call. = FALSE)
  }
  if (!length(sc)) { lg(sprintf("%s: no selected method applies; skipped", tid)); return(invisible(NULL)) }
  S <- do.call(cbind, c(list(`RFind-sc` = as.numeric(rfv)), sc))
  PC[[tid]] <<- data.table(id = tid, cell_id = cell_id, positive = pos, group = grp, S)
  nav <- colSums(is.na(S))
  if (nav[["RFind-sc"]] > 0) stop(sprintf("%s: RFind-sc has %d NA cells", tid, nav[["RFind-sc"]]), call. = FALSE)
  oth <- setdiff(names(nav)[nav > 0], c("RFind-sc", "Beyondcell"))
  if (length(oth)) {
    msg <- sprintf("%s: NA scores in %s", tid, paste(sprintf("%s (%d cells)", oth, nav[oth]), collapse = ", "))
    if (OTHER_NA == "stop") stop(msg, " (OTHER_NA = 'stop')", call. = FALSE)
    lg("  WARNING ", msg, " -> ranked lowest (OTHER_NA = 'as_min')"); for (m in oth) S[is.na(S[, m]), m] <- -Inf
  }
  keep <- rep(TRUE, nrow(S)); note_txt <- ""
  if ("Beyondcell" %in% colnames(S)) {
    raw <- S[, "Beyondcell"]; na <- is.na(raw)
    BCN[[tid]] <<- bc_rule_table(tid, raw, pos, S[, "RFind-sc"], mouse_of(cell_id))
    lg(sprintf("  Beyondcell %s: NaN %d of %d cells (%.2f%%) -> %s", tid, sum(na), length(na), 100 * mean(na), BEYONDCELL_NAN))
    if (any(na)) {
      if (BEYONDCELL_NAN == "as_zero") { S[na, "Beyondcell"] <- 0
      } else if (BEYONDCELL_NAN == "as_min") { S[na, "Beyondcell"] <- -Inf
      } else if (BEYONDCELL_NAN == "worst") { S[na & pos, "Beyondcell"] <- -Inf; S[na & !pos, "Beyondcell"] <- Inf }
    }
  }
  S <- S[keep, , drop = FALSE]; p <- pos[keep]
  r <- boot_task(tid, task, kind, refname, S, p, mouse_of(cell_id)[keep])
  r[, c("n_na", "L", "ref_up", "ref_down", "note") := list(unname(nav[r$method]), LARG, n_up, n_dn, note_txt)]
  if (length(na_note)) r <- rbind(r, data.table(id = tid, task = task, definition = kind, reference = refname, method = names(na_note),
                                                n_pos = sum(pos), n_neg = sum(!pos), n_na = NA_real_, L = LARG, ref_up = n_up, ref_down = n_dn,
                                                note = unname(na_note)), fill = TRUE)
  for (m in intersect(names(m_note), r$method)) r[method == m, note := fifelse(note == "", m_note[[m]], paste(note, m_note[[m]], sep = "; "))]
  rr <- r[r$method == "RFind-sc", ]
  if (!TEST && tid %in% CHANGED) CHK <<- c(CHK, sprintf("  %s RFind-sc: new or changed task in the redo (no stored result to reproduce); AUC %.6f [%.6f, %.6f], %d / %d cells",
                                                       tid, rr$AUC, rr$lo, rr$hi, rr$n_pos, rr$n_neg))
  if (!TEST && all(keep) && !tid %in% CHANGED) {                  # data loading must reproduce the stored AUC and cell numbers (the stored
    ref0 <- if (FAIR) SW[SW$id == tid & SW$L == LARG, ] else BM[BM$id == tid & BM$method == "RFind-sc", ]   # CIs come from a bootstrap over cells: not compared)
    tol <- if (FAIR) 1e-6 else 1e-9
    d <- if (nrow(ref0) == 1L) max(abs(c(rr$AUC - ref0$AUC, rr$n_pos - ref0$n_pos, rr$n_neg - ref0$n_neg))) else Inf
    CHK <<- c(CHK, sprintf("  %s RFind-sc AUC %.6f vs %s: max |diff| %.2g -> %s", tid, rr$AUC,
                           if (FAIR) sprintf("results/fig5_reflength_sweep_auc.csv (L = %s)", LARG) else "results/fig5_final_auc.csv", d, if (d < tol) "reproduced" else "NOT reproduced"))
    if (!(d < tol)) stop(sprintf("%s: RFind-sc does not reproduce the stored AUC (max |diff| %.3g)", tid, d), call. = FALSE)
  }
  RES[[tid]] <<- r
  lg(sprintf("%s (%d cells, %d positive; reference %d / %d): %s", tid, length(pos), sum(pos), n_up, n_dn,
             paste(sprintf("%s %.3f", r$method, r$AUC), collapse = ", ")))
  if (!TEST) fwrite(rbindlist(RES, fill = TRUE), f_out("fig5_redo_auc"))             # partial table after every task
  writeLines(LOG, TXT)
  invisible(r)
}

## ============================================ task definitions (fig5_final.R) ============================================
sig <- readRDS("data/Hammond2019/output/percell_signature_top500.rds"); cells <- sig$cells; univ <- sig$genes; rm(sig); invisible(gc())
meta <- fread("data/Hammond2019/raw/cell_meta_hammond.csv")[match(cells, cell_id)]
rf   <- fread("data/Hammond2019/output/scores_hammond_9axes.csv")[match(cells, cell_id)]
panel <- readRDS("data/Hammond2019/output/references/references_panel9.rds")
FCT <- list(DAM = fread("data/Hammond2019/output/references/dam_KerenShaul.csv"), Maturity = fread("data/Hammond2019/output/references/maturity_Matcovitch_microglia.csv"))
mg <- meta$paper_cluster != "Mono/Mac"
sets <- function(key) list(up = intersect(panel[[key]]$up, univ), dn = intersect(panel[[key]]$down, univ))
TH <- list(
  H1 = list(key = "DAM", rfc = "DAM", ref = "DAM", kind = "condition", task = "Hammond: DAM, LPC vs saline white matter (microglia)",
            ix = which(mg & meta$age %in% c("P100_Injured_WM", "P100_Control_WM")), pos = function(ix) meta$age[ix] == "P100_Injured_WM"),
  H2 = list(key = "Interferon", rfc = "Interferon", ref = "Interferon", kind = "cluster", task = "Hammond: Interferon, cluster 9 vs other microglia",
            ix = which(mg), pos = function(ix) meta$paper_cluster[ix] == "9"),
  H3 = list(key = "Interferon", rfc = "Interferon", ref = "Interferon", kind = "condition", task = "Hammond: Interferon, P540 vs P100",
            ix = which(mg & meta$age %in% c("P540", "P100")), pos = function(ix) meta$age[ix] == "P540"),
  H4 = list(key = "Maturity", rfc = "Maturity", ref = "Maturity", kind = "condition", task = "Hammond: Maturity, P100 vs E14 + P4/P5",
            ix = which(mg & meta$age %in% c("P100", "E14", "P4/P5")), pos = function(ix) meta$age[ix] == "P100"),
  H5 = list(key = "G2M", rfc = "G2M", ref = "G2M", kind = "cluster", task = "Hammond: G2M, cluster 2b vs other microglia",
            ix = which(mg), pos = function(ix) meta$paper_cluster[ix] == "2b"))
if (FAIR) TH <- TH[c("H1", "H4")]                                 # gene sets have no length to vary
for (id in names(TH)) TH[[id]]$ix <- pick(TH[[id]]$ix, TH[[id]]$pos)
dm <- fread("data/DeDonno2026/output/cell_meta_qc.csv")
ds <- fread("data/DeDonno2026/output/scores_dedonno.csv")[match(dm$cell_id, cell_id)]
ref <- readRDS("data/DeDonno2026/output/references_dedonno.rds")$CORT
DT <- list(D1 = "Oligodendrocytes", D2 = "Astro", D3 = "Endothelial")   # the three cell types of the Glucocorticoid tasks
TIX <- lapply(DT, function(ct) pick(which(dm$cell_type == ct & dm$line == "Ctrl"), function(i) dm$condition[i] == "Stress"))
IDS <- c(names(TH), names(DT))

## ========================= fairness mode: RFind-sc at L (fig5_reflength_sweep.R; before any method package is loaded) =========================
RFL <- list()
if (FAIR) {
  lg(sprintf("RFind-sc recomputed at L = %s as in fig5_reflength_sweep.R (cell lists 500; %d cores)", LARG, NCORES))
  dp <- readRDS("data/Hammond2019/output/percell_signature_deep.rds"); stopifnot(identical(dp$cells, cells))
  deep <- list(up = dp$up, down = dp$down, n_deep = dp$n_deep, col = seq_along(cells)); rm(dp); invisible(gc())
  cf0 <- readRDS("data/Hammond2019/output/counts_filtered.rds"); basis <- readRDS("data/Hammond2019/output/delta_P2_basis.rds")
  stopifnot(identical(colnames(cf0), cells), identical(rownames(cf0), basis$genes), identical(names(basis$cs), cells))
  share <- as.numeric(basis$share); csH <- basis$cs; rm(basis)
  expH <- function(idx) outer(share, csH[idx])
  KEYH0 <- vapply(TH, `[[`, "", "key")
  for (key in unique(KEYH0)) {                                    # once per reference, on the union of its task cells
    ids <- names(TH)[KEYH0 == key]; uix <- sort(unique(unlist(lapply(TH[ids], `[[`, "ix")))); e <- panel[[key]]
    check500(key, lists_hammond(FCT[[key]], 500L), e)
    RL <- make_ref_L(e, lists_hammond(FCT[[key]], LVAL), rownames(cf0))
    sw <- sweep_cells(uix, cf0, expH, deep, tolower(rownames(cf0)), list(RL))
    if (sw$mis != 0L) stop(sprintf("%s: recomputed ranking differs from the stored deep ranking in %d lists", key, sw$mis), call. = FALSE)
    for (id in ids) RFL[[id]] <- sw$sc[match(TH[[id]]$ix, uix), 1]
    say(sprintf("RFind-sc at L = %s: %s done (%s; %d cells)", LARG, key, paste(ids, collapse = ", "), length(uix)))
  }
  rm(deep, cf0, share, csH, expH, sw); invisible(gc())
  X <- readRDS("data/DeDonno2026/output/counts_universe.rds"); stopifnot(identical(X$cells, dm$cell_id))
  cu0 <- X$counts; csD <- X$cs; rm(X); invisible(gc())
  keepc <- sort(unique(unlist(TIX)))
  dp <- readRDS("data/DeDonno2026/output/percell_signature_celltype_deep.rds")
  stopifnot(identical(dp$cells, dm$cell_id), identical(dp$baseline, "per_celltype_all_cells"))
  col <- rep(NA_integer_, nrow(dm)); col[keepc] <- seq_along(keepc)
  deep <- list(up = dp$up[, keepc, drop = FALSE], down = dp$down[, keepc, drop = FALSE], n_deep = dp$n_deep, col = col); rm(dp); invisible(gc())
  check500("Glucocorticoid (CORT)", lists_cort(500L), ref)
  RL <- make_ref_L(ref, lists_cort(LVAL), rownames(cu0)); glow0 <- tolower(rownames(cu0))
  for (id in names(DT)) {
    bix <- which(dm$cell_type == DT[[id]])
    base <- Matrix::rowMeans(cu0[, bix, drop = FALSE] %*% Diagonal(x = 1 / csD[bix]))          # per-cell-type baseline share (02_signature.R)
    expD <- function(idx) outer(base, csD[idx])
    sw <- sweep_cells(TIX[[id]], cu0, expD, deep, glow0, list(RL))
    if (sw$mis != 0L) stop(sprintf("%s: recomputed ranking differs from the stored deep ranking in %d lists", id, sw$mis), call. = FALSE)
    RFL[[id]] <- sw$sc[, 1]; say(sprintf("RFind-sc at L = %s: %s done (%d cells)", LARG, id, length(TIX[[id]])))
  }
  rm(deep, cu0, csD, base, expD, sw); invisible(gc())
}

## ================================ Hammond ================================
libsize <- Matrix::colSums(readRDS("data/Hammond2019/output/counts_qc.rds")); invisible(gc())
cf <- readRDS("data/Hammond2019/output/counts_filtered.rds")
stopifnot(identical(colnames(cf), cells), identical(names(libsize), cells), identical(rownames(cf), univ))
ln <- cf %*% Diagonal(x = 1e4 / libsize); ln@x <- log1p(ln@x); dimnames(ln) <- dimnames(cf)
say(sprintf("Hammond: %d genes x %d cells", nrow(cf), ncol(cf)))
ref_lists_H <- function(key) {                                   # the genes every method gets (= fig5_final.R at L = 500)
  s <- sets(key)
  if (!length(s$dn)) return(s)                                    # gene set: membership, fixed
  e <- panel[[key]]; U <- pair_universe(univ, e$eligible)
  l500 <- lists_hammond(FCT[[key]], 500L); check500(key, l500, e)
  if (!identical(list(up = map_genes(l500$up, univ, U), dn = map_genes(l500$down, univ, U)), s))
    stop(key, ": the L = 500 rule does not give the genes of fig5_final.R", call. = FALSE)
  if (!FAIR) return(s)
  l <- lists_hammond(FCT[[key]], LVAL); list(up = map_genes(l$up, univ, U), dn = map_genes(l$down, univ, U))
}
REFH <- lapply(setNames(nm = unique(vapply(TH, `[[`, "", "key"))), ref_lists_H)
for (id in names(TH)) if (length(REFH[[TH[[id]]$key]]$dn)) chk_len(id, REFH[[TH[[id]]$key]])
lg(sprintf("Hammond references (up / down given to every method): %s", paste(sprintf("%s %d / %d", names(REFH), lengths(lapply(REFH, `[[`, "up")),
                                                                                         lengths(lapply(REFH, `[[`, "dn"))), collapse = "; ")))
if ("VISION" %in% SEL) {                                          # before the AUCell rankings (memory)
  VH <- vision_run(cf, libsize, ln, REFH); vision_log("Hammond (all QC cells)", VH, REFH); say("Hammond: VISION done")
}
## AUCell ranks tied genes at random (colRanks(ties.method = "random") in AUCell_buildRankings): the seed is fixed so that the
## scores are reproducible (without the seed, repeated runs differ by up to 0.06 per cell and 0.0008 in AUC)
if ("AUCell" %in% SEL) { set.seed(1); rkH <- AUCell::AUCell_buildRankings(cf, plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE); say("Hammond: rankings built") }
KEYH <- vapply(TH, `[[`, "", "key")
for (key in unique(KEYH)) {                                       # each reference is scored once, on the union of the cells of its tasks
  ids <- names(TH)[KEYH == key]; s <- REFH[[key]]; signed <- length(s$dn) > 0
  if (!signed && all(SEL %in% SIGNED_ONLY)) { lg(sprintf("%s: gene set; no selected method applies; skipped", paste(ids, collapse = ", "))); next }
  ix <- sort(unique(unlist(lapply(TH[ids], `[[`, "ix"))))
  ckw <- CHK_LOAD && any(!ids %in% CHANGED)                       # the data-loading check is needed for a task not in CHANGED only
  sc <- list(); na_note <- character(0); ck <- list(); m_note <- character(0)
  if (signed && "UCellUp" %in% SEL) {                            # UCell with the up genes only (fig5_ucell_up.R); a gene set = UCell
    u <- ucell_scores(cf[, ix], s$up, character(0))
    if (is.null(u)) na_note["UCell (up genes)"] <- sprintf("not applicable: %d signature genes > UCell maxRank 1500 (UCELL_LONG = skip)", length(s$up)) else {
      sc[["UCell (up genes)"]] <- u; if (!is.null(attr(u, "note"))) m_note["UCell (up genes)"] <- attr(u, "note") }
  }
  if ("UCell" %in% SEL) {
    u <- ucell_scores(cf[, ix], s$up, s$dn)
    if (is.null(u)) na_note["UCell"] <- sprintf("not applicable: %d signature genes > UCell maxRank 1500 (UCELL_LONG = skip)", length(s$up) + length(s$dn)) else {
      sc$UCell <- u; if (!is.null(attr(u, "note"))) m_note["UCell"] <- attr(u, "note") }
  }
  if ("AUCell" %in% SEL) sc$AUCell <- aucell_scores(rkH, s$up, cells[ix])
  if ("AddModuleScore" %in% SEL || ckw) {                         # always at L = 500: data-loading check (cheap, seeded)
    a <- ams_scores(ln, s$up, ix); if ("AddModuleScore" %in% SEL) sc$AddModuleScore <- a; if (ckw) ck$AddModuleScore <- a
  }
  if (any(c("ULM", "WMEAN", "singscore") %in% SEL)) { d <- dc_ss(cf[, ix], libsize[ix], s$up, s$dn); for (m in intersect(c("ULM", "WMEAN", "singscore"), SEL)) sc[[m]] <- d[, m] }
  if (signed && any(names(FCKEYS) %in% SEL)) {                    # decoupleR with the log2FC weights
    g <- c(s$up, s$dn); d <- dc_fc(cf[, ix], libsize[ix], g, fcw(FCT[[key]], g)); for (m in FCKEYS[intersect(names(FCKEYS), SEL)]) sc[[m]] <- d[, m]
  }
  if ("ORA" %in% SEL) sc$ORA <- ora_scores(cf[, ix], libsize[ix], s$up, key)   # up genes of a signed reference; a gene set as it is
  if (signed && "scFOCAL" %in% SEL) {                              # scFOCAL-type, centred over the atlas
    g <- c(s$up, s$dn); cc <- cos_focal(ln, g, fcw(FCT[[key]], g), ix, focal = TRUE)
    sc[["scFOCAL-type"]] <- cc[, 3]
  }
  if ("VISION" %in% SEL) sc$VISION <- VH$scores[cells[ix], key]
  if ("Beyondcell" %in% SEL) { bb <- bc_scores(ln, s$up, s$dn, ix); sc$Beyondcell <- bb$score }
  if ("ssGSEA" %in% SEL) sc$ssGSEA <- ssgsea_scores(cf[, ix], libsize[ix], s$up, s$dn, key)
  say(sprintf("%s scored on %d cells (tasks %s)", key, length(ix), paste(ids, collapse = ", ")))
  for (id in ids) {
    Tk <- TH[[id]]; k <- match(Tk$ix, ix); stopifnot(!anyNA(k))
    rfv <- if (FAIR) RFL[[id]] else rf[[Tk$rfc]][Tk$ix]
    finish_task(id, Tk$task, Tk$kind, Tk$ref, cells[Tk$ix], rfv, Tk$pos(Tk$ix), lapply(sc, function(v) v[k]), na_note, length(s$up), length(s$dn),
                ck = lapply(ck, function(v) v[k]), m_note = m_note, grp = meta$age[Tk$ix])
    say(sprintf("%s done (%d cells)", id, length(Tk$ix)))
  }
  rm(sc, ck); invisible(gc())
}
rm(cf, ln); if (exists("rkH")) rm(rkH); if (exists("VH")) rm(VH); invisible(gc())

## ================================ De Donno ================================
X <- readRDS("data/DeDonno2026/output/counts_universe.rds"); cu <- X$counts; colnames(cu) <- X$cells; lib <- X$cs[match(dm$cell_id, X$cells)]; rm(X); invisible(gc())
s0 <- list(up = intersect(ref$up, rownames(cu)), dn = intersect(ref$down, rownames(cu)))          # fig5_final.R
U <- pair_universe(rownames(cu), ref$eligible)
l500 <- lists_cort(500L); check500("Glucocorticoid (CORT)", l500, ref)
if (!identical(list(up = map_genes(l500$up, rownames(cu), U), dn = map_genes(l500$down, rownames(cu), U)), s0))
  stop("Glucocorticoid: the L = 500 rule does not give the genes of fig5_final.R", call. = FALSE)
SD <- if (FAIR) { l <- lists_cort(LVAL); list(up = map_genes(l$up, rownames(cu), U), dn = map_genes(l$down, rownames(cu), U)) } else s0
up <- SD$up; dn <- SD$dn; for (id in names(DT)) chk_len(id, SD)
grt <- fread("data/DeDonno2026/output/DeDonno_GR_glucocorticoid.csv"); wD <- fcw(grt, c(up, dn))
dupg <- unique(rownames(cu)[duplicated(rownames(cu))])
lg(sprintf("De Donno: %d genes x %d cells; Glucocorticoid reference (given to every method) up %d / down %d; duplicated gene names %s (reference genes among them: %s)",
           nrow(cu), ncol(cu), length(up), length(dn), paste(dupg, collapse = ", "), if (length(intersect(dupg, c(up, dn)))) paste(intersect(dupg, c(up, dn)), collapse = ", ") else "none"))
for (id in names(DT)) {
  bix <- which(dm$cell_type == DT[[id]]); tix <- TIX[[id]]; kk <- match(tix, bix); tc <- dm$cell_id[tix]   # kk = task cells within the cell type
  stopifnot(!anyNA(kk))
  lnb <- cu[, dm$cell_id[bix]] %*% Diagonal(x = 1e4 / lib[bix]); lnb@x <- log1p(lnb@x); dimnames(lnb) <- list(rownames(cu), dm$cell_id[bix])
  sc <- list(); na_note <- character(0); ck <- list(); m_note <- character(0)
  if ("UCellUp" %in% SEL) {                                       # UCell with the up genes only (fig5_ucell_up.R)
    u <- ucell_scores(cu[, tc], up, character(0))
    if (is.null(u)) na_note["UCell (up genes)"] <- sprintf("not applicable: %d signature genes > UCell maxRank 1500 (UCELL_LONG = skip)", length(up)) else {
      sc[["UCell (up genes)"]] <- u; if (!is.null(attr(u, "note"))) m_note["UCell (up genes)"] <- attr(u, "note") }
  }
  if ("UCell" %in% SEL) {
    u <- ucell_scores(cu[, tc], up, dn)
    if (is.null(u)) na_note["UCell"] <- sprintf("not applicable: %d signature genes > UCell maxRank 1500 (UCELL_LONG = skip)", length(up) + length(dn)) else {
      sc$UCell <- u; if (!is.null(attr(u, "note"))) m_note["UCell"] <- attr(u, "note") }
  }
  if ("AUCell" %in% SEL) { set.seed(1); rk <- AUCell::AUCell_buildRankings(cu[, tc], plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE); sc$AUCell <- aucell_scores(rk, up, tc); rm(rk) }
  if ("AddModuleScore" %in% SEL || CHK_LOAD) {                    # bins over the cell type; always at L = 500 (check)
    a <- ams_scores(lnb, up, kk); if ("AddModuleScore" %in% SEL) sc$AddModuleScore <- a; if (CHK_LOAD) ck$AddModuleScore <- a
  }
  if (any(c("ULM", "WMEAN", "singscore") %in% SEL)) { d <- dc_ss(cu[, tc], lib[tix], up, dn); for (m in intersect(c("ULM", "WMEAN", "singscore"), SEL)) sc[[m]] <- d[, m] }
  if (any(names(FCKEYS) %in% SEL)) { d <- dc_fc(cu[, tc], lib[tix], c(up, dn), wD); for (m in FCKEYS[intersect(names(FCKEYS), SEL)]) sc[[m]] <- d[, m] }
  if ("ORA" %in% SEL) sc$ORA <- ora_scores(cu[, tc], lib[tix], up, id)          # the up genes of the Glucocorticoid reference
  if ("scFOCAL" %in% SEL) {                                         # scFOCAL-type, centred over the cell type
    cc <- cos_focal(lnb, c(up, dn), wD, kk, focal = TRUE)
    sc[["scFOCAL-type"]] <- cc[, 3]
  }
  if ("VISION" %in% SEL) {                                        # one object per cell type (both conditions; wild-type mice)
    sg <- list(Glucocorticoid = list(up = up, dn = dn))
    VD <- vision_run(cu[, dm$cell_id[bix]], lib[bix], lnb, sg); vision_log(sprintf("%s (%s, %d cells)", id, DT[[id]], length(bix)), VD, sg)
    sc$VISION <- VD$scores[tc, "Glucocorticoid"]; rm(VD)
  }
  if ("Beyondcell" %in% SEL) { bb <- bc_scores(lnb, up, dn, kk); sc$Beyondcell <- bb$score }
  if ("ssGSEA" %in% SEL) sc$ssGSEA <- ssgsea_scores(cu[, tc], lib[tix], up, dn, id)
  rfv <- if (FAIR) RFL[[id]] else ds$CORT[tix]
  finish_task(id, sprintf("De Donno: Glucocorticoid, Stress vs No Stress, %s (control line)", DT[[id]]), "condition", "Glucocorticoid",
              tc, rfv, dm$condition[tix] == "Stress", sc, na_note, length(up), length(dn), ck = ck, m_note = m_note, grp = dm$condition[tix])
  say(sprintf("%s %s: %d cells", id, DT[[id]], length(tix))); rm(lnb); invisible(gc())
}

## ================================ tables and summary ================================
res <- rbindlist(RES, fill = TRUE)
setcolorder(res, intersect(c(names(BM), "n_na", "L", "ref_up", "ref_down", "note"), names(res)))
res <- res[order(match(id, IDS), match(method, c("RFind-sc", METH_ALL)))]
fwrite(res, f_out("fig5_redo_auc"))
fwrite(rbindlist(PC, fill = TRUE), f_out("fig5_redo_percell", "csv.gz"))
if (length(BCN) && !FAIR) fwrite(rbindlist(BCN, fill = TRUE), f_out("fig5_redo_beyondcell"))   # L = 500 only (the tables and Source Data read it)

lgs("", "Checks")
if (length(CHK)) lgs(CHK) else lg("  RFind-sc reproduction: not checked (test run)")
if (length(SSG_CHK)) lgs(SSG_CHK)
if (length(ORA_CHK)) lgs(ORA_CHK)
if (length(UCELL_NOTE)) lg("  UCell: ", paste(unique(UCELL_NOTE), collapse = "; "))
if (!TEST) {                                                      # every method vs the stored table results/fig5_final_auc.csv (L = 500), tasks
                                                                  # not in CHANGED only (log only); no such table is stored for the other lengths
  v1f <- if (LARG == "500") "figs/Fig5_benchmark/output/results/fig5_final_auc.csv" else character(0)
  v1f <- v1f[file.exists(v1f)]
  if (length(v1f)) {
    V1 <- unique(rbindlist(lapply(v1f, function(f) fread(f)[, .(id, method, AUC0 = AUC, lo0 = lo, hi0 = hi)])), by = c("id", "method"))
    cmp <- merge(res[!is.na(AUC) & !id %in% CHANGED, .(id, method, AUC, lo, hi)], V1, by = c("id", "method"))
    if (nrow(cmp)) {
      cmp[, dd := pmax(abs(AUC - AUC0), abs(lo - lo0), abs(hi - hi0))]
      lg(sprintf("  vs the v1 tables of L = %s (%s; unchanged tasks):", LARG, paste(basename(v1f), collapse = ", ")))
      for (m in intersect(c("RFind-sc", METH_ALL), unique(cmp$method)))
        lg(sprintf("    %-22s max |diff| of AUC / CI over %d tasks %.2g -> %s", m, sum(cmp$method == m), max(cmp[method == m]$dd),
                   if (max(cmp[method == m]$dd) < 1e-6) "reproduced" else "DIFFERS"))
    }
  } else lg(sprintf("  vs v1: no v1 table at L = %s", LARG))
}
if (length(BCN)) {
  lgs("", sprintf("Beyondcell NaN (expr.thres %g) and the AUC with the NaN cells set to 0 (main table: %s):", BC_EXPR_THRES, BEYONDCELL_NAN))
  for (x in BCN) lg(sprintf("  %-3s NaN %6d (%5.2f%%)  AUC %.3f", x$id[1], x$n_nan[1], 100 * x$frac_nan[1], if (is.null(x$AUC)) NA_real_ else x$AUC[1]))
}
ab <- c("RFind-sc" = "RFind", "UCell (up genes)" = "UCellUp", UCell = "UCell", AUCell = "AUCell", AddModuleScore = "AMS", ULM = "ULM", WMEAN = "WMEAN",
        "ULM (log2FC)" = "ULMfc", "WMEAN (log2FC)" = "WMEANfc",
        ORA = "ORA", singscore = "singsc", VISION = "VISION", Beyondcell = "Beyond", ssGSEA = "ssGSEA", "scFOCAL-type" = "scFOCAL")
ms <- c("RFind-sc", METH_ALL)[c("RFind-sc", METH_ALL) %in% res$method]
lgs("", sprintf("AUC (* = RFind-sc higher, x = RFind-sc lower, 95%% bootstrap CI of the difference excludes 0; - = not scored) at L = %s:", LARG),
   sprintf("%-4s %s", "id", paste(sprintf("%-9s", ab[ms]), collapse = "")))
for (i in IDS) {
  r <- res[res$id == i]; if (!nrow(r)) next
  v <- vapply(ms, function(m) { x <- r[r$method == m]
    if (!nrow(x) || is.na(x$AUC)) sprintf("%-9s", "-") else sprintf("%-9s", paste0(sprintf("%.3f", x$AUC), if (x$method == "RFind-sc") "" else if (isTRUE(x$d_lo > 0)) "*" else if (isTRUE(x$d_hi < 0)) "x" else "")) }, "")
  lg(sprintf("%-4s %s", i, paste(v, collapse = "")))
}
lgs("", sprintf("Output: %s, %s%s", f_out("fig5_redo_auc"), f_out("fig5_redo_percell", "csv.gz"), if (length(BCN) && !FAIR) paste0(", ", f_out("fig5_redo_beyondcell")) else ""),
   sprintf("Run time %.1f min.", as.numeric(Sys.time() - t0, units = "mins")))
writeLines(LOG, TXT); say("done; summary in ", TXT)
