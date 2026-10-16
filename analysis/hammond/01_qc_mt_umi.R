## Hammond step 1 — QC: mtDNA% and total UMI (MAD-based)
## Input: cell metadata (cell_meta_hammond.csv) + mt_frac_hammond.csv (per-cell mt% from the raw GEO DGE files)
##   (the mt genes were dropped from the annotated h5ad when it was converted, so mt% is taken from the recomputed values)
## Output: data/Hammond2019/output/qc_flags_mt_umi.csv (cell_id + metrics + flags)
suppressPackageStartupMessages(library(data.table))
OUT <- "data/Hammond2019/output"
NMAD <- 5
MT_CAP <- 5   # upper limit (%) of the mtDNA% threshold, the same rule as for the data of De Donno et al. (analysis/dedonno/01_qc_embed.R)

meta <- fread("data/Hammond2019/raw/cell_meta_hammond.csv")
mt   <- fread("data/Hammond2019/raw/mt_frac_hammond.csv")[, .(cell_id, mt_frac)]
d <- merge(meta[, .(cell_id, cell_type_paper, nCount_RNA, nFeature_RNA)], mt, by="cell_id", all.x=TRUE)
stopifnot(nrow(d) == 76149, !any(is.na(d$mt_frac)))

## --- helper: MAD outlier thresholds (mad() uses constant = 1.4826, i.e. SD scale; the usual nmads convention) ---
mad_thr <- function(x, nmad=NMAD) { m <- median(x); s <- mad(x); c(lo=m-nmad*s, hi=m+nmad*s) }

## (1) mtDNA% — upper side only (high mt = dying cells); threshold = the lower of the 5-MAD upper limit and MT_CAP
mt_pct <- d$mt_frac * 100
mt_t <- mad_thr(mt_pct)
mt_hi <- min(mt_t[["hi"]], MT_CAP)
d$flag_mt <- mt_pct > mt_hi

## (2) total UMI — MAD of log10(UMI); low, high and both sides reported
luec <- log10(d$nCount_RNA)
umi_t <- mad_thr(luec)
d$flag_umi_lo <- luec < umi_t["lo"]   # low UMI = low quality / dying cells
d$flag_umi_hi <- luec > umi_t["hi"]   # high UMI = doublet candidates (overlaps with Scrublet)

cat(sprintf("=== QC: mtDNA%% (upper %d-MAD) ===\n", NMAD))
cat(sprintf("  median=%.2f%%  MAD=%.2f%%  threshold(hi)=%.2f%% (5-MAD %.2f%%, cap %.0f%%)\n",
            median(mt_pct), mad(mt_pct), mt_hi, mt_t["hi"], MT_CAP))
cat(sprintf("  flagged (mt > %.2f%%): %d / %d (%.2f%%)\n\n",
            mt_hi, sum(d$flag_mt), nrow(d), 100*mean(d$flag_mt)))

cat(sprintf("=== QC: log10(total UMI) (%d-MAD) ===\n", NMAD))
cat(sprintf("  median(log10)=%.3f  MAD=%.3f  -> UMI lo=%.0f  hi=%.0f\n",
            median(luec), mad(luec), 10^umi_t["lo"], 10^umi_t["hi"]))
cat(sprintf("  flag LOW  (UMI < %.0f): %d (%.2f%%)\n", 10^umi_t["lo"], sum(d$flag_umi_lo), 100*mean(d$flag_umi_lo)))
cat(sprintf("  flag HIGH (UMI > %.0f): %d (%.2f%%)\n", 10^umi_t["hi"], sum(d$flag_umi_hi), 100*mean(d$flag_umi_hi)))
cat(sprintf("  union of both sides: %d (%.2f%%)\n\n", sum(d$flag_umi_lo|d$flag_umi_hi), 100*mean(d$flag_umi_lo|d$flag_umi_hi)))

## summary (mt + UMI, both sides): cells passing at this stage
d$pass_mt_umi <- !(d$flag_mt | d$flag_umi_lo | d$flag_umi_hi)
cat(sprintf("=== pass after mt + UMI (both sides): %d / %d (%.2f%% kept) ===\n",
            sum(d$pass_mt_umi), nrow(d), 100*mean(d$pass_mt_umi)))

fwrite(d, file.path(OUT, "qc_flags_mt_umi.csv"))
cat("\nwrote", file.path(OUT, "qc_flags_mt_umi.csv"), "\n")
