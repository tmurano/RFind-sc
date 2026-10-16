## Hammond: RF scores on the 4 axes DAM, Maturity, Interferon and G2M (references_panel9.rds x per-cell P2 signatures).
## The file names (references_panel9.rds, scores_hammond_9axes.*) are kept for compatibility with the scripts that read them.
## score = 4-term Running Fisher (rfind_score_cell; smallest p over every G x K cut-off pair), clip ±300.
## Background: universe of each target-reference pair = Hammond gene universe ∩ genes the reference could
## contain (analysis/running_fisher.R: pair_universe); per-cell lists are restricted to it before the top 500 are taken (cut_in_universe).
## References: ranked up and down lists (type A of the paper: DAM, Maturity; stored with type = "B" in the panel) and gene sets
## (type B of the paper: Interferon, G2M; stored with type = "C"; up = the set, down = empty).
## Output: data/Hammond2019/output/scores_hammond_9axes.csv (cell_id, cluster, 4 axes) + .rds
suppressPackageStartupMessages({ library(data.table); library(parallel) })
OUTP <- "data/Hammond2019/output"; REFP <- "data/Hammond2019/output/references"   # paths relative to the project root
source("analysis/running_fisher.R")
NCORES <- max(1L, detectCores()-2L); CLIP <- 300; NTOP <- 500L

sig <- readRDS(file.path(OUTP,"percell_signature_deep.rds"))
genes <- sig$genes; glow <- tolower(genes); NDEEP <- sig$n_deep
N <- ncol(sig$up)
cat(sprintf("P2 signatures (deep %d): target universe %d genes × %d cells\n", NDEEP, length(unique(glow)), N))

panel <- readRDS(file.path(REFP,"references_panel9.rds"))
axes <- names(panel)
## gene sets (type = "C") have no order: ranked = FALSE, so the set is used whole
UA  <- lapply(panel, function(e) pair_universe(genes, e$eligible))           # pair universe per reference
INU <- lapply(UA, function(u) glow %in% u)
P   <- vapply(UA, length, integer(1))
b2 <- lapply(axes, function(a) { e <- panel[[a]]
  list(up=intersect(tolower(e$up),UA[[a]]), down=intersect(tolower(e$down),UA[[a]]), ranked=!identical(e$type,"C")) })
names(b2) <- axes
cat(sprintf("Reference axes: %s\n", paste(axes,collapse=", ")))
for(a in axes) cat(sprintf("  [%-11s] type%s pair universe P=%d up=%d down=%d ranked=%s\n", a, panel[[a]]$type, P[[a]], length(b2[[a]]$up), length(b2[[a]]$down), b2[[a]]$ranked))

up <- sig$up; dn <- sig$down
score_chunk <- function(idx){
  m <- matrix(0, length(idx), length(axes), dimnames=list(NULL,axes)); short <- 0L
  for(r in seq_along(idx)){
    j <- idx[r]
    for(a in axes){
      uj <- cut_in_universe(up[,j], INU[[a]], NTOP); dj <- cut_in_universe(dn[,j], INU[[a]], NTOP)
      ## shorter than NTOP is fine when the cell has no more genes of that sign; a truncation only if the deep ranking ran out
      short <- short + (length(uj) < NTOP && sig$n_pos[j] > NDEEP) + (length(dj) < NTOP && sig$n_neg[j] > NDEEP)
      m[r,a] <- rfind_score_cell(glow[uj],glow[dj],b2[[a]],universe_size=P[[a]],clip=CLIP)
    }
  }
  attr(m,"short") <- short; m
}
chunks <- split(seq_len(N), ceiling(seq_len(N)/250))
cat(sprintf("scoring %d cells × all axes (%d cores)...\n", N, NCORES))
t0 <- Sys.time()
res <- mclapply(chunks, score_chunk, mc.cores=NCORES)
n_short <- sum(vapply(res, function(x) attr(x,"short"), integer(1)))
cat(sprintf("lists cut short by the deep ranking (NDEEP=%d): %d  (must be 0)\n", NDEEP, n_short)); stopifnot(n_short == 0L)
M <- do.call(rbind,res); rownames(M) <- sig$cells
cat(sprintf("done %.0fs\n", as.numeric(Sys.time()-t0,units="secs")))

meta <- fread("data/Hammond2019/raw/cell_meta_hammond.csv"); setkey(meta,cell_id)
meta <- meta[sig$cells, .(cell_id, cluster=cell_type_paper)]
out <- cbind(meta, as.data.table(M))
fwrite(out, file.path(OUTP,"scores_hammond_9axes.csv"))
saveRDS(list(scores=M, axes=axes, cells=sig$cells, universe_size=P), file.path(OUTP,"scores_hammond_9axes.rds"))

## median per cluster (sanity check)
ord <- c("1","2a","2b","2c","3","4","5","6","7a","7b","7c","8","9","Mono/Mac")
med <- out[,lapply(.SD,function(x)round(median(x),1)),by=cluster,.SDcols=axes][match(ord,cluster)]
cat("\n=== median per cluster ===\n"); print(med)
cat("\nwrote scores_hammond_9axes.csv + .rds\n")
