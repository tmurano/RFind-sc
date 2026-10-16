## Writes the counts used by fig5_runtime_score.R in a form the Python run time script reads without R (fig5_runtime.sh):
## the 75,572 cells of the data of Hammond et al. after quality control (counts_qc.rds), as the three arrays of the sparse matrix
## (genes x cells; x as float64, i and p as int32), the gene and cell names and the random cell order of fig5_runtime_score.R (seed 1).
## Output folder: $RFIND_RUNTIME_EXPORT (a temporary folder made by fig5_runtime.sh; nothing is written to the project).
## Run (from the project root): RFIND_RUNTIME_EXPORT=<folder> R --no-save --no-restore -f figs/Fig5_benchmark/analysis/fig5_runtime_export.R
suppressPackageStartupMessages(library(Matrix))
D <- Sys.getenv("RFIND_RUNTIME_EXPORT"); stopifnot(nzchar(D), dir.exists(D))
X <- readRDS("data/Hammond2019/output/counts_qc.rds")
set.seed(1); ord <- sample(ncol(X))                                  # the same order as fig5_runtime_score.R
writeBin(as.double(X@x), file.path(D, "x.bin")); writeBin(as.integer(X@i), file.path(D, "i.bin")); writeBin(as.integer(X@p), file.path(D, "p.bin"))
writeLines(rownames(X), file.path(D, "genes.txt")); writeLines(colnames(X), file.path(D, "cells.txt"))
writeLines(as.character(ord), file.path(D, "order.txt")); writeLines(as.character(dim(X)), file.path(D, "dims.txt"))
cat(sprintf("wrote %d genes x %d cells to %s\n", nrow(X), ncol(X), D))
