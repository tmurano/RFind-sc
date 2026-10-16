#!/bin/bash
# Packages of the Fig. 5 benchmark (2026-09-28): VISION v3.0.2, GSVA (ssGSEA), DMwR 0.4.1 and Beyondcell v2.2.2 for the CRAN
# framework build of R 4.5.2 (/Library/Frameworks/R.framework/Versions/4.5-arm64), into its own library, where the other methods of
# the benchmark (Seurat, UCell, AUCell, singscore, decoupleR) are already installed. The framework's 'Current' points to 4.6, and
# the 4.5 bin/R wrapper and etc/Makeconf hard-code /Library/Frameworks/R.framework/Resources (= Current): a source build started
# from 4.5.2 would run R CMD INSTALL under 4.6, with the 4.6 headers, and link to the 4.6 libR. So every install runs from a private
# mirror of the 4.5 R_HOME (symlinks to the real directories; bin/R and etc/Makeconf are patched copies) in a path without spaces,
# in a clean environment (no conda). Nothing under /Library is changed except the new packages. CRAN / Bioconductor packages are
# macOS arm64 binaries for R 4.5 when available; packages already installed are never upgraded (update = FALSE, upgrade = "never"),
# so the methods already installed keep their versions. VISION and Beyondcell come from the GitHub tags, DMwR (needed by Beyondcell)
# from the CRAN archive.
# Usage (project root): bash figs/Fig5_benchmark/analysis/install_r452_bench_work.sh   [R452_MIRROR=<dir without spaces>]
set -euo pipefail
R45=/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources
M="${R452_MIRROR:-${TMPDIR:-/tmp}/r452home}"
case "$M" in *" "*) echo "R452_MIRROR must not contain spaces: $M" >&2; exit 1 ;; esac

# 1. Mirror of the 4.5 R_HOME
rm -rf "$M"; mkdir -p "$M/bin" "$M/etc"
for f in "$R45"/*; do b=$(basename "$f"); case "$b" in bin|etc|R|Rscript) ;; *) ln -s "$f" "$M/$b" ;; esac; done
for f in "$R45"/bin/*; do b=$(basename "$f"); [ "$b" = R ] || ln -s "$f" "$M/bin/$b"; done
sed "s#/Library/Frameworks/R.framework/Resources#$M#g" "$R45/bin/R" > "$M/bin/R"; chmod 755 "$M/bin/R"
cp "$R45"/etc/* "$M/etc/"
sed -i '' 's#^LIBR = .*#LIBR = -L"$(R_HOME)/lib$(R_ARCH)" -lR#' "$M/etc/Makeconf"
grep -q '^LIBR = -L' "$M/etc/Makeconf" || { echo "Makeconf: LIBR was not patched" >&2; exit 1; }

TMPD="$M.tmp"; rm -rf "$TMPD"; mkdir -p "$TMPD"
RUN=(env -i HOME="$HOME" USER="$USER" LANG=en_US.UTF-8 TMPDIR="$TMPD" PATH="$M/bin:/usr/bin:/bin:/usr/sbin:/sbin" R_HOME="$M" RHOME="$M")
"${RUN[@]}" "$M/bin/R" --vanilla --slave -e 'cat(R.version.string, "| R_HOME", R.home(), "| include", Sys.getenv("R_INCLUDE_DIR"), "\n")'

# 2. Preflight: a toy C library built from the mirror must link to the 4.5 libR and not to 4.6
T="$TMPD/shlib_check"; mkdir -p "$T"
printf '#include <R.h>\n#include <Rinternals.h>\nSEXP f(void) { return ScalarInteger(1); }\n' > "$T/f.c"
(cd "$T" && "${RUN[@]}" "$M/bin/R" CMD SHLIB f.c > shlib.log 2>&1) || { cat "$T/shlib.log"; exit 1; }
otool -L "$T/f.so"
otool -L "$T/f.so" | grep -q 'Versions/4.5-arm64/Resources/lib/libR.dylib' || { echo "toy library does not link to the 4.5 libR" >&2; exit 1; }
if otool -L "$T/f.so" | grep -q 'Versions/4.6'; then echo "toy library links to 4.6" >&2; exit 1; fi
echo "preflight OK: builds link to the 4.5 libR"

# 3. Packages
"${RUN[@]}" "$M/bin/R" --no-save --no-restore --slave <<'REOF'
L45 <- "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"
stopifnot(getRversion() == "4.5.2")
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4L, install.packages.compile.from.source = "never")
have <- function(p) p %in% rownames(installed.packages(lib.loc = L45))
ip0 <- installed.packages(lib.loc = L45)[, "Version"]
cran <- c("remotes", "see", "viridis",
          # VISION v3.0.2 Imports (DESCRIPTION of the tag)
          "dplyr", "fastICA", "igraph", "irlba", "jsonlite", "logging", "pbmcapply", "Matrix", "matrixStats", "mclust", "mime",
          "plumber", "RANN", "Rcpp", "rsvd", "Rtsne", "vegan", "ape", "reticulate", "wordspace", "phytools",
          # Beyondcell v2.2.2 Imports other than DMwR, and DMwR 0.4.1 Imports
          "cowplot", "ggplot2", "ggrepel", "patchwork", "plyr", "RColorBrewer", "scales", "stringr",
          "xts", "quantmod", "zoo", "abind", "ROCR")
need <- cran[!vapply(cran, have, TRUE)]
cat("CRAN packages to install:", if (length(need)) need else "none", "\n")
if (length(need)) install.packages(need, lib = L45)
if (!have("GSVA")) BiocManager::install("GSVA", lib = L45, update = FALSE, ask = FALSE)
if (!have("DMwR")) install.packages("https://cran.r-project.org/src/contrib/Archive/DMwR/DMwR_0.4.1.tar.gz", lib = L45,
                                    repos = NULL, type = "source")
if (!have("VISION")) remotes::install_github("YosefLab/VISION@v3.0.2", lib = L45, upgrade = "never")
if (!have("beyondcell")) remotes::install_github("cnio-bu/beyondcell@v2.2.2", lib = L45, dependencies = FALSE, upgrade = "never")

ip1 <- installed.packages(lib.loc = L45)[, "Version"]
common <- intersect(names(ip0), names(ip1)); up <- common[ip0[common] != ip1[common]]
cat("\nPackages already installed whose version changed:", if (length(up)) paste0(up, " ", ip0[up], " -> ", ip1[up], collapse = "; ") else "none", "\n")
new <- setdiff(names(ip1), names(ip0))
cat("Newly installed (", length(new), "): ", paste0(new, " ", ip1[new], collapse = ", "), "\n\n", sep = "")
for (p in c("Seurat", "SeuratObject", "UCell", "AUCell", "singscore", "GSEABase", "decoupleR", "msigdbr",
            "VISION", "GSVA", "DMwR", "beyondcell", "see", "viridis", "remotes")) {
  d <- tryCatch(packageDescription(p, lib.loc = L45), warning = function(w) NULL, error = function(e) NULL)
  cat(sprintf("%-12s %-10s %s\n", p, if (is.null(d)) "MISSING" else d$Version, if (is.null(d)) "" else d$Built))
}
for (p in c("VISION", "GSVA", "beyondcell"))
  cat(p, "loads:", tryCatch({ loadNamespace(p, lib.loc = L45); TRUE }, error = function(e) conditionMessage(e)), "\n")
REOF

# 4. The compiled VISION library must link to the 4.5 libR
if [ -f "$R45/library/VISION/libs/VISION.so" ]; then
  otool -L "$R45/library/VISION/libs/VISION.so"
  if otool -L "$R45/library/VISION/libs/VISION.so" | grep -q 'Versions/4.6'; then echo "VISION.so links to 4.6" >&2; exit 1; fi
fi
echo "done"
