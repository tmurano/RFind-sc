#!/bin/bash
## Hammond pipeline: raw h5ad -> QC -> gene universe -> Delta(P2) -> references -> RFind scores.
## Run from project root:  bash analysis/hammond/run_all.sh   (logs -> logs/runs/hammond_*.log)
## Inputs (README, "Input data"): data/Hammond2019/raw/Hammond2019_microglia_AllSamples_annotated.h5ad,
## data/Matcovitch2016/raw/GSE79812_Timeline.xlsx and the fixed inputs in data/Hammond2019/provided/ and
## data/KerenShaul2017/provided/; step 06 also downloads the MSigDB Hallmark gene sets (msigdbr).
set -euo pipefail
cd "$(dirname "$0")/../.."
export R_HOME="${R_HOME:-/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources}"   # R 4.5 with the analysis packages
R="$R_HOME/bin/exec/R --no-save --no-restore -f"
PY="${PYTHON:-python3}"   # needs scrublet, h5py, scipy (e.g. PYTHON=/path/to/python)
mkdir -p logs/runs data/Hammond2019/raw data/Hammond2019/output data/Hammond2019/output/references
step(){ echo "=== [$(date +%H:%M:%S)] $1"; }
## Per-cell metadata and mitochondrial fraction: fixed inputs, derived from GEO GSE121654 and the annotated object (the
## mitochondrial fraction was computed from the GEO DGE files, as the annotated object has no mitochondrial genes). The
## scripts read them from data/Hammond2019/raw/, so they are copied there when missing.
for f in cell_meta_hammond.csv mt_frac_hammond.csv; do
  if [ ! -f "data/Hammond2019/raw/$f" ]; then
    step "copy data/Hammond2019/provided/$f -> data/Hammond2019/raw/"; cp -p "data/Hammond2019/provided/$f" "data/Hammond2019/raw/$f"
  fi
done
step 01_qc_mt_umi;         $R analysis/hammond/01_qc_mt_umi.R         > logs/runs/hammond_01.log 2>&1
## Doublets: Scrublet's scores depend on the numpy/scipy/scikit-learn/annoy versions (on another computer 145 of the 75,572 kept
## cells changed, 2026-09-15), so the paper uses the calls of 2026-09-14, kept as a fixed input. Delete that file to re-run Scrublet.
DBL=data/Hammond2019/provided/scrublet_doublets_260914.csv
if [ -f "$DBL" ]; then
  step "02_scrublet (fixed calls: $DBL)"; cp "$DBL" data/Hammond2019/output/scrublet_doublets.csv
else
  step 02_scrublet;        $PY analysis/hammond/02_scrublet.py        > logs/runs/hammond_02.log 2>&1
fi
step 03_extract_qc;        $R analysis/hammond/03_extract_qc.R        > logs/runs/hammond_03.log 2>&1
step 04_gene_filter;       $R analysis/hammond/04_gene_filter_relSE.R > logs/runs/hammond_04.log 2>&1
step 05_delta_P2;          $R analysis/hammond/05_delta_P2.R          > logs/runs/hammond_05.log 2>&1
step 06_build_references;  $R analysis/hammond/06_build_references.R  > logs/runs/hammond_06.log 2>&1
step 07_score_axes;        $R analysis/hammond/07_score_axes.R        > logs/runs/hammond_07.log 2>&1
step DONE
