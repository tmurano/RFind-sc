#!/bin/bash
## De Donno pipeline: raw h5ad -> QC/UMAP -> per-cell-type Delta(P2) -> 4 references -> RFind scores.
## Run from project root:  bash analysis/dedonno/run_all.sh   (logs -> logs/runs/dedonno_*.log)
set -euo pipefail
cd "$(dirname "$0")/../.."
export R_HOME="${R_HOME:-/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources}"   # R 4.5 with the analysis packages
R="$R_HOME/bin/exec/R --no-save --no-restore -f"
mkdir -p logs/runs data/DeDonno2026/output
step(){ echo "=== [$(date +%H:%M:%S)] $1"; }
step 01_qc_embed;   $R analysis/dedonno/01_qc_embed.R   > logs/runs/dedonno_01.log 2>&1
step 02_signature;  $R analysis/dedonno/02_signature.R  > logs/runs/dedonno_02.log 2>&1
step 03_references; $R analysis/dedonno/03_references.R > logs/runs/dedonno_03.log 2>&1
step 04_score;      $R analysis/dedonno/04_score.R      > logs/runs/dedonno_04.log 2>&1
step DONE
