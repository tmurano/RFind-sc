#!/bin/bash
## Kang pipeline (Supplementary Fig. S5): raw 10x matrices -> cells and genes -> per-cell-type Delta(P2) -> 4 references -> RFind scores.
## Run from project root:  bash analysis/kang/run_all.sh   (logs -> logs/runs/kang_*.log)
set -euo pipefail
cd "$(dirname "$0")/../.."
export R_HOME="${R_HOME:-/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources}"   # R 4.5 with the analysis packages
R="$R_HOME/bin/exec/R --no-save --no-restore -f"
mkdir -p logs/runs data/Kang2018/output
step(){ echo "=== [$(date +%H:%M:%S)] $1"; }
step 01_qc;         $R analysis/kang/01_qc.R         > logs/runs/kang_01.log 2>&1
step 02_signature;  $R analysis/kang/02_signature.R  > logs/runs/kang_02.log 2>&1
step 03_references; $R analysis/kang/03_references.R > logs/runs/kang_03.log 2>&1
step 04_score;      $R analysis/kang/04_score.R      > logs/runs/kang_04.log 2>&1
step DONE
