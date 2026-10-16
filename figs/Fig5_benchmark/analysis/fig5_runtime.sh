#!/bin/bash
## Run time and peak memory (Supplementary Table 3, sheet "Run time and memory"): scoring the cells of the data of
## Hammond et al. against the DAM reference with RFind-sc (R package RFindsc and Python package rfindsc), ULM of decoupleR and UCell
## (fig5_runtime_score.R / .py), for 1,000 to 75,572 cells, references of 500 genes per direction or the whole list, and 1 or 6 cores;
## and the knock-out x knock-out matching of Supplementary Fig. S3 (figS3_rfind_allpairs.py, 2,048 profiles of TH Prkcd Grin2c Glut,
## 7 workers). Each run is timed under /usr/bin/time -l (macOS; maximum resident set size of the main process, so the memory of the
## worker processes of multi-core runs is not included: the table gives the memory of the single-core runs).
## Run on an otherwise idle computer (the times are meaningless while other jobs run), from the project root (about 1-2 h):
##   bash figs/Fig5_benchmark/analysis/fig5_runtime.sh
## Output: figs/Fig5_benchmark/output/results/fig5_runtime.csv (one row per run) and fig5_runtime_env.tsv (the computer and the
## versions of the run; both go into the sheet "Run time and memory" of fig5_redo_tables.py).
set -u
cd "$(dirname "$0")/../../.." || exit 1
export R_HOME="${R_HOME:-/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources}"
R="$R_HOME/bin/exec/R --no-save --no-restore -f"
PY="${PYTHON:-python3}"; export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="${RFINDSC_PY_SRC:-tools/rfindsc-py/src}${PYTHONPATH:+:$PYTHONPATH}"   # rfindsc from its source, unless installed
RES=figs/Fig5_benchmark/output/results; RAW=$RES/fig5_runtime_raw.tsv; MEM=$RES/fig5_runtime_mem.tsv
rm -f "$RAW" "$MEM"; printf "label\texit\tpeak_rss_bytes\twall_s\n" > "$MEM"
{ printf "item\tvalue\n"; printf "date\t%s\n" "$(date +%Y-%m-%d)"
  printf "computer\t%s (%s, %s cores, %s GB memory)\n" "$(sysctl -n hw.model)" "$(sysctl -n machdep.cpu.brand_string)" \
    "$(sysctl -n hw.ncpu)" "$(( $(sysctl -n hw.memsize) / 1073741824 ))"
  printf "operating system\t%s %s\n" "$(sw_vers -productName)" "$(sw_vers -productVersion)"
  "$R_HOME/bin/exec/R" --no-save --no-restore -s -e 'cat(sprintf("R\t%s.%s\n", R.version$major, R.version$minor));
    for (p in c("RFindsc", "decoupleR", "UCell", "Matrix")) cat(sprintf("%s\t%s\n", p, as.character(packageVersion(p))))'
  "$PY" -c 'import sys, numpy, scipy, pandas, rfindsc; print(f"Python\t{sys.version.split()[0]}")
for m in (rfindsc, numpy, scipy, pandas): print(f"{m.__name__}\t{m.__version__}")'
} > "$RES/fig5_runtime_env.tsv"
TMPD=$(mktemp -d); export RFIND_RUNTIME_EXPORT="$TMPD"; trap 'rm -rf "$TMPD"' EXIT
measure() {  # measure <label> <command...>
  local label=$1; shift; local t0=$(date +%s)
  /usr/bin/time -l "$@" > "$TMPD/run.log" 2>&1; local rc=$?
  local rss=$(awk '/maximum resident set size/ {print $1}' "$TMPD/run.log")
  printf "%s\t%s\t%s\t%s\n" "$label" "$rc" "${rss:-NA}" "$(( $(date +%s) - t0 ))" | tee -a "$MEM"
  [ "$rc" = 0 ] || tail -5 "$TMPD/run.log"
}
$R figs/Fig5_benchmark/analysis/fig5_runtime_export.R > "$TMPD/export.log" 2>&1 || { cat "$TMPD/export.log"; exit 1; }
for spec in "1000 500 1" "10000 500 1" "all 500 1" "all 500 6" "10000 all 1" "all all 6"; do
  set -- $spec
  measure "rfindsc_R|$1|$2|$3"  $R figs/Fig5_benchmark/analysis/fig5_runtime_score.R --args method=rfindsc n=$1 L=$2 cores=$3
  measure "rfindsc_py|$1|$2|$3" $PY figs/Fig5_benchmark/analysis/fig5_runtime_score.py n=$1 L=$2 cores=$3
done
for spec in "all 500 1" "all all 1"; do
  set -- $spec; measure "ulm|$1|$2|$3" $R figs/Fig5_benchmark/analysis/fig5_runtime_score.R --args method=ulm n=$1 L=$2 cores=$3
done
for spec in "all 500 1" "all 500 6" "all all 6"; do
  set -- $spec; measure "ucell|$1|$2|$3" $R figs/Fig5_benchmark/analysis/fig5_runtime_score.R --args method=ucell n=$1 L=$2 cores=$3
done
measure "allpairs_g151|2048|500|7" $PY figs/Supp_S3_Perturbseq/analysis/figS3_rfind_allpairs.py 7 delta g151
$PY - <<'EOF'
import pandas as pd
RES = 'figs/Fig5_benchmark/output/results'
raw = pd.read_csv(f'{RES}/fig5_runtime_raw.tsv', sep='\t'); mem = pd.read_csv(f'{RES}/fig5_runtime_mem.tsv', sep='\t')
mem[['method', 'n', 'L', 'cores']] = mem.label.str.split('|', expand=True)
name = {'rfindsc': 'RFind-sc (R package RFindsc)', 'rfindsc_py': 'RFind-sc (Python package rfindsc)', 'ulm': 'ULM (decoupleR)', 'ucell': 'UCell'}
rows = []
for _, m in mem.iterrows():
    key = 'rfindsc' if m.method == 'rfindsc_R' else m.method
    if m.method == 'allpairs_g151':
        rows.append({'Task': 'Knock-out x knock-out matching (Supplementary Fig. S3)', 'Method': 'RFind-sc (analysis code, Python)',
                     'Cells or profiles': 2048, 'Reference genes (up / down)': '500 / 500', 'Cores': 7, 'Time (s)': m.wall_s,
                     'Peak memory (GB)': None, 'Exit code': m.exit})
        continue
    r = raw[(raw.method == key) & (raw.cores == int(m.cores)) & (raw.L.astype(str) == m.L)]
    r = r[r.cells == (r.cells.max() if m.n == 'all' else int(m.n))]
    rows.append({'Task': 'Scoring cells against the DAM reference (data of Hammond et al.)', 'Method': name[key],
                 'Cells or profiles': int(r.cells.iloc[0]) if len(r) else None,
                 'Reference genes (up / down)': f'{int(r.up.iloc[0])} / {int(r.down.iloc[0])}' if len(r) else None,
                 'Cores': int(m.cores), 'Time (s)': round(float(r.seconds.iloc[0]), 1) if len(r) else None,
                 'Peak memory (GB)': round(float(m.peak_rss_bytes) / 1e9, 2) if m.cores == '1' and m.peak_rss_bytes != 'NA' else None,
                 'Exit code': m.exit})
pd.DataFrame(rows).to_csv(f'{RES}/fig5_runtime.csv', index=False)
print(pd.DataFrame(rows).to_string(index=False))
EOF
