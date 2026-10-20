#!/bin/bash
## Fig. 5 benchmark — the per-cell scores of every method at every reference length, in one command
## (run by figs/run_all_figs.sh, which then runs the other Fig. 5 scripts in order; see README "Which script makes what"). From any folder:
##   bash figs/Fig5_benchmark/analysis/run_fig5_redo.sh          the test step, then everything (about 2 h on Apple silicon)
##   bash figs/Fig5_benchmark/analysis/run_fig5_redo.sh test     the test step only (about 5 min)
## Steps (logs and run_summary.txt in logs/runs/fig5_redo_<yymmdd>/, <yymmdd> = the day of the start; a finished step leaves <step>.ok and
## is skipped when the script is started again with the same log folder, so an interrupted run resumes where it stopped: on a later
## day, set FIG5_LOGD=logs/runs/fig5_redo_<yymmdd of the first start>; FORCE=1 re-runs every step):
##   test    fig5_bench_redo.R test: all methods, 600 cells per task, outputs in FIG5_TEST_DIR (default logs/runs/fig5_test/)
##           -> stops here on any error (packages, the ORA check against decoupleR::run_ora, ...)
##   L500    fig5_bench_redo.R L=500 (the main results; checks against the stored results) -> stops on any error
##   L1000, L2000, Lall, L250, L100   the other reference lengths -> an error is logged and the other steps go on
## Supplementary Table 3 (fig5_redo_tables.py) and the figure panels are NOT made here: they also need the added type A tasks
## (fig5_extra_refs.R, fig5_extra_score.R), fig5_ucell_zero.R, fig5_best_length.py, fig5_best_length_8tasks.py and
## fig5_balanced_accuracy.py, which figs/run_all_figs.sh runs after this script.
## Time (the run of 2026-10-07 on an Apple M3 with 16 GB, R 4.5.2): test 2 min, L500 39 min, L1000 22 min, L2000 25 min,
## Lall 47 min, L250 15 min, L100 15 min.
## Memory: about 10 GB per benchmark step (the steps run one after another; 16 GB RAM or more). On macOS the computer is kept awake.
## Do not edit the analysis scripts while this runs: the copies in logs/runs/fig5_redo_<yymmdd>/scripts/ are what runs.
## Environment: FIG5_R_HOME or R_HOME = R 4.5 with the benchmark packages (README, "Requirements"); FIG5_LOGD = the log folder
##   (default logs/runs/fig5_redo_<yymmdd of today>; set it to resume a run started on an earlier day); FIG5_TEST_DIR = the output folder of
##   the test step (default logs/runs/fig5_test/; also read by the R scripts it runs).
{   # the whole script is one block, read completely before it runs (bash also reads a script while it runs)
set -u
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"; cd "$ROOT" || exit 1
[ -f README.md ] && [ -d figs ] && [ -f analysis/running_fisher.R ] || { echo "not the RFind-sc project: $ROOT"; exit 1; }
R_HOME="${FIG5_R_HOME:-${R_HOME:-/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources}}"   # R 4.5.2 in the paper
RBIN="$R_HOME/bin/exec/R"
[ -x "$RBIN" ] || { echo "R not found: $RBIN (set R_HOME or FIG5_R_HOME)"; exit 1; }
MODE="${1:-all}"
case "$MODE" in all|test|after) ;; *) echo "usage: $0 [test|after]"; exit 1 ;; esac
D=$(date +%y%m%d); LOGD="${FIG5_LOGD:-logs/runs/fig5_redo_$D}"; mkdir -p "$LOGD"   # FIG5_LOGD: resume a run started on an earlier day
if command -v caffeinate >/dev/null 2>&1; then caffeinate -ims -w $$ & fi                  # macOS: no sleep while this script runs
SUM="$LOGD/run_summary.txt"; FAILED=""
say() { echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$SUM"; }
step() {                                                                 # step <name> <critical 1/0> <command...>
  local name=$1 crit=$2; shift 2
  if [ -z "${FORCE:-}" ] && [ -f "$LOGD/$name.ok" ]; then say "skip  $name (finished before)"; return 0; fi
  say "start $name"
  local t0; t0=$(date +%s)
  "$@" > "$LOGD/$name.log" 2>&1
  local rc=$?; local t=$(( ($(date +%s) - t0) / 60 ))
  if [ $rc -eq 0 ]; then touch "$LOGD/$name.ok"; say "done  $name ($t min)"; return 0; fi
  say "FAILED $name (exit $rc after $t min); the last lines of $LOGD/$name.log:"; tail -15 "$LOGD/$name.log" | tee -a "$SUM"
  FAILED="$FAILED $name"
  if [ "$crit" = 1 ]; then say "stopped: a step the others depend on failed"; exit $rc; fi
  return 0
}
rrun() { env R_HOME="$R_HOME" "$RBIN" --no-save --no-restore -f "$@"; }
## R reads a script given with -f while it runs, so a script edited during the run (for example through a synced folder) breaks
## the run. The scripts are therefore copied first and the copies are run; they also record the version used.
A="$LOGD/scripts"; mkdir -p "$A"
SCRIPTS="fig5_bench_redo.R"
for s in $SCRIPTS; do cp -p "figs/Fig5_benchmark/analysis/$s" "$A/$s" || exit 1; done
md5s() { if command -v md5 >/dev/null 2>&1; then md5 -q "$@"; else md5sum "$@" | cut -d" " -f1; fi; }
say "Fig. 5 benchmark, mode $MODE, project $ROOT, logs $LOGD; scripts copied to $A (md5: $(cd "$A" && md5s $SCRIPTS | cut -c1-8 | tr '\n' ' '))"
if [ "$MODE" != after ]; then
  step test 1 rrun "$A/fig5_bench_redo.R" --args test
  if [ "$MODE" = test ]; then say "test only: finished (outputs in ${FIG5_TEST_DIR:-logs/runs/fig5_test}/)"; exit 0; fi
  step L500 1 rrun "$A/fig5_bench_redo.R" --args L=500 logdir="$LOGD"
  for L in 1000 2000 all 250 100; do step "L$L" 0 rrun "$A/fig5_bench_redo.R" --args L=$L logdir="$LOGD"; done
fi
if [ -n "$FAILED" ]; then say "finished with failures:$FAILED (see the logs; the finished steps are skipped when the script is started again with FIG5_LOGD=$LOGD)"; exit 1; fi
say "all steps finished; summaries in $LOGD/fig5_redo_L<L>.txt; the tables and panels follow in figs/run_all_figs.sh"
exit 0
}
