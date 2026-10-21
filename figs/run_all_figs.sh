#!/bin/bash
## Every figure, Supplementary Tables 1-3 and the Source Data workbooks from the pipeline outputs, in dependency order.
## Needs the pipelines to have been run (analysis/hammond, analysis/dedonno, analysis/kang, analysis/perturbseq).
## Run from the project root:  bash figs/run_all_figs.sh        (logs -> logs/runs/<step>.log, summary -> logs/runs/run_all_figs_summary.tsv)
## Environment: R_HOME = R 4.5 with the analysis packages; PYTHON = Python with numpy, scipy, pandas, matplotlib, openpyxl
## (e.g. PYTHON=/path/to/python; when unset, the first of python3 and python on the PATH that imports these packages, otherwise
## python3). Continues after a failed step; check the summary.
## Time: about 4 h in all (run of 2026-10-19 on an Apple M3 with 16 GB): the Fig. 5 benchmark about 2 h, the rest of Fig. 5 about
## 1 h, Supplementary Fig. S3 about 1 h, everything else a few minutes.
## The Fig. 5 benchmark needs about 10 GB of memory (16 GB RAM or more).
cd "$(dirname "$0")/.." || exit 1
[ -f README.md ] && [ -d figs ] && [ -f analysis/running_fisher.R ] || { echo "not the RFind-sc project: $(pwd)"; exit 1; }
export R_HOME="${R_HOME:-/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources}"
R="$R_HOME/bin/exec/R --no-save --no-restore -f"
if [ -z "${PYTHON:-}" ]; then
  PYTHON=python3
  for c in python3 python; do
    if command -v "$c" > /dev/null 2>&1 && "$c" -c 'import numpy, scipy, pandas, matplotlib, openpyxl' > /dev/null 2>&1; then
      PYTHON=$c; break
    fi
  done
fi
export PYTHON; PY="$PYTHON"
LOGS=logs/runs; mkdir -p "$LOGS"; SUM="$LOGS/run_all_figs_summary.tsv"; : > "$SUM"   # run logs in logs/runs/ (logs/ itself keeps only notes)
run() {  # $1 = log name, rest = command
  local name=$1; shift; local t0=$(date +%s)
  echo "=== [$(date +%H:%M:%S)] $name"
  "$@" > "$LOGS/$name.log" 2>&1; local rc=$?
  printf "%s\t%s\t%ss\n" "$name" "$rc" "$(( $(date +%s) - t0 ))" | tee -a "$SUM"
}

## Supplementary Fig. S2 (Delta): data (R), then the panels and the whole figure (Python)
run figS2_delta_benchmark  $R figs/Supp_S2_Delta/analysis/figS2_delta_benchmark.R
run figS2_panels           $PY figs/Supp_S2_Delta/analysis/figS2_panels.py
## Fig. 1c (assembled on the slide; only the Delta formula is an image, 7 pt)
run fig1c_formula          $PY figs/Fig1_overview/analysis/fig1c_formula.py
## Fig. 2 (Hammond)
for f in fig2_tsne_original fig2_violins fig2_cluster_table fig2_scatter_pairs; do run "fig2_$f" $R "figs/Fig2_Hammond/analysis/$f.R"; done
## Fig. 3 (De Donno): violins first (the UMAP and the table read its list of cell types, Fig3_celltype_ncells.csv)
for f in fig3_violins fig3_umap fig3_celltype_table fig3_scatter; do run "fig3_$f" $R "figs/Fig3_DeDonno/analysis/$f.R"; done
## Supplementary Fig. S5 (Kang, human PBMC): the markers of b from the 2.7k PBMC tutorial, then a and b, the violins (c, d), the
## scatter maps (e, f) and the whole figure (about 1 min)
for f in figS5_pbmc3k_markers figS5_ab figS5_violins figS5_scatter; do run "$f" $R "figs/Supp_S5_Kang/analysis/$f.R"; done
run figS5_figure           $PY figs/Supp_S5_Kang/analysis/figS5_compose.py
## Fig. 1
## a and b are drawn from simulated data
run fig1ab_simulated       $R figs/Fig1_overview/analysis/fig1ab_simulated.R
run fig1b_3d               $PY figs/Fig1_overview/analysis/fig1b_3d.py
run fig1d                  $R figs/Fig1_overview/analysis/fig1d_matching.R
## Supplementary Fig. S1 (Running Fisher schematic)
run figS1_running_fisher   $PY figs/Supp_S1_RunningFisher/analysis/figS1_running_fisher.py
## Fig. 4 (Perturb-seq)
run fig4a_schematic       $PY figs/Fig4_Perturbseq/analysis/fig4a_schematic.py      # a (schematic)
run fig4b_schematic       $PY figs/Fig4_Perturbseq/analysis/fig4b_schematic.py      # b (whole panel; needs the De Donno QC metadata)
run fig4c_CA1              $PY figs/Fig4_Perturbseq/analysis/fig4c_rank_tables.py CA1 CA1
run fig4d_GR-complex       $R figs/Fig4_Perturbseq/analysis/fig4d_go_set.R --args GR-complex   # 4d: all knock-outs ranked by RFind score, the five GR-complex knock-outs marked
run fig4e_expression       $PY figs/Fig4_Perturbseq/analysis/fig4e_expression_vs_score.py   # 4e: Nr3c1 and Ncoa1 expression (De Donno No Stress) vs the RFind score of their knock-outs
run supp_table_ko_rankings $PY figs/Fig4_Perturbseq/analysis/supp_table_ko_rankings.py   # needs openpyxl
run fig4_ranking_checks     $PY figs/Fig4_Perturbseq/analysis/fig4_ranking_checks.py   # the check of Methods (lower limit of 5 nuclei)
run fig4_random_sets      $PY figs/Fig4_Perturbseq/analysis/fig4_random_set_null.py     # Fig. 4d: the five knock-outs against random sets of the library (P in Results; Source Data)
## Fig. 5, Supplementary Fig. S4 and Supplementary Table 2 (capability benchmark), in the order of their inputs:
##  1  the benchmark: per-cell scores of every method at every reference length (fig5_redo_{auc,percell}_L<L>.*,
##     results/fig5_redo_beyondcell_L500.csv);
##     about 2 h and ~10 GB RAM; logs and a run summary in logs/runs/fig5_redo_<yymmdd>/ (<yymmdd> = the day of the start); a finished
##     step is skipped when started again on the same day (on a later day, set FIG5_LOGD=logs/runs/fig5_redo_<yymmdd of the first start>).
##     The per-cell files are not in the repository, so this step is always needed on a fresh copy.
run fig5_benchmark           bash figs/Fig5_benchmark/analysis/run_fig5_redo.sh
##  1b the two added type A references of Supplementary Fig. S4 (Forced swim, GSE169505; PS2APP, E-GEOD-75431; downloads in
##     README, Input data): extra_refs/fig5_extra_refs.rds and the fold-change tables in data/GSE169505/output/, data/GSE75431/output/
##     (a few seconds)
run fig5_extra_refs          $R figs/Fig5_benchmark/analysis/fig5_extra_refs.R
##  1c the four added type A tasks (DA1, AS1-AS3) at every reference length: extra_refs/fig5_extra_{auc,percell}_L<L>.*
##     (about 45 min in all and ~10 GB RAM; summaries in logs/runs/fig5_extra/; the per-cell files are not in the repository)
for L in 500 100 250 1000 2000 all; do
  run "fig5_extra_score_L$L"  $R figs/Fig5_benchmark/analysis/fig5_extra_score.R --args L=$L
done
##  2  UCell's signed score split into U+ and U- (results/fig5_ucell_zero.csv; the shares in the Source Data of Fig. 5; a few min)
run fig5_ucell_zero          $R figs/Fig5_benchmark/analysis/fig5_ucell_zero.R
##  3  every method at its best reference length (results/fig5_best_length.csv and results/fig5_mean_auc_diff.csv; reads the
##     per-cell files of step 1)
run fig5_best_length         $PY figs/Fig5_benchmark/analysis/fig5_best_length.py
##  3b the same over the eight type A tasks (results/fig5_best_length_8tasks.csv and results/fig5_mean_auc_diff_8tasks.csv; reads the per-cell files of
##     steps 1 and 1c and checks the outputs of step 3 first; about 3 min)
run fig5_best_length_8tasks  $PY figs/Fig5_benchmark/analysis/fig5_best_length_8tasks.py
##  3c the mean over the tasks of each type (results/fig5_task_means.csv; Fig. 5e and the sheet Mean over tasks of Supplementary Table 2;
##     reads the per-cell files of steps 1 and 1c at L = 500 and checks the type A differences against step 3b; about 3 min)
run fig5_task_means          $PY figs/Fig5_benchmark/analysis/fig5_task_means.py
##  4  output/data/Fig5c_embedding.png and Fig5d_auc.png (the t-SNE and the DAM-index distributions that panel b of Fig. 5 is cut
##     from; reads step 1 only). Run here, before
##     steps 5, 6 and 7: fig5_rebuild_figure_draft.py (step 7, imported by step 6) reads these two images when it is loaded.
run fig5_redo_figure         $PY figs/Fig5_benchmark/analysis/fig5_redo_figure.py
##  5  Supplementary Table 2 (SuppTable2.xlsx; reads steps 1, 1c, 3 and 3b)
run fig5_supp_table3         $PY figs/Fig5_benchmark/analysis/fig5_redo_tables.py
##  6  Supplementary Fig. S4: all tasks and the type A tasks not in Fig. 5c (output/FigS4.png and the panels FigS4a_tasks.png and
##     FigS4b_typeA_tasks.png, directly in output/); reads steps 1c, 5 [SuppTable2.xlsx, sheet Tasks] and 4
##     [Fig5c_embedding.png and Fig5d_auc.png, read when fig5_rebuild_figure_draft.py is imported] and imports the script of step 7
##     (which reads step 1)
run figS4_capability_map     $PY figs/Supp_S4_Benchmark/analysis/figS4_capability_map.py
##  7  Fig. 5: the whole figure (output/Fig5.png) and the panels pasted in the slide (output/Fig5a_methods_table.png ... Fig5e_task_means.png, Fig5f_reference_length.png and Fig5_methods_key.png); e and f average
##     the eight type A tasks (reads steps 1, 1c and 4)
run fig5_figure              $PY figs/Fig5_benchmark/analysis/fig5_rebuild_figure_draft.py
## The scripts of an earlier stored run of the Fig. 5 benchmark are not in the repository; its two result tables kept in
## figs/Fig5_benchmark/output/ (results/fig5_final_auc.csv, results/fig5_reflength_sweep_auc.csv; the tasks not in CHANGED of
## fig5_bench_redo.R) are read by fig5_bench_redo.R only as a check.
## Supplementary Fig. S3 (Perturb-seq knock-outs as references; baseline = the summed counts of all columns, as Fig. 4; b-d in the thalamic
## glutamatergic neurons TH Prkcd Grin2c Glut = g151, the cell type Shi et al. name in their main-text comparison of knock-outs;
## e = the thalamus and the cell types of Fig. 4)
run figS3_allpairs_g151    $PY figs/Supp_S3_Perturbseq/analysis/figS3_rfind_allpairs.py 7 delta g151  # ~14 min
run figS3_allpairs_GABA    $PY figs/Supp_S3_Perturbseq/analysis/figS3_rfind_allpairs.py 7 delta GABA  # ~15 min
run figS3_allpairs_HIP     $PY figs/Supp_S3_Perturbseq/analysis/figS3_rfind_allpairs.py 7 delta HIP   # e, hippocampus: DG + CA1 + CA3 (Delta summed per column; ~27 min)
run figS3_panels           $PY figs/Supp_S3_Perturbseq/analysis/figS3_rfind_panels.py                 # b, c, d (thalamus)
run figS3_celltypes        $PY figs/Supp_S3_Perturbseq/analysis/figS3e_celltypes.py                   # e
run figS3_schematic        $R figs/Supp_S3_Perturbseq/analysis/figS3_schematic.R
## Supplementary Table 3 (the references and their numbers of genes; reads the reference files of the Hammond, De Donno and Kang
## pipelines and of step 1b of Fig. 5): the numbers (R), then the workbook
run supp_table3_sizes      $R figs/Supp_Table3_References/analysis/supp_table_references.R
run supp_table3_references $PY figs/Supp_Table3_References/analysis/supp_table_references.py  # needs openpyxl
## Source Data (one workbook per figure; read the outputs above, Fig. 5 and S4 also those of step 1c): Figs. 2-5 and
## Supplementary Fig. S4, then Supplementary Figs. S2, S3 and S5 (S1 is a made-up example: no Source Data)
run source_data              $PY figs/make_source_data.py
run source_data_supp         $PY figs/make_source_data_supp.py
echo "=== [$(date +%H:%M:%S)] DONE (summary: $SUM)"
