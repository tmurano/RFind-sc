# -*- coding: utf-8 -*-
"""Supplementary Table 1 — every knock-out ranked by its RFind score against the stress response of
the same cell type, for the four hippocampal cell types of Fig. 4 (the CA1, CA3 and GABA lists are given as a supplementary
table, not as figure panels; DG is included so that the table is complete).

  <python> figs/Fig4_Perturbseq/analysis/supp_table_ko_rankings.py        (from the project root)

Input : analysis/perturbseq/results/axis_{DG_vs_g017_DG, CA1_vs_CA1, CA3_vs_g017_CA3, GABA_vs_GABA3}.csv  (analysis/perturbseq/03_axis_vs_ko.py)
Output: figs/Fig4_Perturbseq/output/SuppTable_KO_rankings.xlsx (README + "All" + one sheet per cell type; no CSV copies:
        fig4_random_set_null.py and the Source Data of Fig. 4c read the sheets of the workbook)
"All" = every list in one table: one row per knock-out, rank and score in each cell type side by side.
Only neuronal cell types exist: the Perturb-seq nuclei are all neuronal, so astrocytes, microglia, oligodendrocytes and
endothelial cells (Fig. 3) have no knock-out references.
The lists are descriptive: no significance columns (the rule for list panels).
Rank: counted from the knock-out most opposed to the
stress response (rank 1 = most opposed = n - rank_among_KO + 1), as in Fig. 4c-e; the RFind scores are unchanged, and the
pipeline's rank_among_KO (1 = most aligned) is only turned round for display.
"""
import os
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

RES, OUT = "analysis/perturbseq/results", "figs/Fig4_Perturbseq/output"
SHEETS = [  # sheet, csv, Perturb-seq subclass(es) of Shi et al. 2026, De Donno et al. cell type of the stress response
    ("DG",   "axis_DG_vs_g017_DG.csv",  "037 DG Glut (dentate gyrus)", "DG Glut"),
    ("CA1",  "axis_CA1_vs_CA1.csv",     "016 CA1-ProS Glut (CA1)", "CA1 Glut"),
    ("CA3",  "axis_CA3_vs_g017_CA3.csv", "017 CA3 Glut (CA3)",     "CA3 Glut"),
    ("GABA", "axis_GABA_vs_GABA3.csv",  "052 Pvalb + 053 Sst + 049 Lamp5 Gaba, merged", "GABA"),
]
NAME = {"DG": "DG Glut", "CA1": "CA1 Glut", "CA3": "CA3 Glut", "GABA": "GABA"}
COLS = [("rank_from_opposed", "Rank"), ("target", "Knock-out"), ("rfind", "RFind score"), ("z_vs_KO", "z (among knock-outs)"),
        ("n_cells", "Nuclei"), ("target_cpm_in_non", "Target CPM (non-targeting nuclei)")]
FMT = {"RFind score": "0.00", "z (among knock-outs)": "0.00", "Target CPM (non-targeting nuclei)": "0.0"}
F, FB = Font(name="Arial", size=10), Font(name="Arial", size=10, bold=True)
HEAD = PatternFill("solid", fgColor="DDE3EA")

tables = {}
for sheet, csv, _, _ in SHEETS:
    d = pd.read_csv(os.path.join(RES, csv))
    d = d[d["kind"] == "KO"].copy()
    d["rank_among_KO"] = d["rank_among_KO"].astype(int)                      # pipeline: 1 = most aligned
    assert sorted(d["rank_among_KO"]) == list(range(1, len(d) + 1)), f"{sheet}: pipeline ranks are not 1..n"
    d["rank_from_opposed"] = len(d) - d["rank_among_KO"] + 1                   # shown: 1 = most opposed
    d = d.sort_values("rank_from_opposed")[[c for c, _ in COLS]].rename(columns=dict(COLS))
    assert d["Rank"].tolist() == list(range(1, len(d) + 1)), f"{sheet}: ranks are not 1..n"
    assert d["RFind score"].is_monotonic_increasing, f"{sheet}: scores do not follow the ranks"
    tables[sheet] = d

ALL = None
for sheet, *_ in SHEETS:
    t = tables[sheet][["Knock-out", "Rank", "RFind score"]].rename(columns={"Rank": f"{NAME[sheet]} rank", "RFind score": f"{NAME[sheet]} score"})
    ALL = t if ALL is None else ALL.merge(t, on="Knock-out", how="outer")
ALL = ALL.sort_values("Knock-out", key=lambda k: k.str.lower()).reset_index(drop=True)
tables_all = {"All": ALL}

wb = Workbook(); ws = wb.active; ws.title = "README"
readme = [
    ("Supplementary Table 1 | Rankings of the Perturb-seq knock-outs against the stress response of each hippocampal cell type", FB),
    ("", F),
    ("Each cell-type sheet lists every knock-out with at least 10 single-guide nuclei in that cell type, ranked by its RFind score against the "
     "stress response of the same cell type in the wild-type mice of De Donno et al. (rank 1 = most opposed to the stress response; the last rank = most aligned).", F),
    ("The lists are descriptive; no significance is assigned to individual knock-outs.", F),
    ("", F),
    ("Sheet", FB),
    (f"All: every knock-out ({len(ALL):,}) with its rank and RFind score in each of the four cell types side by side (blank = fewer than 10 nuclei in that cell type).", F),
    *[(f"{NAME[s]}: Perturb-seq {p}; stress response of De Donno et al. {c} (wild-type mice); {len(tables[s]):,} knock-outs", F) for s, _, p, c in SHEETS],
    ("Only neuronal cell types are included because the Perturb-seq nuclei are all neuronal.", F),
    ("", F),
    ("Columns", FB),
    ("Rank — position among the knock-outs of the cell type, counted from the knock-out most opposed to the stress response (rank 1 = most opposed; the last rank = most aligned), as in Fig. 4c–e. "
     "The direction of counting is chosen for display only; the RFind scores are not affected by it.", F),
    ("Knock-out — target gene.", F),
    ("RFind score — signed similarity score (S_total) of the knock-out's effect vector to the stress response; positive, same direction; negative, opposite.", F),
    ("z (among knock-outs) — the RFind score standardized among all knock-outs of the cell type.", F),
    ("Nuclei — number of single-guide nuclei of the knock-out whose UMI counts were summed into its pseudo-bulk profile.", F),
    ("Target CPM (non-targeting nuclei) — expression of the target gene in the non-targeting control nuclei of the cell type (counts per million).", F),
    ("", F),
    ("Method", FB),
    ("Perturb-seq: in vivo Perturb-seq of the mouse brain (doi:10.64898/2026.03.16.711480; data HuggingFace perturbai/wholebrain_crispr_atlas, CC-BY-4.0), nuclei carrying exactly one guide, the authors' cell-type labels. "
     "Knock-out effect vector: raw counts summed per knock-out (pseudo-bulk), Delta (P2 residual) against the summed counts of all knock-outs and controls of the cell type (columns with ≥ 10 nuclei). "
     "Stress response: De Donno et al., control-line cells of the cell type summed per condition, Delta of Stress against the No Stress profile. "
     "RFind-sc: top and bottom 500 genes of both vectors, four-term Running Fisher (smallest p over every top-G × top-K cut-off pair).", F),
    (f"Source files: analysis/perturbseq/results/axis_*.csv (analysis/perturbseq/03_axis_vs_ko.py); table made by "
     f"figs/Fig4_Perturbseq/analysis/supp_table_ko_rankings.py.", F),
]
for i, (text, font) in enumerate(readme, start=1):
    c = ws.cell(row=i, column=1, value=text); c.font = font; c.alignment = Alignment(wrap_text=True, vertical="top")
ws.column_dimensions["A"].width = 140

for sheet, d in [("All", ALL)] + [(sh, tables[sh]) for sh, *_ in SHEETS]:
    s = wb.create_sheet(NAME.get(sheet, sheet))   # sheet names = the names of De Donno et al.; CSV names keep the keys
    for j, name in enumerate(d.columns, start=1):
        c = s.cell(row=1, column=j, value=name); c.font = FB; c.fill = HEAD
        c.alignment = Alignment(wrap_text=True, horizontal="center", vertical="center")
    for i, row in enumerate(d.itertuples(index=False), start=2):
        for j, (name, v) in enumerate(zip(d.columns, row), start=1):
            c = s.cell(row=i, column=j, value=(None if pd.isna(v) else (v.item() if hasattr(v, "item") else v)))
            c.font = F
            if name in FMT or name.endswith(" score"): c.number_format = FMT.get(name, "0.00")
    for j, w in enumerate([14] + [10, 11] * 4 if sheet == "All" else [8, 14, 15, 13, 9, 20], start=1): s.column_dimensions[get_column_letter(j)].width = w
    s.row_dimensions[1].height = 30; s.freeze_panes = "A2"

path = os.path.join(OUT, "SuppTable_KO_rankings.xlsx"); wb.save(path)
print("wrote", path, "|", ", ".join(f"{s} {len(tables[s]):,} KOs" for s, *_ in SHEETS))
