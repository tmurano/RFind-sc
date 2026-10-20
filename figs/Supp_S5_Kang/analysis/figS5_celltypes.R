## Shared cell-type definitions for Supplementary Fig. S5 (Kang et al., 2018), as Fig. 3a of the original study: the eight author
## cell types in the order of the original key (myeloid: cM, ncM, DC, Mkc; lymphoid: TH, B, Tc, NK), the abbreviations of the
## original figure (CT_LAB, plotmath strings; CT_TXT, plain text) and its colours (RColorBrewer "Set2" in the alphabetical order
## of the authors' labels); conditions as in the original (Paired blue and red).
CT_AUTH <- c("B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells")
CT_PAL <- setNames(RColorBrewer::brewer.pal(8, "Set2"), CT_AUTH)
CT_LV <- c("CD14+ Monocytes", "FCGR3A+ Monocytes", "Dendritic cells", "Megakaryocytes", "CD4 T cells", "B cells", "CD8 T cells", "NK cells")
CT_LAB <- c(`CD14+ Monocytes` = "cM", `FCGR3A+ Monocytes` = "ncM", `Dendritic cells` = "DC", Megakaryocytes = "Mkc",
            `CD4 T cells` = "T[H]", `B cells` = "B", `CD8 T cells` = "Tc", `NK cells` = "NK")
CT_TXT <- sub("T\\[H\\]", "TH", CT_LAB)
CT_MYELOID <- c("CD14+ Monocytes", "FCGR3A+ Monocytes", "Dendritic cells")   # myeloid leukocytes (Megakaryocytes apart)
COND_LV <- c("ctrl", "stim"); COND_LAB <- c(ctrl = "Control", stim = "IFN-β")   # internal labels
COND_SHOW <- c(Control = "Unstimulated", `IFN-β` = "IFN-β, 6 h")    # labels shown: say the stimulus (as S5a)
COND_PAL <- c(Control = RColorBrewer::brewer.pal(12, "Paired")[2], `IFN-β` = RColorBrewer::brewer.pal(12, "Paired")[6])
FIG <- "figs/Supp_S5_Kang/output"; dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
