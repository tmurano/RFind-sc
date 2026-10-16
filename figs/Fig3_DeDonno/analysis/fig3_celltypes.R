## Shared cell-type definitions for Fig3 (De Donno): the 22 author types (author `cell_type_short` with the
## number suffix removed), their order, display labels, colours and neuron class. In the wild-type cells analysed, 21 of
## them remain (Choroid Plexus has 11 cells and is dropped by the 20-cell rule of 01_qc_embed.R).
## Sourced by fig3_umap.R, fig3_violins.R, fig3_scatter.R and fig3_celltype_table.R so all panels agree.

CT_LV <- c("DG Glut", "CA1 Glut", "CA2 Glut", "CA3 Glut", "Mossy Glut", "GABA", "Cajal-Retzius",
           "Granule Neuroblasts", "NPCs",
           "Astro", "Microglia", "Macrophages", "Oligodendrocytes", "MFOLs", "OPCs", "COPs",
           "Endothelial", "Vascular", "Pericytes", "Meningeal", "Ependymal", "Choroid Plexus")

## labels used on the figures = the names of the original study;
## the descriptions (e.g. "Dentate gyrus granule neurons") are in the Fig. 3b table
CT_LAB <- setNames(CT_LV, CT_LV)

## colours: neurons blues/purples, progenitors pinks, astro green, microglia/macrophage reds,
## oligodendrocyte lineage oranges/yellows, vascular and other non-neural cells blues/violet/teal/brown
CT_PAL <- c(`DG Glut` = "#1F6FB4", `CA1 Glut` = "#3F9FD6", `CA2 Glut` = "#A6CEE3", `CA3 Glut` = "#6FC1E6",
            `Mossy Glut` = "#23395D", GABA = "#7E57C2", `Cajal-Retzius` = "#B39DDB",
            `Granule Neuroblasts` = "#E879B8", NPCs = "#AD1457",
            Astro = "#2E9E5B", Microglia = "#E0322B", Macrophages = "#F4A09C",
            Oligodendrocytes = "#F08A24", MFOLs = "#B35806", OPCs = "#F6C445", COPs = "#D9C27A",
            Endothelial = "#4363D8", Vascular = "#8E9CE8", Pericytes = "#911EB4", Meningeal = "#469990",
            Ependymal = "#8C564B", `Choroid Plexus` = "#9E9D24")

CT_NEURON <- c("DG Glut", "CA1 Glut", "CA2 Glut", "CA3 Glut", "Mossy Glut", "GABA", "Cajal-Retzius",
               "Granule Neuroblasts", "NPCs")   # neurons and neuronal progenitors (left of the dotted line in violins)

stopifnot(setequal(names(CT_LAB), CT_LV), setequal(names(CT_PAL), CT_LV), all(CT_NEURON %in% CT_LV))
