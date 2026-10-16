"""Fig. 1b -- illustration: cells placed on three RFind-sc indices, 3D scatter coloured by cell type.
a and b are drawn from SIMULATED (virtual) data, made by fig1ab_simulated.R: a count matrix with seven cell types,
processed as real data (Delta with a dataset-wide baseline, per-cell UP/DOWN lists, 4-term Running Fisher against three
references). The axes carry generic names (Reference A-C-index, from fig1b_axes.csv), since the values are an example.
Colours as in Fig. 1a; no key. Points outside the 1st-99th percentile of any axis are not drawn.
Input : figs/Fig1_overview/output/{fig1_sim_cells.csv (rows with in_b), fig1b_axes.csv}   (fig1ab_simulated.R)
Output: figs/Fig1_overview/output/Fig1b_3D.png   (pasted size, 49.0 x 38.6 mm, 600 dpi)
Run   : <python> figs/Fig1_overview/analysis/fig1b_3d.py   (from project root; matplotlib)
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox
import pandas as pd

OUT = Path("figs/Fig1_overview/output")
d = pd.read_csv(OUT / "fig1_sim_cells.csv")
d = d[d["in_b"]]                                          # the 8,000 cells drawn in b
axes = pd.read_csv(OUT / "fig1b_axes.csv")               # the three indices and their axis names
AX, LAB = list(axes["column"]), list(axes["label"])
print("axes:", AX)
lo, hi = d[AX].quantile(0.01), d[AX].quantile(0.99)
keep = ((d[AX] >= lo) & (d[AX] <= hi)).all(axis=1)
print(f"cells drawn: {keep.sum()} of {len(d)}")
d = d[keep].sample(frac=1, random_state=1)

## drawn at the size b is pasted into Fig. 1 (49.0 x 38.6 mm, 600 dpi). Arial; axis names 7 pt, tick labels 6 pt; axis
## lines and ticks 0.5 pt. Point sizes and pads are given for a drawing 1 / SHRINK times larger and scaled by SHRINK;
## the points are also scaled by ENLARGE (49.0 / 48.2, the panel width over the 48.2-mm width the point size was set for),
## while the text stays 7 / 6 pt.
W_MM, H_MM, DPI = 49.0, 38.6, 600                       # the size on the slide
MM = 1 / 25.4
SHRINK = 0.30                                            # scale of the point sizes and pads below
ENLARGE = W_MM / 48.2                                    # further scale of the points (panel width / 48.2 mm)
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5})
fig = plt.figure(figsize=(round(W_MM * MM * DPI) / DPI, round(H_MM * MM * DPI) / DPI))   # whole pixels at 600 dpi
ax = fig.add_axes([0, 0, 1, 1], projection="3d")
# colours = the "colour" column written by fig1ab_simulated.R, the palette of Fig. 1a, so a cell type is the
# same colour in both panels
ax.scatter(d[AX[0]], d[AX[1]], d[AX[2]], c=d["colour"], s=6.0 * (SHRINK * ENLARGE) ** 2, alpha=0.7, linewidths=0,
           depthshade=False)   # point size 6.0, scaled to the pasted size
## Unlike the UMAP axes of a, these axes are indices that can be read as amounts, so they carry names and tick numbers
## (at most three ticks per axis; axis names 7 pt, tick labels 6 pt at the pasted size).
## mplot3d places the labels (pad + 21) pt and the tick labels (pad + 8) pt from the axis, so the whole offset is scaled
LABEL_OFF, TICK_OFF = 21, 8                              # mplot3d's fixed offsets (axis3d.py, matplotlib 3.10)
for setlab, lab, pad in zip((ax.set_xlabel, ax.set_ylabel, ax.set_zlabel), LAB, (22, 22, 18)):
    setlab(lab, labelpad=(pad + LABEL_OFF) * SHRINK - LABEL_OFF, fontsize=7)
ax.zaxis.set_rotate_label(False); ax.zaxis.label.set_rotation(90)
for a in (ax.xaxis, ax.yaxis, ax.zaxis):
    a.pane.set_facecolor("white"); a.pane.set_edgecolor("grey"); a.pane.set_linewidth(0.4)
    a._axinfo["grid"].update(color="#DDDDDD", linewidth=0.25)
    a._axinfo["axisline"]["linewidth"] = 0.5
    a._axinfo["tick"]["linewidth"] = {True: 0.5, False: 0.5}
    a.line.set_linewidth(0.5)
## at most three ticks per axis: multiples of the smallest "nice" step that leaves no more than three inside the drawn range
def nice_ticks(lo_, hi_):
    for step in (1, 2, 2.5, 5, 10, 20, 25, 50, 100):
        t = [k * step for k in range(int(lo_ // step) - 1, int(hi_ // step) + 2) if lo_ <= k * step <= hi_]
        if len(t) <= 3:
            return t
    return [0]
for setticks, a in zip((ax.set_xticks, ax.set_yticks, ax.set_zticks), AX):
    setticks(nice_ticks(d[a].min(), d[a].max()))
ax.tick_params(labelsize=6, pad=(6 + TICK_OFF) * SHRINK - TICK_OFF, width=0.5)
ax.view_init(elev=18, azim=-55)
## no key: the colours only separate cell types, and which type is which is not the point of this panel
## The figure is the pasted size, so the axes are placed (not the saved area cropped): the union of the panes, the axis
## lines with their ticks and the three axis labels is fitted to the figure less a small pad, by moving and resizing the
## axes a few times (the text does not scale with the axes, so one step is not exact).
## (the Axes' own tight bbox also holds the empty 2-D axes area above the box, hence the union of the parts)
PAD = 0.3 * MM
def content_bbox():
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = Bbox.union([a.pane.get_window_extent(r) for a in (ax.xaxis, ax.yaxis, ax.zaxis)]
                    + [a.get_tightbbox(r) for a in (ax.xaxis, ax.yaxis, ax.zaxis)]
                    + [a.label.get_window_extent(r) for a in (ax.xaxis, ax.yaxis, ax.zaxis)])
    return bb.transformed(fig.dpi_scale_trans.inverted())            # inches
fw, fh = fig.get_size_inches()
for _ in range(30):
    bb = content_bbox()
    l, b, w, h = ax.get_position().bounds                            # figure fraction
    # grow/shrink the axes by the misfit of the content, and shift it so the content starts at the pad
    sx = (fw - 2 * PAD) / bb.width; sy = (fh - 2 * PAD) / bb.height
    w, h = w * sx, h * sy
    ax.set_position([l, b, w, h])
    bb = content_bbox()
    l, b = l + (PAD - bb.x0) / fw, b + (PAD - bb.y0) / fh
    ax.set_position([l, b, w, h])
bb = content_bbox()
print(f"content {bb.x0 / MM:.2f}-{bb.x1 / MM:.2f} x {bb.y0 / MM:.2f}-{bb.y1 / MM:.2f} mm of {W_MM} x {H_MM} mm")
OUT_PPTX = OUT                                  # the pasted-size panel goes to output/
OUT_PPTX.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PPTX / "Fig1b_3D.png", dpi=DPI, facecolor="white")
print("wrote Fig1b_3D.png")
