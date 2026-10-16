"""Fig. 1c: the formula of Delta as a separate image for the slide (Fig. 1c is assembled on the slide; the formula goes under
the arrow (5)). Arial, 7 pt (Nature: text 5-7 pt), main letters italic and Delta upright, as in Fig. 1c and Supplementary
Fig. S2; transparent margins cut to the drawn content plus 0.3 mm.
Output: figs/Fig1_overview/output/Fig1c_formula.png (paste at 100%)
Run from the project root: ~/miniforge3/bin/python3 figs/Fig1_overview/analysis/fig1c_formula.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                     "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})
DPI, FS = 600, 7
fig = plt.figure(figsize=(40 / 25.4, 15 / 25.4), dpi=DPI)
fig.text(0.5, 0.5, r"$\Delta_{ij} = \dfrac{x_{ij} - e_{ij}}{\sqrt{e_{ij} + 1}}$", ha="center", va="center", fontsize=FS)
out = "figs/Fig1_overview/output/Fig1c_formula.png"
fig.savefig(out, dpi=DPI, facecolor="white")
im = plt.imread(out); ink = (im[..., :3] < 0.9).any(axis=2)
ys, xs = np.nonzero(ink); m = int(round(0.3 / 25.4 * DPI))
im = im[max(ys.min() - m, 0):ys.max() + m + 1, max(xs.min() - m, 0):xs.max() + m + 1]
plt.imsave(out, im, dpi=DPI)
print(f"wrote {out}: {im.shape[1] / DPI * 25.4:.1f} x {im.shape[0] / DPI * 25.4:.1f} mm")
