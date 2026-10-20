"""Supplementary Fig. S5 - all panels on one page at their drawn sizes (600 dpi), with panel letters (Arial bold 8 pt).
Rows: a (the authors' released t-SNE, with full cell-type names in the key) + b (cell types and markers); c + d (violins); e + f (scatter maps).
Run from the project root after the panel scripts: python figs/Supp_S5_Kang/analysis/figS5_compose.py
Output: figs/Supp_S5_Kang/output/FigS5.png
"""
import os
from PIL import Image, ImageDraw, ImageFont
D = 'figs/Supp_S5_Kang/output/'; DPI = 600
px = lambda mm: int(round(mm / 25.4 * DPI))
W = 183.0
ROWS = [  # (file, x mm, letter)
    [('FigS5a_tSNE.png', 0, 'a'), ('FigS5b_celltype_table.png', 131.0, 'b')],   # x of b printed by figS5_ab.R
    [('FigS5c_violin_IFNb_PGE2.png', 0, 'c'), ('FigS5d_violin_Interferon_Chemokine.png', 91.0, 'd')],
    [('FigS5e_scatter_IFNb_vs_PGE2.png', 0, 'e'), ('FigS5f_scatter_Interferon_vs_Chemokine.png', 91.0, 'f')],
]
GAP, TOP = 4.0, 3.5
imgs = [[(Image.open(D + f).convert('RGB'), x, l) for f, x, l in r] for r in ROWS]
heights = [max(im.size[1] for im, _, _ in r) for r in imgs]
Htot = px(TOP) + sum(heights) + px(GAP) * (len(imgs) - 1) + px(1)
canvas = Image.new('RGB', (px(W), Htot), 'white'); dr = ImageDraw.Draw(canvas)
font = None
for f in ['/System/Library/Fonts/Supplemental/Arial Bold.ttf', '/Library/Fonts/Arial Bold.ttf']:
    if os.path.exists(f): font = ImageFont.truetype(f, int(8 / 72 * DPI)); break
y = px(TOP)
for r, h in zip(imgs, heights):
    for im, x, l in r:
        canvas.paste(im, (px(x) + px(2.5), y))
        if l: dr.text((px(x), y - px(TOP) + px(0.3)), l, fill='black', font=font)
    y += h + px(GAP)
canvas.save(D + 'FigS5.png', dpi=(DPI, DPI))
print('wrote FigS5.png', canvas.size, 'px =', round(canvas.size[0] / DPI * 25.4, 1), 'x', round(canvas.size[1] / DPI * 25.4, 1), 'mm')
