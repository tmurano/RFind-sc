"""Mouse-level bootstrap for the Fig. 5 scripts in Python (fig5_best_length.py, fig5_best_length_8tasks.py),
with the rule of boot_task in fig5_bench_redo.R.
Mouse = the sample of a Hammond cell (the cell name before its barcode, e.g. P100_M_SALINE_A3) or the donor of a De Donno cell
(donor_id of data/DeDonno2026/output/cell_meta_qc.csv; De Donno cell names have no '_').
A resample draws mice with replacement and takes all cells of every drawn mouse: within the positive and the negative group when every
mouse of the task is in one group, else from all mice of the task (the tasks defined by a cluster, whose mice have cells in both
groups). Tasks with the same mice and grouping share the draws: the draws of a set of mice are made once, from the generator given,
when a task first asks for them, so the order in which the tasks ask fixes the draws.
"""
import numpy as np
import pandas as pd

DONOR_FILE = 'data/DeDonno2026/output/cell_meta_qc.csv'
_DONOR = None


def mouse_of(cells):
    """the mouse of each cell (array of str)"""
    global _DONOR
    c = pd.Series(np.asarray(cells, dtype=str))
    m = c.str.rsplit('_', n=1).str[0]
    dd = ~c.str.contains('_', regex=False)
    if dd.any():
        if _DONOR is None:
            _DONOR = pd.read_csv(DONOR_FILE, usecols=['cell_id', 'donor_id']).set_index('cell_id').donor_id
        m[dd] = c[dd].map(_DONOR)
    assert not m.isna().any(), 'cells without a mouse'
    return m.to_numpy(dtype=str)


class MouseDraws:
    """B resamples of mice per set of mice, shared by the tasks with the same mice and grouping"""

    def __init__(self, rng, B):
        self.rng, self.B, self.draws = rng, B, {}

    def cells(self, mouse, pos):
        """yields the B resamples of a task as arrays of cell indices"""
        mouse, pos = np.asarray(mouse, dtype=str), np.asarray(pos, dtype=bool)
        mice = np.unique(mouse)
        idx = {k: np.flatnonzero(mouse == k) for k in mice}
        if all(pos[i].all() or not pos[i].any() for i in idx.values()):          # every mouse in one group: draw within the groups
            mp, mn = np.unique(mouse[pos]), np.unique(mouse[~pos])
            key = ('within groups', tuple(mp), tuple(mn))
            if key not in self.draws:
                self.draws[key] = [np.concatenate([self.rng.choice(mp, len(mp)), self.rng.choice(mn, len(mn))]) for _ in range(self.B)]
        else:                                                                       # mice in both groups: draw from all mice
            key = ('all mice', tuple(mice))
            if key not in self.draws:
                self.draws[key] = [self.rng.choice(mice, len(mice)) for _ in range(self.B)]
        for d in self.draws[key]:
            i = np.concatenate([idx[k] for k in d])
            assert pos[i].any() and not pos[i].all(), 'a resample without one of the groups'
            yield i
