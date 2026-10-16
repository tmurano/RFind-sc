"""Running Fisher kernel for the Python scripts of this repository: the same computation as running_fisher() in
analysis/running_fisher.R.

Copied from rfindsc-py 0.3.8, file src/rfindsc/running_fisher.py
(https://github.com/tmurano/rfindsc-py). Below the imports, this file is lines 12-166 of that file,
unchanged: running_fisher() and the functions it calls (_corners, _set_corners, _running_fisher_many,
_running_fisher_one_universe, _positions). Only this docstring and the imports are new:
log_upper_tail() comes from analysis/_hyper.py (the copy of src/rfindsc/_hyper.py of the same version) instead of the
package module rfindsc._hyper, and pandas and Mapping, used only by the parts not copied, are not imported. Not copied
(not used here): bid_running_fisher(), to_b2(), rfind_score_cell() and their helpers.

running_fisher(g1, g2, universe_size, clip=300.0, g2_ranked=True) = -log10 of the smallest hypergeometric upper-tail p,
P(X >= M), over every pair of cut-offs (the top G genes of g1 and the top K genes of g2, M = genes shared), universe
universe_size, clipped to [0, clip]; with g2_ranked=False, g2 is a gene set used whole (K = its number of distinct genes).

Used by analysis/perturbseq/03_axis_vs_ko.py, figs/Supp_S3_Perturbseq/analysis/figS3_rfind_allpairs.py and demo/run_demo.py
(BASE = the project root):
    sys.path.insert(0, os.path.join(BASE, 'analysis'))
    from running_fisher import running_fisher
"""
from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

from _hyper import log_upper_tail

_LN10 = math.log(10.0)


def _corners(pos2: np.ndarray):
    """Cut-off pairs to evaluate for one pair of lists.

    ``pos2[i]`` is the 0-based position in the second list of element ``i`` of
    the first list (-1 when absent; first occurrence, as R's ``match``).
    Returns ``None`` when the lists share nothing, else ``(M, G, K, top)``:
    the overlap ``M`` of the top ``G`` of list 1 and the top ``K`` of list 2 at
    the evaluated corners, and the largest ``G`` / ``K`` over all corners.

    The corners are the pairs of ranks of shared genes (R evaluates all
    ``m x m`` of them). A corner whose overlap equals that of the corner one
    shared gene earlier in either list has a larger ``G`` or ``K`` at the same
    ``M`` and so a p that is not smaller; those corners are skipped, which
    leaves the minimum unchanged.
    """
    sel = np.flatnonzero(pos2 >= 0)
    if sel.size == 0:
        return None
    r1 = sel + 1                             # ranks in list 1 (increasing)
    r2 = pos2[sel] + 1                       # their ranks in list 2
    K = np.sort(r2)
    M = np.cumsum(r2[:, None] <= K[None, :], axis=0)   # M[a, b]: top r1[a] x top K[b]
    keep = ((np.diff(M, axis=0, prepend=-1) != 0)
            & (np.diff(M, axis=1, prepend=-1) != 0))
    a, b = np.nonzero(keep)
    return M[a, b], r1[a], K[b], max(int(r1[-1]), int(K[-1]))


def _set_corners(pos2: np.ndarray, set_size: int):
    """Cut-offs to evaluate when the second list is an unordered gene set.

    The set is used whole (``K`` = ``set_size``, the number of distinct genes
    in it) and only the top-``G`` cut-offs of list 1 are scanned: for ``G`` at
    the rank of the ``j``-th shared gene the overlap is ``M = j`` (R:
    ``running_fisher(g2_ranked = FALSE)``). Same return value as ``_corners``.
    """
    sel = np.flatnonzero(pos2 >= 0)
    if sel.size == 0:
        return None
    set_size = int(set_size)
    M = np.arange(1, sel.size + 1, dtype=np.int64)
    return M, sel + 1, np.full(sel.size, set_size, dtype=np.int64), max(int(sel[-1]) + 1, set_size)


def _running_fisher_many(pos2_list, universe_size, clip, set_sizes=None) -> np.ndarray:
    """Running Fisher scores of many list pairs at once.

    Each element of ``pos2_list`` is a ``pos2`` array (see ``_corners``) or
    ``None`` for an empty list. ``universe_size`` is one size for all pairs or
    one per pair (each panel has its own pair universe since 0.3.5; pairs are
    grouped by size). ``set_sizes``, if given, has one entry per pair: ``None``
    when the second list is ranked, or the number of distinct genes of the
    second list when it is an unordered gene set (see ``_set_corners``). All
    tail probabilities of one universe size are evaluated in one vectorised
    call.
    """
    if np.ndim(universe_size) == 0:
        return _running_fisher_one_universe(pos2_list, universe_size, clip, set_sizes)
    Ps = np.rint(np.asarray(universe_size, dtype=np.float64)).astype(np.int64)
    if Ps.size != len(pos2_list):
        raise ValueError("`universe_size` must be a scalar or have one value per pair.")
    out = np.zeros(len(pos2_list))
    for P in np.unique(Ps):
        idx = np.flatnonzero(Ps == P)
        out[idx] = _running_fisher_one_universe(
            [pos2_list[t] for t in idx], int(P), clip,
            None if set_sizes is None else [set_sizes[t] for t in idx])
    return out


def _running_fisher_one_universe(pos2_list, universe_size, clip, set_sizes=None) -> np.ndarray:
    P = int(round(universe_size))
    out = np.zeros(len(pos2_list))
    Ms, Gs, Ks, owner = [], [], [], []
    for t, pos2 in enumerate(pos2_list):
        if pos2 is None or len(pos2) == 0:
            continue
        k_set = None if set_sizes is None else set_sizes[t]
        c = _corners(pos2) if k_set is None else _set_corners(pos2, k_set)
        if c is None:
            continue
        M, G, K, top = c
        if top > P:
            # A list is longer than the universe: R's phyper() returns NaN, and
            # a non-finite score is returned as `clip`.
            out[t] = clip
            continue
        Ms.append(M)
        Gs.append(G)
        Ks.append(K)
        owner.append(t)
    if owner:
        sizes = np.fromiter((m.size for m in Ms), dtype=np.int64, count=len(Ms))
        starts = np.concatenate(([0], np.cumsum(sizes)[:-1]))
        lp = log_upper_tail(np.concatenate(Ms), np.concatenate(Gs),
                            np.concatenate(Ks), P)
        s = -np.minimum.reduceat(lp, starts) / _LN10
        out[owner] = np.where(np.isfinite(s), np.minimum(s, clip), clip)
    return out


def _positions(g1: Sequence, g2: Sequence) -> np.ndarray:
    """0-based position in ``g2`` (first occurrence) of each element of ``g1``."""
    first = {}
    for i, g in enumerate(g2):
        first.setdefault(g, i)
    return np.fromiter((first.get(g, -1) for g in g1), dtype=np.int64,
                       count=len(g1))


def running_fisher(g1, g2, universe_size, clip=300.0, g2_ranked=True) -> float:
    """Running Fisher score of two ranked gene lists.

    ``-log10`` of the smallest hypergeometric upper-tail probability
    ``P(X >= M)`` over every pair of cut-offs: the top ``G`` genes of ``g1``
    and the top ``K`` genes of ``g2`` (``M`` = genes shared by the two
    prefixes; Kupershmidt et al. 2010, G x K form). The probabilities are
    intermediate quantities of the similarity score, computed on the log scale.
    Only cut-offs at ranks of shared genes are evaluated, which gives the same
    minimum as the full G x K scan.

    A membership gene set (e.g. an MSigDB Hallmark set) has no order. With
    ``g2_ranked=False`` it is used whole, as a biogroup in the original Running
    Fisher: ``K`` is fixed at the number of distinct genes in ``g2`` and only
    the top-``G`` cut-offs of ``g1`` are scanned, so the score does not depend
    on the order of ``g2`` (added in 0.3.4, as in RFindsc 0.3.4).

    Parameters
    ----------
    g1 : sequence of str
        Ranked gene symbols, most extreme first.
    g2 : sequence of str
        Ranked gene symbols, or a membership gene set when ``g2_ranked=False``.
    universe_size : int
        Size ``P`` of the gene universe (hypergeometric null).
    clip : float
        Score cap (default 300).
    g2_ranked : bool
        ``False`` when ``g2`` is an unordered gene set (default ``True``).

    Returns
    -------
    float
        The score, clipped to ``[0, clip]``; 0 when either list is empty or the
        lists share no gene.
    """
    g1 = list(g1)
    g2 = list(g2)
    if not g1 or not g2:
        return 0.0
    sizes = None if g2_ranked else [len(set(g2))]
    return float(_running_fisher_many([_positions(g1, g2)], universe_size, clip, sizes)[0])
