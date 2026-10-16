#' Running Fisher enrichment between two ranked gene lists (all pairs of cut-offs)
#'
#' For every pair of cut-offs — the top G genes of `g1` and the top K genes of `g2` — the overlap M is tested with
#' Fisher's exact test (hypergeometric upper tail, P(X >= M), universe `universe_size`), and the score of the pair of
#' lists is -log10 of the smallest p over all G x K combinations (Supplementary Fig. S2).
#' The minimum always lies where G and K are ranks of overlapping genes (with M fixed, p only rises as G or K grows),
#' so only those m x m corners are evaluated (m = number of genes the two lists share); the result is identical to the
#' full G x K scan. p is computed on the log scale, so no underflow.
#'
#' @param g1 Character vector: ranked gene symbols (first list, most extreme first).
#' @param g2 Character vector: ranked gene symbols (second list), or a membership gene set when `g2_ranked = FALSE`.
#' @param universe_size Integer: total gene universe `P` used for the hypergeometric null.
#' @param clip Numeric: score cap (default 300).
#' @param g2_ranked Logical: `FALSE` when `g2` is an unordered gene set; then K is fixed at the set size and only the
#'   cut-offs of `g1` are scanned.
#'
#' @return Numeric: -log10(min p), clipped to `[0, clip]`. Returns 0 when the lists share no gene.
#'
#' @examples
#' set.seed(1)
#' pool <- paste0("g", 1:500)
#' g1 <- sample(pool, 50)
#' g2 <- c(g1[1:10], sample(setdiff(pool, g1), 40))
#' running_fisher(g1, g2, universe_size = 500)
#'
#' @export
running_fisher <- function(g1, g2, universe_size, clip = 300, g2_ranked = TRUE) {
  if (!length(g1) || !length(g2)) return(0)
  r1 <- which(g1 %in% g2)                                   # ranks in g1 of the shared genes (increasing)
  m <- length(r1)
  if (!m) return(0)
  if (!g2_ranked) {
    ## g2 is a membership gene set (e.g. MSigDB Hallmark): it has no order, so it is used whole (K = set size) and only
    ## the top-G cut-offs of the ranked list g1 are scanned, as for a biogroup in the original Running Fisher.
    K <- length(unique(g2))
    lp <- stats::phyper(seq_len(m) - 1, K, universe_size - K, r1, lower.tail = FALSE, log.p = TRUE)
    s <- -min(lp) / log(10)
    if (!is.finite(s)) return(clip)
    return(min(s, clip))
  }
  r2 <- match(g1[r1], g2)                                   # their ranks in g2
  K <- sort(r2)
  M <- crossprod(outer(r1, r1, "<="), outer(r2, K, "<="))   # M[a, b] = #{shared genes with rank1 <= r1[a] and rank2 <= K[b]}
  Kv <- rep(K, each = m)
  lp <- stats::phyper(as.vector(M) - 1, Kv, universe_size - Kv, rep(r1, times = m), lower.tail = FALSE, log.p = TRUE)
  s <- -min(lp) / log(10)
  if (!is.finite(s)) return(clip)
  min(s, clip)
}

#' Gene universe shared by a target and a reference: the background P of Running Fisher
#'
#' One rule for every dataset (one method, so one code path for all datasets): the universe of a
#' target-reference pair is the target's gene universe intersected with the genes the reference could contain, i.e. the
#' genes the reference study reported (finite, non-zero change; plus any expression floor used to build the list).
#' A membership gene set (e.g. MSigDB Hallmark) can contain any gene, so `eligible = NULL` gives the target universe.
#' Both lists are restricted to this universe before they are cut, and P = its size. This is the background of the
#' original Running Fisher (the overlap of the two platforms; Illumina Technical Note) and the one of the Perturb-seq
#' comparison (analysis/perturbseq/03_axis_vs_ko.py: genes present in both datasets).
#'
#' @param target_genes Character: the target's gene universe.
#' @param eligible Character or NULL: genes the reference could contain; NULL for a membership gene set.
#' @return Character: lower-case unique gene symbols of the pair universe (P = its length).
#' @export
pair_universe <- function(target_genes, eligible = NULL) {
  tg <- unique(tolower(target_genes))
  if (is.null(eligible)) tg else intersect(tg, unique(tolower(eligible)))
}

#' Cut a ranked target list to its top `n_top` genes inside a pair universe
#'
#' @param idx Integer: gene indices of one ranked list (most extreme first; NA-padded), from a deep per-cell ranking.
#' @param in_u Logical over the target genes: TRUE if the gene is in the pair universe.
#' @param n_top Integer: list length after restriction (500 in this study).
#' @return Integer: at most `n_top` indices, in rank order.
#' @keywords internal
cut_in_universe <- function(idx, in_u, n_top) {
  idx <- idx[!is.na(idx)]
  utils::head(idx[in_u[idx]], n_top)
}

#' Four-term RFind-sc score for a single cell against a single panel
#'
#' Composes the signed RFind-sc score:
#' \deqn{S = RF(b1u, b2u) - RF(b1u, b2d) - RF(b1d, b2u) + RF(b1d, b2d)}
#' where `b1u/b1d` are the per-cell UP/DOWN signatures and `b2u/b2d` are the panel UP/DOWN gene lists (an unsigned
#' gene set is passed as `up` with an empty `down`). Positive = cell aligned with the panel's case direction,
#' negative = aligned with control.
#'
#' @param b1_up,b1_down Character vectors: per-cell UP/DOWN signatures (ranked).
#' @param b2 Panel object in internal form (list with `up`, `down`, `n_up`, `n_down`, and `ranked = FALSE` for a membership
#'   gene set, which is then used whole). Use [to_b2()] to convert a bioset entry.
#' @param universe_size Integer: total gene universe `P`.
#' @param clip Score magnitude cap (symmetric). Default 300.
#'
#' @return Numeric: signed RFind-sc score, clipped to `[-clip, clip]`.
#' @keywords internal
rfind_score_cell <- function(b1_up, b1_down, b2, universe_size, clip = 300) {
  rk <- !isFALSE(b2$ranked)                                 # b2$ranked = FALSE for a membership gene set
  s <- running_fisher(b1_up,   b2$up,   universe_size, clip = clip, g2_ranked = rk) -
       running_fisher(b1_up,   b2$down, universe_size, clip = clip, g2_ranked = rk) -
       running_fisher(b1_down, b2$up,   universe_size, clip = clip, g2_ranked = rk) +
       running_fisher(b1_down, b2$down, universe_size, clip = clip, g2_ranked = rk)
  max(-clip, min(s, clip))
}

#' Convert a bioset panel entry into internal scoring form
#'
#' @param panel A list with `$up` and `$down` data.frames each containing a
#'   `gene` column (as in [default_biosets()]).
#'
#' @return A list with lowercase gene vectors and sizes, ready for
#'   [rfind_score_cell()].
#' @keywords internal
to_b2 <- function(panel) {
  up_sym <- tolower(trimws(panel$up$gene))
  dn_sym <- tolower(trimws(panel$down$gene))
  list(up = up_sym, down = dn_sym,
       n_up = length(up_sym), n_down = length(dn_sym))
}
