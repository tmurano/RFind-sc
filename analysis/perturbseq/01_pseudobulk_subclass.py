# -*- coding: utf-8 -*-
"""Perturb-seq pseudo-bulk per gene_target for single subclasses of group 017 (CA3 + DG + CA2).

  PERTURB_DIR=<folder with the g017 shards> <python> analysis/perturbseq/01_pseudobulk_subclass.py

The g017 shards (Shi et al. 2026, predicted_group '017 CA3 CA2-FC DG Glut', QC-passing cells; 317 parquet files,
g017_000.parquet ... g017_316.parquet, extracted from the atlas by a separate script, not included in this repository;
see README) are READ, never modified, from PERTURB_DIR/*.parquet, or, when that folder holds none, from
PERTURB_DIR/analysis_260814/g017_shards/*.parquet (the folder layout of that extraction). PERTURB_DIR is
required. Outputs go to data/Perturbseq/ in this project:
  pb_g017_DG.npz    cells of subclass '037 DG Glut'  (dentate gyrus granule cells)
  pb_g017_CA3.npz   cells of subclass '017 CA3 Glut'
Aggregation: raw counts summed over single-guide cells per gene_target; nothing normalised, centred or filtered.
One pass over the shards builds both.
"""
import os, glob, time, numpy as np, duckdb

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))                        # project root
PDIR = os.environ.get('PERTURB_DIR', '')
if not PDIR:
    raise SystemExit('PERTURB_DIR is not set: point it to the folder with the g017 parquet shards (README, "How to run")')
SHARDS = (sorted(glob.glob(os.path.join(PDIR, '*.parquet'))) or       # the shard folder itself
          sorted(glob.glob(os.path.join(PDIR, 'analysis_260814', 'g017_shards', '*.parquet'))))   # folder layout of the extraction
if not SHARDS:
    raise SystemExit(f'no g017 shards found: neither {os.path.join(PDIR, "*.parquet")} nor '
                     f'{os.path.join(PDIR, "analysis_260814", "g017_shards", "*.parquet")} matches any file')
OUT = os.path.join(BASE, 'data', 'Perturbseq'); os.makedirs(OUT, exist_ok=True)
GM = os.path.join(OUT, 'gene_metadata.parquet')                      # gene table of the atlas (README, Input data)
SUBSETS = {'DG': ['037 DG Glut'], 'CA3': ['017 CA3 Glut']}
t0 = time.time(); say = lambda *a: print(f'[{time.time()-t0:6.1f}s]', *a, flush=True)
say(f'{len(SHARDS)} shards (read-only) from {os.path.dirname(SHARDS[0])}')

con = duckdb.connect()
con.execute("SET enable_progress_bar=false; SET memory_limit='4GB'; SET threads=4;")
gm = con.execute(f"SELECT gene_token_id, gene_name FROM read_parquet('{GM}') ORDER BY gene_token_id").df()
NG = int(gm.gene_token_id.max()) + 1
gene_names = gm.set_index('gene_token_id').gene_name.reindex(range(NG)).values

meta = con.execute(f"""SELECT predicted_subclass sc, gene_target, count(*) n, sum(num_rna_umi) u
    FROM read_parquet({SHARDS!r}) WHERE num_guides=1 GROUP BY 1,2""").df()
say('single-guide cells per subclass: ' +
    ', '.join(f'{k} {v:,}' for k, v in meta.groupby('sc').n.sum().sort_values(ascending=False).items()))

def subset_meta(scs):
    m = meta if scs is None else meta[meta.sc.isin(scs)]
    m = m.groupby('gene_target', as_index=False).agg(n=('n', 'sum'), u=('u', 'sum'))
    return m.sort_values(['n', 'gene_target'], ascending=[False, True]).reset_index(drop=True)

S = {}
for key, scs in SUBSETS.items():
    m = subset_meta(scs)
    S[key] = dict(scs=scs, meta=m, tix={t: i for i, t in enumerate(m.gene_target)},
                  M=np.zeros((NG, len(m)), np.float64), det=np.zeros(NG, np.int64))

for j, f in enumerate(SHARDS):
    r = con.execute(f"""SELECT predicted_subclass sc, gene_target, g, sum(e)::DOUBLE s, count(*) d FROM
        (SELECT predicted_subclass, gene_target, unnest(genes) g, unnest(expressions) e
         FROM read_parquet('{f}') WHERE num_guides=1) GROUP BY 1,2,3""").fetchnumpy()
    sc, tg, g, s, d = r['sc'], r['gene_target'], r['g'], r['s'], r['d']
    for key, z in S.items():
        keep = np.ones(len(sc), bool) if z['scs'] is None else np.isin(sc, z['scs'])
        if not keep.any():
            continue
        ci = np.array([z['tix'][t] for t in tg[keep]])
        np.add.at(z['M'], (g[keep], ci), s[keep])
        np.add.at(z['det'], g[keep], d[keep])
    if (j + 1) % 50 == 0:
        say(f'  {j+1}/{len(SHARDS)}')

kind = lambda t: 'NTC' if t == 'Non_target' else 'safe' if t.startswith('Safe_target') else 'KO'
for key, z in S.items():
    m = z['meta']; out = os.path.join(OUT, f'pb_g017_{key}.npz')
    np.savez_compressed(out, counts=z['M'].astype(np.float32), genes=gene_names,
                        targets=m.gene_target.values.astype(str), n_cells=m.n.values,
                        tot_umi=m.u.values, det_cells=z['det'])
    k = np.array([kind(t) for t in m.gene_target])
    nko = m.n.values[k == 'KO']
    say(f'wrote {os.path.relpath(out, BASE)}  {z["M"].shape[0]:,} genes x {z["M"].shape[1]:,} targets; '
        f'KO {int((k=="KO").sum()):,} (cells/KO median {int(np.median(nko))}, >=20 cells: {int((nko>=20).sum()):,}), '
        f'safe {int((k=="safe").sum())}, non {int(m.n.values[k=="NTC"].sum()):,} cells')

