# -*- coding: utf-8 -*-
"""Scan the full Shi et al. 2026 atlas once and build per-gene_target pseudo-bulk for several subclasses.

  <python> analysis/perturbseq/01b_scan_subclasses.py

The atlas is read over HTTPS from HuggingFace (perturbai/wholebrain_crispr_atlas; 317 parquet files, 81.3 GB;
URL list = data/Perturbseq/urls.txt, in the repository: one line per parquet file under data/ of the dataset, sorted by
name; the part files below are numbered in its order). One pass takes 2-4 h and is network bound
(parallel reads do not help). For every file, the QC-passing single-guide cells of the listed subclasses are
aggregated at once (raw counts summed per subclass x gene_target x gene) and written as two small part files,
so an interrupted run can simply be restarted (finished files are skipped). When every part exists, the parts
are assembled into data/Perturbseq/pb_<tag>.npz with the layout of 01_pseudobulk_subclass.py
(counts, genes, targets, n_cells, tot_umi, det_cells). Nothing is normalised or filtered here.

Subclasses (Fig4d): CA1 and three GABAergic subclasses (the De Donno GABA type is not split).
"""
import os, time, numpy as np, duckdb, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(BASE, 'data', 'Perturbseq')
## part files total ~10 GB (about 31 MB per atlas file) -> keep them outside cloud-synced folders via SCAN_PARTS
PARTS = os.environ.get('SCAN_PARTS', os.path.join(DATA, 'scan_parts')); os.makedirs(PARTS, exist_ok=True)
SUB = {'CA1': '016 CA1-ProS Glut', 'Pvalb': '052 Pvalb Gaba', 'Sst': '053 Sst Gaba', 'Lamp5': '049 Lamp5 Gaba'}
t0 = time.time(); say = lambda *a: print(f'[{(time.time()-t0)/60:6.1f} min]', *a, flush=True)

urls = [u.strip() for u in open(os.path.join(DATA, 'urls.txt')) if u.strip()]
con = duckdb.connect()
con.execute("INSTALL httpfs; LOAD httpfs; SET preserve_insertion_order=false; SET memory_limit='6GB'; "
            "SET threads=4; SET enable_progress_bar=false;")
insql = ', '.join(f"'{s}'" for s in SUB.values())
part = lambda kind, i: os.path.join(PARTS, f'{kind}_{i:03d}.parquet')

## ---------------- 1. scan (resumable) ----------------
done = sum(os.path.exists(part('expr', i)) and os.path.exists(part('meta', i)) for i in range(len(urls)))
say(f'{len(urls)} files, {done} already done; subclasses: ' + ', '.join(SUB.values()))
for i, u in enumerate(urls):
    if os.path.exists(part('expr', i)) and os.path.exists(part('meta', i)):
        continue
    for attempt in range(4):
        try:
            con.execute(f"""CREATE OR REPLACE TEMP TABLE t AS
                SELECT predicted_subclass sc, gene_target, num_rna_umi, genes, expressions
                FROM read_parquet('{u}')
                WHERE passes_qc AND num_guides = 1 AND predicted_subclass IN ({insql})""")
            con.execute(f"""COPY (SELECT sc, gene_target, count(*) n, sum(num_rna_umi) u FROM t GROUP BY 1, 2)
                TO '{part('meta', i)}.tmp' (FORMAT parquet)""")
            con.execute(f"""COPY (SELECT sc, gene_target, g, sum(e)::DOUBLE s, count(*) d FROM
                (SELECT sc, gene_target, unnest(genes) g, unnest(expressions) e FROM t) GROUP BY 1, 2, 3)
                TO '{part('expr', i)}.tmp' (FORMAT parquet)""")
            os.replace(part('meta', i) + '.tmp', part('meta', i)); os.replace(part('expr', i) + '.tmp', part('expr', i))
            break
        except Exception as e:
            say(f'  retry {attempt + 1} on file {i}: {e}'); time.sleep(10)
    else:
        say(f'  FAILED file {i}: {u}'); continue
    done += 1
    if done % 10 == 0 or done == len(urls):
        el = time.time() - t0
        say(f'  {done}/{len(urls)} files')

missing = [i for i in range(len(urls)) if not (os.path.exists(part('expr', i)) and os.path.exists(part('meta', i)))]
if missing:
    say(f'{len(missing)} files missing (e.g. {missing[:5]}); rerun this script to finish, then it assembles'); raise SystemExit(1)

## ---------------- 2. assemble one pseudo-bulk per subclass ----------------
gm = pd.read_parquet(os.path.join(DATA, 'gene_metadata.parquet'), columns=['gene_token_id', 'gene_name'])
NG = int(gm.gene_token_id.max()) + 1
gene_names = gm.set_index('gene_token_id').gene_name.reindex(range(NG)).values
meta = pd.concat([pd.read_parquet(part('meta', i)) for i in range(len(urls))])
kind = lambda t: 'NTC' if t == 'Non_target' else 'safe' if t.startswith('Safe_target') else 'KO'
for tag, sc in SUB.items():
    m = (meta[meta.sc == sc].groupby('gene_target', as_index=False).agg(n=('n', 'sum'), u=('u', 'sum'))
         .sort_values(['n', 'gene_target'], ascending=[False, True]).reset_index(drop=True))
    tix = {t: j for j, t in enumerate(m.gene_target)}
    M = np.zeros((NG, len(m))); det = np.zeros(NG, np.int64)
    for i in range(len(urls)):
        x = pd.read_parquet(part('expr', i)); x = x[x.sc == sc]
        if len(x):
            np.add.at(M, (x.g.values, x.gene_target.map(tix).values), x.s.values)
            np.add.at(det, x.g.values, x.d.values)
    out = os.path.join(DATA, f'pb_{tag}.npz')
    np.savez_compressed(out, counts=M.astype(np.float32), genes=gene_names, targets=m.gene_target.values.astype(str),
                        n_cells=m.n.values, tot_umi=m.u.values, det_cells=det)
    k = np.array([kind(t) for t in m.gene_target]); nko = m.n.values[k == 'KO']
    say(f'wrote {os.path.relpath(out, BASE)} ({sc}): {int(m.n.sum()):,} cells, KO {int((k=="KO").sum()):,} '
        f'(median {int(np.median(nko))} cells/KO, >=10 cells: {int((nko>=10).sum()):,}), safe {int((k=="safe").sum())}')
