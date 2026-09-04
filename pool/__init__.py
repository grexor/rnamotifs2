"""
Run motif jobs in-process across a multiprocessing pool.

The parent loads the event table, the genomic sequences and the permutation
index once; forked workers inherit them copy-on-write, so a motif job no longer
pays for a fresh interpreter + re-import + 8 MB pickle reload (~1.4 s each,
tens of thousands of times per run).
"""

import multiprocessing as mp
import os
import traceback

import rnamotifs2


def _run(task):
    kind = task[0]
    try:
        if kind == "motif":
            rnamotifs2.motifjob.run_motif(*task[1:], preloaded=True)
        elif kind == "cluster":
            rnamotifs2.motifjob.run_motif_cluster(*task[1:], preloaded=True)
        else:
            raise ValueError("unknown task kind %r" % kind)
        return (task, None)
    except Exception:
        return (task, traceback.format_exc())


def run(tasks, cores):
    """Execute a list of (kind, *args) tasks; raise if any worker failed."""
    tasks = list(tasks)
    if not tasks:
        return
    cores = max(1, min(int(cores), len(tasks)))
    if cores == 1 or os.environ.get("RNAMOTIFS2_NOPOOL") == "1":
        results = [_run(t) for t in tasks]
    else:
        ctx = mp.get_context("fork")
        with ctx.Pool(cores) as pool:
            results = pool.map(_run, tasks, chunksize=1)

    errors = [(t, tb) for (t, tb) in results if tb]
    if errors:
        for t, tb in errors[:5]:
            print("TASK FAILED: %r\n%s" % (t, tb))
        raise RuntimeError("%d/%d motif tasks failed" % (len(errors), len(tasks)))
