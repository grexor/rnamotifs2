"""
Single-motif and single-cluster-motif jobs as importable functions.

These hold the logic that used to live only in ``bin/rnamotifs2.motif`` and
``bin/rnamotifs2.motif.cluster``. The pipeline now runs them in-process through
a multiprocessing pool (see ``rnamotifs2.pool``) instead of spawning one Python
interpreter per motif; the bin scripts remain as thin CLI wrappers.

Set ``RNAMOTIFS2_REFERENCE=1`` to force the original (slow) ``search.v17``
instead of the vectorised ``fastsearch.v17`` — used by the golden tests.
"""

import os
import pickle

import rnamotifs2


def _search_fn(sf):
    if sf == "v17" and os.environ.get("RNAMOTIFS2_REFERENCE") != "1":
        return rnamotifs2.fastsearch.v17
    return getattr(rnamotifs2.search, sf)


def _ensure_loaded(comps, genome, motif):
    if getattr(rnamotifs2.data, "data", None) in (None, []):
        rnamotifs2.data.read(comps)
    if getattr(rnamotifs2.sequence, "sequence", None) is None:
        rnamotifs2.sequence.load(comps)
    rnamotifs2.perm.compute(comps, genome, motif)


def _read_base_motif(comps, region, cn):
    rfile = os.path.join(rnamotifs2.path.comps_folder, comps, region, "results%s.tab" % cn)
    if not os.path.exists(rfile):
        return None, None
    with open(rfile, "rt") as f:
        f.readline()
        r = f.readline().replace("\n", "").replace("\r", "").split("\t")
    return r[0], int(r[1])


def _prev_tree_rfilter(comps, region, cn, pickle_folder):
    """rfilter carried over from the previous cluster tree (cn>0)."""
    if cn <= 0:
        return {}
    tree_path = os.path.join(rnamotifs2.path.comps_folder, comps, region, "tree%s.tab" % (cn - 1))
    with open(tree_path, "rt") as f:
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        m1 = m2 = None
        c = 0
        for line in f:
            row = line.replace("\r", "").replace("\n", "").split("\t")
            c += 1
            data = dict(zip(header, row))
            if float(data["fisher"]) > 0.01 and c == 1:
                m1 = data["motif"].split("_")
                m2 = data["cmotif"].split("_")
                break
            if float(data["fisher"]) > 0.01 and c > 1:
                break
            m1 = data["motif"].split("_")
            m2 = data["cmotif"].split("_")
    if m2 != ['']:
        pf = os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn - 1, m1[0], "_".join(sorted(m2))))
    else:
        pf = os.path.join(pickle_folder, "c%s.%s.pickle" % (cn - 1, m1[0]))
    _, _, _, rfilter, _, _, _ = pickle.load(open(pf, "rb"))
    return rfilter


def run_motif(comps, genome, region, motif, pth, cn, sf, preloaded=False):
    """One raw (unfiltered) motif / motif-cluster search + Fisher test."""
    motif = motif if isinstance(motif, list) else motif.split("_")
    pth, cn = float(pth), int(cn)
    pickle_folder = os.path.join(rnamotifs2.path.comps_folder, comps, region, "pickle")
    os.makedirs(pickle_folder, exist_ok=True)

    if not preloaded:
        _ensure_loaded(comps, genome, motif)

    print("%s.%s.%s: processing; pth=%s" % (comps, genome, motif, pth))
    rfilter = _prev_tree_rfilter(comps, region, cn, pickle_folder)

    results = _search_fn(sf)(comps, genome, region=region, motif=motif,
                             rfilter=rfilter, step=0, pth=pth)
    if results is None:
        return None
    area, rcounts, h, rfilter, nums, present = results
    test_result = rnamotifs2.compute.rtest(comps, genome, motif, rcounts, nums)
    out = os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(sorted(motif))))
    pickle.dump((area, test_result, h, rfilter, nums, rcounts, present),
                open(out, "wb"), protocol=2)
    return out


def run_motif_cluster(comps, genome, region, motif, cn, cmotif, pth, sf, preloaded=False):
    """One filtered cluster-growth search (motif added on top of cmotif)."""
    motif = motif if isinstance(motif, list) else motif.split("_")
    cmotif = cmotif if isinstance(cmotif, list) else cmotif.split("_")
    pth, cn = float(pth), int(cn)
    search_motif = motif + cmotif
    pickle_folder = os.path.join(rnamotifs2.path.comps_folder, comps, region, "pickle")
    os.makedirs(pickle_folder, exist_ok=True)

    if not preloaded:
        _ensure_loaded(comps, genome, motif)

    print("%s.%s.%s: processing; cmotif=%s" % (comps, genome, motif, cmotif))

    if len(cmotif) == 1:
        src = os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(cmotif)))
    else:
        src = os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, cmotif[0], "_".join(sorted(cmotif[1:]))))
    _, _, _, rfilter, _, _, _ = pickle.load(open(src, "rb"))

    step = len(cmotif)
    base_motif, base_h = _read_base_motif(comps, region, cn)
    print("base_motif=%s, motif_h=%s" % (base_motif, base_h))

    results = _search_fn(sf)(comps, genome, region=region, motif=search_motif,
                             rfilter=rfilter, step=step, pth=pth,
                             base_motif=base_motif, base_h=base_h)
    if results is None:
        return None
    area, rcounts, h, rfilter, nums, present = results
    test_result = rnamotifs2.compute.rtest(comps, genome, search_motif, rcounts, nums)
    out = os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, "_".join(motif), "_".join(sorted(cmotif))))
    pickle.dump((area, test_result, h, rfilter, nums, rcounts, present),
                open(out, "wb"), protocol=2)
    return out
