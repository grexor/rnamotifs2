import os
import sys
import rnamotifs2
import operator
import pickle
import shutil
import pybio

max_steps = 4

def start_motifs(comps, region, cn, rank=0):
    """rank=0 is the single best base motif (what plain greedy search uses);
    beam search calls this with rank=1,2,... to also grow chains starting
    from the 2nd-best, 3rd-best, ... base motif."""
    region_folder = os.path.join(rnamotifs2.path.comps_folder, comps, region)
    motifs = []
    filename = os.path.join(region_folder, "results%s.tab" % cn)
    f = open(filename, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        motif = data["motif"]
        h = data["h"]
        fisher = float(data["fisher"])
        ig = float(data["ig"])
        p_emp = float(data["p_emp"])
        motifs.append((motif, h, fisher, ig, ig, p_emp))
        r = f.readline()
    f.close()
    motifs = sorted(motifs, key=operator.itemgetter(2))
    cmotif, h, fisher, ig, raw_ig, p_emp = motifs[rank]
    motifs = [[x] for x in rnamotifs2.results.get_motifs(comps, region) if x!=cmotif]
    return motifs, cmotif.split("_"), fisher, ig, raw_ig, p_emp, h

def next_motifs(comps, region, cn, step, tag=""):
    region_folder = os.path.join(rnamotifs2.path.comps_folder, comps, region)
    motifs = []
    filename = os.path.join(region_folder, "%sc%s.temp%s.tab" % (tag, cn, step-1))
    f = open(filename, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        cmotif = data["motif"].split("_") + data["cmotif"].split("_")
        h = int(data["h"])
        fisher = float(data["fisher"])
        ig = float(data["ig"])
        raw_ig = float(data["raw_ig"])
        p_emp = float(data["p_emp"])
        motifs.append((data["motif"].split("_"), cmotif, h, fisher, ig, raw_ig, p_emp))
        r = f.readline()
    f.close()
    motifs = sorted(motifs, key=operator.itemgetter(3))
    _, cmotif, h, fisher, ig, raw_ig, p_emp = motifs[0]
    del motifs[0]
    motifs = [x[0] for x in motifs]
    return motifs, cmotif, fisher, ig, raw_ig, p_emp, h

def next_cluster(comps, genome, region, cn, pth=0.5, sf="r", rank=0, tag=""):
    """Greedy growth of a single cluster chain, starting from the `rank`-th
    best base motif in results{cn}.tab (rank=0, the default, is what a plain
    run always used). `tag`, a filename prefix, keeps a beam-search caller's
    parallel chains (see next_cluster_beam) from overwriting each other's
    c{cn}.temp{step}.tab / tree{cn}.tab; motif-set-keyed pickles are shared
    across chains/beams as before."""
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    region_folder = os.path.join(comps_folder, region)
    pickle_folder = os.path.join(region_folder, "pickle")
    rnamotifs2.data.read(comps)
    rnamotifs2.sequence.save(comps, genome)
    rnamotifs2.sequence.load(comps)
    rnamotifs2.perm.compute(comps, genome, ["_base_"])

    for step in range(0, max_steps):
        if step==0:
            motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.start_motifs(comps, region, cn, rank=rank)
        else:
            motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.next_motifs(comps, region, cn, step, tag=tag)
            if fisher>rnamotifs2.data.cluster_stop_thr:
                draw(comps, genome, region, cn, steps=step, rank=rank, tag=tag)
                return # stop tree construction

        tasks = []
        for motif in motifs:
            cluster = motif + cmotif
            pickle_file = os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, "_".join(sorted(motif)), "_".join(sorted(cmotif))))
            if not os.path.exists(pickle_file):
                tasks.append(("cluster", comps, genome, region, "_".join(motif), cn, "_".join(cmotif), pth, sf))
            # raw results, without filtering
            pickle_file = os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(sorted(cluster))))
            if not os.path.exists(pickle_file):
                tasks.append(("motif", comps, genome, region, "_".join(cluster), pth, 0, sf))
        rnamotifs2.pool.run(tasks, rnamotifs2.data.cores)

        assemble(comps, genome, region, cn, step, rank=rank, tag=tag)

    draw(comps, genome, region, cn, steps=max_steps, rank=rank, tag=tag)

def next_cluster_beam(comps, genome, region, cn, pth=0.5, sf="r", beam_width=5):
    """Beam search over the starting (base) motif: instead of only ever
    growing a chain from the single best base motif (next_cluster's rank=0),
    independently grow one chain per each of the top `beam_width` base
    motifs in results{cn}.tab, using the exact same greedy step-by-step
    growth for each. The chain with the best final fisher is promoted to the
    canonical tree{cn}.tab (so everything downstream - get_motifs(), the
    next tree's rfilter carryover, bin/rnamotifs2.draw - is unaffected); every
    beam's result is also written to beams{cn}.tab for inspection.

    This diversifies the highest-leverage decision (choosing among ~320
    candidates at step 0) without the cost of re-ranking all chains against
    each other at every subsequent step (full recursive beam search) - a
    weaker individual base motif that turns out to combine better than the
    single best one is no longer invisible, but a promising *combination*
    that only becomes visible a few steps into some other chain still is.
    """
    region_folder = os.path.join(rnamotifs2.path.comps_folder, comps, region)
    results_file = os.path.join(region_folder, "results%s.tab" % cn)
    with open(results_file, "rt") as f:
        f.readline()
        rows = [l.replace("\r", "").replace("\n", "").split("\t") for l in f if l.strip()]

    width = min(beam_width, len(rows))
    print("%s.%s.%s: beam search, width=%d" % (comps, genome, region, width))

    beams = []
    for rank in range(width):
        if float(rows[rank][2]) > rnamotifs2.data.base_motif_thr:
            break  # results{cn}.tab is sorted by fisher; later ranks are only worse
        tag = "beam%d." % rank
        next_cluster(comps, genome, region, cn, pth=pth, sf=sf, rank=rank, tag=tag)
        tree_file = os.path.join(region_folder, "%stree%s.tab" % (tag, cn))
        with open(tree_file, "rt") as f:
            header = f.readline().replace("\r", "").replace("\n", "").split("\t")
            last = None
            for line in f:
                if line.strip():
                    last = dict(zip(header, line.replace("\r", "").replace("\n", "").split("\t")))
        beams.append((rank, tag, float(last["fisher"]), last))

    if not beams:
        return

    beams.sort(key=lambda x: x[2])  # best (lowest) fisher first
    winner_rank, winner_tag, winner_fisher, _ = beams[0]
    shutil.copyfile(os.path.join(region_folder, "%stree%s.tab" % (winner_tag, cn)),
                    os.path.join(region_folder, "tree%s.tab" % cn))

    with open(os.path.join(region_folder, "beams%s.tab" % cn), "wt") as f:
        f.write("rank\tfisher\tmotif\tcmotif\twinner\n")
        for rank, tag, fisher, last in beams:
            f.write("%d\t%s\t%s\t%s\t%s\n" % (rank, fisher, last["motif"], last["cmotif"], rank == winner_rank))
    print("%s.%s.%s: beam search winner = rank %d (fisher=%s) of %d beams" %
          (comps, genome, region, winner_rank, winner_fisher, len(beams)))

def _tree_row(comps, region, cn, step, motif, cmotif, fisher, ig, raw_ig, h):
    """One row of a tree{cn}.tab-style file: nums/rcounts/rfilter bookkeeping
    for a known (motif, cmotif) state at a given step. Pickle-only (keyed by
    motif set, same convention everywhere else in this module) - unlike
    draw()'s per-step loop this does not need to re-derive (motif, cmotif)
    from start_motifs/next_motifs, so it works for a chain reconstructed from
    next_cluster_beam_recursive's beam archive, which has no single
    c{cn}.temp{step}.tab per step (each step mixes candidates from several
    beams)."""
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    region_folder = os.path.join(comps_folder, region)
    pickle_folder = os.path.join(region_folder, "pickle")
    cluster_region, control_region = "t", "c"

    if cmotif != []:
        _, _, _, _, nums, rcounts, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, "_".join(sorted(motif)), "_".join(sorted(cmotif)))), "rb"))
    else:
        _, _, _, _, nums, rcounts, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(sorted(motif)))), "rb"))

    if cn > 0:
        f = open(os.path.join(region_folder, "tree%s.tab" % (cn - 1)), "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        m1, m2, c = None, None, 0
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            c += 1
            data = dict(zip(header, r))
            if float(data["fisher"]) > 0.01 and c == 1:
                m1 = data["motif"].split("_"); m2 = data["cmotif"].split("_"); break
            if float(data["fisher"]) > 0.01 and c > 1:
                break
            m1 = data["motif"].split("_"); m2 = data["cmotif"].split("_")
            r = f.readline()
        f.close()
        if m2 != ['']:
            pf = os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn - 1, m1[0], "_".join(sorted(m2))))
        else:
            pf = os.path.join(pickle_folder, "c%s.%s.pickle" % (cn - 1, m1[0]))
        _, _, _, rfilter, _, _, _ = pickle.load(open(pf, "rb"))
    else:
        rfilter = {}

    if len(cmotif) > 1:
        _, _, _, rfilter, _, _, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, cmotif[0], "_".join(sorted(cmotif[1:])))), "rb"))
    elif len(cmotif) == 1:
        _, _, _, rfilter, _, _, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(cmotif))), "rb"))

    removed_region = sum(1 for k in rfilter if k.startswith("t"))
    removed_control = sum(1 for k in rfilter if k.startswith("c"))

    return [step, fisher, ig, raw_ig, h, removed_region, removed_control,
            rcounts[cluster_region], rcounts[control_region],
            nums[cluster_region] - rcounts[cluster_region],
            nums[control_region] - rcounts[control_region],
            "compute_specificity", "_".join(motif), "_".join(cmotif)]

def next_cluster_beam_recursive(comps, genome, region, cn, pth=0.5, sf="r", beam_width=5):
    """Full recursive beam search: unlike next_cluster_beam (which only
    diversifies the choice of base motif and then grows each chain greedily
    on its own), this globally re-ranks every active beam's candidate
    extensions at *every* step and keeps only the best beam_width overall -
    so a promising combination that only becomes visible a few steps into
    some chain can still take over compute from a chain that looked better
    early on.

    A beam that loses that global competition, or that has no candidate left
    clearing cluster_stop_thr, stops being extended - but its current state
    is kept in an archive of completed candidate answers (mirroring
    next_cluster's own "stop and report wherever growth stopped" semantics,
    just applied per lineage instead of to a single chain). The single best
    archived state overall is promoted to tree{cn}.tab; every archived state
    is written to beams{cn}.tab.

    If use_FDR, the fisher column is corrected per step across ALL active
    beams' candidates pooled together (a larger, honest test count vs.
    next_cluster's single-chain-sized correction).

    Costs up to beam_width times more work per step than next_cluster_beam
    (which only pays that cost once, for step 0); with the vectorised search
    this is still minutes on comps/paper.bh-sized data, not hours.
    """
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    region_folder = os.path.join(comps_folder, region)
    pickle_folder = os.path.join(region_folder, "pickle")
    rnamotifs2.data.read(comps)
    rnamotifs2.sequence.save(comps, genome)
    rnamotifs2.sequence.load(comps)
    rnamotifs2.perm.compute(comps, genome, ["_base_"])

    all_motifs = [[x] for x in rnamotifs2.results.get_motifs(comps, region)]

    results_file = os.path.join(region_folder, "results%s.tab" % cn)
    with open(results_file, "rt") as f:
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        rows = [dict(zip(header, l.replace("\r", "").replace("\n", "").split("\t")))
                for l in f if l.strip()]
    rows.sort(key=lambda r: float(r["fisher"]))

    beams = []
    for row in rows[:beam_width]:
        fisher = float(row["fisher"])
        if fisher > rnamotifs2.data.base_motif_thr:
            break  # rows are pre-sorted by fisher; later ranks are only worse
        beams.append({"cmotif": row["motif"].split("_"), "fisher": fisher,
                      "ig": float(row["ig"]), "raw_ig": float(row["ig"]),
                      "p_emp": float(row["p_emp"]), "h": row["h"]})
    if not beams:
        return
    print("%s.%s.%s: recursive beam search, width=%d" % (comps, genome, region, len(beams)))

    finished = list(beams)  # every beam state that ever existed - a lineage
                            # that stops being extended is "done", not lost

    for step in range(1, max_steps):
        tasks = []
        pending = []  # (beam_index, candidate_motif)
        for bi, b in enumerate(beams):
            for m in all_motifs:
                if m[0] in b["cmotif"]:
                    continue
                pending.append((bi, m))
                cluster = m + b["cmotif"]
                filt = os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, "_".join(sorted(m)), "_".join(sorted(b["cmotif"]))))
                if not os.path.exists(filt):
                    tasks.append(("cluster", comps, genome, region, "_".join(m), cn, "_".join(b["cmotif"]), pth, sf))
                raw = os.path.join(pickle_folder, "c0.%s.pickle" % "_".join(sorted(cluster)))
                if not os.path.exists(raw):
                    tasks.append(("motif", comps, genome, region, "_".join(cluster), pth, 0, sf))
        rnamotifs2.pool.run(tasks, rnamotifs2.data.cores)

        scored = []  # [fisher, bi, m, ig, raw_ig, p_emp, h] (fisher mutable -> list)
        for bi, m in pending:
            b = beams[bi]
            k_string = "_".join(sorted(m + b["cmotif"]))
            filt = os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, "_".join(sorted(m)), "_".join(sorted(b["cmotif"]))))
            if not os.path.exists(filt):
                continue  # candidate skipped (e.g. h-search found nothing)
            _, test_result, h, _, _, _, _ = pickle.load(open(filt, "rb"))
            fisher, p_emp_list, ig = test_result
            p_emp = sum(1 if x <= fisher else 0 for x in p_emp_list)
            p_emp = (1 + p_emp) / float(1 + rnamotifs2.config.perms)
            raw_ig = ig
            raw = os.path.join(pickle_folder, "c0.%s.pickle" % k_string)
            if os.path.exists(raw):
                _, rtest_result, _, _, _, _, _ = pickle.load(open(raw, "rb"))
                _, _, raw_ig = rtest_result
            scored.append([fisher, bi, m, ig, raw_ig, p_emp, h])

        if not scored:
            break  # nothing, anywhere, could be extended this step

        if rnamotifs2.data.use_FDR:
            corrected = pybio.utils.FDR([s[0] for s in scored])
            for s, c in zip(scored, corrected):
                s[0] = c

        scored.sort(key=lambda s: s[0])

        new_beams = []
        seen_sets = set()
        for fisher, bi, m, ig, raw_ig, p_emp, h in scored:
            if fisher > rnamotifs2.data.cluster_stop_thr:
                break  # sorted ascending - nothing further qualifies either
            new_cmotif = m + beams[bi]["cmotif"]
            setkey = tuple(sorted(new_cmotif))
            if setkey in seen_sets:
                continue  # two beams reached the identical motif set; keep one
            seen_sets.add(setkey)
            nb = {"cmotif": new_cmotif, "fisher": fisher, "ig": ig,
                  "raw_ig": raw_ig, "p_emp": p_emp, "h": h}
            new_beams.append(nb)
            finished.append(nb)
            if len(new_beams) >= beam_width:
                break

        beams = new_beams
        if not beams:
            break

    winner = min(finished, key=lambda b: b["fisher"])

    by_set = {tuple(sorted(b["cmotif"])): b for b in finished}
    chain = winner["cmotif"]  # newest-first
    L = len(chain)
    rows_out = []
    for step in range(L):
        motif = [chain[L - 1 - step]]
        cmotif = chain[L - step:]
        state = by_set[tuple(sorted(motif + cmotif))]
        rows_out.append(_tree_row(comps, region, cn, step, motif, cmotif,
                                  state["fisher"], state["ig"], state["raw_ig"], state["h"]))

    num_region = rnamotifs2.data.dist[region[-1]]
    num_control = rnamotifs2.data.dist["c"]
    header = ["step", "fisher", "ig", "raw_ig", "h",
              "removed.t [%s]" % num_region, "removed.c [%s]" % num_control,
              "found.t", "found.c", "absent.t", "absent.c", "specificity",
              "motif", "cmotif"]
    with open(os.path.join(region_folder, "tree%s.tab" % cn), "wt") as f:
        f.write("\t".join(header) + "\n")
        for row in rows_out:
            row[11] = num_control * row[7] / max(1, float(num_region * row[8]))  # specificity
            f.write("\t".join(str(x) for x in row) + "\n")

    with open(os.path.join(region_folder, "beams%s.tab" % cn), "wt") as f:
        f.write("fisher\tmotif\tcmotif\twinner\n")
        for b in sorted(finished, key=lambda b: b["fisher"]):
            f.write("%s\t%s\t%s\t%s\n" % (b["fisher"], b["cmotif"][0],
                                          "_".join(b["cmotif"][1:]), b is winner))
    print("%s.%s.%s: recursive beam search winner fisher=%s, cluster=%s (%d states archived)" %
          (comps, genome, region, winner["fisher"], "+".join(winner["cmotif"]), len(finished)))

def assemble(comps, genome, region, cn, step, rank=0, tag=""):
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    region_folder = os.path.join(comps_folder, region)
    pickle_folder = os.path.join(region_folder, "pickle")

    if step==0:
        motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.start_motifs(comps, region, cn, rank=rank)
    else:
        motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.next_motifs(comps, region, cn, step, tag=tag)

    area = {}
    h = {}
    test_results = {}
    rtest_results = {}
    index = 0
    for motif in motifs:
        k_string = "_".join(sorted(motif + cmotif))
        pickle_filename = os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, "_".join(sorted(motif)), "_".join(sorted(cmotif))))
        if os.path.exists(pickle_filename):
            index += 1
            print("%s.%s: loading %s (%s)" % (comps, genome, k_string, index))
            _, test_results[k_string], h[k_string], _, _, _, _ = pickle.load(open(pickle_filename, "rb"))
        pickle_filename = os.path.join(pickle_folder, "c0.%s.pickle" % k_string)
        if os.path.exists(pickle_filename):
            _, rtest_results[k_string], _, _, _, _, _ = pickle.load(open(pickle_filename, "rb"))

    data = []
    index = 0
    for motif in motifs:
        k_string = "_".join(sorted(motif + cmotif))
        if test_results.get(k_string, None)==None:
            continue
        fisher, p_emp, ig = test_results[k_string]
        p_emp = [1 if x<=fisher else 0 for x in p_emp]
        p_emp = sum(p_emp)
        p_emp = (1+p_emp)/float(1+rnamotifs2.config.perms)
        if rtest_results.get(k_string, None)!=None:
            _, _, raw_ig = rtest_results[k_string]
        else:
            row_id = ""
        row = ["_".join(motif), "_".join(cmotif), h[k_string], fisher, ig, raw_ig, p_emp]
        data.append(row)

    data = sorted(data, key=operator.itemgetter(3))
    out_filename = os.path.join(region_folder, "%sc%s.temp%s.tab" % (tag, cn, step))
    f = open(out_filename, "wt")
    header = ["motif", "cmotif", "h", "fisher", "ig", "raw_ig", "p_emp"]
    f.write("\t".join(header)+"\n")
    for row in data:
        f.write("\t".join(str(x) for x in row)+"\n")
    f.close()

    # correct for fdr?
    rnamotifs2.data.read_config(comps)
    if rnamotifs2.data.use_FDR:
        pybio.utils.FDR_tab(out_filename, "fisher")

def draw(comps, genome, region, cn, steps=4, rank=0, tag=""):
    cluster_region = "t"
    control_region = "c"

    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    region_folder = os.path.join(comps_folder, region)
    pickle_folder = os.path.join(region_folder, "pickle")

    # read data
    rnamotifs2.data.read(comps)
    rnamotifs2.sequence.load(comps)
    motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.start_motifs(comps, region, cn, rank=rank)

    fout = open(os.path.join(region_folder, "%stree%s.tab" % (tag, cn)), "wt")
    header = ["step", "fisher", "ig", "raw_ig", "h"]
    num_region = rnamotifs2.data.dist[region[-1]] # dist contains only keys: s, e, c
    num_control = rnamotifs2.data.dist[control_region]
    header.append("removed.%s [%s]" % (cluster_region, num_region))
    header.append("removed.%s [%s]" % (control_region, num_control))
    header.append("found.t")
    header.append("found.c")
    header.append("absent.t")
    header.append("absent.c")
    header.append("specificity")
    header.append("motif")
    header.append("cmotif")
    fout.write("\t".join(header) + "\n")
    selection = []
    data_all = []

    follow = True
    for step in range(0, steps):
        if step==0:
            motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.start_motifs(comps, region, cn, rank=rank)
            motif, cmotif = cmotif, []
        else:
            motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.next_motifs(comps, region, cn, step, tag=tag)
            motif, cmotif = [cmotif[0]], cmotif[1:]
            cluster = motif+cmotif

        removed_region = 0
        removed_control = 0

        if cmotif!=[]:
            _, _, _, _, nums, rcounts, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, "_".join(sorted(motif)), "_".join(sorted(cmotif)))), "rb"))
        else:
            _, _, _, _, nums, rcounts, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(sorted(motif)))), "rb"))

        # find out where the previous tree was cut
        # and get filtered exons
        if cn>0:
            f = open(os.path.join(comps_folder, region, "tree%s.tab" % (cn-1)), "rt")
            header = f.readline().replace("\r", "").replace("\n", "").split("\t")
            r = f.readline()
            m1, m2, c = None, None, 0
            while r:
                r = r.replace("\r", "").replace("\n", "").split("\t")
                c+=1
                data = dict(zip(header, r))
                if float(data["fisher"])>0.01 and c==1:
                    m1 = data["motif"].split("_")
                    m2 = data["cmotif"].split("_")
                    break
                if float(data["fisher"])>0.01 and c>1:
                    break
                m1 = data["motif"].split("_")
                m2 = data["cmotif"].split("_")
                r = f.readline()
            f.close()

            if m2!=['']:
                pickle_filename = os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn-1, m1[0], "_".join(sorted(m2))))
            else:
                pickle_filename = os.path.join(pickle_folder, "c%s.%s.pickle" % (cn-1, m1[0]))
            _, _, _, rfilter, _, _, _ = pickle.load(open(pickle_filename, "rb"))
        else:
            rfilter = {}

        if len(cmotif)>1:
            _, _, _, rfilter, _, _, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, cmotif[0], "_".join(sorted(cmotif[1:])))), "rb"))
        elif len(cmotif)==1:
            _, _, _, rfilter, _, _, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(cmotif))), "rb"))

        for k in rfilter.keys():
            if k.startswith("t"):
                removed_region += 1
            if k.startswith("c"):
                removed_control += 1
        row = [step, fisher, ig, raw_ig, h, removed_region, removed_control, rcounts[cluster_region], rcounts[control_region], nums[cluster_region]-rcounts[cluster_region], nums[control_region]-rcounts[control_region], "compute_specificity", "_".join(motif), "_".join(cmotif)]
        data_all.append(row)

        # add * to show where the list was cut
        #if follow and fisher>0.01:
        #    if len(data_all)>1:
        #        data_all[-2][0] = "*" + str(data_all[-2][0])
        #        follow = False

    for row in data_all:
        dif_reg = row[7]
        dif_con = row[8]
        last_reg = row[7]
        last_con = row[8]
        row[11] = num_control*dif_reg / max(1, float(num_region*dif_con)) # specificity
        fout.write("\t".join(str(x) for x in row)+"\n")
    fout.close()

    """
    motif = "_".join(selection[-1][-2].split("_") + selection[-1][-1].split("_"))
    row = [cn, region, motif, selection[-1][0], selection[-1][4]]
    f.write("\t".join(str(x) for x in row) + "\n")
    f.close()

    # write removed regions
    motif = "_".join(sorted(motif.split("_")))
    filename = os.path.join(pickle_folder, "c%s.removed.pickle" % cluster_number)
    _, _, _, rfilter, _, _, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.pickle" % (cluster_number, motif))))
    # filter out not cluster regions (e.g. keep only r1.s and r1.c)
    for k in rfilter.keys():
        if not k.startswith(cluster_region) and not k.startswith(control_region):
            del rfilter[k]
    pickle.dump(rfilter, open(filename, "wb"), protocol=2)
    """

def redraw(comps, genome):
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)

    cluster_file = os.path.join(comps_folder, "clusters.tab")
    if os.path.exists(cluster_file):
        os.remove(cluster_file)

    for cn in range(1, 20):
        last_file = os.path.join(comps_folder, "c%s.%s.tab" % (cn, max_steps))
        if os.path.exists(last_file):
            draw(comps, genome, cn)
