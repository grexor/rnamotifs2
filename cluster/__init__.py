import os
import sys
import rnamotifs2
import operator
import cPickle as pickle
from Queue import Queue
from threading import Thread
import pybio

max_steps = 4

def start_motifs(comps, region, cn):
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
    cmotif, h, fisher, ig, raw_ig, p_emp = motifs[0]
    motifs = [[x] for x in rnamotifs2.results.get_motifs(comps, region) if x!=cmotif]
    return motifs, cmotif.split("_"), fisher, ig, raw_ig, p_emp, h

def next_motifs(comps, region, cn, step):
    region_folder = os.path.join(rnamotifs2.path.comps_folder, comps, region)
    motifs = []
    filename = os.path.join(region_folder, "c%s.temp%s.tab" % (cn, step-1))
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

def next_cluster(comps, genome, region, cn, pth=0.5, sf="r"):
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    region_folder = os.path.join(comps_folder, region)
    pickle_folder = os.path.join(region_folder, "pickle")
    rnamotifs2.data.read(comps)
    rnamotifs2.sequence.save(comps, genome)

    for step in range(0, max_steps):
        if step==0:
            motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.start_motifs(comps, region, cn)
        else:
            motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.next_motifs(comps, region, cn, step)
            if fisher>rnamotifs2.data.cluster_stop_thr:
                draw(comps, genome, region, cn, steps=step)
                return # stop tree construction

        num_worker_threads = rnamotifs2.data.cores
        q = Queue()
        def worker():
            while True:
                task = q.get()
                os.system(task)
                q.task_done()
        tasks = []
        for motif in motifs:
            cluster = motif + cmotif
            pickle_file = os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, "_".join(sorted(motif)), "_".join(sorted(cmotif))))
            if not os.path.exists(pickle_file):
                command = "rnamotifs2.motif.cluster %s %s %s %s %s %s %s %s" % (comps, genome, region, "_".join(motif), cn, "_".join(cmotif), pth, sf)
                tasks.append(command)
            # raw results, without filtering
            pickle_file = os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(sorted(cluster))))
            if not os.path.exists(pickle_file):
                command = "rnamotifs2.motif %s %s %s %s %s %s %s" % (comps, genome, region, "_".join(cluster), pth, 0, sf)
                tasks.append(command)
        for i in range(num_worker_threads):
             t = Thread(target=worker)
             t.daemon = True
             t.start()
        for task in tasks:
            q.put(task)
        q.join()

        assemble(comps, genome, region, cn, step)

    draw(comps, genome, region, cn, steps=max_steps)

def assemble(comps, genome, region, cn, step):
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    region_folder = os.path.join(comps_folder, region)
    pickle_folder = os.path.join(region_folder, "pickle")

    if step==0:
        motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.start_motifs(comps, region, cn)
    else:
        motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.next_motifs(comps, region, cn, step)

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
            print "%s.%s: loading %s (%s)" % (comps, genome, k_string, index)
            _, test_results[k_string], h[k_string], _, _, _, _ = pickle.load(open(pickle_filename))
        pickle_filename = os.path.join(pickle_folder, "c0.%s.pickle" % k_string)
        if os.path.exists(pickle_filename):
            _, rtest_results[k_string], _, _, _, _, _ = pickle.load(open(pickle_filename))

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
    f = open(os.path.join(region_folder, "c%s.temp%s.tab" % (cn, step)), "wt")
    header = ["motif", "cmotif", "h", "fisher", "ig", "raw_ig", "p_emp"]
    f.write("\t".join(header)+"\n")
    for row in data:
        f.write("\t".join(str(x) for x in row)+"\n")
    f.close()

    # correct for fdr?
    rnamotifs2.data.read_config(comps)
    if rnamotifs2.data.use_FDR:
        pybio.utils.FDR_tab(os.path.join(region_folder, "c%s.temp%s.tab" % (cn, step)), "fisher")

def draw(comps, genome, region, cn, steps=4):
    cluster_region = "t"
    control_region = "c"

    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    region_folder = os.path.join(comps_folder, region)
    pickle_folder = os.path.join(region_folder, "pickle")

    # read data
    rnamotifs2.data.read(comps)
    rnamotifs2.sequence.load(comps)
    motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.start_motifs(comps, region, cn)

    fout = open(os.path.join(region_folder, "tree%s.tab" % cn), "wt")
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
            motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.start_motifs(comps, region, cn)
            motif, cmotif = cmotif, []
        else:
            motifs, cmotif, fisher, ig, raw_ig, p_emp, h = rnamotifs2.cluster.next_motifs(comps, region, cn, step)
            motif, cmotif = [cmotif[0]], cmotif[1:]
            cluster = motif+cmotif

        removed_region = 0
        removed_control = 0

        if cmotif!=[]:
            _, _, _, _, nums, rcounts, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, "_".join(sorted(motif)), "_".join(sorted(cmotif))))))
        else:
            _, _, _, _, nums, rcounts, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(sorted(motif))))))

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
            _, _, _, rfilter, _, _, _ = pickle.load(open(pickle_filename))
        else:
            rfilter = {}

        if len(cmotif)>1:
            _, _, _, rfilter, _, _, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, cmotif[0], "_".join(sorted(cmotif[1:]))))))
        elif len(cmotif)==1:
            _, _, _, rfilter, _, _, _ = pickle.load(open(os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(cmotif)))))

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
