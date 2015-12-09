"""
RNAmotifs2
"""

import rnamotifs2
import os
import pybio
import data
import search
import compute
import draw
import path
import config
import perm
import sequence
import cluster
import results
from Queue import Queue
from threading import Thread
import random
import operator
import pickle
import shutil
random.seed(42)

path.init()

def get_motifs():
    return pybio.genomes.make_motifs_nr(4)+pybio.genomes.make_motifs_nr(3)

def start(comps, region, cn, pth):
    rnamotifs2.data.read_config(comps)
    if rnamotifs2.data.data_type=="apa":
        sf = "v17_apa"
    if rnamotifs2.data.data_type=="splice":
        sf = "v17"
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    region_folder = os.path.join(comps_folder, region)
    comps_filename = os.path.join(comps_folder, "%s.tab" % comps)
    pickle_folder = os.path.join(region_folder, "pickle")

    if os.path.exists(region_folder):
        shutil.rmtree(region_folder)
    os.makedirs(os.path.join(region_folder))

    if not os.path.exists(os.path.join(pickle_folder)):
        os.makedirs(os.path.join(pickle_folder))

    start_cluster(comps, rnamotifs2.data.genome, region, cn, pth, sf)

def start_cluster(comps, genome, region, cn, pth, sf):
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    region_folder = os.path.join(comps_folder, region)
    comps_filename = os.path.join(comps_folder, "%s.tab" % comps)
    pickle_folder = os.path.join(region_folder, "pickle")

    # read data
    rnamotifs2.data.read(comps)
    motifs = rnamotifs2.results.get_motifs(comps, region)

    # make sequences
    print "%s.%s: saving sequences to pickle" % (comps, genome)
    rnamotifs2.sequence.save(comps, genome)

    num_worker_threads = 40
    q = Queue()
    def worker():
        while True:
            task = q.get()
            os.system(task)
            q.task_done()
    tasks = []
    for motif in motifs:
        pickle_file = os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(sorted(motif.split("_")))))
        if not os.path.exists(pickle_file):
            command = "rnamotifs2.motif %s %s %s %s %s %s %s" % (comps, genome, region, "_".join(motif.split("_")), pth, cn, sf)
            print "COMMAND=%s" % command
            tasks.append(command)
    for i in range(num_worker_threads):
         t = Thread(target=worker)
         t.daemon = True
         t.start()
    for task in tasks:
        q.put(task)
    q.join()

    base_motif_fisher = assemble_results(comps, genome, region, cn)
    if base_motif_fisher<=0.01:
        continue_cluster(comps, genome, region, cn, pth, sf)
    return base_motif_fisher

def continue_cluster(comps, genome, region, cn, pth, sf):
    rnamotifs2.cluster.next_cluster(comps, genome, region, cn, pth=pth, sf=sf)
    # try to start new cluster, if base_motif will have fisher < thr, it will stop processing
    start_cluster(comps, genome, region, cn+1, pth, sf)

def assemble_results(comps, genome, region, cn):

    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    region_folder = os.path.join(comps_folder, region)
    comps_filename = os.path.join(comps_folder, "%s.tab" % comps)
    pickle_folder = os.path.join(region_folder, "pickle")

    motifs = rnamotifs2.get_motifs()

    test_results = {}
    rtest_results = {}
    h = {}
    index = 0
    for motif in motifs:
        pickle_filename = os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(sorted(motif.split("_")))))
        if not os.path.exists(pickle_filename):
            continue
        index += 1
        print "%s.%s: loading %s (%s)" % (comps, genome, "_".join(motif.split("_")), index)
        _, test_results[motif], h[motif], _, _, _, _ = pickle.load(open(pickle_filename))

    # FDR
    """
    print "%s.%s: correct p-values" % (comps, genome)
    for rt in ["r1", "r2", "r3"]:
        for event_class in ["s", "e"]:
            # fdr correct p-value
            temp = []
            for motif in motifs_all:
                motif = "_".join(motif)
                temp.append(test_results[motif]["%s.%s" % (rt, event_class)][0])
            temp = rnamotifs2.fdr(temp)
            for index, motif in enumerate(motifs_all):
                motif = "_".join(motif)
                test_results[motif]["%s.%s" % (rt, event_class)][0] = temp[index]

            # fdr correct bootstrapped p-values
            temp = []
            for p in range(0, rnamotifs2.config.perms):
                temp.append([])
            for p in range(0, rnamotifs2.config.perms):
                for motif in motifs_all:
                    motif = "_".join(motif)
                    temp[p].append(test_results[motif]["%s.%s" % (rt, event_class)][1][p])
            for p in range(0, rnamotifs2.config.perms):
                temp[p] = rnamotifs2.fdr(temp[p])
                for index, motif in enumerate(motifs_all):
                    motif = "_".join(motif)
                    test_results[motif]["%s.%s" % (rt, event_class)][1][p] = temp[p][index]
    """

    data = []
    for motif in motifs:
        # some motifs were not considered
        if test_results.get(motif, None)==None:
            continue
        row = [motif, h[motif]]
        fisher, p_emp, ig = test_results[motif]
        row.append(fisher)
        row.append(ig)
        row.append(ig)
        p_emp = [1 if x<=fisher else 0 for x in p_emp]
        p_emp = sum(p_emp)
        p_emp = (1+p_emp)/float(1+rnamotifs2.config.perms)
        row.append(p_emp)
        data.append(row)

    data = sorted(data, key = operator.itemgetter(2))

    f = open(os.path.join(region_folder, "results%s.tab" % cn), "wt")
    header = ["motif", "h", "fisher", "ig", "raw_ig", "p_emp"]
    f.write("\t".join(header)+"\n")
    for row in data:
        f.write("\t".join(str(x) for x in row)+"\n")
    f.close()

    base_motif_fisher, _, _ = test_results[data[0][0]] # get fisher value of top (base) motif
    return base_motif_fisher

def fdr(pvalues, correction_type = "Benjamini-Hochberg"):
    """
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in xrange(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues
