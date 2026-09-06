import os
import sys
import rnamotifs2
import bisect
import copy
import numpy as np
np.random.seed(42)

ec_perm = {}
ec_dist = {}
generation = 0  # bumped on every compute() call; lets fastsearch cache the
                # per-region permutation label matrix without going stale if
                # compute() is ever called again for the same event table

def compute(comps, genome, motif):
    global generation
    generation += 1
    ec_perm.clear()
    ec_dist.clear()
    print("%s.%s.%s: computing permutations" % (comps, genome, motif))
    if rnamotifs2.config.perms == 0:
        return
    lim_s = rnamotifs2.data.data_class.count("s")
    lim_e = lim_s + rnamotifs2.data.data_class.count("e")
    lim_c = lim_e + rnamotifs2.data.data_class.count("c")
    for p in range(0, rnamotifs2.config.perms):
        x = np.random.randint(low=0, high=len(rnamotifs2.data.data_class), size=len(rnamotifs2.data.data_class))
        ec_perm[p] = copy.copy(x)
        x.sort()
        pnum_s = bisect.bisect_left(x, lim_s)
        pnum_e = bisect.bisect_left(x, lim_e) - pnum_s
        pnum_c = bisect.bisect_left(x, lim_c) - bisect.bisect_left(x, lim_e)
        ec_dist[p] = {"s":pnum_s, "e":pnum_e, "c":pnum_c}
