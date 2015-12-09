"""
Search functions for motifs around splice events
"""

import os
import sys
import rnamotifs2
from os.path import join as pjoin
import numpy as np
import pybio
import math
from fisher import pvalue
import random
from collections import Counter
import pickle
import operator

cache = {}

def dist(x):
    y = np.bincount(x)
    f = {}
    for i in range(1, len(y)): # 1: not interested in the frequency of 0
        f[i] = y[i]
    return f

def filter(v, thr):
    v = np.array(v)
    v[v<thr] = 0 # first filter <
    v[v>=thr] = 1 # then filter >=
    return v

def choose_h(distances, pth, debug=False):
    if debug:
        for r in distances:
            print r
    nd = [(dist, h) for (h, perc, dist) in distances if 3<=perc<=7] # perc between 3 and 7 %
    nd = sorted(nd, key=operator.itemgetter(0, 1)) # search for closest to 4% pth
    
    if len(nd)>0:
        return nd[0][1]
    else:
        return None

def v17(comps, genome, region="r1s", motif="YCAY", hw=15, pth=4, rfilter={}, step=0, base_motif="TTTT", base_h=None):
    nums = Counter() # frequencies of s/e/c, key = r1.s, r1.e, ...
    area_size = 0
    vectors = {}
    vectorst = {}
    seqs = {}
    consider_exons = set()
    for (eid, chr, strand, skip_start, in_start, in_stop, skip_stop, event_class) in rnamotifs2.data.data:
        if event_class not in [region[-1], "c"]:
            continue
        event_class = {"c":"c", "e":"t", "s":"t"}[event_class]

        nums[event_class+".all"] += 1

        if rfilter.get("%s.%s" % (event_class, eid), None)!=None:
            continue

        nums[event_class] += 1

        (a1_exon, a1_intron, a1_start, a1_stop), (a2_exon, a2_intron, a2_start, a2_stop), (a3_exon, a3_intron, a3_start, a3_stop), (a4_exon, a4_intron, a4_start, a4_stop) = rnamotifs2.sequence.coords(strand, skip_start, in_start, in_stop, skip_stop)
        a2_len = (a2_stop-a2_start+1)
        a3_len = (a3_stop-a3_start+1)
               
        a1_core, a2_core, a3_core, a4_core = rnamotifs2.sequence.sequence[eid]
        a2_core = coverage(a2_core, a2_len, hw, motif)
        a3_core = coverage(a3_core, a3_len, hw, motif)
        assert(len(a2_core)==a2_len)
        assert(len(a3_core)==a3_len)
        a2_core = [0] * (200-a2_intron) + a2_core + [0] * (50-a2_exon)
        a3_core = [0] * (50-a3_exon) + a3_core + [0] * (200-a3_intron)

        a1_seq, a2_seq, a3_seq, a4_seq = rnamotifs2.sequence.sequence[eid]
        a2_seq = a2_seq[hw:-hw]
        a3_seq = a3_seq[hw:-hw]
        a2_seq = "-" * (200-a2_intron) + a2_seq + "-" * (50-a2_exon)
        a3_seq = "-" * (50-a3_exon) + a3_seq + "-" * (200-a3_intron)

        assert(len(a2_core)==251)
        assert(len(a3_core)==251)

        if region in ["r1s", "r1e"]:
            r = a2_core[-95:-55+1]
            rseq = a2_seq[-95-15:-55+1+15]
        if region in ["r2s", "r2e"]:
            r = a2_core[-50:-20+1]+a3_core[20:50+1]
            rseq = a2_seq[-50-15:-20+1]+a3_seq[20:50+1+15]
        if region in ["r3s", "r3e"]:
            r = a3_core[60:100+1]
            rseq = a3_seq[60-15:100+1+15]
        vectors[(eid, event_class)] = r
        
        if step==0:
            consider_exons.add(eid)
        else:
            basemotif_hmin = max(4, base_h/2)
            if len(rseq)>30:
                # motif 0 is the new motif of the cluster (see rnamotifs2.motif.cluster)
                rseq = coverage(rseq, len(rseq)-30, hw, [base_motif, motif[0]], strict=True)
                #assert(len(rseq) in [31, 41, 62]) # r1 = 31, r2 = 62, r3 = 41
                if len(rseq)>0:
                    if max(rseq)>=basemotif_hmin:
                        consider_exons.add(eid)
                        if event_class=="t":
                            nums["t1"] += 1
                        else:
                            nums["c1"] += 1
    
    if step==0:
        print "%s.%s.%s: looking for h closest to threshold" % (comps, genome, motif)
        distances = []
        for h in range(4, 32):
            exons_present = 0
            for (eid, event_class), r in vectors.items():
                if event_class!="t":
                    continue
                if rfilter.get("%s.%s" % (event_class, eid), None)==None:
                    if max(r)>=h:
                        exons_present += 1
            perc = exons_present/float(nums["t"]) * 100
            distances.append((h, perc, abs(perc-pth)))
        choosen_h = choose_h(distances, pth, debug=False)
        
        # don't consider this motif if there is no 3%<=pth<=7%
        if choosen_h==None:
            return None

    if step>0:
        thr = max(2, int(0.04 * nums["t"]))
        thr50 = 0.5 * nums["t1"]
        print nums["t"], 0.04*nums["t"]
        distances = []
        for h in range(4, 32):
            exons_present = 0
            for (eid, event_class), r in vectors.items():
                if eid not in consider_exons:
                    continue
                if event_class!="t":
                    continue
                assert(rfilter.get("%s.%s" % (event_class, eid), None)==None)
                if max(r)>=h:
                    exons_present += 1
            d = abs(exons_present - thr)
            if exons_present<thr50:
                distances.append((d, h, exons_present, thr))
        distances.sort()
        if len(distances)>0:
            choosen_h = max(4, distances[0][1])
        else:
            return None
        
    print "%s.%s.%s: h=%s" % (comps, genome, motif, choosen_h)
    
    fvectors = {}
    # filter vectors and compute sums
    print "%s.%s.%s: filter" % (comps, genome, motif)
    for (eid, event_class), r in vectors.items():
        if eid in consider_exons:
            fvectors[(eid, event_class)] = filter(r, choosen_h)
    
    print "%s.%s.%s: sum and count" % (comps, genome, motif)
    # sum vectors
    vectors_sum = {}
    rcounts = Counter()
    present = {}
    for index, ((eid, event_class), r) in enumerate(fvectors.items()):
        rs = vectors_sum.get(event_class, [0]*len(r))
        rs = np.add(r, rs)
        vectors_sum[event_class] = rs
        
        # count if motif present in specific regions (r1, r2, r3) across s/e/c
        sum_r = 1 if sum(r)>0 else 0
        
        rcounts[event_class] += sum_r
        if sum_r==1:
            present["%s.%s" % (event_class, eid)] = 1
            if max(vectors[(eid, event_class)])>=14:
                rfilter["%s.%s" % (event_class, eid)] = 1
        
        for p in range(0, len(rnamotifs2.perm.ec_perm)):
            ec = rnamotifs2.data.data_class[rnamotifs2.perm.ec_perm[p][index]]
            rcounts["%s.p%s" % (ec, p)] += sum_r
    
    #print nums, rcounts
    return vectors_sum, rcounts, choosen_h, rfilter, nums, present

def v17_apa(comps, genome, region="r1s", motif="YCAY", hw=15, pth=4, rfilter={}, step=0, base_motif="TTTT", base_h=None):
    nums = Counter() # frequencies of s/e/c, key = r1.s, r1.e, ...
    area_size = 0
    vectors = {}
    vectorst = {}
    seqs = {}
    consider_exons = set()
    for (eid, chr, strand, pos, event_class) in rnamotifs2.data.data:

        if event_class not in [region[-1], "c"]:
            continue

        event_class = {"c":"c", "e":"t", "s":"t"}[event_class]
        
        nums[event_class+".all"] += 1

        if rfilter.get("%s.%s" % (event_class, eid), None)!=None:
            continue

        nums[event_class] += 1

        a1_core = rnamotifs2.sequence.sequence[eid]
        a1_core = coverage(a1_core, 201, hw, motif)

        a1_seq = rnamotifs2.sequence.sequence[eid]
        #a1_seq = a1_seq[hw:-hw]

        if region in ["r1s", "r1e"]:
            r = a1_core[0:60+1]
            rseq = a1_seq[0:60+30+1]
        if region in ["r2s", "r2e"]:
            r = a1_core[60:120+1]
            rseq = a1_seq[60:120+30+1]
        if region in ["r3s", "r3e"]:
            r = a1_core[120:180+1]
            rseq = a1_seq[120:180+30+1]
        vectors[(eid, event_class)] = r
        
        if step==0:
            consider_exons.add(eid)
        else:
            basemotif_hmin = max(4, base_h/2)
            if len(rseq)>30:
                # motif 0 is the new motif of the cluster (see rnamotifs2.motif.cluster)
                rseq = coverage(rseq, len(rseq)-30, hw, [base_motif, motif[0]], strict=True)
                if len(rseq)>0:
                    if max(rseq)>=basemotif_hmin:
                        consider_exons.add(eid)
                        if event_class=="t":
                            nums["t1"] += 1
                        else:
                            nums["c1"] += 1
    
    if step==0:
        print "%s.%s.%s: looking for h closest to threshold" % (comps, genome, motif)
        distances = []
        for h in range(4, 32):
            exons_present = 0
            for (eid, event_class), r in vectors.items():
                if event_class!="t":
                    continue
                if rfilter.get("%s.%s" % (event_class, eid), None)==None:
                    if max(r)>=h:
                        exons_present += 1
            perc = exons_present/float(nums["t"]) * 100
            distances.append((h, perc, abs(perc-pth)))
        #print motif, region
        #print distances
        #print
        choosen_h = choose_h(distances, pth, debug=False)
        
        # don't consider this motif if there is no 3%<=pth<=7%
        if choosen_h==None:
            return None

    if step>0:
        thr = max(2, int(0.04 * nums["t"]))
        thr50 = 0.5 * nums["t1"]
        print nums["t"], 0.04*nums["t"]
        distances = []
        for h in range(4, 32):
            exons_present = 0
            for (eid, event_class), r in vectors.items():
                if eid not in consider_exons:
                    continue
                if event_class!="t":
                    continue
                assert(rfilter.get("%s.%s" % (event_class, eid), None)==None)
                if max(r)>=h:
                    exons_present += 1
            d = abs(exons_present - thr)
            if exons_present<thr50:
                distances.append((d, h, exons_present, thr))
        distances.sort()
        if len(distances)>0:
            choosen_h = max(4, distances[0][1])
        else:
            return None
        
    print "%s.%s.%s: h=%s" % (comps, genome, motif, choosen_h)
    
    fvectors = {}
    # filter vectors and compute sums
    print "%s.%s.%s: filter" % (comps, genome, motif)
    for (eid, event_class), r in vectors.items():
        if eid in consider_exons:
            fvectors[(eid, event_class)] = filter(r, choosen_h)
    
    print "%s.%s.%s: sum and count" % (comps, genome, motif)
    # sum vectors
    vectors_sum = {}
    rcounts = Counter()
    present = {}
    for index, ((eid, event_class), r) in enumerate(fvectors.items()):
        rs = vectors_sum.get(event_class, [0]*len(r))
        rs = np.add(r, rs)
        vectors_sum[event_class] = rs
        
        # count if motif present in specific regions (r1, r2, r3) across s/e/c
        sum_r = 1 if sum(r)>0 else 0
        
        rcounts[event_class] += sum_r
        if sum_r==1:
            present["%s.%s" % (event_class, eid)] = 1
            if max(vectors[(eid, event_class)])>=14:
                rfilter["%s.%s" % (event_class, eid)] = 1
        
        for p in range(0, len(rnamotifs2.perm.ec_perm)):
            ec = rnamotifs2.data.data_class[rnamotifs2.perm.ec_perm[p][index]]
            rcounts["%s.p%s" % (ec, p)] += sum_r
    
    return vectors_sum, rcounts, choosen_h, rfilter, nums, present

def areas(comps, motif="YCAY", hw=15, h=None):
    nums = Counter() # frequencies of s/e/c, key = r1.s, r1.e, ...
    area = 0
    vectors = {}
    seqs = {}
    for (eid, chr, strand, skip_start, in_start, in_stop, skip_stop, event_class) in rnamotifs2.data.data:      
        (r1_exon, r1_intron, r1_start, r1_stop), (r2_exon, r2_intron, r2_start, r2_stop), (r3_exon, r3_intron, r3_start, r3_stop), (r4_exon, r4_intron, r4_start, r4_stop) = rnamotifs2.sequence.coords(strand, skip_start, in_start, in_stop, skip_stop)
        r1_len = (r1_stop-r1_start+1)
        r2_len = (r2_stop-r2_start+1)
        r3_len = (r3_stop-r3_start+1)
        r4_len = (r4_stop-r4_start+1)
        r1_core, r2_core, r3_core, r4_core = rnamotifs2.sequence.sequence[eid]
        r1_core = coverage(r1_core, r1_len, hw, motif)
        r2_core = coverage(r2_core, r2_len, hw, motif)
        r3_core = coverage(r3_core, r3_len, hw, motif)
        r4_core = coverage(r4_core, r4_len, hw, motif)
        assert(len(r1_core)==r1_len)
        assert(len(r2_core)==r2_len)
        assert(len(r3_core)==r3_len)
        assert(len(r4_core)==r4_len)
        r1_core = [0] * (50-r1_exon) + r1_core + [0] * (200-r1_intron)
        r2_core = [0] * (200-r2_intron) + r2_core + [0] * (50-r2_exon)
        r3_core = [0] * (50-r3_exon) + r3_core + [0] * (200-r3_intron)
        r4_core = [0] * (200-r4_intron) + r4_core + [0] * (50-r4_exon)

        assert(len(r1_core)==251)
        assert(len(r2_core)==251)
        assert(len(r3_core)==251)
        assert(len(r4_core)==251)
        vectors[(eid, event_class)] = (r1_core, r2_core, r3_core, r4_core)

    print "%s.%s: h=%s" % (comps, motif, h)
    
    stats = Counter()
    # filter vectors and compute sums
    print "%s.%s: filter" % (comps, motif)
    for (eid, event_class), (v1, v2, v3, v4) in vectors.items():
        v1 = filter(v1, h)
        v2 = filter(v2, h)
        v3 = filter(v3, h)
        v4 = filter(v4, h)
        vectors[(eid, event_class)] = (v1, v2, v3, v4)

        stats[event_class] += 1
        
        temp = v2[-95:-55+1]
        temp = 1 if sum(temp)>0 else 0
        stats["r1%s" % event_class] += temp
        
        temp = v2[-50:-20+1]+v3[20:50+1]
        temp = 1 if sum(temp)>0 else 0
        stats["r2%s" % event_class] += temp

        temp = v3[60:100+1]
        temp = 1 if sum(temp)>0 else 0
        stats["r3%s" % event_class] += temp
    
    print "%s.%s: sum and count" % (comps, motif)
    # sum vectors
    vectors_sum = {}
    for index, ((eid, event_class), (v1, v2, v3, v4)) in enumerate(vectors.items()):
        (v1s, v2s, v3s, v4s) = vectors_sum.get(event_class, ([0]*251, [0]*251, [0]*251, [0]*251))
        v1s = np.add(v1, v1s)
        v2s = np.add(v2, v2s)
        v3s = np.add(v3, v3s)
        v4s = np.add(v4, v4s)
        vectors_sum[event_class] = (v1s, v2s, v3s, v4s)
    return vectors_sum, h, stats

def areas_apa(comps, motif="YCAY", hw=15, h=None, pth=4):
    nums = Counter() # frequencies of s/e/c, key = r1.s, r1.e, ...
    area = 0
    vectors = {}
    seqs = {}
    for (eid, chr, strand, pos, event_class) in rnamotifs2.data.data:      
        r1_core = rnamotifs2.sequence.sequence[eid]
        r1_core = coverage(r1_core, 201, hw, motif)
        vectors[(eid, event_class)] = (r1_core)
        nums[event_class] += 1
    
    if h==None:
        print "%s.%s: looking for h closest to threshold" % (comps, motif)
        distances = []
        for h in range(4, 32):
            exons_present = 0
            for (eid, event_class), r in vectors.items():
                if event_class=="c":
                    continue
                if max(r)>=h:
                    exons_present += 1
            perc = exons_present/float(nums["e"]+nums["s"]) * 100
            distances.append((h, perc, abs(perc-pth)))
        h = choose_h(distances, pth, debug=False)

    print "%s.%s: h=%s" % (comps, motif, h)
    
    # filter vectors and compute sums
    stats = Counter()
    print "%s.%s: filter" % (comps, motif)
    for (eid, event_class), (v1) in vectors.items():
        v1 = filter(v1, h)
        vectors[(eid, event_class)] = (v1)
        
        stats[event_class] += 1
        
        temp = v1[0:60+1]
        temp = 1 if sum(temp)>0 else 0
        stats["r1%s" % event_class] += temp
        
        temp = v1[60:120+1]
        temp = 1 if sum(temp)>0 else 0
        stats["r2%s" % event_class] += temp

        temp = v1[120:180+1]
        temp = 1 if sum(temp)>0 else 0
        stats["r3%s" % event_class] += temp
    
    print "%s.%s: sum and count" % (comps, motif)
    # sum vectors
    vectors_sum = {}
    for index, ((eid, event_class), (v1)) in enumerate(vectors.items()):
        (v1s) = vectors_sum.get(event_class, ([0]*201))
        v1s = np.add(v1, v1s)
        vectors_sum[event_class] = (v1s)
    return vectors_sum, h, stats

def coverage(seq, seq_len, hw, motif, strict=False):
    if seq=="":
        return [0] * seq_len
    _, vector1 = pybio.sequence.search(seq, motif, strict=strict) # strict=True : require all motifs in list to be found
    vector2 = np.convolve(vector1, [1]*(hw*2+1), "same")
    vector = np.multiply(vector1, vector2)
    return list(vector[hw:-hw])

