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
import pickle

PAS_hexamers = [
    'AATAAA',
    'ATTAAA',
    'AGTAAA',
    'TATAAA',
    'CATAAA',
    'GATAAA',
    'AATATA',
    'AATACA',
    'AATAGA',
    'ACTAAA',
    'AAGAAA',
    'AATGAA'
]

def coords(strand, skip_start, in_start, in_stop, skip_stop, max_intron=200, max_exon=50):
    upintron_len = in_start-skip_start+1
    exon_len = in_stop - in_start + 1
    downintron_len = skip_stop - in_stop + 1
    if strand=="+":
        r1_exon = max_exon
        r1_intron = min(max_intron, upintron_len/2)
        r1_start = skip_start-r1_exon
        r1_stop = skip_start + r1_intron

        r2_intron = min(max_intron, upintron_len/2)
        r2_exon = min(max_exon, exon_len/2)
        r2_start = in_start - r1_intron
        r2_stop = in_start + r2_exon

        r3_intron = min(max_intron, downintron_len/2)
        r3_exon = min(50, exon_len/2)
        r3_start = in_stop - r3_exon
        r3_stop = in_stop + r3_intron

        r4_exon = max_exon
        r4_intron = min(max_intron, downintron_len/2)
        r4_start = skip_stop - r4_intron
        r4_stop = skip_stop + r4_exon
    else:
        r1_exon = max_exon
        r1_intron = min(max_intron, downintron_len/2)
        r1_start = skip_stop - r1_intron
        r1_stop = skip_stop + r1_exon

        r2_intron = min(max_intron, downintron_len/2)
        r2_exon = min(max_exon, exon_len/2)
        r2_start = in_stop - r2_exon
        r2_stop = in_stop + r2_intron

        r3_exon = min(max_exon, exon_len/2)
        r3_intron = min(max_intron, upintron_len/2)
        r3_start = in_start - r3_intron
        r3_stop = in_start + r3_exon

        r4_exon = max_exon
        r4_intron = min(max_intron, upintron_len/2)
        r4_start = skip_start - r4_exon
        r4_stop = skip_start + r4_intron
    return (r1_exon, r1_intron, r1_start, r1_stop), (r2_exon, r2_intron, r2_start, r2_stop), (r3_exon, r3_intron, r3_start, r3_stop), (r4_exon, r4_intron, r4_start, r4_stop)

def load(comps):
    print "%s: loading sequences" % comps
    pickle_folder = os.path.join(rnamotifs2.path.comps_folder, comps, "pickle")
    pickle_filename = os.path.join(pickle_folder, "sequence.pickle")
    rnamotifs2.sequence.sequence = pickle.load(open(pickle_filename))

def save(comps, genome, hw=15):
    if rnamotifs2.data.data_type=="apa":
        save_apa(comps, genome, hw=hw)
    else:
        save_splice(comps, genome, hw=hw)

def save_splice(comps, genome, hw=15):
    pickle_folder = os.path.join(rnamotifs2.path.comps_folder, comps, "pickle")
    if not os.path.exists(os.path.join(pickle_folder)):
        os.makedirs(os.path.join(pickle_folder))
    pickle_filename = os.path.join(pickle_folder, "sequence.pickle")
    if os.path.exists(pickle_filename):
        return
    sequence = {}
    for (eid, chr, strand, skip_start, in_start, in_stop, skip_stop, event_class) in rnamotifs2.data.data:
        upintron_len = in_start-skip_start+1
        exon_len = in_stop - in_start + 1
        downintron_len = skip_stop - in_stop + 1
        (r1_exon, r1_intron, r1_start, r1_stop), (r2_exon, r2_intron, r2_start, r2_stop), (r3_exon, r3_intron, r3_start, r3_stop), (r4_exon, r4_intron, r4_start, r4_stop) = coords(strand, skip_start, in_start, in_stop, skip_stop)
        seq1 = pybio.genomes.seq_direct(genome, chr.replace("chr", ""), strand, r1_start-hw, r1_stop+hw)
        seq2 = pybio.genomes.seq_direct(genome, chr.replace("chr", ""), strand, r2_start-hw, r2_stop+hw)
        seq3 = pybio.genomes.seq_direct(genome, chr.replace("chr", ""), strand, r3_start-hw, r3_stop+hw)
        seq4 = pybio.genomes.seq_direct(genome, chr.replace("chr", ""), strand, r4_start-hw, r4_stop+hw)

        seq1_len = len(seq1)
        seq2_len = len(seq2)
        seq3_len = len(seq3)
        seq4_len = len(seq4)

        # this can happen because of chrUn_*
        if seq1_len==0 or seq2_len==0 or seq3_len==0 or seq4_len==0:
            continue

        # mask exon-intron and intron-exon junctions
        # 5' splice site: -3, +6 mask (CAG + GTAAGT)
        # 3' splice site: -3, +2 mask (CAG + GT)

        start = 15+r1_exon-min(3, r1_exon) # here we don't do -1 since we are at the first intron nt
        stop = 15+r1_exon+min(6, r1_intron)-1
        seq1 = seq1[:start] + "N"*(stop-start+1) + seq1[stop+1:]

        start = 15+r2_intron-min(3, r2_intron)
        stop = 15+r2_intron+min(2, r2_exon)-1
        seq2 = seq2[:start] + "N"*(stop-start+1) + seq2[stop+1:]

        start = 15+r3_exon-min(3, r3_exon)
        stop = 15+r3_exon+min(6, r3_intron)-1
        seq3 = seq3[:start] + "N"*(stop-start+1) + seq3[stop+1:]

        start = 15+r4_intron-min(3, r4_intron)
        stop = 15+r4_intron+min(2, r4_exon)-1
        seq4 = seq4[:start] + "N"*(stop-start+1) + seq4[stop+1:]

        assert(seq1_len==len(seq1))
        assert(seq2_len==len(seq2))
        assert(seq3_len==len(seq3))
        assert(seq4_len==len(seq4))

        sequence[eid] = (seq1, seq2, seq3, seq4)
    pickle.dump(sequence, open(pickle_filename, "wb"))

def save_apa(comps, genome, hw=15):
    pickle_folder = os.path.join(rnamotifs2.path.comps_folder, comps, "pickle")
    pickle_folder = os.path.join(rnamotifs2.path.comps_folder, comps, "pickle")
    if not os.path.exists(os.path.join(pickle_folder)):
        os.makedirs(os.path.join(pickle_folder))
    pickle_filename = os.path.join(pickle_folder, "sequence.pickle")
    sequence = {}
    for (eid, chr, strand, pos, event_class) in rnamotifs2.data.data:
        seq1 = pybio.genomes.seq_direct(genome, chr.replace("chr", ""), strand, pos-rnamotifs2.data.flanking-hw, pos+rnamotifs2.data.flanking+hw)
        # mask poly-A signal
        for h in PAS_hexamers:
            if seq1.find(h)!=-1:
                seq1 = seq1.replace(h, "NNNNNN")
                break
        sequence[eid] = (seq1)
    pickle.dump(sequence, open(pickle_filename, "wb"))
