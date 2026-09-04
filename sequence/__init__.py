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
        r1_intron = min(max_intron, upintron_len//2)
        r1_start = skip_start-r1_exon
        r1_stop = skip_start + r1_intron

        r2_intron = min(max_intron, upintron_len//2)
        r2_exon = min(max_exon, exon_len//2)
        r2_start = in_start - r1_intron
        r2_stop = in_start + r2_exon

        r3_intron = min(max_intron, downintron_len//2)
        r3_exon = min(50, exon_len//2)
        r3_start = in_stop - r3_exon
        r3_stop = in_stop + r3_intron

        r4_exon = max_exon
        r4_intron = min(max_intron, downintron_len//2)
        r4_start = skip_stop - r4_intron
        r4_stop = skip_stop + r4_exon
    else:
        r1_exon = max_exon
        r1_intron = min(max_intron, downintron_len//2)
        r1_start = skip_stop - r1_intron
        r1_stop = skip_stop + r1_exon

        r2_intron = min(max_intron, downintron_len//2)
        r2_exon = min(max_exon, exon_len//2)
        r2_start = in_stop - r2_exon
        r2_stop = in_stop + r2_intron

        r3_exon = min(max_exon, exon_len//2)
        r3_intron = min(max_intron, upintron_len//2)
        r3_start = in_start - r3_intron
        r3_stop = in_start + r3_exon

        r4_exon = max_exon
        r4_intron = min(max_intron, upintron_len//2)
        r4_start = skip_start - r4_exon
        r4_stop = skip_start + r4_intron
    return (r1_exon, r1_intron, r1_start, r1_stop), (r2_exon, r2_intron, r2_start, r2_stop), (r3_exon, r3_intron, r3_start, r3_stop), (r4_exon, r4_intron, r4_start, r4_stop)

def load(comps):
    print("%s: loading sequences" % comps)
    pickle_folder = os.path.join(rnamotifs2.path.comps_folder, comps, "pickle")
    pickle_filename = os.path.join(pickle_folder, "sequence.pickle")
    rnamotifs2.sequence.sequence = pickle.load(open(pickle_filename, "rb"))

def save(comps, genome, hw=15):
    if rnamotifs2.data.data_type=="apa":
        save_apa(comps, genome, hw=hw)
    else:
        save_splice(comps, genome, hw=hw)


def _cache_signature(comps, genome, hw):
    """Hash of everything the extracted sequences depend on."""
    import hashlib
    comps_filename = os.path.join(rnamotifs2.path.comps_folder, comps, "%s.tab" % comps)
    h = hashlib.sha1()
    h.update(("%s|%s|%s\n" % (genome, hw, rnamotifs2.data.data_type)).encode())
    with open(comps_filename, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def _cache_valid(pickle_filename, sig):
    sig_filename = pickle_filename + ".sig"
    return (os.path.exists(pickle_filename) and os.path.exists(sig_filename)
            and open(sig_filename).read().strip() == sig)


class _ChrReader:
    """Reads slices from pybio's ``<chr>.string`` files, reusing file handles.

    Replaces one open()+seek() per event with one open() per chromosome.
    Falls back to pybio.core.genomes.seq_direct for anything unusual
    (missing .string file, negative coordinates).
    """

    def __init__(self, genome):
        import pybio
        self.species, self.gv = rnamotifs2.genomes.resolve(genome)
        self.genome = genome
        self.folder = os.path.join(
            pybio.config.genomes_folder,
            "%s.assembly.%s" % (self.species, self.gv))
        self._fh = {}

    def _handle(self, chr):
        if chr not in self._fh:
            path = os.path.join(self.folder, "%s.string" % chr)
            self._fh[chr] = open(path, "rt") if os.path.exists(path) else None
        return self._fh[chr]

    def get(self, chr, strand, start, stop):
        import pybio
        if start > stop:
            start, stop = stop, start
        if start < 0:
            return pybio.core.genomes.seq_direct(
                self.species, chr, strand, start, stop, genome_version=self.gv)
        fh = self._handle(chr)
        if fh is None:
            return ""
        fh.seek(start, 0)
        seq = fh.read(stop - start + 1)
        seq += "N" * ((stop - start + 1) - len(seq))
        if strand == "-":
            return pybio.sequence.reverse_complement(seq)
        return seq

    def close(self):
        for fh in self._fh.values():
            if fh is not None:
                fh.close()
        self._fh.clear()


def save_splice(comps, genome, hw=15):
    pickle_folder = os.path.join(rnamotifs2.path.comps_folder, comps, "pickle")
    if not os.path.exists(os.path.join(pickle_folder)):
        os.makedirs(os.path.join(pickle_folder))
    pickle_filename = os.path.join(pickle_folder, "sequence.pickle")
    sig = _cache_signature(comps, genome, hw)
    if _cache_valid(pickle_filename, sig):
        return
    reader = _ChrReader(genome)
    sequence = {}
    for (eid, chr, strand, skip_start, in_start, in_stop, skip_stop, event_class) in rnamotifs2.data.data:
        upintron_len = in_start-skip_start+1
        exon_len = in_stop - in_start + 1
        downintron_len = skip_stop - in_stop + 1
        (r1_exon, r1_intron, r1_start, r1_stop), (r2_exon, r2_intron, r2_start, r2_stop), (r3_exon, r3_intron, r3_start, r3_stop), (r4_exon, r4_intron, r4_start, r4_stop) = coords(strand, skip_start, in_start, in_stop, skip_stop)
        _chr = chr.replace("chr", "")
        seq1 = reader.get(_chr, strand, r1_start-hw, r1_stop+hw)
        seq2 = reader.get(_chr, strand, r2_start-hw, r2_stop+hw)
        seq3 = reader.get(_chr, strand, r3_start-hw, r3_stop+hw)
        seq4 = reader.get(_chr, strand, r4_start-hw, r4_stop+hw)

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
    reader.close()
    pickle.dump(sequence, open(pickle_filename, "wb"), protocol=2)
    with open(pickle_filename + ".sig", "wt") as f:
        f.write(sig)

def save_apa(comps, genome, hw=15):
    pickle_folder = os.path.join(rnamotifs2.path.comps_folder, comps, "pickle")
    if not os.path.exists(os.path.join(pickle_folder)):
        os.makedirs(os.path.join(pickle_folder))
    pickle_filename = os.path.join(pickle_folder, "sequence.pickle")
    sig = _cache_signature(comps, genome, hw)
    if _cache_valid(pickle_filename, sig):
        return
    reader = _ChrReader(genome)
    sequence = {}
    for (eid, chr, strand, pos, event_class) in rnamotifs2.data.data:
        seq1 = reader.get(chr.replace("chr", ""), strand, pos-rnamotifs2.data.flanking-hw, pos+rnamotifs2.data.flanking+hw)
        # mask poly-A signal
        for h in PAS_hexamers:
            if seq1.find(h)!=-1:
                seq1 = seq1.replace(h, "NNNNNN")
                break
        sequence[eid] = (seq1)
    reader.close()
    pickle.dump(sequence, open(pickle_filename, "wb"), protocol=2)
    with open(pickle_filename + ".sig", "wt") as f:
        f.write(sig)
