import os
import sys
from collections import Counter
import rnamotifs2
import sys

m = sys.modules[__name__]

data = []
data_class = []
data_type = None
genome = None

def read_config(comps):
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    config_file = os.path.join(comps_folder, "%s.config" % comps)
    if not os.path.exists(config_file):
        return None
    f = open(config_file, "rt")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        if r[0].startswith("#"):
            r = f.readline()
            continue
        pars = r[0].split("=")
        if len(pars)==2:
            try:
                setattr(m, pars[0], eval(pars[1]))
            except:
                try:
                    setattr(m, pars[0], int(pars[1]))
                except:
                    setattr(m, pars[0], pars[1])
        r = f.readline()
        continue
    f.close()

def read(comps):
    read_config(comps)
    if data_type=="apa":
        read_apa(comps)
    else:
        read_splice(comps)

def read_splice(comps):
    dist = Counter()
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    comps_filename = os.path.join(comps_folder, "%s.tab" % comps)
    temp = []
    f = open(comps_filename, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        temp.append((r[0], r[1], r[2], int(r[3]), int(r[4]), int(r[5]), int(r[6]), r[7]))
        data_class.append(r[7])
        dist[r[7]] += 1
        r = f.readline()
    f.close()
    # [s,s,s,s,e,e,e,e,c,c,c,c,c,c,c] : we want to have this sorted because of permutations and counting frequencies of s/e/c with bisection
    rnamotifs2.data.data = sorted(temp, key=lambda x: x[-1], reverse=True) # reverse is cosmetic, order: s/e/c
    rnamotifs2.data.dist = dist
    data_class.sort(reverse=True) # reverse is cosmetic, order: s/e/c

def read_apa(comps):
    dist = Counter()
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    comps_filename = os.path.join(comps_folder, "%s.tab" % comps)
    temp = []
    f = open(comps_filename, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        temp.append((r[0], r[1], r[2], int(r[3]), r[4]))
        data_class.append(r[4])
        dist[r[4]] = dist.setdefault(r[4], 0) + 1
        r = f.readline()
    f.close()
    # [s,s,s,s,e,e,e,e,c,c,c,c,c,c,c] : we want to have this sorted because of permutations and counting frequencies of s/e/c with bisection
    rnamotifs2.data.data = sorted(temp, key=lambda x: x[-1], reverse=True) # reverse is cosmetic, order: s/e/c
    data_class.sort(reverse=True) # reverse is cosmetic, order: s/e/c
    rnamotifs2.data.dist = dist
