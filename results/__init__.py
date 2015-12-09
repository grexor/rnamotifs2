import os
import sys
import rnamotifs2
import pybio
import glob

def get_motifs(comps, region):
    motifs_all = rnamotifs2.get_motifs()
    motifs_remove = []
    # remove all significant motifs from previous trees!
    region_folder = os.path.join(rnamotifs2.path.comps_folder, comps, region)
    trees = glob.glob(region_folder+"/tree*.tab")
    for filename in trees:
        f = open(filename, "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            if float(data["fisher"])<0.1:
                motifs_remove.append(data["motif"])
            r = f.readline()
    motifs_all = [m for m in motifs_all if m not in motifs_remove]
    return motifs_all
