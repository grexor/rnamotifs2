import os
import sys
import rnamotifs2
import math
import numpy as np
import pybio

def highlight_region(mpatches, ax, x0, width, edgecolor):
    """Yellow fill for the analyzed window, with a colored border showing
    which event class (blue=silenced, red=enhanced) it was detected in.

    The rectangle is sized to the axes' actual y-limits (rather than a fixed
    height) so its top and bottom edges land inside the visible plot area —
    a fixed height taller than the data range clips the top border out of
    view entirely. Because the top/bottom edges then sit exactly on the axes'
    clip boundary, the border patch has clipping disabled — otherwise only
    the inner half of its stroke survives there, making those two edges look
    thinner than the (strictly interior) left/right edges.
    """
    y0, y1 = ax.get_ylim()
    ax.add_patch(mpatches.Rectangle([x0, y0], width, y1 - y0, color="#FFFF00", alpha=0.1, zorder=2))
    ax.add_patch(mpatches.Rectangle([x0, y0], width, y1 - y0, facecolor='none',
                                    edgecolor=edgecolor, linewidth=4, zorder=3, clip_on=False))

def read_tree(filename):
    motif = []
    rows = []
    h = None
    f = open(filename, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        rows.append(r)
        r = f.readline()
    f.close()
    # fisher, h, motif_cluster
    if int(rows[-1][0])>0: #step > 0, combine the cluster motifs
        return float(rows[-1][1]), int(rows[-1][4]), rows[-1][-2].split("_")+rows[-1][-1].split("_")
    else: # step = 0, just return one single motif
        return float(rows[-1][1]), int(rows[-1][4]), rows[-1][-2].split("_")

def read_tree_full(filename):
    """Every column of a tree{cn}.tab's last (final) row, as a dict, plus
    p_emp - the permutation empirical p-value, which tree{cn}.tab itself
    never carries (only results{cn}.tab / c{cn}.temp{step}.tab do) - read
    back from the same motif-set-keyed pickle everything else in this
    package uses, so it needs no schema change to tree{cn}.tab.

    Derives region/cn from `filename` (.../<comps>/<region>/tree<cn>.tab).
    """
    import pickle as _pickle
    region_folder, base = os.path.split(filename)
    region = os.path.basename(region_folder)
    cn = int(base[len("tree"):-len(".tab")])

    f = open(filename, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    rows = [l.replace("\r", "").replace("\n", "").split("\t") for l in f if l.strip()]
    f.close()
    last = dict(zip(header, rows[-1]))

    motif = last["motif"].split("_")
    cmotif = last["cmotif"].split("_") if last["cmotif"] else []
    pickle_folder = os.path.join(region_folder, "pickle")
    if cmotif:
        pf = os.path.join(pickle_folder, "c%s.%s.filter.%s.pickle" % (cn, "_".join(sorted(motif)), "_".join(sorted(cmotif))))
    else:
        pf = os.path.join(pickle_folder, "c%s.%s.pickle" % (cn, "_".join(sorted(motif))))

    p_emp = None
    if os.path.exists(pf):
        _, test_result, _, _, _, _, _ = _pickle.load(open(pf, "rb"))
        fisher_raw, p_emp_list, _ = test_result
        if p_emp_list:
            k = sum(1 for x in p_emp_list if x <= fisher_raw)
            p_emp = (1 + k) / float(1 + len(p_emp_list))

    last["region"] = region
    last["cn"] = cn
    last["cluster"] = motif + cmotif
    last["p_emp"] = p_emp
    last["n_steps"] = int(last["step"]) + 1
    return last

def area(motif, s, e, filename, area=None, region=None, limy=None, stats=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import math
    import gzip
    from matplotlib import cm as CM
    import numpy
    import matplotlib.patches as mpatches
    import matplotlib.ticker as mticker
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.colors as mcolors
    c = mcolors.ColorConverter().to_rgb

    # styling
    matplotlib.rcParams['axes.labelsize'] = 30
    matplotlib.rcParams['axes.titlesize'] = 22
    matplotlib.rcParams['xtick.labelsize'] = 30
    matplotlib.rcParams['ytick.labelsize'] = 30
    matplotlib.rcParams['legend.fontsize'] = 30
    matplotlib.rc('axes', edgecolor='gray')
    matplotlib.rcParams['axes.linewidth'] = 0.3
    matplotlib.rcParams['legend.frameon'] = 'False'

    v = [x+y for x,y in zip(s,e)]

    fig = plt.figure(figsize=(20, 5))
    a = plt.axes([0.1, 0.3, 0.8, 0.6])
    a.grid(color="lightgray", linestyle="--", linewidth=0.4)
    a.set_axisbelow(True)

    plt.xlabel("distance [nt]")
    #plt.ylabel("-2 * log(p_value)")
    plt.ylabel("enrichment score")

    s = pybio.utils.smooth(s)
    e = pybio.utils.smooth(e)

    plt.fill_between(range(0, len(s)), 0, s, facecolor='blue', alpha=0.5, interpolate=True)
    plt.plot(range(0, len(s)), s, color='blue', alpha=1)

    plt.fill_between(range(0, len(e)), 0, e, facecolor='red', alpha=0.5, interpolate=True)
    plt.plot(range(0, len(e)), e, color='red', alpha=1)

    a.set_xlim(0, 250)
    a.set_ylim(0, math.ceil(limy))

    if area==0:
        p = mpatches.Rectangle([50, 0], 0.01, 200, facecolor='none', edgecolor=(0.7, 0.7, 0.7))
        plt.gca().add_patch(p)
        p = mpatches.Rectangle([0, 0], 50, 200, color=(0,0,0), alpha=0.1)
        plt.gca().add_patch(p)
        plt.xticks([0,25,50,75,100,125,150,175,200,225,250], [-50, -25, "skip.start", 25, 50, 75, 100, 125, 150, 175, 200])
    # the analyzed window is highlighted in yellow, with a border in the color
    # of the event class it was detected in: blue for silenced (region ends
    # in "s"), red for enhanced (region ends in "e")
    highlight = "blue" if region[-1]=="s" else "red"

    if area==1:
        if region.startswith("r1"):
            highlight_region(mpatches, plt.gca(), 155, 40, highlight)
        if region.startswith("r2"):
            highlight_region(mpatches, plt.gca(), 200, 30, highlight)
        p = mpatches.Rectangle([200, 0], 0.01, 200, facecolor='none', edgecolor=(0.7, 0.7, 0.7))
        plt.gca().add_patch(p)
        p = mpatches.Rectangle([200, 0], 50, 200, color=(0,0,0), alpha=0.1)
        plt.gca().add_patch(p)
        plt.xticks([0,25,50,75,100,125,150,175,200,225,250], [-200, -175, -150, -125, -100, -75, -50, -25, "in.start", 25, 50])
    if area==2:
        if region.startswith("r2"):
            highlight_region(mpatches, plt.gca(), 20, 30, highlight)
        if region.startswith("r3"):
            highlight_region(mpatches, plt.gca(), 55, 40, highlight)
        p = mpatches.Rectangle([50, 0], 0.01, 200, facecolor='none', edgecolor=(0.7, 0.7, 0.7))
        plt.gca().add_patch(p)
        p = mpatches.Rectangle([0, 0], 50, 200, color=(0,0,0), alpha=0.1)
        plt.gca().add_patch(p)
        plt.xticks([0,25,50,75,100,125,150,175,200,225,250], [-50, -25, "in.stop", 25, 50, 75, 100, 125, 150, 175, 200])
    if area==3:
        p = mpatches.Rectangle([200, 0], 0.01, 200, facecolor='none', edgecolor=(0.7, 0.7, 0.7))
        plt.gca().add_patch(p)
        p = mpatches.Rectangle([200, 0], 50, 200, color=(0,0,0), alpha=0.1)
        plt.gca().add_patch(p)
        plt.xticks([0,25,50,75,100,125,150,175,200,225,250], [-200, -175, -150, -125, -100, -75, -50, -25, "skip.stop", 25, 50])

    # only display statistics on areas with region info (yellow plots)
    if (area==1 and region[:2] in ["r1", "r2"]) or (area==2 and region[:2] in ["r2", "r3"]):
        r = region[0]+region[1]
        title = "%s : %s=%s (of %s=%.2f%%), %s=%s (of %s=%.2f%%)" % (motif, region, stats[region], stats[region[-1]], stats[region]*100/float(stats[region[-1]]), r+"c", stats[r+"c"], stats["c"], stats[r+"c"]*100/float(stats["c"]))
    else:
        title = motif

    plt.title(title)
    plt.savefig(filename+".png", dpi=80, transparent=True)
    plt.savefig(filename+".pdf")
    plt.close()
    return sum(v)

def area_apa(motif, s, e, filename, area=None, region=None, limy=None, fisher=None):
    print(motif)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import math
    import gzip
    from matplotlib import cm as CM
    import numpy
    import matplotlib.patches as mpatches
    import matplotlib.ticker as mticker
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.colors as mcolors
    c = mcolors.ColorConverter().to_rgb

    # styling
    matplotlib.rcParams['axes.labelsize'] = 30
    matplotlib.rcParams['axes.titlesize'] = 22
    matplotlib.rcParams['xtick.labelsize'] = 30
    matplotlib.rcParams['ytick.labelsize'] = 30
    matplotlib.rcParams['legend.fontsize'] = 30
    matplotlib.rc('axes', edgecolor='gray')
    matplotlib.rcParams['axes.linewidth'] = 0.3
    matplotlib.rcParams['legend.frameon'] = 'False'

    v = [x+y for x,y in zip(s,e)]

    fig = plt.figure(figsize=(20, 4))
    a = plt.axes([0.04, 0.1, 0.92, 0.75])

    a.grid(color="lightgray", linestyle="--", linewidth=0.4)
    a.set_axisbelow(True)

    s = pybio.utils.smooth(s)
    e = pybio.utils.smooth(e)

    plt.fill_between(range(0, len(s)), 0, s, facecolor='blue', alpha=0.5, interpolate=True)
    plt.plot(range(0, len(s)), s, color='blue', alpha=1)

    plt.fill_between(range(0, len(e)), 0, e, facecolor='red', alpha=0.5, interpolate=True)
    plt.plot(range(0, len(e)), e, color='red', alpha=1)

    a.set_xlim(0, 200)
    a.set_ylim(0, math.ceil(limy))

    p = mpatches.Rectangle([100, 0], 0.01, 100, facecolor='none', edgecolor=(0.7, 0.7, 0.7))
    plt.gca().add_patch(p)
    plt.xticks([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200], [-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100])

    # the analyzed window is highlighted in yellow, with a border in the color
    # of the event class it was detected in: blue for silenced (region ends
    # in "s"), red for enhanced (region ends in "e")
    highlight = "blue" if region[-1]=="s" else "red"

    if region.startswith("r1"):
        plt.text(20, limy-2, "R1 [-100, -40]", fontsize=25)
        highlight_region(mpatches, plt.gca(), 0, 60, highlight)

    if region.startswith("r2"):
        plt.text(80, limy-2, "R2 [-40, 20]", fontsize=25)
        highlight_region(mpatches, plt.gca(), 60, 60, highlight)

    if region.startswith("r3"):
        plt.text(140, limy-2, "R3 [20, 80]", fontsize=25)
        highlight_region(mpatches, plt.gca(), 120, 60, highlight)

    if fisher!=None:
        plt.title("%s, p-value = %.2e" % (motif, fisher), y=1.06)
    else:
        plt.title("%s, p-value = cluster manually added" % (motif))

    #plt.tight_layout()
    plt.savefig(filename+".png", dpi=80, transparent=True)
    plt.savefig(filename+".pdf")
    plt.close()
    return sum(v)
