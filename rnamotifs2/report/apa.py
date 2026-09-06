"""
RNA-map report for alternative-polyadenylation (`data_type=apa`) comparisons.

This is the older table/highslide layout - the splice report
(`rnamotifs2.report.build_splice`) was rewritten, this one was left as-is
because there is no bundled apa example to validate a rewrite against. It is
kept importable and callable so `rnamotifs2.report.build()` can dispatch to
it.
"""
import os
import time
import datetime
import glob

import rnamotifs2
import pybio

WEB_FOLDER = "/rnamotifs2"  # expressRNA static-asset prefix


def build_apa(comps):
    nocache = time.mktime(datetime.datetime.today().timetuple())
    comps_folder = os.path.join(rnamotifs2.path.comps_folder, comps)
    rnamotifs2.data.read(comps)
    rnamotifs2.sequence.load(comps)

    print(rnamotifs2.data.dist.get("s", 0), rnamotifs2.data.dist.get("e", 0),
          rnamotifs2.data.dist.get("c", 0))

    configs = {}
    for region in ["r1s", "r1e", "r2s", "r2e", "r3s", "r3e"]:
        for index, tree_fname in enumerate(sorted(glob.glob(os.path.join(comps_folder, region, "tree*.tab")))):
            fisher, h, motif_cluster = rnamotifs2.draw.read_tree(tree_fname)
            areas, h, stats = rnamotifs2.search.areas_apa(comps, motif=motif_cluster, h=h)
            configs.setdefault(region, []).append((fisher, h, motif_cluster, areas))

    os.makedirs(os.path.join(comps_folder, "rnamap"), exist_ok=True)

    f = open(os.path.join(comps_folder, "rnamap", "index.html"), "wt")
    f.write("<html>\n")
    head = """<head>
<script type="text/javascript" src="%s/highslide/highslide/highslide.js"></script>
<link rel="stylesheet" type="text/css" href="%s/highslide/highslide/highslide.css" />
<script type="text/javascript">
    hs.graphicsDir = '%s/highslide/highslide/graphics/';
    hs.showCredits = false;
</script>
<style>
.highslide img { border: 0px; outline: none; }
a { text-decoration: none; }
</style>
</head>""" % (WEB_FOLDER, WEB_FOLDER, WEB_FOLDER)
    f.write(head + "\n")
    f.write("<body>\n<center>")
    f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 12px;'>"
            "<tr><td width=15px></td><td width=15px></td><td width=15px></td><td></td></tr>\n")

    max_es = 0
    for region, tree_list in configs.items():
        for fisher, h, motif_list, area in tree_list:
            area.setdefault("s", [0] * 201)
            area.setdefault("e", [0] * 201)
            logs, loge = rnamotifs2.compute.es_apa(area["s"], area["e"], area["c"])
            max_es = max(max_es, max(pybio.utils.smooth(logs[0])), max(pybio.utils.smooth(loge[0])))

    for region, tree_list in configs.items():
        for index, (fisher, h, motif_list, area) in enumerate(tree_list):
            area.setdefault("s", [0] * 201)
            area.setdefault("e", [0] * 201)
            logs, loge = rnamotifs2.compute.es_apa(area["s"], area["e"], area["c"])
            image_filename = "%s_tree%s_area%s" % (region, index, 1)
            rnamotifs2.draw.area_apa("+".join(motif_list), logs[0], loge[0],
                                     os.path.join(comps_folder, "rnamap", image_filename),
                                     area=0, region=region, limy=max_es, fisher=fisher)

    for region in sorted(configs.keys()):
        for index, (fisher, h, motif_list, area) in enumerate(configs[region]):
            tree_file = "%s/%s/%s/tree%s.tab" % (WEB_FOLDER, comps, region, index)
            f.write("<tr>")
            for r in ["r1", "r2", "r3"]:
                if region[:-1] == r:
                    color = "#ff0000" if region[-1] == "e" else "#0000ff"
                else:
                    color = "#aaaaaa"
                f.write("<td align=center width=15px>")
                if fisher is not None:
                    f.write("<svg xmlns='http://www.w3.org/2000/svg' width='18px' height='18px'>"
                            "<circle cx=9 cy=9 r=7 stroke=#ffffff stroke-width=1 fill='%s'/></svg>" % color)
                f.write("</td>")
            img = "%s/%s/rnamap/%s_tree%s_area1.png?nocache=%s" % (WEB_FOLDER, comps, region, index, nocache)
            f.write("<td align=center>%s (<a href=%s target=_new>tree=%s</a>, h=%s)<br>"
                    "<a href=%s class='highslide' onclick='return hs.expand(this)'>"
                    "<img src=%s.png width=350></a></td>" % (
                        "+".join(motif_list), tree_file, index, h, img, img))
            f.write("</tr><tr style='height:10px;'></tr>")
    f.write("\n</table></body></html>\n")
    f.close()
