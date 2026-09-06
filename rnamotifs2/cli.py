"""
Command-line interface. Installed as the ``rnamotifs2`` console script.

    rnamotifs2 run   -comps <name>          # full analysis + report
    rnamotifs2 draw  -comps <name>          # (re)build the RNA-map report only
    rnamotifs2 motif <comps> <genome> <region> <motif> <pth> <cn> <sf>
    rnamotifs2 motif-cluster <comps> <genome> <region> <motif> <cn> <cmotif> <pth> <sf>

`rnamotifs2 run` reads `comps/<name>/<name>.config`, wipes stale per-region
output (keeping the signature-checked sequence cache), runs every region and
then regenerates the report. The comps directory is
`$RNAMOTIFS2_COMPS`, else `./comps` - or set it in code with
`rnamotifs2.path.set_comps_folder(...)` (what expressRNA should do).
"""
import argparse
import glob
import os
import shutil
import sys

import rnamotifs2


def run(comps):
    rnamotifs2.data.read_config(comps)
    rnamotifs2.data.read(comps)

    comps_dir = os.path.join(rnamotifs2.path.comps_folder, comps)
    for f in glob.glob(os.path.join(comps_dir, "*")):
        if os.path.isdir(f) and os.path.basename(f) != "pickle":
            shutil.rmtree(f)

    for region in ["r1s", "r1e", "r2s", "r2e", "r3s", "r3e"]:
        if rnamotifs2.data.dist.get("s", None) is None and region[-1] == "s":
            continue
        if rnamotifs2.data.dist.get("e", None) is None and region[-1] == "e":
            continue
        rnamotifs2.start(comps, region, 0, rnamotifs2.data.pth)

    rnamotifs2.report.build(comps)


def draw(comps):
    rnamotifs2.report.build(comps)


def _motif(argv):
    comps, genome, region, motif, pth, cn, sf = argv[:7]
    rnamotifs2.motifjob.run_motif(comps, genome, region, motif.split("_"),
                                  float(pth), int(cn), sf)


def _motif_cluster(argv):
    comps, genome, region, motif, cn, cmotif, pth, sf = argv[:8]
    rnamotifs2.motifjob.run_motif_cluster(comps, genome, region, motif.split("_"),
                                          int(cn), cmotif.split("_"), float(pth), sf)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)

    parser = argparse.ArgumentParser(prog="rnamotifs2")
    parser.add_argument("--comps-folder", default=None,
                        help="directory holding the per-comparison subdirs "
                             "(default: $RNAMOTIFS2_COMPS or ./comps)")
    sub = parser.add_subparsers(dest="cmd")

    p_run = sub.add_parser("run", help="full analysis + report")
    p_run.add_argument("-comps", required=True)

    p_draw = sub.add_parser("draw", help="(re)build the RNA-map report only")
    p_draw.add_argument("-comps", required=True)

    sub.add_parser("motif", add_help=False)
    sub.add_parser("motif-cluster", add_help=False)

    # back-compat: `rnamotifs2 -comps <name>` with no subcommand == `run`
    if argv and argv[0] not in ("run", "draw", "motif", "motif-cluster", "-h", "--help") \
            and not argv[0].startswith("--comps-folder"):
        argv = ["run"] + argv

    if argv and argv[0] in ("motif", "motif-cluster"):
        rest = argv[1:]
        (_motif if argv[0] == "motif" else _motif_cluster)(rest)
        return

    args, _ = parser.parse_known_args(argv)
    if args.comps_folder:
        rnamotifs2.path.set_comps_folder(args.comps_folder)

    if args.cmd == "run":
        run(args.comps)
    elif args.cmd == "draw":
        draw(args.comps)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
