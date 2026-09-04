#!/usr/bin/env python3
"""
Lift comps/paper.bh/paper.bh.tab from hg19 (GRCh37) to GRCh38 / Ensembl 115.

The original table stores four genomic coordinates per splice event
(skip_start, in_start, in_stop, skip_stop) on the UCSC hg19 assembly.
This script converts them to GRCh38 with pyliftover (UCSC hg19 -> hg38 chain,
auto-downloaded on first run) and writes:

  paper.bh.tab                 <- lifted table (genome: homo_sapiens.ensembl115)
  paper.bh.hg19.tab            <- verbatim backup of the original
  paper.bh.lift_dropped.tab    <- events that could not be lifted cleanly

An event is kept only if all four coordinates lift to the same chromosome and
strand and preserve their original ordering (so the region geometry in
rnamotifs2.sequence.coords stays valid).
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "paper.bh.tab")
BACKUP = os.path.join(HERE, "paper.bh.hg19.tab")
DROPPED = os.path.join(HERE, "paper.bh.lift_dropped.tab")

COORD_COLS = ["skip_start", "in_start", "in_stop", "skip_stop"]


def main():
    from pyliftover import LiftOver

    if os.path.exists(BACKUP):
        src = BACKUP  # already lifted once; always re-lift from the pristine hg19 table
    else:
        src = SRC

    with open(src, "rt") as f:
        lines = [l.rstrip("\n") for l in f]

    # keep leading comment lines untouched
    comments = []
    i = 0
    while i < len(lines) and lines[i].startswith("#"):
        comments.append(lines[i])
        i += 1
    header = lines[i].split("\t")
    rows = [l.split("\t") for l in lines[i + 1:] if l.strip()]

    col = {name: idx for idx, name in enumerate(header)}
    for c in ["chr", "strand"] + COORD_COLS:
        if c not in col:
            sys.exit("column %r not found in %s" % (c, src))

    print("loading hg19 -> hg38 chain (downloads on first run) ...")
    lo = LiftOver("hg19", "hg38")

    kept, dropped = [], []
    for r in rows:
        chrom = r[col["chr"]]
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        strand = r[col["strand"]]

        newpos = {}
        ok = True
        reason = ""
        for c in COORD_COLS:
            pos = int(r[col[c]])
            conv = lo.convert_coordinate(chrom, pos, strand)
            if not conv:
                ok, reason = False, "%s unmapped" % c
                break
            nchrom, npos, nstrand = conv[0][0], conv[0][1], conv[0][2]
            if nchrom != chrom:
                ok, reason = False, "%s -> %s" % (c, nchrom)
                break
            if nstrand != strand:
                ok, reason = False, "%s strand flip" % c
                break
            newpos[c] = npos

        if ok:
            vals = [newpos[c] for c in COORD_COLS]
            if vals != sorted(vals):
                ok, reason = False, "coordinate order broken after lift"

        if not ok:
            dropped.append(r + [reason])
            continue

        nr = list(r)
        for c in COORD_COLS:
            nr[col[c]] = str(newpos[c])
        kept.append(nr)

    # write pristine backup once
    if not os.path.exists(BACKUP):
        with open(BACKUP, "wt") as f:
            f.write("\n".join(lines) + "\n")

    with open(SRC, "wt") as f:
        for c in comments:
            f.write(c + "\n")
        f.write("\t".join(header) + "\n")
        for r in kept:
            f.write("\t".join(r) + "\n")

    with open(DROPPED, "wt") as f:
        f.write("\t".join(header + ["reason"]) + "\n")
        for r in dropped:
            f.write("\t".join(r) + "\n")

    print("kept    : %d" % len(kept))
    print("dropped : %d  (see %s)" % (len(dropped), os.path.basename(DROPPED)))
    from collections import Counter
    cls = Counter(r[col["event_class"]] for r in kept) if "event_class" in col else {}
    print("classes : %s" % dict(cls))


if __name__ == "__main__":
    main()
