# RNAmotifs2: cluster motif analysis

RNAmotifs2 is a Python software for identifying clusters of motifs underlying regulatory principles of alternative splicing and alternative polyadenylation. RNAmotifs2 can be used as a standalone software (requires [pybio](https://github.com/grexor/pybio)), however it is also integrated with [expressRNA](http://expressRNA.org).

## Short description

Using high-through sequencing and bioinformatics analysis, we obtain sets of enhanced, repressed and control features (e.g. exons or polyA sites) comparing control and test experimental conditions. The research question is then of how to identify sequence short-motifs in regions surrounding the regulated features.

Initially, we search all features for all possible short motifs (3, 4, 5-mers). The detected motif signals are convoluted with a short sliding window (15nt half-window). This accounts for the fact that RNA-protein binding affinity is influenced by several closely-spaced short motifs. After the signal is filtered, the Fisher test is performed on each of the two comparisons (enhanced vs. control, repressed vs. control) and each of the motif signals (AAA, AAT, AAC, etc.) separately. At this step (after FDR), we identified the strongest motif best separating the sequences, i.e. one motif for the enhanced vs. control, another for the repressed vs. control comparison.

However, since several proteins can regulate pre-mRNA processing by binding simultaneously around regulated features, the features with the strongest identified motif signal are removed. This is then compensated by searching the remaining feature space with the already identified motif (or cluster) paired with all other possible short-motifs. The search is reiterated until we reach a cluster of max. 4 motifs.

Finally, we compute an enrichment score (ES) on the super-imposed sequences of all the features (exons, polyA sites) and draw a motif regulatory RNA-map.

## Installation and running

RNAmotifs2 runs on Python 3 and depends on [pybio](https://github.com/grexor/pybio)
(genome handling), plus `numpy`, `scipy`, `matplotlib`, `fisher` and
`pyliftover` (only for lifting example coordinates).

```bash
micromamba create -n rnamotifs2 -c conda-forge -c bioconda python=3.12 \
    numpy scipy matplotlib pandas pysam psutil beautifulsoup4 requests pip
micromamba run -n rnamotifs2 pip install fisher pyliftover pybio
```

A comparison lives in `comps/<name>/` and needs two files:

* `<name>.tab` — tab-separated events with columns
  `id chr strand skip_start in_start in_stop skip_stop event_class`, where
  `event_class` is `s` (silenced), `e` (enhanced) or `c` (control)
* `<name>.config` — `data_type=splice`, `genome=<species>.<version>`
  (e.g. `homo_sapiens.ensembl115`), `hw=15`, `cores=<n>`

Run the whole analysis (motif search, cluster growth, RNA maps) with:

```bash
./run_example.sh <name>        # defaults to paper.bh
```

Output lands in `comps/<name>/`: per-region `results*.tab` / `tree*.tab` and
`rnamap/index.html` with the RNA maps.

The motif search is vectorised (`rnamotifs2.fastsearch`) and runs in a
multiprocessing pool. Set `RNAMOTIFS2_REFERENCE=1` to fall back to the original
pure-Python `search.v17` (much slower); `RNAMOTIFS2_NOPOOL=1` runs the pool
serially. `tests/golden_r1s.py` checks the fast path reproduces the reference
outputs bit-for-bit.

### Example: `comps/paper.bh`

The bundled brain/heart splicing example ships in hg19 coordinates. It has
been lifted to GRCh38 / Ensembl 115 with
`comps/paper.bh/lift_hg19_to_ensembl115.py` (original kept as
`paper.bh.hg19.tab`). Just run `./run_example.sh`.

### Multiple testing: `comps/paper.bh.strict`

Each cluster-growth step scans up to ~320 candidate motifs and keeps only the
single best (minimum p-value) one — a real single comparison of `paper.bh`
runs ~26,000 individual Fisher tests, of which ~70 become "the" reported
motif for some tree/step. The minimum of many tests is not itself a valid
p-value (it's biased low), so two corrections are available:

* `perms=<n>` — permutation/bootstrap empirical p-value (`p_emp` column):
  for each motif, `n` random relabelings of the event classes give a null
  distribution to compare its own real signal against. Was silently
  non-functional (a module mismatch meant `perms=` in a `.config` file never
  reached the code that reads it, and the empirical-p computation itself had
  a variable-name bug) — both are now fixed, and the per-motif accumulation
  is vectorised so `perms=200` costs no measurable time per motif.
* `use_FDR=True` — Benjamini-Hochberg-corrects the `fisher` column within
  each step's ranking before anything downstream (including the
  significance thresholds that decide whether a tree keeps growing) reads it.

Neither one corrects for "picked the best of ~320 candidates" on its own —
that would need permuting labels and rescoring the *entire* candidate set per
permutation (Westfall-Young / max-T), not yet implemented. They do correct
the two most misleading things: the reported p understating how many motifs
were tried (FDR), and whether a motif's own association could plausibly arise
by chance relabeling (permutation).

`comps/paper.bh.strict` is the same lifted table as `paper.bh` with both
turned on (`perms=200`, `use_FDR=True`) — run it with
`./run_example.sh paper.bh.strict` and compare its `results0.tab` `fisher` /
`p_emp` columns to `paper.bh`'s.

## Authors

[RNAmotifs2](https://github.com/grexor/rnamotifs2) is maintained by [Gregor Rot](https://grexor.github.io) in collaboration with several research laboratories worldwide.

The development started in 2010 when Matteo Cereda wrote and published the first version of [RNAmotifs](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-1-r20). In 2014, Gregor Rot refactored the RNAmotifs code to Python and created a new branch of the software adding cluster analysis and integrating the software with [expressRNA](http://expressRNA.org). This cluster branch of the software is now called RNAmotifs2.

## Citing RNAmotifs2

[High-resolution RNA maps suggest common principles of splicing and polyadenylation regulation by TDP-43](http://www.cell.com/cell-reports/abstract/S2211-1247(17)30522-3)<br />
Rot, G., Wang, Z., Huppertz, I., Modic, M., Lenče, T., Hallegger, M., Haberman, N., Curk, T., von Mering, C., Ule, J.<br />
Cell Reports , Volume 19 , Issue 5 , 1056 - 1067

## Reporting problems

Use the [issues page](https://github.com/grexor/rnamotifs2/issues) to report issues and leave suggestions.
