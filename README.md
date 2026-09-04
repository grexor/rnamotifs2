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
`rnamap/index.html` with the RNA maps. The report is fully self-contained
(no external JS/CSS — the old jquery/highslide paths only ever resolved when
deployed under expressRNA's own asset layout, so standalone they silently
did nothing):

* a **run settings** panel (grouped into Input / Events / Search strategy /
  Statistics / Compute) — data file, genome, permutation count, FDR on/off,
  search mode, both cluster-growth thresholds
* per motif cluster, a **stat card**: FDR q-value or raw p-value (labeled
  correctly either way), permutation p-value, cluster size, `h`, specificity
  and information gain (all in scientific notation), support counts
* **R1/R2/R3 dots** beside the cluster name marking which of the three areas
  along the splice junction it was found in, and a small **dot on the
  specific area plot(s)** that actually carry the highlighted (regulated)
  window — both blue for silenced/repressed, red for enhanced
* click any plot for an in-page **zoom/pan lightbox** (wheel to zoom, drag to
  pan, double-click/Escape to reset) — plain JS, nothing external or bundled

See `comps/paper.bh.strict.beam.fullrecursive/rnamap/index.html` (checked
into the repo — open it directly, no need to run anything) for what it looks
like end to end.

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

### Greedy search: `beam_width` and `comps/paper.bh.strict.beam`

FDR and permutations fix the reported *p-value*; they don't touch a separate
problem in the cluster-growth *search itself*. It's a greedy algorithm: at
each step it commits to the single best-scoring motif, removes its support,
and only ever builds on top of that one choice. A motif that scores worse
individually can still combine into a much stronger cluster than the "best"
individual motif ever reaches — pure greedy search never looks, because it
locked onto the top-ranked motif at step 0 and never backtracks.

`beam_width=<k>` (default 1 = plain greedy, `cluster.next_cluster`) grows `k`
independent chains in parallel, one per each of the top-`k` base motifs
(`cluster.next_cluster_beam`), and keeps whichever chain ends up with the
best final score — every beam's result is written to `beams<n>.tab` for
inspection. This diversifies the highest-leverage decision (choosing among
~320 candidates at step 0); it does not re-rank chains against each other at
every later step too (full recursive beam search), which would cost
`beam_width` times more at every step instead of just the first.

On `comps/paper.bh.strict.beam` (`beam_width=5` on top of `paper.bh.strict`'s
`perms=200`/`use_FDR=True`), region r1s's rank-2 starting motif (`TTCA`,
individually weaker than the greedy #1 pick `TCAT`) won: its final cluster
(`TGT`+`ATT`+`TTCA`) reached fisher `5.4e-10`, about five orders of magnitude
better than the plain-greedy `TCAT` chain's `3.0e-05` — exactly the failure
mode beam search exists to catch.

### Full recursive beam search: `beam_recursive` and `...fullrecursive`

`beam_width` alone (`cluster.next_cluster_beam`) only diversifies the step-0
choice; each of the `k` chains is then grown by its own independent greedy
search, so a promising *combination* that only becomes visible a few steps
into some chain is still invisible to the others. `beam_recursive=True`
(`cluster.next_cluster_beam_recursive`) removes that limit: at *every* step
it scores every active beam's candidate extensions, globally re-ranks all of
them together, and keeps only the best `beam_width` overall — which can mean
one strong beam produces several of the survivors while a weaker one drops
out, or a beam that looked mediocre a step ago suddenly pulls ahead. A beam
that drops out (no viable next candidate, or loses the global competition)
isn't discarded — its current state joins an archive of completed candidate
answers, and the single best state in that whole archive (not necessarily
the longest-lived beam) is promoted to `tree<n>.tab`; the full archive is
written to `beams<n>.tab`. If `use_FDR`, the FDR correction is applied across
*all* active beams' pooled candidates for that step, so a step exploring more
combinations is held to a correspondingly stricter bar.

This costs up to `beam_width` times more work at every step (vs. paying that
cost only once, at step 0, for `next_cluster_beam`) — 23 min end to end on
`comps/paper.bh`, vs. `paper.bh.strict.beam`'s 15 min and plain `paper.bh.strict`'s
3.5 min (this run also carries `perms=1000` vs. the other two's `200`).
`comps/paper.bh.strict.beam.fullrecursive` (`paper.bh.strict` +
`beam_width=5`, `beam_recursive=True`, `perms=1000`) shows all three outcomes
a step-wise re-ranking search can produce, on the same three regions that
found anything significant at all:

| region | greedy (`paper.bh.strict`) | top-`k`-restart (`...beam`) | full recursive |
|---|---|---|---|
| r1s | `ATTT+TGTG+TCAT`, fisher 1.4e-05 | same as greedy | **`TGTG+TCAT`** (2 motifs), fisher 1.6e-05 |
| r1e | `TCA+TGCT+TCT`, fisher 3.3e-05 | `TAAC+TCT+TGC+TTC`, fisher 9.5e-08 | `CTGT+TCT+CTT`, fisher 1.2e-07 |
| r3e | `TCTC+TGT+TCAT+CAT`, fisher 9.4e-04 | same as greedy | **`CATT+TGTG+CATC+CAT`**, fisher 4.5e-06 |

r3e is the clean win recursive search exists for: neither greedy nor
restarting from a different base motif ever found this combination — only
re-ranking every beam's candidates jointly, every step, surfaced it, ~200x
better than the other two methods' shared answer. r1e shows recursive search
matching (not exceeding) the improvement `next_cluster_beam` already found by
a different route. r1s is the interesting one: the *raw* 3-motif extension
`ATTT+TGTG+TCAT` was tried here too (it's in `beams0.tab`) and still cleared
`cluster_stop_thr`, so a chain kept extending past it — but its FDR-corrected
q-value came out worse than the 2-motif state one step earlier (`TGTG+TCAT`),
because that q-value is computed by pooling *all 5 beams'* candidates that
step (~1,575 tests) rather than one chain's ~315, a correspondingly stricter
bar. The archive mechanism catches exactly this: growing further isn't always
better, and only tracking "the best state ever seen" (not "wherever growth
stopped") reports it correctly.

The rendered report for this run is checked into the repo at
[`comps/paper.bh.strict.beam.fullrecursive/rnamap/index.html`](comps/paper.bh.strict.beam.fullrecursive/rnamap/index.html) —
open it directly (no need to run anything) to see the settings panel,
per-region stat cards, R1/R2/R3 + per-area regulation dots, and the
click-to-zoom plots described above.

### A note on `perms`

Bumped the example configs' `perms` from `200` to `1000` going forward: the
empirical p-value floor is `1/(perms+1)`, so `1000` gives a floor
(`~0.001`) that lines up with `cluster_stop_thr`, at a modest few-minutes
cost increase over `200` (per-motif cost scales roughly with `perms`, but
most of a run's total time is fixed per-region/per-step overhead, not
per-motif). `10000`+ only pays off if you need to resolve significance
*below* `0.001`, not just detect that something clears it.

## Authors

[RNAmotifs2](https://github.com/grexor/rnamotifs2) is maintained by [Gregor Rot](https://grexor.github.io) in collaboration with several research laboratories worldwide.

The development started in 2010 when Matteo Cereda wrote and published the first version of [RNAmotifs](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-1-r20). In 2014, Gregor Rot refactored the RNAmotifs code to Python and created a new branch of the software adding cluster analysis and integrating the software with [expressRNA](http://expressRNA.org). This cluster branch of the software is now called RNAmotifs2.

## Citing RNAmotifs2

[High-resolution RNA maps suggest common principles of splicing and polyadenylation regulation by TDP-43](http://www.cell.com/cell-reports/abstract/S2211-1247(17)30522-3)<br />
Rot, G., Wang, Z., Huppertz, I., Modic, M., Lenče, T., Hallegger, M., Haberman, N., Curk, T., von Mering, C., Ule, J.<br />
Cell Reports , Volume 19 , Issue 5 , 1056 - 1067

## Reporting problems

Use the [issues page](https://github.com/grexor/rnamotifs2/issues) to report issues and leave suggestions.
