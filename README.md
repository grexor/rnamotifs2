# RNAmotifs2: cluster motif analysis

RNAmotifs2 is a Python software for identifying clusters of motifs underlying regulatory principles of alternative splicing and alternative polyadenylation. RNAmotifs2 can be used as a standalone software (requires [pybio](https://github.com/grexor/pybio)), however it is also integrated with [expressRNA](http://expressRNA.org).

## Short description

Using high-through sequencing and bioinformatics analysis, we obtain sets of enhanced, repressed and control features (e.g. exons or polyA sites) comparing control and test experimental conditions. The research question is then of how to identify sequence short-motifs in regions surrounding the regulated features.

Initially, we search all features for all possible short motifs (3, 4, 5-mers). The detected motif signals are convoluted with a short sliding window (15nt half-window). This accounts for the fact that RNA-protein binding affinity is influenced by several closely-spaced short motifs. After the signal is filtered, the Fisher test is performed on each of the two comparisons (enhanced vs. control, repressed vs. control) and each of the motif signals (AAA, AAT, AAC, etc.) separately. At this step (after FDR), we identified the strongest motif best separating the sequences, i.e. one motif for the enhanced vs. control, another for the repressed vs. control comparison.

However, since several proteins can regulate pre-mRNA processing by binding simultaneously around regulated features, the features with the strongest identified motif signal are removed. This is then compensated by searching the remaining feature space with the already identified motif (or cluster) paired with all other possible short-motifs. The search is reiterated until we reach a cluster of max. 4 motifs.

Finally, we compute an enrichment score (ES) on the super-imposed sequences of all the features (exons, polyA sites) and draw a motif regulatory RNA-map.

## Authors

[RNAmotifs2](https://github.com/grexor/rnamotifs2) is maintained by [Gregor Rot](https://grexor.github.io) in collaboration with several research laboratories worldwide.

The development started in 2010 when Matteo Cereda wrote and published the first version of [RNAmotifs](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-1-r20). In 2014, Gregor Rot refactored the RNAmotifs code to Python and make another branch of the software adding cluster analysis and integrating the software with [expressRNA](http://expressRNA.org). This cluster branch of the software is now called RNAmotifs2.

## Citing RNAmotifs2

[High-resolution RNA maps suggest common principles of splicing and polyadenylation regulation by TDP-43](http://www.cell.com/cell-reports/abstract/S2211-1247(17)30522-3)<br />
Rot, G., Wang, Z., Huppertz, I., Modic, M., Lenče, T., Hallegger, M., Haberman, N., Curk, T., von Mering, C., Ule, J.<br />
Cell Reports , Volume 19 , Issue 5 , 1056 - 1067

## Reporting problems

Use the [issues page](https://github.com/grexor/rnamotifs2/issues) to report issues and leave suggestions.
