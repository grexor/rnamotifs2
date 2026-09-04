"""
Compatibility layer over pybio's genome API.

Keeps rnamotifs2 working across pybio versions without modifying pybio:

  * older pybio exposed ``pybio.genomes``; current pybio exposes
    ``pybio.core.genomes``
  * ``seq_direct()`` now takes ``(species, chr, strand, start, stop,
    genome_version=...)`` instead of a single ``genome`` string such as
    ``"hg19"``

A comparison's ``genome=`` config value is resolved here into the
``(species, genome_version)`` pair pybio expects.
"""

import pybio

try:
    _g = pybio.genomes
except AttributeError:
    _g = pybio.core.genomes

# legacy genome-string aliases -> (species, genome_version)
_ALIASES = {
    "hg19": ("homo_sapiens", "ensembl115"),
    "hg38": ("homo_sapiens", "ensembl115"),
    "grch38": ("homo_sapiens", "ensembl115"),
    "mm10": ("mus_musculus", "ensembl115"),
    "mm39": ("mus_musculus", "ensembl115"),
}


def resolve(genome):
    """Map a config ``genome=`` value to ``(species, genome_version)``.

    Accepts ``"species.genome_version"`` (e.g. ``"homo_sapiens.ensembl115"``),
    a bare species (``"homo_sapiens"`` -> pybio default version), or a legacy
    alias (``"hg19"``).
    """
    key = str(genome).strip()
    if key.lower() in _ALIASES:
        return _ALIASES[key.lower()]
    if "." in key:
        species, genome_version = key.split(".", 1)
        return species, genome_version
    return key, None


def make_motifs_nr(motif_size):
    return _g.make_motifs_nr(motif_size)


def seq_direct(genome, chr, strand, start, stop, flank="N"):
    species, genome_version = resolve(genome)
    return _g.seq_direct(species, chr, strand, start, stop, flank=flank,
                         genome_version=genome_version)
