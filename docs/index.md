# `eris` 🧬🧞‍♀🔮️
Uncovering IS-mediated discord in bacterial genomes

[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/eris.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/eris)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/eris?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/eris)
[![Wheel](https://img.shields.io/pypi/wheel/eris.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/eris/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/eris.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/eris/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/eris.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/eris/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/tomdstanton/eris/)
[![Issues](https://img.shields.io/github/issues/tomdstanton/eris.svg?style=flat-square&maxAge=600)](https://github.com/tomdstanton/eris/issues)
[![Docs](https://img.shields.io/readthedocs/eris/latest?style=flat-square&maxAge=600)](https://eris.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/tomdstanton/eris/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/eris?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/eris)

## Introduction 🌐
`eris` is a Python package for finding IS elements in bacterial genomes and quantifying their effect on other genes.
IS elements are known to move, disrupt and even promote genes and whilst there are many tools to find IS elements in
genomes, few attempt to report the resulting effects. 

Like many bioinformatics tools `eris` is designed to work from 
the command-line, but is built on top of a robust API with few dependencies, and can be easily installed and
incorporated into other programs, scripts and pipelines.

## Installation ⚙️

### Requires 🧰
```
python >=3.9
minimap2 >=2.18
pyrodigal >=3.5.0 (for ORF prediction only)
```
**NOTE: eris is not yet on PyPI or Bioconda, please install from source until it is released**

### From source:
```shell
# First clone the repo
git clone https://github.com/tomdstanton/eris.git && cd eris
# Then install with pip
pip install .  # -e for editable, developers only!
# or install with pixi
pixi install
```

**NOTE: For [Pyrodigal](https://pyrodigal.readthedocs.io/en/stable/), you should install the `orf` or `dev`
environments with `pip` or [`pixi`](https://pixi.sh/dev/)**

```shell
# First clone the repo
git clone https://github.com/tomdstanton/eris.git && cd eris
# Then install with pip
pip install .[orf]  # -e for editable, developers only!
# or install with pixi
pixi install -e orf
```

## Usage 🧑‍💻
The information below explains how to use the `eris` CLI. 
For API usage, please refer to the [reference documentation](https://tomdstanton.github.io/eris/reference/eris/).

### Scan 🔍

#### Quickstart
`eris scan *.{fasta,gfa,gb} > results.tsv`

#### Arguments
```shell
usage: eris scan <genome> <genome...> [options]

========================|> eris |>========================
             Scan for IS in bacterial genomes             

Inputs:
  
  Note, input file(s) may be compressed
  Note, Genome(s) in FASTA/GFA format can paired up with GFA/BED
  annotation files with the same prefix.

  <genome>              Genome(s) in FASTA, GFA or Genbank format;
                        reads from stdin by default.
  -a [ ...], --annotations [ ...]
                        Optional genome annotations in GFF3/BED format;
                        These will be matched up to input genomes (FASTA/GFA) 
                        with corresponding filenames

Outputs:
  
  Note, text outputs accept "-" or "stdout" for stdout
  If a directory is passed, individual files will be written per input genome

  --tsv []              Path to output tabular results (default: stdout)
  --ffn []              Path to output Feature DNA sequences in FASTA format;
                        defaults to "./[genome]_eris_results.ffn" when passed without arguments.
  --faa []              Path to output Feature Amino acid sequences in FASTA format;
                        defaults to "./[genome]_eris_results.faa" when passed without arguments.
  --no-tsv-header       Suppress header in TSV output

Other options:

  --progress            Show progress bar
  -v, --version         Show version number and exit
  -h, --help            Show this help message and exit

For more help, visit: eris.readthedocs.io
```

#### The algorithm 
1. Given a bacterial genome as an assembly (FASTA), assembly-graph (GFA) or annotation file (Genbank), the `scan` pipeline
will align IS element nucleotide sequences from the [ISFinder](https://isfinder.biotoul.fr/) database against the 
assembly contigs using [minimap2](https://lh3.github.io/minimap2/).
1. These alignments are then sorted by their target contig, and culled such that each region aligned contains the highest 
scoring query.
1. Each Element alignment is then considered to be a "mobile-element" Feature, and added to the list
of Features on the respective contig.
1. If the genome is from a sequence file (FASTA/GFA), ORFs are predicted with 
[Pyrodigal](https://pyrodigal.readthedocs.io/en/stable) and CDS Features are added to each contig.
1. The genome is then converted into a **Feature graph**, whereby Features on each contig, sorted by their respective
start coordinates, are connected to their flanking Features; and if the contig is connected to other contigs 
(GFA input), Features on the termini of connected contigs are also connected to each other.
1. Promoters are searched for in each Element Feature using a regular expression. 
1. For each Element Feature, the **[Breadth-first search](https://en.wikipedia.org/wiki/Breadth-first_search) (BFS)**
algorithm traverses the **Feature graph** to find CDS that either **overlap** (part of the element) or **flank**
the element.
1. The relative effect of the Element on each flanking CDS Feature is predicted.

#### Performance 
`eris scan` is very fast, especially when providing annotations or once the 
[`pyrodigal.GeneFinder`](https://pyrodigal.readthedocs.io/en/stable/api/gene_finder.html#genefinder) instance has been
trained (this occurs on the first input genome); it should only take <1 second per assembly 🚀

### Outputs

The main output of `eris scan` is the TSV tabular result which are written to `stdout` by default. Sequence information
can be written to separate files via the respective CLI flags.

#### TSV tabular output

The TSV tabular output reports one line per **Feature** of interest, which can either be the Element itself from the
resulting alignments, the CDS _inside_ the Element, or the CDS _flanking_ the Element. If the context is the
Element, information about the element from [ISFinder](https://isfinder.biotoul.fr/) will be reported. If the context is a CDS, information about
the ORF/translation will be reported.

The TSV columns are as follows:

1. **Genome**: The name of the input genome.
1. **Feature**: The unique identifier of the Feature in question.
1. **Type**: The [Feature type](https://www.insdc.org/submitting-standards/feature-table/), currently only CDS are 
supported if annotations are provided; Elements are annotated as "`mobile_element`" and promoters are annotated as
"`regulatory`".
1. **Contig**: The name of the contig the Element is on.
1. **Start**: The start coordinate (0-based) of the Feature.
1. **End**: The end coordinate of the Feature.
1. **Strand**: The strand of the Feature (1 or -1).
1. **Partial**: Whether the Feature overlaps with the start or end the contig.
1. **Element**: The unique identifier of the Element from the current context.
1. **Element_distance**: Signifies the distance of the Feature from the Element from the current context.
1. **Element_location**: Signifies the relative location of the Element from the current context.
1. **Element_strand**: Signifies the relative strand of the Element from the current context.
1. **Element_effect**: Signifies the relative effect of the Element from the current context on the CDS.
1. **Percent_identity**: The percent identity of the Element.
1. **Percent_coverage**: The percent coverage of the Element.
1. **Name**: The name of the Element if known (from ISFinder).
1. **Family**: The family of the Element if known (from ISFinder).
1. **Group**: The group of the Element if known (from ISFinder).
1. **Synonyms**: The synonyms of the Element if known (from ISFinder).
1. **Origin**: The origin of the Element if known (from ISFinder).
1. **IR**: The relative location of the inverted repeat of the Element if known (from ISFinder).
1. **DR**: The relative location of the direct repeat of the Element if known (from ISFinder).

#### Example output
This is a single Element context in TSV format from a _K. pneumoniae_ genome, supplied from a GFA and BED file.
You can see that the Element (`0db94675-30ef-4b28-a212-f26cbad7409e`) contains one CDS, one promoter and is of the
IS1380 family. On the same strand downstream is one CDS, annotated as a CTX-M-15 gene, known to confer ESBL resistance.
As the Element contains a promoter and is on the same strand downstream, and in close proximity (48bp),
we predict this gene is being upregulated.

| Genome     | Feature                                         | Type           | Contig | Start | End  | Strand | Partial | Element                              | Element_distance | Element_location | Element_strand | Element_effect | Percent_identity | Percent_coverage | Name                                              | Family | Group | Synonyms | Origin           | IR    | DR |
| ---------- | ----------------------------------------------- | -------------- | ------ | ----- | ---- | ------ | ------- | ------------------------------------ | ---------------- | ---------------- | -------------- | -------------- | ---------------- | ---------------- | ------------------------------------------------- | ------ | ----- | -------- | ---------------- | ----- | -- |
| ERR4920392 | KNDCPA_05188                                    | CDS            | 65     | 834   | 1710 | \-1    | FALSE   | 0db94675-30ef-4b28-a212-f26cbad7409e | 48bp             | downstream       | same strand    | upregulated    | \-               | \-               | extended-spectrum class A beta-lactamase CTX-M-15 | \-     | \-    | \-       | \-               | \-    | \- |
| ERR4920392 | 0db94675-30ef-4b28-a212-f26cbad7409e            | mobile_element | 65     | 1758  | 2157 | \-1    | TRUE    | 0db94675-30ef-4b28-a212-f26cbad7409e | \-               | \-               | \-             | \-             | 100              | 24.0942029       | ISEc9                                             | IS1380 |       | None     | Escherichia coli | 13/22 | NA |
| ERR4920392 | KNDCPA_05189                                    | CDS            | 65     | 1862  | 1985 | \-1    | FALSE   | 0db94675-30ef-4b28-a212-f26cbad7409e | \-               | inside           | same strand    | \-             | \-               | \-               | MobQ family relaxase                              | \-     | \-    | \-       | \-               | \-    | \- |
| ERR4920392 | 0db94675-30ef-4b28-a212-f26cbad7409e_promoter_1 | regulatory     | 65     | 2053  | 2083 | \-1    | FALSE   | 0db94675-30ef-4b28-a212-f26cbad7409e | \-               | inside           | \-             | \-             | \-               | \-               | \-                                                | \-     | \-    | \-       | \-               | \-    | \- |

### Pan 🦘
`eris pan` quantifies the effect of IS-mediated events in pan-genome graphs.

Coming soon!

### Map 🗺️
Re-implementation of [ISMapper](https://github.com/jhawkey/IS_mapper), with reference-based and reference-free options.

Coming soon!

![eris](https://static.wikia.nocookie.net/dreamworks/images/6/6f/Sinbad-disneyscreencaps.com-1100.jpg/revision/latest?cb=20240311205845)
