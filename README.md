# `eris` 🧬🧞‍♀🔮️
Uncovering IS-mediated discord in bacterial genomes

## Introduction
`eris` is a Python package for finding IS elements in bacterial genomes and quantifying their effect on other genes.
IS elements are known to move, disrupt and even promote genes and whilst there are many tools to find IS elements in
genomes, few attempt to report the resulting effects. 

Like many bioinformatics tools `eris` is designed to work from 
the command-line, but is built on top of a robust API with few dependencies, and can be easily installed and
incorporated into other programs, scripts and pipelines.

## Installation

### Requires
```
python >=3.9
minimap2 >=2.18
pyrodigal >=3.5.0 (for ORF prediction only)
```
**NOTE: eris is not yet on PyPI or Bioconda, please install from source until it is released**

### From source:
```shell
git clone https://github.com/tomdstanton/eris.git
cd eris
pip install -e .
# OR
pixi install
```

## Usage
The information below explains how to use the `eris` CLI. For API usage, please refer to the `docs/`.

### Scan

#### Quickstart
`eris scan *.{fasta,gfa,gb} > results.tsv`

#### Arguments
```shell
usage: eris scan <genome> <genome...> [options]

========================|> eris |>========================
             Scan for IS in bacterial genomes             

Inputs:
  
  Note, input file(s) may be compressed.

  <genome>         Genome(s) in FASTA, GFA or Genbank format; reads from stdin by default.
                   Genome(s) in FASTA/GFA format can paired up with GFA/BED
                   annotation files with the same prefix.

Outputs:
  
  Note, text outputs accept "-" or "stdout" for stdout
  If a directory is passed, individual files will be written per input genome

  --tsv []         Path to output tabular results (default: stdout)
  --ffn []         Path to output feature DNA sequences in FASTA format
                   Defaults to "./[genome]_eris_results.ffn" when passed without arguments
  --faa []         Path to output feature Amino acid sequences in FASTA format
                   Defaults to "./[genome]_eris_results.faa" when passed without arguments
  --no-tsv-header  Suppress header in TSV output

Other options:

  -v, --version    Show version number and exit
  -h, --help       Show this help message and exit

For more help, visit: eris.readthedocs.io
```

#### The algorithm 
- Given a bacterial genome as an assembly (FASTA), assembly-graph (GFA) or annotation file (Genbank), the `scan` pipeline
will align IS element nucleotide sequences from the ISFinder database against the assembly contigs using minimap2.
- These alignments are then sorted by their target contig, and culled such that each region aligned contains the highest 
scoring query.
- Each IS element alignment is then considered to be a "mobile-element" feature, and added to the list
of features on the respective contig.
- If the genome is from a sequence file (FASTA/GFA), ORFs are predicted with Pyrodigal and CDS features are added to
each contig.
- The genome is then converted into a **feature graph**, whereby features on each contig, sorted by their respective
start coordinates, are connected to their flanking features; and if the contig is connected to other contigs 
(GFA input), features on the termini of connected contigs are also connected to each other.
- For each IS element feature, the **[Breadth-first search](https://en.wikipedia.org/wiki/Breadth-first_search) (BFS)**
algorithm traverses the **feature graph** to find CDS that either **overlap** (part of the element) or **flank**
the element.

#### Performance 
`eris scan` is very fast, especially when providing annotations or once the `pyrodigal.GeneFinder` instance has been
trained (this occurs on the first input genome); it should only take <1 second per assembly 🚀

### Outputs

The main output of `eris scan` is the TSV tabular result which are written to `stdout` by default. Sequence information
can be written to separate files via the respective CLI flags.

#### TSV tabular output

The TSV tabular output reports one line per **feature** of interest, which can either be the IS element itself from the
resulting alignments, the CDS _inside_ the IS element, or the CDS _flanking_ the IS element. If the context is the
IS element, information about the element from ISFinder will be reported. If the context is a CDS, information about
the ORF/translation will be reported.

The TSV columns are as follows:

1. **Genome**: The name of the input genome.
1. **IS_element**: The unique identifier of the IS element.
1. **Context**: Signifies the IS element from alignment, an CDS inside or an CDS outside the IS element.
1. **Feature**: The unique identifier of the feature in question.
1. **Type**: The [feature type](https://www.insdc.org/submitting-standards/feature-table/), currently only CDS are 
supported if annotations are provided; IS elements are annotated as "`mobile_element`".
1. **Contig**: The name of the contig the IS element is on.
1. **Start**: The start coordinate (0-based) of the feature.
1. **End**: The end coordinate of the feature.
1. **Strand**: The strand of the feature (1 or -1).
1. **Partial_start**: Whether the CDS overlaps with the start of the sequence.
1. **Partial_end**: Whether the CDS overlaps with the end of the sequence.
1. **Translation_start**: The first amino acid residue of the CDS translation (should be `M`).
1. **Translation_end**: The last amino acid residue of the CDS translation (should be `*`).
1. **Percent_identity**: The percent identity of the IS element.
1. **Percent_coverage**: The percent coverage of the IS element.
1. **Name**: The name of the IS element if known (from ISFinder).
1. **Family**: The family of the IS element if known (from ISFinder).
1. **Group**: The group of the IS element if known (from ISFinder).
1. **Synonyms**: The synonyms of the IS element if known (from ISFinder).
1. **Origin**: The origin of the IS element if known (from ISFinder).
1. **IR**: The relative location of the inverted repeat of the IS element if known (from ISFinder).
1. **DR**: The relative location of the direct repeat of the IS element if known (from ISFinder).

#### Caveats and things to do
 - Exact promoter region quantification for flanking genes to add more evidence of gene promotion.
 - If providing multiple sequence file inputs from different species, the `pyrodigal.GeneFinder` instance
will be trained on the first genome, which may impact how well it finds ORFs in other species.
 - Consider implementing [FragGeneScanRs](github.com/unipept/FragGeneScanRs) instead of Pyrodigal 
(faster + works on reads)


### Map
Re-implementation of [ISMapper](https://github.com/jhawkey/IS_mapper), with reference-based and reference-free options.

Coming soon!


![eris](https://static.wikia.nocookie.net/dreamworks/images/6/6f/Sinbad-disneyscreencaps.com-1100.jpg/revision/latest?cb=20240311205845)
