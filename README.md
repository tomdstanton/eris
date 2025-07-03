# insta :camera::next_track_button::dna:
Annotate IS elements and interrupted genes in bacterial genomes


## Background
Insertion sequence (IS) elements facilitate the means to rapidly alter bacterial phenotypes, without
the aquisition of new genes.
Studies have shown that IS elements can insert themselves within the open-reading frames (ORFs) of genes,
resulting in the transcriptional/translational loss of the gene product.
For example, in Klebsiella pneumoniae, insertions disrupting mgrB can confer colistin resistance, while insertions 
disrupting wcaJ can "switch off" the capsule.
Inversley, IS elements can harbour strong promoters, which when interted upstream of an ORF, can result in the
over-expression of the gene product. This phenomenon has been described in Acinetobacter baumannii whereby
ISAba1 insertions upstream of blaOXA-51 confer carbapenem resistance.
The repetitive nature of IS elements make them notoriously

## `scan`

### What it does
1. Takes bacterial assemblies in FASTA, GFA or GenBank format.
1. Finds IS elements.
1. Predicts the effect each IS element is having.

### How it does it
- If assembly graphs are supplied,



### Interpreting the results
The main output is tabular, however, 



## `insert`
<img src="https://i.kym-cdn.com/entries/icons/original/000/041/943/1aa1blank.png" width="150">

### What it does
1. Takes bacterial assemblies in FASTA, GFA or GenBank format.
1. Finds IS elements.
1. Predicts the effect each IS element is having.

### How it does it
- If assembly graphs are supplied,

### Outputs
**All** sequence based outputs are 

IS elements with no nucleotide alignments to the IS database, but with HMM domain hits are suffixed with
'novel', in case you wish to submit them later.

## `remove`
<img src="https://freepngimg.com/thumb/walter_white/134277-white-walter-free-transparent-image-hq.png" width="100">

### What it does
1. Takes bacterial assemblies in FASTA, GFA or GenBank format.
1. Finds IS elements.
1. Predicts the effect each IS element is having.

### How it does it
- If assembly graphs are supplied,

## Notes
Classes are built on top of Biopython's `SeqRecord` and `SeqFeature` objects
to allow the seamless conversion of results between different biological
text formats.

## Collaborate!
If there is an IS that hasn't been detected, please let us know so we
can develop better detection methods!

Similarly, if there is an IS-related event that would be good to
detect and report, please reach out!
