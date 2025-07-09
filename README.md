# eris
Uncovering IS-mediated discord in bacterial genomes


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
1. Finds IS elements with Minimap2.
1. Predicts ORFs if needed
1. Quantifies the relative relationship between each IS element and flanking/overlappign ORFs.

