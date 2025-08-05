"""
Module for running the scan pipeline on bacterial genomes.
"""
from _operator import attrgetter
from typing import Literal, Union, Generator, IO
from warnings import warn
from pathlib import Path
from dataclasses import asdict
from collections import deque

from eris import ErisWarning, require
from eris.db import Database
from eris.io import Genome, SeqFile, GeneFinderConfig
from eris.seq import Feature, Record, Qualifier
from eris.alignment import cull_all
from eris.external import Minimap2, Minimap2AlignConfig, Minimap2IndexConfig
from eris.utils import grouper


# Constants ------------------------------------------------------------------------------------------------------------
TSV_HEADER = (
    'Genome\tFeature\tType\tContig\tStart\tEnd\tStrand\tPartial\t'
    'Element\tElement_distance\tElement_location\tElement_strand\tElement_effect\t'
    'Percent_identity\tPercent_coverage\tName\tFamily\tGroup\tSynonyms\tOrigin\tIR\tDR\n'
)

# Classes --------------------------------------------------------------------------------------------------------------
class ScannerResult:
    def __init__(self, genome_id: str):
        self.genome_id = genome_id
        self.features: list[Feature] = []
        self.edges: list['Edge'] = []
        self.contexts: dict[str, Context] = {}

    def __repr__(self):
        return f'ScannerResult({self.genome_id})'

    def __iter__(self):
        return iter(self.features)

    def __len__(self):
        return len(self.features)
    #
    def __format__(self, __format_spec: Literal['tsv', 'fna', 'bed', 'faa', 'ffn'] = '') -> str:
        if __format_spec == '':
            return self.__str__()
        elif __format_spec in {'fna', 'faa', 'ffn', 'bed', 'gff'}:
            return ''.join(format(i, __format_spec) for i in self)
        elif __format_spec == 'tsv':
            lines = []
            for f in self.features:
                line = f'{self.genome_id}\t{f:tsv}\t{f.partial()}\t{self.contexts[f.id]}\t'
                if f.kind == 'mobile_element':
                    line += (f"{f['identity']}\t{f['query_coverage']}\t{f['Name']}\t{f['Family']}\t"
                             f"{f['Group']}\t{f['Synonyms']}\t{f['Origin']}\t{f['IR']}\t{f['DR']}\n")
                elif f.kind == 'CDS':
                    line += f"-\t-\t{f['Name']}\t-\t-\t-\t-\t-\t-\n"
                elif f.kind == 'regulatory':
                    line += f"-\t-\t-\t-\t-\t-\t-\t-\t-\n"
                lines.append(line)
            return ''.join(lines)
        else:
            raise NotImplementedError(f'Invalid format: {__format_spec}')


class Context:
    def __init__(self, id_: str, distance: int = None, location: str = '-', strand: str = '-', effect: str = '-'):
        self.id = id_
        self.distance = distance
        self.location = location
        self.strand = strand
        self.effect = effect

    def __str__(self):
        return f'{self.id}\t{'-' if self.distance is None else f"{self.distance}bp"}\t{self.location}\t{self.strand}\t{self.effect}'


class ScannerWarning(ErisWarning):
    pass


class ScannerError(Exception):
    pass


class Scanner:
    """
    A class that runs the eris scan pipeline on bacterial genomes. The scan pipeline takes a single bacterial genome,
    predicts ORFs if necessary, and uses Minimap2 to find IS elements with nucleotide-nucleotide alignment against
    contigs. Then, ORFs belonging to- or flanked by- the IS elements are identified
    (including potential graph traversal), and are returned in a result object for later reporting.

    Attributes:
        db (Database): The ISFinder sequence eris.Database object.
        align_config (Minimap2AlignConfig): Minimap2 alignment configuration.
        index_config (Minimap2IndexConfig): Minimap2 index configuration.
        gene_finder_config (GeneFinderConfig): Configuration for the pyrodigal.GeneFinder instance if needed.
        min_overlap (float): Minimum %overlap for a CDS to be part of an element.
        min_proximity (int): Minimum distance in bp for an element to have an effect on a CDS.

    Examples:
        >>> from eris.scanner import Scanner
        ... from pathlib import Path
        ... with Scanner() as scanner:
        ...     results = list(scanner.scan(*Path.iterdir('genomes/')))
    """
    def __init__(self, align_config: Minimap2AlignConfig = None, index_config: Minimap2IndexConfig = None,
                 gene_finder_config: GeneFinderConfig = None, pool: 'Executor' = None, min_overlap: float = 0.8,
                 min_proximity: int = 150):
        self.db: Database = Database()
        self.align_config: Minimap2AlignConfig = align_config or Minimap2AlignConfig(c=True)
        self.index_config: Minimap2IndexConfig = index_config or Minimap2IndexConfig()
        self.gene_finder_config: GeneFinderConfig = gene_finder_config or GeneFinderConfig()
        self.min_overlap: float = min_overlap
        self.min_proximity: int = min_proximity
        self._pool: 'Executor' = pool  # Only init the pool if we need it for the gene finder
        self._gene_finder: 'pyrodigal.GeneFinder' = None

    @property
    @require('pyrodigal')
    def gene_finder(self):
        if not self._gene_finder:
            from pyrodigal import GeneFinder
            self._gene_finder = GeneFinder(**asdict(self.gene_finder_config))
        return self._gene_finder

    def cleanup(self):
        self._gene_finder = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()

    def __del__(self):
        self.cleanup()

    def __call__(self, *args, **kwargs):
        return self._pipeline(*args, **kwargs)

    def scan(self, *genomes: Union[str, Path, IO, SeqFile, Genome]) -> Generator[ScannerResult, None, None]:
        """
        Runs the pipeline on input genomes

        Arguments:
            genomes: Input genomes to scan, can be strings or Paths to files, SeqFile instances or Genome instances

        Yields:
            ScannerResult instances for each genome processed or None
        """
        yield from map(self._pipeline, genomes)

    def _pipeline(self, genome: Union[str, Path, IO, SeqFile, Genome]) -> Union[ScannerResult, None]:
        if not isinstance(genome, Genome):
            try:
                genome = Genome.from_file(genome)
            except Exception as e:
                raise ScannerError(f"Could not parse {genome} as Genome") from e

        contigs_with_is: list[Record] = []
        with Minimap2(targets=genome, index_config=self.index_config) as aligner:  # Align IS to contigs
            for contig, alignments in grouper(aligner.align(self.db, config=self.align_config), 'target'):
                (contig := genome[contig]).add_features(*(i.as_feature('mobile_element') for i in cull_all(alignments)))
                contigs_with_is.append(contig)  # Add contig name to list, also indicator that we have alignments

        if not contigs_with_is:  # No alignments found, warn and exit early
            return warn('No alignments found', ScannerWarning)

        if not genome.is_annotated:  # Get ORFs if necessary, only do this if we have alignments to exit early
            if self.gene_finder is None:
                raise ScannerError(f"Could not predict ORFs in {genome}; is Pyrodigal installed?")
            genes = {gene.id: gene for gene in genome.find_genes(self.gene_finder, self._pool)}
        else:  # Get the existing annotations from the genome, TODO: Should we support other feature types?
            genes = {gene.id: gene for contig in genome.contigs.values() for gene in contig.features if gene.kind == 'CDS'}

        graph = genome.as_feature_graph()  # Get the gene graph for traversing connected genes
        result = ScannerResult(genome.id)  # Init the result container

        for contig in contigs_with_is:  # Iterate over contigs with IS element alignments
            for element in filter(lambda i: i.kind == 'mobile_element', contig):  # type: Feature

                # This block performs a graph traversal (BFS) starting from the neighbors of the IS element
                # to find all genes that are part of the IS element vs. those that are merely flanked by it.
                element.extract(contig, store_seq=True)  # Extract and store nucleotide sequence of element
                element.qualifiers += self.db[element['name']].qualifiers  # Pull over the qualifiers from the database
                element.qualifiers.append(Qualifier('mobile_element_type', 'insertion sequence'))
                result.features.append(element)  # Add element to result features
                result.contexts[element.id] = Context(element.id)

                promoter = None
                for promoter in element.find_promoters(contig):
                    result.features.append(promoter)  # Add element promoter to result features
                    result.contexts[promoter.id] = Context(element.id, location='inside')  # Add promoter to result context

                # A queue for edges to visit. Start with the direct neighbors of the IS element.
                # We use a deque for efficient appends and pops from the left (FIFO) and  a set to keep track of
                # edges we've already processed to avoid cycles and redundant work.
                processed_edges = set(edges_to_visit := deque(graph.get_neighbors(element.id)))

                while edges_to_visit:
                    edge = edges_to_visit.popleft()
                    if not (gene := genes.get(edge.to)) or gene.id in result.contexts:
                        continue  # Use result to check processed genes

                    result.edges.append(edge)  # Add edge to result for gfa output
                    result.features.append(gene)  # Add gene to result features
                    gene_context = Context(
                        element.id, location='flanking',
                        strand='same strand' if element.location.strand == gene.location.strand else 'opposite strand'
                    )
                    translation = gene.translate(parent=genome[gene.location.parent_id])  # Extract translation
                    # Calculate the context of the gene relative to the element.
                    # A gene can only be part of the IS element if it's on the same contig.
                    # Otherwise, it's a flanking gene by definition.
                    # TODO: Add stitching logic for elements that span multiple contigs
                    if gene.location.parent_id == element.location.parent_id:  # On the same contig
                        if element.overlap(gene) / len(gene) >= self.min_overlap:  # This gene is part of the IS element
                            gene_context.location = 'inside'  # Since it's part of the IS, we need to explore its
                            # neighbors as well to find the full extent of the element.
                            for next_edge in graph.get_neighbors(gene.id):
                                if next_edge not in processed_edges:
                                    processed_edges.add(next_edge)
                                    edges_to_visit.append(next_edge)

                        else:  # Non-overlapping gene on same contig
                            if element.location.start > gene.location.start: # Element comes after gene
                                gene_context.distance = element.location.start - gene.location.end
                                gene_context.location = 'upstream' if element.location.strand == 1 else 'downstream'
                            else:
                                gene_context.distance = gene.location.start - element.location.end
                                gene_context.location = 'downstream' if element.location.strand == 1 else 'upstream'

                            if gene_context.distance <= self.min_proximity:  # Predict effects
                                if promoter and gene_context.location == 'downstream' and gene_context.strand == 'same strand':
                                    gene_context.effect = 'upregulated'
                                elif not gene.partial() and (translation[0] not in {'M', 'V'} or translation[-1] != '*'):
                                    gene_context.effect = 'disrupted'

                    result.contexts[gene.id] = gene_context

        result.features.sort(key=attrgetter('location.parent_id', 'location.start'))
        return result
