"""
Module for running the scan pipeline on bacterial genomes.
"""
from typing import Literal, Union, Generator, IO
from warnings import warn
from concurrent.futures import Executor, ThreadPoolExecutor
from pathlib import Path
from dataclasses import asdict
from collections import deque
from itertools import chain

from eris import ErisWarning, RESOURCES, require
from eris.db import Database
from eris.io import Genome, SeqFile, GeneFinderConfig
from eris.seq import Record, Feature
from eris.alignment import cull_all, group_alignments
from eris.external import Minimap2, Minimap2AlignConfig, Minimap2IndexConfig
from eris.graph import Edge


# Constants ------------------------------------------------------------------------------------------------------------
TSV_HEADER = ('Genome\tIS_element\tContext\tFeature\tType\tContig\tStart\tEnd\tStrand\tPartial_start\t'
                                 'Partial_end\tTranslation_start\tTranslation_end\tPercent_identity\tPercent_coverage\t'
                                 'Name\tFamily\tGroup\tSynonyms\tOrigin\tIR\tDR\n')

# Classes --------------------------------------------------------------------------------------------------------------
class InsertionSequence:
    """
    Class to store an insertion sequence feature from an alignment and relevant CDS features.

    Attributes:
        feature: Feature instance representing the IS element
        CDS_in_element: set[Feature] CDS contained within the element
        CDS_flanking_element: set[Feature] CDS flanking the element on the same or connected contigs
        edges: set[Edge] edges in the feature graph connected to the IS element feature
    """
    def __init__(self, genome_id: str, feature: Feature, contig: Record):
        self.genome_id = genome_id
        self.feature: Feature = feature
        self.CDS_in_element: set[Feature] = set()
        self.CDS_flanking_element: set[Feature] = set()
        self.edges: set[Edge] = set()

    def __repr__(self):
        return self.feature.__repr__()

    def __str__(self):
        return self.feature.__str__()

    def __len__(self):
        return len(self.feature)

    def __format__(self, __format_spec: Literal['tsv', 'fna', 'bed', 'faa', 'ffn'] = '') -> str:
        if __format_spec == '':
            return self.__str__()
        elif __format_spec in {'fna', 'bed'}:
            return ''.join(
                format(i, __format_spec) for i in chain([self.feature], self.CDS_in_element, self.CDS_flanking_element))
        elif __format_spec in {'faa', 'ffn'}:
            return ''.join(
                format(i, __format_spec) for i in chain(self.CDS_in_element, self.CDS_flanking_element))
        elif __format_spec == 'tsv':
            results = [
                f'{self.genome_id}\t{self.feature.id}\tAlignment\t{self.feature:tsv}\t'
                f'{self.feature.location.partial_start}\t{self.feature.location.partial_end}\t-\t-\t'
                f'{self.feature['identity']}\t{self.feature['query_coverage']}\t'
                f'{self.feature['Name']}\t{self.feature['Family']}\t{self.feature['Group']}\t'
                f'{self.feature['Synonyms']}\t{self.feature['Origin']}\t{self.feature['IR']}\t{self.feature['DR']}\n'
            ]
            for feature in self.CDS_in_element:
                translation = feature.translate()  # Feature should aready have translation Qualifier
                results.append(
                    f'{self.genome_id}\t{self.feature.id}\tCDS_in_element\t{feature:tsv}\t'
                    f'{feature.location.partial_start}\t{feature.location.partial_end}\t'
                    f'{translation[0]}\t{translation[-1]}\t'
                    f'-\t-\t-\t-\t-\t-\t-\t-\t-\n'
                )
            for feature in self.CDS_flanking_element:
                translation = feature.translate()  # Feature should aready have translation Qualifier
                results.append(
                    f'{self.genome_id}\t{self.feature.id}\tCDS_flanking_element\t{feature:tsv}\t'
                    f'{feature.location.partial_start}\t{feature.location.partial_end}\t'
                    f'{translation[0]}\t{translation[-1]}\t'
                    f'-\t-\t-\t-\t-\t-\t-\t-\t-\n'
                )
            return ''.join(results)
        else:
            raise NotImplementedError(f'Invalid format: {__format_spec}')


class ScannerResult:
    def __init__(self, genome_id: str):
        self.genome_id = genome_id
        self.insertion_sequences: list[InsertionSequence] = []

    def __repr__(self):
        return f'ScannerResult({self.genome_id})'

    def __iter__(self):
        return iter(self.insertion_sequences)

    def __len__(self):
        return len(self.insertion_sequences)

    def __format__(self, __format_spec: Literal['tsv', 'fna', 'bed', 'faa', 'ffn'] = '') -> str:
        if __format_spec == '':
            return self.__str__()
        elif __format_spec in {'tsv', 'fna', 'bed', 'faa', 'ffn'}:
            return ''.join(format(i, __format_spec) for i in self.insertion_sequences)
        else:
            raise NotImplementedError(f'Invalid format: {__format_spec}')


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
        db: The ISFinder sequence eris.Database object
        align_config: Minimap2 alignment configuration
        index_config: Minimap2 index configuration
        gene_finder_config: Configuration for the pyrodigal.GeneFinder instance if needed
        pool: A ThreadPoolExecutor or ProcessPoolExecutor instance if needed

    Examples:
        >>> from eris.scanner import Scanner
        ... from pathlib import Path
        ... with Scanner() as scanner:
        ...     results = list(scanner.scan(*Path.iterdir('genomes/')))
    """
    def __init__(self, align_config: Minimap2AlignConfig = None, index_config: Minimap2IndexConfig = None,
                 gene_finder_config: GeneFinderConfig = None, pool: Executor = None):
        self.db: Database = Database()
        self.align_config: Minimap2AlignConfig = align_config or Minimap2AlignConfig(c=True)
        self.index_config: Minimap2IndexConfig = index_config or Minimap2IndexConfig()
        self.gene_finder_config: GeneFinderConfig = gene_finder_config or GeneFinderConfig()
        self.pool: Executor = pool  # Only init the pool if we need it for the gene finder
        self._gene_finder: 'pyrodigal.GeneFinder' = None

    @property
    @require('pyrodigal')
    def gene_finder(self):
        if not self._gene_finder:
            from pyrodigal import GeneFinder
            self._gene_finder = GeneFinder(**asdict(self.gene_finder_config))
        if self.pool is None:
            self.pool = ThreadPoolExecutor(min(32, RESOURCES.available_cpus + 4))
        return self._gene_finder

    def cleanup(self):
        if self.pool:
            self.pool.shutdown(wait=False, cancel_futures=True)
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

        contigs_with_is = []
        with Minimap2(targets=genome, index_config=self.index_config) as aligner:  # Align IS to contigs
            for contig, alignments in group_alignments(aligner.align(self.db, config=self.align_config), 'target'):
                genome[contig].add_features(*(i.as_feature('mobile_element') for i in cull_all(alignments)))
                contigs_with_is.append(contig)  # Add contig name to list, also indicator that we have alignments

        if not contigs_with_is:  # No alignments found, warn and exit early
            return warn('No alignments found', ScannerWarning)

        if not genome.is_annotated:  # Get ORFs if neccessary, only do this if we have alignments to exit early
            genes = {gene.id: gene for gene in genome.find_genes(self.gene_finder, self.pool)}
        else:  # Get the existing annotations from the genome
            # TODO: Should we support other feature types?
            genes = {gene.id: gene for contig in genome.contigs.values() for gene in contig.features if gene.kind == 'CDS'}

        graph = genome.as_feature_graph()  # Get the gene graph for traversing connected genes
        result = ScannerResult(genome.id)  # Init the result container

        for contig in contigs_with_is:  # Iterate over contigs with IS element alignments
            for feature in (contig := genome[contig]):  # type: Feature
                if feature.kind == 'mobile_element':  # This is an IS element, let's see if it overlaps any ORFs
                    feature.qualifiers += self.db[feature['name']].qualifiers  # Pull over the qualifiers from the database
                    result.insertion_sequences.append(is_element := InsertionSequence(genome.id, feature, contig))
                    # This block performs a graph traversal (BFS) starting from the neighbors of the IS element
                    # to find all genes that are part of the IS element vs. those that are merely flanked by it.

                    # A queue for edges to visit. Start with the direct neighbors of the IS element.
                    # We use a deque for efficient appends and pops from the left (FIFO).
                    edges_to_visit = deque(graph.get_neighbors(feature.id))
                    # A set to keep track of edges we've already processed to avoid cycles and redundant work.
                    processed_edges = set(edges_to_visit)
                    processed_genes = set()  # <-- Add this to track visited genes

                    while edges_to_visit:
                        edge = edges_to_visit.popleft()
                        is_element.edges.add(edge)  # Use .add() for the new set
                        if not (gene := genes.get(edge.to)):
                            continue

                        if gene.id in processed_genes:  # <-- Check if we've already processed this gene
                            continue
                        processed_genes.add(gene.id)  # <-- Mark this gene as processed

                        # Check if the gene is substantially part of the IS element.
                        # A gene can only be part of the IS element if it's on the same contig.
                        # Otherwise, it's a flanking gene by definition.
                        if (gene.location.parent_id == feature.location.parent_id and
                                feature.overlap(gene) / len(gene) >= 0.8):
                            # This gene is part of the IS. Add it to the 'overlaps' set.
                            is_element.CDS_in_element.add(gene)  # Use .add()
                            # Since it's part of the IS, we need to explore its neighbors as well
                            # to find the full extent of the element.
                            for next_edge in graph.get_neighbors(gene.id):
                                if next_edge not in processed_edges:
                                    processed_edges.add(next_edge)
                                    edges_to_visit.append(next_edge)
                        else:
                            # This gene is not part of the IS (it's on another contig or doesn't overlap
                            # sufficiently), so it's a flanking gene.
                            is_element.CDS_flanking_element.add(gene)  # Use .add()
                            # We stop the traversal down this path by not adding its neighbors to the queue.
        return result
