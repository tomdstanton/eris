"""
Copyright 2025 Tom Stanton (tomdstanton@gmail.com)
https://github.com/tomdstanton/eris

This file is part of eris. eris is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. eris is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with eris.
If not, see <https://www.gnu.org/licenses/>.
"""
from typing import Literal, Union, IO
from os import fstat
from io import IOBase
from warnings import warn
from concurrent.futures import Executor, ThreadPoolExecutor
from pathlib import Path
from dataclasses import asdict
from collections import deque

from eris import ErisWarning, Requires, RESOURCES
from eris.db import Database
from eris.io import Genome, SeqFile, GeneFinderConfig
from eris.seq import Record, Feature
from eris.alignment import cull_all, group_alignments
from eris.external import Minimap2, Minimap2AlignConfig, Minimap2IndexConfig
from eris.graph import Edge


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
        # # self.contig_gc: float = contig.seq.GC()
        # self.contg_copy_number: float = next((v for k, v in contig.qualifiers if k == 'dp'), 0)
        self.CDS_in_element: set[Feature] = set()
        self.CDS_flanking_element: set[Feature] = set()
        self.edges: set[Edge] = set()

    def __repr__(self):
        return self.feature.__repr__()

    def __str__(self):
        return self.feature.__str__()

    def __len__(self):
        return len(self.feature)

    def __format__(self, __format_spec: Literal['tsv'] = '') -> str:
        if __format_spec == '':
            return self.__str__()
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

    def __format__(self, __format_spec: Literal['tsv'] = '') -> str:
        if __format_spec == '':
            return self.__str__()
        elif __format_spec == 'tsv':
            return ''.join(format(i, __format_spec) for i in self.insertion_sequences)
        else:
            raise NotImplementedError(f'Invalid format: {__format_spec}')


class ScannerWarning(ErisWarning):
    pass


class ScannerError(Exception):
    pass

class ScannerResultWriter:
    def __init__(
                self,
                tsv: Union[str, Path, IO] = None, bed: Union[str, Path, IO] = None, gfa: Union[str, Path, IO] = None,
                fna: Union[str, Path, IO] = None, ffn: Union[str, Path, IO] = None, faa: Union[str, Path, IO] = None,
                png: Union[str, Path, IO] = None, suffix: str = '_eris_results', no_tsv_header: bool = False
        ):
        self._files = {}
        self._handles = {}
        self._suffix: str = suffix
        self._tsv_header: str = ('Genome\tIS_element\tContext\tFeature\tType\tContig\tStart\tEnd\tStrand\tPartial_start\t'
                                 'Partial_end\tTranslation_start\tTranslation_end\tPercent_identity\tPercent_coverage\t'
                                 'Name\tFamily\tGroup\tSynonyms\tOrigin\tIR\tDR\n')
        self._no_tsv_header: bool = no_tsv_header
        self._results_written: int = 0
        for k, v in (('tsv', tsv), ('bed', bed), ('gfa', gfa), ('fna', fna), ('ffn', ffn), ('faa', faa),
                     ('png', png)):
            if v is not None:
                if isinstance(v, str):
                    self._files[k] = Path(v)
                elif isinstance(v, Path):
                    self._files[k] = v
                elif isinstance(v, IOBase):
                    self._handles[k] = v
                else:
                    raise TypeError(f"Argument {k} must be a string, Path or IO, not {type(v)}")

    def __len__(self):
        return self._results_written

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def open(self, *args, **kwargs):
        """Opens the handles that aren't stdout or stderr"""
        for handle in self._handles.values():
            if handle.name not in {'<stdout>', '<stderr>'}:
                handle.open(*args, **kwargs)

    def close(self):
        """Close the handles that aren't stdout or stderr"""
        for handle in self._handles.values():
            if handle.name not in {'<stdout>', '<stderr>'}:
                handle.close()

    def write(self, result: Union[ScannerResult, None]):
        """Write the typing result to files or file handles."""
        if result is None:
            return None

        for fmt, file in self._files.items():  # Write to the file outputs
            with open(file / f'{result.genome_id}{self._suffix}.{fmt}', 'wt') as handle:
                handle.write(format(result, fmt))

        for fmt, handle in self._handles.items():  # Write to the handle outputs
            # Logic for non-repeating TSV header, this limits TSV to a handle but should be the case
            if fmt == 'tsv' and not self._no_tsv_header and self._results_written == 0:
                handle.write(self._tsv_header)
            handle.write(format(result, fmt))

        self._results_written += 1  # Incremment the counter


class Scanner:
    def __init__(self, align_config: Minimap2AlignConfig = None, index_config: Minimap2IndexConfig = None,
                 gene_finder_config: GeneFinderConfig = None, pool: Executor = None):
        self.db = Database()
        self.align_config = align_config or Minimap2AlignConfig(c=True)
        self.index_config = index_config or Minimap2IndexConfig()
        self.gene_finder_config = gene_finder_config or GeneFinderConfig()
        self.pool = pool  # Only init the pool if we need it for the gene finder
        self._gene_finder: 'pyrodigal.GeneFinder' = None

    @property
    @Requires(requires='pyrodigal')
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

    def __call__(self, *genomes: Genome):
        yield from map(self._pipeline, genomes)

    def _pipeline(self, genome: Union[str, Path, SeqFile, Genome]) -> Union[ScannerResult, None]:
        """
        The scanner pipeline takes a single bacterial genome, predicts ORFs if necessary, and uses Minimap2 to find IS
        elements with nucleotide-nucleotide alignment against contigs. Then, ORFs belonging to- or flanked by- the IS
        elements are identified (including potential graph traversal), and are returned in a result object for later
        reporting.
        """
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
            return warn('No alignments found', ISScannerWarning)

        if not genome.is_annotated:  # Get ORFs if neccessary, only do this if we have alignments to exit early
            genes = {gene.id: gene for gene in genome.find_genes(self.gene_finder, self.pool)}
        else:
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
