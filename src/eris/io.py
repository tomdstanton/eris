"""
Module for parsing and managing bacterial sequence files and data.
"""
from functools import partial
from warnings import warn
from re import compile
from pathlib import Path
from typing import Union, Generator, IO, Literal, Iterable, TextIO, get_args, Callable
from concurrent.futures import Executor, ThreadPoolExecutor
from itertools import chain
from io import BufferedIOBase, RawIOBase
from shutil import copyfileobj
from tempfile import NamedTemporaryFile
from random import Random
from uuid import uuid4
from dataclasses import dataclass, asdict
from sys import stdin

from eris import ErisWarning, require, RESOURCES
from eris.alphabet import DNA
from eris.seq import Record, Feature, Seq, Qualifier, Location
from eris.graph import Graph, Edge
from eris.utils import xopen, Config, grouper

# Constants ------------------------------------------------------------------------------------------------------------
_SUPPORTED_FORMATS = Literal['fasta', 'gfa', 'genbank', 'fastq', 'gff']
_TAG2TYPE = {'f': float, 'i': int, 'Z' : str}  # See: https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#optional-fields
_FEATURE_KINDS = {'CDS'}  # 'gene', 'misc_feature', 'ncRNA', 'rRNA', 'regulatory', 'tRNA', 'tmRNA'}
_SEQUENCE_FILE_REGEX = compile(
    r'\.('
    r'(?P<fasta>f(asta|a|na|fn|as|aa))|'
    r'(?P<fastq>f(ast)?q)|'
    r'(?P<gfa>gfa)|'
    r'(?P<genbank>g(b|bff|bk|enbank))|'
    r'(?P<gff>gff(3)?)|'
    r'(?P<bed>bed)'
    r')\.?(?P<compression>(gz|bz2|xz|zst))?$'
)
# Regex for Illumina BaseSpace read files as specified here:
# https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
_ILLUMINA_READ_REGEX = compile(
    r'(?P<sample_name>.+)_S(?P<sample_number>\d+)_L(?P<lane_number>\d+)_R(?P<read_number>\d+)_(?P<index_number>\d+)$')
_SHORT_READ_REGEX = compile(r'(?P<sample_name>.+)_R?(?P<read_number>[12])$')
_TOPOLOGY_REGEX = compile(r'(?i)(\bcircular\b|\bcircular\s*=\s*true\b)')
# _COPY_NUMBER_REGEX = compile(r"depth=(\d+\.\d+)x")
_GENBANK_LOCATION_REGEX = compile(r'(?P<partial_start><)?(?P<start>[0-9]+)\.\.(?P<partial_end>>)?(?P<end>[0-9]+)')

# Classes --------------------------------------------------------------------------------------------------------------
@dataclass
class GeneFinderConfig(Config):
    meta: bool = False
    closed: bool = False
    mask: bool = False
    min_mask: int = 50
    min_gene: int = 90
    min_edge_gene: int = 60
    max_overlap: int = 60
    backend: str = "detect"


class ParserWarning(ErisWarning):
    pass

class SeqFileWarning(ErisWarning):
    pass


class SeqFileError(Exception):
    pass


class SeqFile:
    """
    Class for handling a (possibly compressed) file or stream of biological formats.

    :param file: Path to file (str/Path) or an IO stream (binary or text).
    :param format_: Format of the file. If None, the format will be guessed from the
                    filename extension or the stream content.
    :return: A SeqFile instance
    """
    def __init__(self, file: Union[str, Path, IO], format_: _SUPPORTED_FORMATS = None, temp_prefix='seqfile_'):
        self.id: str = "unknown"
        self.path: Union[Path, None] = None
        self.format: Union[str, None] = format_
        self._handle: Union[IO, None] = None
        self._from_stream: bool = False
        self._open_func: Union[None, Callable] = None

        if file == '-':  # Handle stdin symbol, the rest of the logic should deal with the stream
            file = stdin

        if hasattr(file, 'read') and not isinstance(file, (str, Path)):  # --- Handle IO Stream Input ---
            if not isinstance(file, (BufferedIOBase, RawIOBase)):
                if buffer := getattr(file, 'buffer', None):
                    file = buffer
                else:
                    raise SeqFileError("Input stream is not binary and lacks an accessible binary buffer (.buffer).")

            # The existing, clever solution: dump the stream to a temp file to make it seekable.
            with NamedTemporaryFile(prefix=temp_prefix, suffix='.tmp', delete=False, mode='wb') as temp_f_handle:
                copyfileobj(file, temp_f_handle)
                self.path = Path(temp_f_handle.name)
                self.id = self.path.stem
                self._from_stream = True

            # Now, we can safely open and guess from the temp file.
            if self.format is None:
                try:
                    # xopen handles decompression automatically based on magic numbers
                    with xopen(self.path, mode='rt') as f:
                        self.format = _guess_format_from_handle(f)
                except Exception as e:
                    # Clean up the temp file on failure
                    self.path.unlink(missing_ok=True)
                    raise e

            # The suffix for the temp file doesn't matter since we use xopen's 'magic' method
            self._open_func = partial(xopen, self.path, method='magic', mode='rt')

        elif isinstance(file, (str, Path)):  # --- Handle File Path Input ---
            self.path = file if isinstance(file, Path) else Path(file)
            if m := _SEQUENCE_FILE_REGEX.search(self.path.name):
                self.id = self.path.name.rstrip(m.group())
                # If format is not provided, guess from extension.
                self.format = format_ or next(fmt for fmt in get_args(_SUPPORTED_FORMATS) if m[fmt])
                self._open_func = partial(xopen, self.path, method=m['compression'] or 'uncompressed', mode='rt')
            else:
                # If no valid extension, try guessing from content
                if self.format is None:
                    try:
                        with xopen(self.path, mode='rt') as f:
                            self.format = _guess_format_from_handle(f)
                        self._open_func = partial(xopen, self.path, method='magic', mode='rt')
                    except Exception as e:
                         raise SeqFileError(f'Unsupported file extension and could not guess format for: {self.path.name}') from e
                else: # Format was provided, but extension is weird. Trust the user.
                    self._open_func = partial(xopen, self.path, method='magic', mode='rt')
        else:
            raise TypeError(f"Input must be a file path (str or Path) or an IO stream, not {type(file)}")

    def __repr__(self):
        return f'SeqFile({self.id}, format={self.format})'

    def __str__(self):
        return self.id

    def __enter__(self):
        self._ensure_handle_open()  # Open the handle when entering context, using the appropriate _open_func
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Ensure the primary handle (_handle) is closed
        self.close()
        # Clean up the temporary file *only* if it was created from a stream
        if self._from_stream:
            try:  # Use missing_ok=True for robustness against race conditions etc.
                self.path.unlink(missing_ok=True)
            except OSError as e:  # Log or warn about failure to delete temp file
                 warn(f"Could not delete temporary file {self.path}: {e}", SeqFileWarning)

    def __del__(self):  # Minimal cleanup: try to close the handle and delete temp file if applicable
        self.close()  # Subject to usual __del__ caveats (timing, interpreter state)
        if getattr(self, '_from_stream', False) and hasattr(self, 'path'):
             try:  # Avoid print/warn in __del__ if possible, or use sys.stderr
                 self.path.unlink(missing_ok=True)
             except Exception:
                 pass # Ignore errors in __del__

    def _ensure_handle_open(self):
        """Internal helper to open/reopen the file via _open_func."""
        if self._handle is None or self._handle.closed:
             self.close()  # Close any potentially lingering closed handle first
             self._handle = self._open_func()  # Open fresh using the configured function (xopen for paths/streams)

    def close(self):
        """Closes the associated file handle (_handle), if open."""
        if self._handle and not self._handle.closed:
            try:
                self._handle.close()
            except Exception as e:
                 pass   # Ignore potential errors during close, especially in __del__ context
        self._handle = None  # Mark as closed

    def read(self) -> str:
        """Reads the entire content of the file (decompressed)."""
        self._ensure_handle_open()  # Ensure we have a fresh handle from the start
        content = self._handle.read()  # Read from the current handle (which was just opened/rewound)
        # Close handle after full read? Optional, depends on expected usage.
        # self.close()
        return content

    def rewind(self):
        """Resets the file to be read from the beginning."""
        # The simplest way to rewind when using xopen dynamically is to close
        # the current handle and let _ensure_handle_open get a new one.
        self.close()
        # Next call to read() or __iter__() will automatically reopen via _ensure_handle_open()

    def __iter__(self) -> Generator[Union['Record', 'Edge', 'Feature'], None, None]:
        """Returns an iterator (generator) over records in the file."""
        self._ensure_handle_open() # Ensure handle is open/reopened
        if not self.format:
             raise SeqFileError("Cannot parse file: format is unknown.")
        yield from parse(self._handle, self.format)
        # Note: Iteration consumes the handle. Re-iteration requires rewind()/re-opening.

    def peek(self) -> str:
        """Returns the first line and resets to the beginning."""
        self._ensure_handle_open()
        first_line = self._handle.readline()
        self.rewind()  # Rewind by closing and letting the next operation reopen
        return first_line


class ReadFileError(Exception):
    pass


class ReadFile(SeqFile):
    """
    A subclass of :class:`~eris.core.io.SeqFile` that specifically handles a single file of sequencing reads
    in FASTQ format.
    This class extracts and stores read-specific information from the filename, and is designed to be used with the
    :class:`~eris.core.io.ReadSet` class.

    Attributes:
        sample_name: The name of the sample from which the reads were derived.
        read_type: The type of read ('short', 'illumina', 'long', 'pb', 'ont', or 'unknown').
        sample_number: The sample number (for Illumina reads).
        lane_number: The lane number (for Illumina reads).
        read_number: The read number (1 or 2 for paired-end reads).
        index_number: The index number (for Illumina reads).
    """
    def __init__(self, file: Union[str, Path]):
        """
        :param file: Path to file (str/Path) or an IO stream (binary or text).
        """
        super().__init__(file)
        self.sample_name = self.id
        self._read_type: Literal['unknown', 'short', 'illumina', 'long', 'pb', 'ont'] = 'unknown'
        self._sample_number: int = 0
        self._lane_number: int = 0
        self._read_number: int = 0
        self._index_number: int = 0
        if m := _ILLUMINA_READ_REGEX.match(self.id):
            self.sample_name = m.group('sample_name')
            self._read_type = "illumina"
            self._sample_number = int(m.group('sample_number'))
            self._lane_number = int(m.group('lane_number'))
            self._read_number = int(m.group('read_number'))
            self._index_number = int(m.group('index_number'))
        elif m := _SHORT_READ_REGEX.match(self.id):
            self.sample_name = m.group('sample_name')
            self._read_type = "short"
            self._read_number = int(m.group('read_number'))

    @property
    def read_type(self):
        return self._read_type

    @property
    def sample_number(self):
        return self._sample_number

    @property
    def lane_number(self):
        return self._lane_number

    @property
    def read_number(self):
        return self._read_number

    @property
    def index_number(self):
        return self._index_number


class ReadSetError(Exception):
    pass


class ReadSet:
    """
    Whilst still representing genomes, these are different from Genome instances as they will not hold
    sequence information in memory and may consist of multiple files.
    """
    def __init__(self, *files: Union[str, Path, ReadFile], id_: str = None):
        self.id = id_
        self._files = []
        read_types = set()
        for file in files:
            if not isinstance(file, ReadFile):
                file = ReadFile(file)
            if self.id is None:
                self.id = file.sample_name
            else:
                if self.id != file.sample_name:
                    raise ReadSetError(f'{file.id} sample name {file.sample_name} does not match ReadSet {self.id}')
            self._files.append(file)
            read_types.add(file.read_type)
        if len(read_types) > 1:
            raise ReadSetError(f'Hybrid read set {self.id} not yet supported')
        else:
            self._set_type = read_types.pop()

    @property
    def set_type(self):
        return self._set_type

    def __iter__(self):
        return iter(self._files)

    def __repr__(self):
        return self.id

    def __str__(self):
        return ' '.join([str(read) for read in self._files])

    def __len__(self):
        return len(self._files)


class GenomeError(Exception):
    pass

class Genome:
    """
    A class representing a single genome assembly in memory with contigs and potentially edges

    Attributes:
        id: The ID of the genome.
        contigs: A dictionary of contig IDs as keys and Record objects as values.
        edges: A list of Edges connecting the contigs.
        is_annotated: True if the genome has been loaded with annotations (e.g. from Genbank / GFF / BED).
    """
    def __init__(self, id_: str, contigs: dict[str: 'Record'] = None, edges: list[Edge] = None):
        """Represents a single bacterial genome to be loaded into memory from a file"""
        self.id: str = id_
        self.contigs: dict[str: 'Record'] = contigs or {}
        self.edges: list[Edge] = edges or []
        self.is_annotated: bool = False

    def __len__(self):
        return sum(len(i) for i in self.contigs.values())

    def __iter__(self):
        return iter(self.contigs.values())

    def __str__(self):
        return self.id

    def __getitem__(self, item: str) -> 'Record':
        return self.contigs[item]

    def __format__(self, __format_spec: Literal['fasta', 'fna', 'ffn', 'faa', 'bed', 'gfa'] = ''):
        if __format_spec == '':
            return self.__str__()
        elif __format_spec in {'fasta', 'fna', 'ffn', 'faa', 'bed'}:
            return ''.join(format(i, __format_spec) for i in self.contigs.values())
        elif __format_spec == 'gfa':
            return ''.join(format(i, __format_spec) for i in chain(self.contigs.values(), self.edges))
        else:
            raise NotImplementedError(f'Invalid format: {__format_spec}')

    @classmethod
    def from_file(cls, file: Union[str, Path, IO, SeqFile], annotations: Union[str, Path, SeqFile] = None):
        """
        Loads a genome from a file with optional annotations, e.g. a FASTA with a BED/GFF file.
        :param file: Path to file (str/Path) or a SeqFile instance.
        :param annotations: Optional annotations as a file path (str/Path) or SeqFile instance.
        :return: A Genome instance representing the genome in memory
        """
        file = SeqFile(file) if not isinstance(file, SeqFile) else file
        self = cls(file.id)
        for record in file:
            if isinstance(record, Record):
                self.contigs[record.id] = record
                if not file.format == 'gfa' and next((v for k, v in record.qualifiers if k == 'topology'),
                                                     'linear') == 'circular':
                    self.edges.append(Edge(record.id, record.id, 1, 1))  # Add self loop for circular genomes
                    # We assume a GFA file already has this edge but this may not be the case
            elif isinstance(record, Edge):
                self.edges.append(record)
        if annotations:
            if file.format not in {'fasta', 'gfa'}:
                raise GenomeError(f'Can only provide annotations to FASTA and GFA files, not {file.format}')
            annotations = SeqFile(annotations) if not isinstance(annotations, SeqFile) else annotations
            if not annotations.format in {'gff', 'bed'}:
                raise GenomeError(f'Annotations must be in GFF or BED format, not {annotations.format}')
            for contig_id, features in grouper(annotations, 'location.parent_id'):  # Sort by contig id
                if contig := self.contigs.get(contig_id):  # type: Record
                    contig.add_features(*features)
                    self.is_annotated = True
        return self

    @classmethod
    def random(
            cls, id_: str = None, genome: Record = None, rng: Random = RESOURCES.rng, n_contigs: int = None,
            min_contigs: int = 1, max_contigs: int = 1000, gc: float = 0.5, length: int = None, min_len: int = 10,
            max_len: int = 5000000
    ):
        """
        Generates a random genome assembly for testing purposes.

        :param id_: ID of the genome, if not provided a random one will be generated
        :param genome: Initial genome record; if not provided a random one will be generated
        :param rng: Random number generator.
        :param n_contigs: Number of contigs in the assembly. If not provided, a random number of contigs will be generated.
        :param min_contigs: Minimum number of contigs if n_contigs is not specified.
        :param max_contigs: Maximum number of contigs if n_contigs is not specified.
        :param gc: GC content of the genome.
        :param length: Total length of the genome. If not provided, a random length will be generated.
        :param min_len: Minimum length of the genome if length is not specified.
        :param max_len: Maximum length of the genome if length is not specified.
        :return: A Genome instance representing the random assembly.
        """
        n_contigs = n_contigs or rng.randint(min_contigs, max_contigs)
        genome = genome or Record.random(None, DNA, rng, gc, length, min_len, max_len)
        return cls(id_ or str(uuid4()), {i.id: i for i in genome.shred(rng, n_contigs)})

    def as_assembly_graph(self, directed: bool = False) -> Graph:
        """
        Returns the assembly as a graph where contigs are nodes
        """
        return Graph(*self.edges, directed=directed)

    def as_feature_graph(self, directed: bool = False) -> Graph:
        """
        Returns the assembly as a graph where features are nodes.
        Adjacent features on the same contig are connected. Features at the termini
        of connected contigs are also connected.
        """
        feature_graph = Graph(directed=directed)

        # 1. Add intra-contig edges (connecting adjacent features on the same contig)
        for contig in self.contigs.values():
            if contig.features:
                # This assumes contig.as_edge_list() correctly connects adjacent features
                for edge in contig.as_edge_list():
                    feature_graph.add_edge(edge)

        # 2. Add inter-contig edges by iterating through all assembly connections
        for edge in self.as_assembly_graph(directed=directed).edges:
            # Get the original 'from' and 'to' contigs for this edge
            fr_contig_orig, to_contig_orig = self.contigs.get(edge.fr), self.contigs.get(edge.to)

            # Skip if either contig is missing or has no features
            if not (fr_contig_orig and fr_contig_orig.features) or not (to_contig_orig and to_contig_orig.features):
                continue

            # Determine the correct orientation based on the assembly edge attributes
            fr_contig = fr_contig_orig.reverse_complement() if edge.fr_attribute == -1 else fr_contig_orig
            to_contig = to_contig_orig.reverse_complement() if edge.to_attribute == -1 else to_contig_orig

            # Identify the terminal features on the correctly oriented contigs
            fr_feature = fr_contig.features[-1]  # The last feature on the 'from' contig
            to_feature = to_contig.features[0]  # The first feature on the 'to' contig

            # Calculate the gap distance between the features
            distance = (len(fr_contig) - fr_feature.location.end) + to_feature.location.start

            feature_graph.add_edge(  # Add the new edge connecting the two features
                Edge(
                    fr_feature.id,
                    to_feature.id,
                    fr_feature.location.strand,
                    to_feature.location.strand,
                    distance
                )
            )

        return feature_graph

    @require('pyrodigal')
    def find_genes(self, gene_finder: 'pyrodigal.GeneFinder' = None, pool: Executor = None,
                   config: GeneFinderConfig = None) -> Generator[Feature, None, None]:
        """
        Predicts ORFs in the Genome using Pyrodigal.

        Parameters:
            gene_finder: An optional pre-trained `pyrodigal.GeneFinder` instance. If not provided, one will be
                         initialized using the `config` parameter.
            pool: An optional `Executor` (e.g., `ThreadPoolExecutor` or `ProcessPoolExecutor`)
                  to parallelize gene finding across contigs. If `None`, a `ThreadPoolExecutor` will be created.
            config: A `GeneFinderConfig` object to configure the Pyrodigal GeneFinder. Only used if `gene_finder`
                    is `None`.

        Yields:
            `Feature` objects representing the predicted CDS (Coding Sequence) genes.

        Notes:
            - Requires `pyrodigal` to be installed.
        """
        from pyrodigal import __version__ as _pyrodigal_version
        from pyrodigal import Gene

        if gene_finder is None:
            from pyrodigal import GeneFinder
            gene_finder = GeneFinder(**asdict(config or GeneFinderConfig()))

        if pool is None:
            pool = ThreadPoolExecutor(max_workers=min(32, RESOURCES.available_cpus + 4))

        if gene_finder.training_info is None:
            gene_finder.train(*(bytes(contig.seq) for contig in self.contigs.values()))

        n = 0  # Gene counter
        for contig, genes in zip(self.contigs.values(), pool.map(lambda contig: gene_finder.find_genes(bytes(contig.seq)), self.contigs.values())):
            features = []
            for gene in genes:  # type: Gene
                features.append(cds := Feature(
                    id_=f'{contig.id}_{n + 1:05d}', kind='CDS', seq=Seq(gene.sequence(), 'DNA'),
                    location=Location(gene.begin - 1, gene.end, gene.strand, gene.partial_begin, gene.partial_end,
                                      contig.id),
                    qualifiers=[
                        Qualifier('translation', Seq(gene.translate(), 'Amino')),
                        Qualifier('transl_table', gene.translation_table),
                        Qualifier('inference', f'ab initio prediction:Pyrodigal:{_pyrodigal_version}'),
                        Qualifier('GC', gene.gc_cont),
                        Qualifier('rbs_motif', gene.rbs_motif),
                        Qualifier('rbs_spacer', gene.rbs_spacer)
                    ]
                ))
                n += 1  # Increment the counter
                yield cds  # Yield the CDS for potential use

            contig.add_features(*features)  # Use the add_features method to sort with existing features

        if n > 0:  # If genes were found, mark the genome as annotated
            self.is_annotated = True


# Functions ------------------------------------------------------------------------------------------------------------
def _guess_format_from_handle(handle: TextIO) -> _SUPPORTED_FORMATS:
    """
    Guesses the file format by peeking at the first line of a text handle.
    Assumes the handle is seekable.
    """
    try:
        first_line = handle.readline().strip()
        handle.seek(0)  # Rewind the handle for the actual parser

        if not first_line:
            raise SeqFileError("Cannot guess format from an empty file.")

        if first_line.startswith('>'):
            return 'fasta'
        elif first_line.startswith('@'):
            return 'fastq'
        elif first_line.startswith('LOCUS'):
            return 'genbank'
        elif first_line.startswith('S\t'):
            return 'gfa'
        elif first_line.startswith('##gff-version 3'):
            return 'gff'
        else:
            raise SeqFileError(f"Could not guess file format from first line: '{first_line[:100]}...'")

    except Exception as e:
        raise SeqFileError("Failed to read from stream to guess format.") from e


def parse(handle: TextIO, format_: _SUPPORTED_FORMATS = 'guess') -> Generator[Union[Record, Edge], None, None]:
    """
    Simple parser for fasta, gfa and genbank formats,
    similar to `Biopython <https://biopython.org/docs/latest/api/Bio.SeqIO.html#Bio.SeqIO.parse>`_

    :param handle: A file handle opened in text-mode / text stream
    :param format_: The format of the file; must be one of the supported formats or guess by default
    :returns: A Generator of Record objects
    """
    if format_ == 'guess':
        if m := _SEQUENCE_FILE_REGEX.search(handle.name):
            format_ = next(fmt for fmt in get_args(_SUPPORTED_FORMATS) if m[fmt])
        else:
            raise SeqFileError(f'Unsupported SeqFile format or extension: {handle.name}')
    if parser := {'fasta': _parse_fasta, 'gfa': _parse_gfa, 'genbank': _parse_genbank, 'fastq': _parse_fastq,
                  'gff': _parse_gff}.get(format_):
        yield from parser(handle)
    else:
        raise NotImplementedError(f'Format "{format_}" not supported')


def _parse_fasta(handle: TextIO) -> Generator[Record, None, None]:
    """
    Simple FASTA parser

    :param handle: A file handle opened in text-mode / text stream
    :returns: A Generator of Record objects
    """
    header, seq = '', []
    for line in handle:  # Loop over lines in chunk
        if not (line := line.strip()):
            continue
        if line.startswith('>'):
            if header and seq:
                name, desc = header.split(' ', 1) if ' ' in header else (header, '')
                yield Record(name, Seq(''.join(seq)), desc, qualifiers=[
                    Qualifier('topology', 'circular' if _TOPOLOGY_REGEX.search(desc) else 'linear')])
            header, seq = line[1:], []
        else:
            seq.append(line)
    if header and seq:
        name, desc = header.split(' ', 1) if ' ' in header else (header, '')
        _TOPOLOGY_REGEX.search(desc)
        yield Record(name, ''.join(seq), desc, qualifiers=[
                    Qualifier('topology', 'circular' if _TOPOLOGY_REGEX.search(desc) else 'linear')])
    else:
        warn('No records parsed', ParserWarning)


def _parse_fastq(handle: TextIO) -> Generator[Record, None, None]:
    """
    Simple FASTQ parser

    :param handle: A file handle opened in text-mode / text stream
    :returns: A Generator of Record objects
    """
    while True:
        if not (header := handle.readline().strip()):  # 1. Read the header line
            break  # End of file
        if not header.startswith('@'):
            raise ValueError(f"Expected FASTQ record header starting with '@', but got: {header}")
        if not (seq := handle.readline().strip()):  # 2. Read the sequence line
            raise ValueError(f"Unexpected end of file after reading header: {header}")
        if not (sep := handle.readline().strip()):  # 3. Read the separator line
            raise ValueError(f"Unexpected end of file after reading sequence for header: {header}")
        if not sep.startswith('+'):
            raise ValueError(f"Expected FASTQ separator line starting with '+', but got: {sep}")
        if not (qual := handle.readline().strip()):   # 4. Read the quality line
            raise ValueError(f"Unexpected end of file after reading separator for header: {header}")
        name, desc = header.split(' ', 1) if ' ' in header else (header, '')
        yield Record(name[1:], Seq(seq, 'DNA'), desc, qualifiers=[Qualifier('quality', qual)])


def _parse_gfa(handle: TextIO) -> Generator[Union[Record, Edge], None, None]:
    """
    Simple GFA parser

    :param handle: A file handle opened in text-mode / text stream
    :returns: A Generator of Record objects
    """
    records = False
    for line in handle:  # Iterate over file lines
        if line.startswith('S\t'):  # Segment contains contig info
            parts = line[2:].strip().split('\t', 2)  # Split into 3 parts: name, sequence and description
            if len(parts[1]) >= 1:  # Check sequence is at least 1bp
                yield Record(parts[0], parts[1], qualifiers=list(_parse_tags(parts[2])) if len(parts) > 2 else [])
                records = True
        elif line.startswith('L\t'):  # Add links once all contigs are added
            fr, fr_strand, to, to_strand = line[2:].strip().split('\t')[:4]
            yield Edge(fr, to, 1 if fr_strand == '+' else -1, 1 if to_strand == '+' else -1)
    if not records:
        warn('No records parsed', ParserWarning)


def _parse_bed(handle: TextIO) -> Generator[Feature, None, None]:
    """
    Simple BED parser
    
    :param handle: A file handle opened in text-mode / text stream
    :returns: A Generator of Feature objects
    """
    raise NotImplementedError
    # TODO: Implement this
    # See: https://samtools.github.io/hts-specs/BEDv1.pdf
    # bed_format = 0
    # while True:
    #     if not (line := handle.readline().strip()):  # 1. Read the header line
    #         break  # End of file
    #     if n_columns := len(line := line.split('\t')) < 3:
    #         raise ValueError(f"Expected at least 3 columns but got {n_columns}")
    #     if not bed_format:
    #         bed_format = n_columns  # Use header to determine the format of the file, e.g. BED12 or BED6
    #     elif n_columns != bed_format:
    #         raise ValueError(f"Expected {bed_format} BED columns, but got {n_columns}")
    #
    #     location = Location(int(line[1]), int(line[2]), -1 if (bed_format > 3 and line[5] == '-') else 1, ref=line[0])


def _parse_gff(handle: TextIO, feature_kinds: set[str] = frozenset({'CDS'})) -> Generator[Feature, None, None]:
    """
    Simple GFF3 parser, see https://ensembl.org/info/website/upload/gff3.html

    :param handle: A file handle opened in text-mode / text stream
    :param feature_kinds: Set of feature kinds to parse; note the 'source' feature will populate the record's qualifiers
    :returns: A Generator of Feature objects
    """
    feature = None
    for line in handle:
        if line.startswith('#'):
            continue
        if len(parts := line.strip().split('\t')) >= 9 and parts[2] in feature_kinds:
            feature = Feature(
                Location(int(parts[3]) - 1, int(parts[4]), -1 if parts[6] == '-' else 1, parent_id=parts[0]),
                kind=parts[2],
                qualifiers=[Qualifier(*i.split('=', 1)) for i in parts[8].split(';')] + [
                    Qualifier('source', parts[1])]
            )
            feature.id = feature['ID'] or 'unknown'  # Update the feature.id using the ID qualifier
            yield feature
    if feature is None:
        warn('No features parsed', ParserWarning)

def _parse_genbank(handle: TextIO, feature_kinds: set[str] = frozenset({'CDS'})) -> Generator[Record, None, None]:
    """
    Simple genbank parser

    :param handle: A file handle opened in text-mode / text stream
    :param feature_kinds: Set of feature kinds to parse; note the 'source' feature will populate the record's qualifiers
    :returns: A Generator of Record objects
    """
    record = []
    for line in handle:  # Loop over lines in chunk
        if line := line.strip():
            if line.startswith('LOCUS'):  # This is the beginning of the new record
                if len(record) > 1:  # Records must all consist of more than 1 line
                    yield _parse_genbank_record(record, feature_kinds)  # Yield the previous record
                    record = []
            if not line.startswith('//'):  # Add lines until the end of the record, signified by "//"
                record.append(line)
    if len(record) > 1:  # Records must all consist of more than 1 line
        yield _parse_genbank_record(record, feature_kinds)
    else:
        warn('No records parsed', ParserWarning)


def _parse_genbank_record(record: list[str], feature_kinds: set[str] = frozenset({'CDS'})) -> Record:
    """Parser for a single record in Genbank format"""
    # attributes, features = '\n'.join(record).split('FEATURES', 1)
    # features, origin = features.split('\nORIGIN\n', 1)
    features, origin = '\n'.join(record).split('FEATURES', 1)[1].split('\nORIGIN\n', 1)
    record = Record(
        id_=record[0].split()[1], desc=record[1].split(maxsplit=1)[1],
        seq=parse_genbank_origin(origin.strip().split('\n')),
        qualifiers=[Qualifier('topology', 'circular' if _TOPOLOGY_REGEX.search(record[0]) else 'linear')]
    )
    current_feature = []
    for line in features.split('\n')[1:]:
        if not line.startswith('/') and len(line.split()) == 2 and '..' in line:  # New feature
            if current_feature:
                if feature := _parse_genbank_feature(current_feature, feature_kinds, record.id):
                    if feature.kind == 'source':
                        record.qualifiers.extend(feature.qualifiers)
                    else:
                        record.features.append(feature)
                current_feature = []
        current_feature.append(line)
    if current_feature:
        if feature := _parse_genbank_feature(current_feature, feature_kinds, record.id):
            if feature.kind == 'source':
                record.qualifiers.extend(feature.qualifiers)
            else:
                record.features.append(feature)
    return record


def parse_genbank_origin(origin: list[str]) -> Seq:
    return Seq(''.join(chain.from_iterable(i.split()[1:] for i in origin if i)), alphabet='DNA')


def _parse_genbank_feature(feature_lines: list[str], feature_kinds: set[str] = frozenset({'CDS'}), parent_id: str = None
                           ) -> Union[Feature, None]:
    feature_kind, location = feature_lines[0].split(maxsplit=1)
    if feature_kind not in feature_kinds and feature_kind != 'source':
        return None
    feature = Feature(kind=feature_kind, location=_parse_genbank_location(location, parent_id))
    if len(feature_lines) > 1:
        for qualifier in _parse_genbank_qualifiers(feature_lines[1:]):
            if qualifier.key == 'locus_tag' and feature.id == 'unknown':  # Use locus tag as ID
                feature.id = qualifier.value  # Don't add as qualifier
            elif qualifier.key == 'translation':  # Turn the translation into a Seq object
                # Add stop codon if not already there
                # qualifier.value = Seq(qualifier.value + ('' if qualifier.value.endswith('*') else '*'), alphabet='Amino')
                qualifier.value = Seq(qualifier.value, alphabet='Amino')
                feature.qualifiers.append(qualifier)
            else:
                feature.qualifiers.append(qualifier)
    return feature


def _parse_genbank_location(location: str, parent_id: str = None) -> Location:
    locations = []
    strand: Literal[1, -1] = -1 if 'complement' in location else 1
    for match in _GENBANK_LOCATION_REGEX.finditer(location):
        locations.append(
            new_location := Location(int(match.group('start')) - 1, int(match.group('end')), strand, parent_id=parent_id))
        if match.group('partial_start'):
            new_location.partial_start = True
        if match.group('partial_end'):
            new_location.partial_end = True
    if not locations:
        raise ValueError(f'Could not parse location: {location}')

    location = locations.pop(0)  # type: Location
    location.joins.extend(locations)
    return location


def _parse_genbank_qualifiers(lines: list[str]) -> Generator[Qualifier, None, None]:
    """Parse the attribute lines of a genbank record"""
    current_qualifier = []
    for line in lines:
        if line.startswith('/'):
            if current_qualifier:
                yield _format_qualifier(current_qualifier)
            current_qualifier = []
        current_qualifier.append(line)
    if current_qualifier:
        yield _format_qualifier(current_qualifier)


def _format_qualifier(qualifier: list[str]) -> Qualifier:
    # So far, I think only multiline translations need to be joined without whitespace
    # TODO: We could add a step to replace multiple spaces with a single space
    qualifier = ('' if qualifier[0].startswith('/translation') else ' ').join(qualifier).lstrip('/')
    if '=' not in qualifier:
        return Qualifier(qualifier, True)  # I prefer having "True" as the value rather than None
    key, value = qualifier.split('=', 1)
    if value.startswith('"'):
        value = value.strip('"')
    elif value.startswith('('):
        value = value.lstrip('(').rstrip(')')
    elif '.' in value:
        value = float(value)
    else:
        value = int(value)
    return Qualifier(key, value)


def _parse_tags(line: str, col_delim: str = '\t', tag_delimiter: str = ':'
                ) -> Generator[Qualifier, None, None]:
    """Parse tag column and yield a tuple of the tag and value in the correct type"""
    for item in line.split(col_delim):
        if item.count(tag_delimiter) == 2:
            tag, typ, val = item.split(tag_delimiter)  # type: str, str, str
            yield Qualifier(tag, _TAG2TYPE.get(typ, str)(val))


def group_genomes(genomes: Iterable[Union[Path, str, IO, SeqFile]]) -> Generator[Genome, None, None]:
    """
    Helper function, mostly for CLI, for grouping together sequence and annotation files for the same genome.
    Useful for adding annotations in BED/GFF format to a genome in FASTA/GFA format.

    Note:
        Whilst this function can accept input streams, everything is coerced into SeqFiles instances and will be
        written to disk, so this should be taken into consideration before processing many streams.

    Arguments:
        genomes: An iterable of files as strings or Path objects, SeqFiles or IO streams.

    Yields:
        Genome instances

    Example:
        >>> from pathlib import Path
        >>> for genome in group_genomes(Path('genomes').iterdir()):
        >>>     print(f'{genome:fasta}', end='')

    """
    for id_, files in grouper(((SeqFile(i) if not isinstance(i, SeqFile) else i) for i in genomes), 'id'):
        if len(files := list(files)) > 2:
            raise GenomeError(f'More than 2 files found for {id_}: {files}')
        seq_formats, annotation_formats = [], []
        for file in files:
            if file.format in {'fasta', 'gfa', 'genbank'}:
                seq_formats.append(file)
            elif file.format in {'bed', 'gff'}:
                annotation_formats.append(file)
                
        if (n_seqfiles := len(seq_formats)) != 1:
            raise GenomeError(f'{n_seqfiles} sequence files found for {id_}: {seq_formats}')
        
        if len(seq_formats) == 1 and len(annotation_formats) == 0:
            yield Genome.from_file(seq_formats[0])
        
        if len(seq_formats) == 1 and len(annotation_formats) == 1:
            yield Genome.from_file(seq_formats[0], annotations=annotation_formats[0])


def group_reads(reads: Iterable[Union[Path, str, IO, ReadFile]]) -> Generator[ReadSet, None, None]:
    """
    Helper function, mostly for CLI, for grouping together read files for the same sample.
    Useful for grouping paired-end or hybrid (long + short) readsets.

    Note:
        Whilst this function can accept input streams, everything is coerced into ReadFile instances and will be
        written to disk, so this should be taken into consideration before processing many streams.

    Arguments:
        reads: An iterable of files as strings or Path objects, ReadFile or IO streams.

    Yields:
        ReadSet instances

    Example:
        >>> from pathlib import Path
        >>> for readset in group_reads(Path('reads').iterdir()):
        >>>     print(f'{readset:fasta}', end='')

    """
    for sample_name, files in grouper(((ReadFile(i) if not isinstance(i, ReadFile) else i) for i in reads), 'sample_name'):
        yield ReadSet(*files, id_=sample_name)
