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
from functools import partial
from operator import attrgetter
from warnings import warn
from re import compile
from pathlib import Path
from typing import Union, Generator, IO, Literal, Iterable, TextIO, get_args, Callable
from concurrent.futures import Executor
from itertools import chain, groupby
from io import BufferedIOBase, RawIOBase, BytesIO
from shutil import copyfileobj
from tempfile import NamedTemporaryFile
from random import Random
from uuid import uuid4

from eris import ErisWarning, Requires, RESOURCES
from eris.alphabet import DNA
from eris.seq import Record, Feature, Seq, Qualifier, Location
from eris.graph import Graph, Edge
from eris.utils import xopen

# Constants ------------------------------------------------------------------------------------------------------------
_SUPPORTED_FORMATS = Literal['fasta', 'gfa', 'genbank', 'fastq', 'bed']
_TAG2TYPE = {'f': float, 'i': int, 'Z' : str}  # See: https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#optional-fields
_FEATURE_KINDS = {'CDS'}  # 'gene', 'misc_feature', 'ncRNA', 'rRNA', 'regulatory', 'tRNA', 'tmRNA'}
_SEQUENCE_FILE_REGEX = compile(
    r'\.((?P<fasta>f(asta|a|na|fn|as|aa))|(?P<fastq>f(ast)?q)|(?P<gfa>gfa)|(?P<genbank>g(b|bff|bk|enbank))|(?P<bed>bed))'
    r'\.?(?P<compression>(gz|bz2|xz|zst))?$'
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
    :param format_: Format of the file (required if input is a stream).
    :return: A SeqFile instance
    """
    def __init__(self, file: Union[str, Path, IO], format_: _SUPPORTED_FORMATS = None):
        self.id: str = "unknown"
        self.path: Union[Path, None] = None
        self._handle: Union[IO, None] = None
        self._from_stream: bool = False
        self._format: Union[str, None] = format_
        self._open_func: Union[None, Callable] = None

        if hasattr(file, 'read') and not isinstance(file, (str, Path)):  # --- Handle IO Stream Input ---
            if format_ is None:
                raise SeqFileError('Input format must be specified when providing an IO stream.')

            if not isinstance(file, (BufferedIOBase, RawIOBase)):
                # Try to get underlying binary buffer (common for TextIOBase like sys.stdin)
                if buffer := getattr(file, 'buffer', None):
                    file = buffer  # Use the buffer instead
                else:  # It's some other kind of stream we can't easily get bytes from
                    raise SeqFileError("Input stream is not binary and lacks an accessible binary buffer (.buffer).")

            with NamedTemporaryFile(prefix='seqfile_', suffix=f'.{self._format}', delete=False, mode='wb') as temp_f_handle:
                copyfileobj(file, temp_f_handle)  # Dump raw bytes from input stream to temp file
                self.path = Path(temp_f_handle.name)  # Store path to the temp file (which contains raw/compressed data)
                self.id = self.path.stem
                self._from_stream = True
                self._open_func = partial(xopen, self.path, method='magic', mode='rt')

        elif isinstance(file, (str, Path)):  # --- Handle File Path Input ---
            self.path = file if isinstance(file, Path) else Path(file)
            if m := _SEQUENCE_FILE_REGEX.search(self.path.name):
                self.id = self.path.name.rstrip(m.group())
                self._format = next(fmt for fmt in get_args(_SUPPORTED_FORMATS) if m[fmt])
                self._open_func = partial(xopen, self.path, method=m['compression'] or 'uncompressed', mode='rt')
            else:
                raise SeqFileError(f'Unsupported SeqFile format or extension: {self.path.suffixes}')
        else:
            raise TypeError(f"Input must be a file path (str or Path) or an IO stream, not {type(file)}")

    @property
    def format(self):
        return self._format

    def __repr__(self):
        return f'SeqFile({self.id}, format={self._format})'

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

    def __iter__(self) -> Generator[Union['Record', 'Edge'], None, None]:
        """Returns an iterator (generator) over records in the file."""
        self._ensure_handle_open() # Ensure handle is open/reopened
        if not self._format:
             raise SeqFileError("Cannot parse file: format is unknown.")
        yield from parse(self._handle, self._format)
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
    def __init__(self, *files: Union[str, Path, ReadFile]):
        self.id = None
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
    """A class representing a single genome assembly with contigs and edges"""
    def __init__(self, id_: str, contigs: dict[str: 'Record'] = None, edges: list[Edge] = None):
        """Represents a single bacterial genome to be loaded into memory from a file"""
        self.id: str = id_
        self.contigs: dict[str: 'Record'] = contigs or {}
        self.edges: list[Edge] = edges or []
        self._is_annotated: bool = False

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
    def from_file(cls, file: Union[str, Path, SeqFile], annotations: Union[str, Path, SeqFile] = None):
        file = SeqFile(file) if not isinstance(file, SeqFile) else file
        self = cls()
        for record in file:
            if isinstance(record, Record):
                self.contigs[record.id] = record
                if not file.format == 'gfa' and next((v for k, v in record.qualifiers if k == 'topology'),
                                                     'linear') == 'circular':
                    self.edges.append(Edge(record.id, record.id, 1, 1))  # Add self loop for circular genomes
                    # We assume a GFA file already has this edge but this may not be the case
            elif isinstance(record, Edge):
                self.edges.append(record)
        return self

    @classmethod
    def random(
            cls, genome: Record = None, rng: Random = RESOURCES.rng, n_contigs: int = None, min_contigs: int = 1,
            max_contigs: int = 1000, gc: float = 0.5, length: int = None, min_len: int = 10, max_len: int = 5000000
    ):
        """
        Generates a random genome assembly for testing purposes.

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
        return cls(BytesIO(''.join(format(i, 'fasta') for i in genome.shred(rng, n_contigs)).encode()), 'fasta')

    def as_assembly_graph(self, directed: bool = False) -> Graph:
        """
        Returns the assembly as a graph where contigs are nodes
        """
        graph = Graph(*self.edges, directed=directed)
        for i in self.contigs:
            graph.add_node(i)
        return graph

    # def as_feature_graph(self, directed: bool = False) -> Graph:
    #     """
    #     Returns the assembly as a graph where features are nodes
    #     Adjacent features are connected to one another by an edge, if the Genome has edges connecting contigs,
    #     features that are on the termini of connected contigs will also be connected to each other.
    #     """
    #     edges = []
    #     for contig in self.contigs.values():
    #         for n, to in enumerate(contig.features[1:], start=1):
    #             edges.append(
    #                 Edge((fr := contig.features[n - 1]).name, to.name, fr.location.strand,
    #                      to.location.strand, to.location.start - fr.location.end)
    #             )
    #     for edge in self.edges:
    #         if (fr_contig := self.contigs.get(edge.fr)) and (to_contig := self.contigs.get(edge.to)):
    #             # TODO: Calculate distance between features across different contigs for weight
    #             # TODO: Deal with multiple locations
    #             edges.append(
    #                 Edge((fr := fr_contig.features[-1 if edge.fr_attribute == 1 else 0]).name,
    #                      (to := to_contig.features[0 if edge.to_attribute == 1 else -1]).name,
    #                      fr.location.strand, to.location.strand)
    #             )
    #     return Graph(*edges, directed=directed)

    @Requires(requires='pyrodigal')
    def find_genes(self, gene_finder: 'GeneFinder', pool: Executor) -> Generator[Feature, None, None]:
        from pyrodigal import __version__ as _pyrodigal_version
        if not gene_finder.training_info:
            gene_finder.train(*(bytes(contig.seq) for contig in self.contigs.values()))
        n = 0
        for contig, genes in zip(self.contigs.values(), pool.map(lambda contig: gene_finder.find_genes(bytes(contig.seq)), self.contigs.values())):
            for gene in genes:
                contig.features.append(cds := Feature(
                    id_=f'{contig.name}_{n + 1:05d}', kind='CDS', seq=Seq(gene.sequence(), 'DNA'),
                    location=Location(gene.begin - 1, gene.end, gene.strand, gene.partial_begin, gene.partial_end),
                    qualifiers=[
                        Qualifier('translation', Seq(gene.translate(), 'Amino')),
                        Qualifier('transl_table', gene.translation_table),
                        Qualifier('inference', f'ab initio prediction:Pyrodigal:{_pyrodigal_version}')
                    ]
                ))
                n += 1
                yield cds

        self._is_annotated = True

# Functions ------------------------------------------------------------------------------------------------------------
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
                  'bed': _parse_bed}.get(format_):
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
                yield Record(parts[0], parts[1], qualifiers=list(_parse_tags(parts[2])))
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


def group_reads(reads: Iterable[Union[Path, str, ReadFile]]) -> Generator[ReadSet, None, None]:
    """
    Groups reads by sample name

    E.g. to load in ReadSets from a directory of reads, use: `group_reads(Path('reads').iterdir())`
    """
    for _, readset in groupby(sorted(((ReadFile(i) if not isinstance(i, ReadFile) else i) for i in reads),
                                     key=attrgetter('sample_name')), key=attrgetter('sample_name')):
        yield ReadSet(*readset)
