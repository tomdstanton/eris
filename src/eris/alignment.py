"""
Module for parsing, viewing and managing sequence alignments.
"""
from typing import Iterable, Generator, Literal, Union
from operator import attrgetter
from re import compile as regex

from eris.seq import HasLocation, Location, Seq, Feature, Qualifier

# Constants ------------------------------------------------------------------------------------------------------------
_CIGAR_OPERATIONS = regex(r'(?P<n>[0-9]+)(?P<operation>[MIDNSHP=X])')
_QUERY_CONSUMING_OPERATIONS = {"M", "I", "S", "=", "X"}
_TARGET_CONSUMING_OPERATIONS = {"M", "D", "N", "=", "X"}

# Classes --------------------------------------------------------------------------------------------------------------
class AlignmentError(Exception):
    pass


class Alignment(HasLocation):
    """
    Class representing an alignment between a query sequence and target sequence.

    Attributes:
        query: Query sequence name
        query_length: Query sequence length
        query_start: Query start coordinate (0-based)
        query_end: Query end coordinate (0-based)
        target: Target sequence name
        target_length: Target sequence length
        length: Number of residues in the alignment including gaps
        cigar: CIGAR string
        score: Alignment score
        E: Alignment entropy
        n_matches: Number of matching residues in the alignment
        quality: Mapping quality (0-255 with 255 for missing)
        tags: {tag: value} pairs
        identity: Percentage of matching residues
        query_coverage: Percentage of query sequence covered by the alignment
        target_coverage: Percentage of target sequence covered by the alignment
        aligned_seqs: Tuple of aligned sequences
    """

    def __init__(self, query: str, query_start: int, query_end: int, target: str, location: Location,
                 query_length: int = 0, target_length: int = 0, length: int = 0, cigar: str = None,
                 n_matches: int = 0, quality: int = 0, tags: dict[str, Union[str, int, float]] = None, score: float = 0,
                 E: float = 0, identity: float = 0,  query_coverage: float = 0, target_coverage: float = 0,
                 aligned_seqs: tuple[Seq, Seq] = None):
        super().__init__(location)
        self.query: str = query  # Query sequence name
        self.query_length: int = query_length  # Query sequence length
        self.query_start: int = query_start  # Query start coordinate (0-based)
        self.query_end: int = query_end  # Query end coordinate (0-based)
        self.target: str = target  # Target sequence name
        self.target_length: int = target_length  # Target sequence length
        self.length: int = length  # Number of residues in the alignment including gaps
        self.cigar: str = cigar
        self.score: float = score
        self.E: float = E
        self.n_matches: int = n_matches  # Number of matching residues in the alignment
        self.quality: int = quality  # Mapping quality (0-255 with 255 for missing)
        self.tags: dict[str, Union[str, int, float]] = tags or {}  # {tag: value} pairs)
        self.identity: float = identity
        self.query_coverage: float = query_coverage
        self.target_coverage: float = target_coverage
        self.aligned_seqs: tuple[Seq, Seq] = aligned_seqs

    def __repr__(self):
        return f'{self.query}:{self.query_start}-{self.query_end} {self.target}:{self.location}'

    def __len__(self):
        return self.length

    @classmethod
    def from_paf(cls, paf_line: str):
        """
        Parse a line in PAF format and return an Alignment object.

        Parameters:
            paf_line: A text string representing a single alignment
        """
        if len(paf_line := paf_line.strip().split('\t')) < 12:
            raise AlignmentError(f"PAF Line has < 12 columns: {paf_line}")
        try:
            self = Alignment(  # Parse standard fields
                query=paf_line[0], query_length=int(paf_line[1]), query_start=int(paf_line[2]),
                query_end=int(paf_line[3]), target=paf_line[5],
                location=Location(int(paf_line[7]), int(paf_line[8]), 1 if paf_line[4] == '+' else -1),
                target_length=int(paf_line[6]), n_matches=int(paf_line[9]), length=int(paf_line[10]),
                quality=int(paf_line[11]), tags={
                    (x := t.split(":", 2))[0]:
                        int(x[2]) if x[1] == "i" else float(x[2]) if x[1] == "f" else x[2] for
                    t in paf_line[12:]
                }
            )
            self.identity = (self.n_matches / self.length) * 100
            self.query_coverage = (self.length / self.query_length) * 100  # TODO: Exclude gaps from these counts
            self.target_coverage = (self.length / self.target_length) * 100  # TODO: Exclude gaps from these counts
            if 'cg' in self.tags:  # Add cigar to attributes
                self.cigar = self.tags.pop('cg')
            if 'AS' in self.tags:  # Add alignment score to attributes
                self.score = self.tags.pop('AS')
            self.location.parent_id = self.target

            is_partial_at_query_start = self.query_start > 0
            is_partial_at_query_end = self.query_end < self.query_length

            if self.location.strand == 1:
                self.location.partial_start = is_partial_at_query_start
                self.location.partial_end = is_partial_at_query_end
            else:
                self.location.partial_start = is_partial_at_query_end
                self.location.partial_end = is_partial_at_query_start

            return self
        except Exception as e:
            raise AlignmentError(f"Error parsing PAF line: {paf_line}: {e}")

    def get_aligned_seqs(self, query: Seq = None, target: Seq = None, gap_character: str = '-') -> tuple[Seq, Seq]:
        """
        Calculates the aligned sequences from the query and target sequences using the cigar.
        Will also set the ``aligned_seqs`` attribute.

        :param query: The query sequence
        :param target: The target sequence
        :param gap_character: The character to use to represent gaps
        :return: A tuple of aligned sequences
        """
        if not self.aligned_seqs:
            assert query, AlignmentError('Need a query sequence to compute traceback')
            assert target, AlignmentError('Need a target sequence to compute traceback')
            assert len(query) == self.query_length, AlignmentError('Query length does not match alignment length')
            assert len(target) == self.target_length, AlignmentError('Target length does not match alignment length')
            query_trimmed = str(query[self.query_start:self.query_end])
            target_trimmed = str(target[self.location.start:self.location.end])
            query_traceback, target_traceback = [], []
            query_position, target_position = 0, 0
            for operation, n in parse_cigar(self.cigar):
                if operation in _QUERY_CONSUMING_OPERATIONS:
                    query_traceback.append(query_trimmed[query_position:query_position + n])
                    query_position += n
                else:
                    query_traceback.append(gap_character * n)
                if operation in _TARGET_CONSUMING_OPERATIONS:
                    target_traceback.append(target_trimmed[target_position:target_position + n])
                    target_position += n
                else:
                    target_traceback.append(gap_character * n)

            self.aligned_seqs = (Seq(''.join(query_traceback), query.alphabet), Seq(''.join(target_traceback), target.alphabet))
        return self.aligned_seqs

    def print(self, wrap: int = 70):
        assert self.aligned_seqs, AlignmentError(f'Aligned sequences have not been computed for {self}')
        query, target = self.aligned_seqs
        matches = ''.join('|' if x == y else '.' for x, y in zip(query, target))
        rjust = len(str(max(self.query_start, self.location.start)))
        for i in range(0, self.length, wrap):
            print(f"target {self.location.start + i:>{rjust}} {target[i:i + wrap]}\n"
                  f"       {i:>{rjust}} {matches[i:i + wrap]}\n"
                  f"query  {self.query_start + i:>{rjust}} {query[i:i + wrap]}")

    def as_feature(self, kind: str = 'misc_feature') -> Feature:
        return Feature(
            location=self.location, kind=kind, qualifiers=[
                Qualifier('name', self.query),
                Qualifier('query_start', self.query_start),
                Qualifier('query_end', self.query_end),
                Qualifier('query_length', self.query_length),
                Qualifier('query_coverage', self.query_coverage),
                Qualifier('target_length', self.target_length),
                Qualifier('target_coverage', self.target_coverage),
                Qualifier('identity', self.identity),
                Qualifier('length', self.length),
                Qualifier('cigar', self.cigar),
                Qualifier('n_matches', self.n_matches),
                Qualifier('quality', self.quality),
                Qualifier('score', self.score),
                Qualifier('E', self.E)
            ]
        )


# Functions ------------------------------------------------------------------------------------------------------------
def cull(keep: Alignment, alignments: Iterable[Alignment], max_overlap_fraction: float = 0.1
         ) -> Generator[Alignment, None, None]:
    """
    Filters alignments by excluding those that overlap a reference alignment beyond a specified fraction.

    This function takes a reference alignment and iterates through a collection of
    alignments, yielding only those alignments whose target does not match the
    target of the reference alignment, or whose overlap ratio with the reference
    alignment (calculated as the overlap divided by the length of the alignment)
    is less than a provided maximum overlap fraction.

    Args:
        keep (Alignment): The reference alignment used for comparison.
        alignments (Iterable[Alignment]): A collection of alignments to be filtered.
        max_overlap_fraction (float): The maximum allowed fraction of overlap
            between the reference alignment and any given alignment. Defaults to 0.1.

    Yields:
        Generator[Alignment, None, None]: A generator that yields alignments
            meeting the specified criteria.
    """
    for a in alignments:
        if a.target != keep.target or (a.overlap(keep) / a.length) < max_overlap_fraction:
            yield a


def cull_all(alignments: Iterable[Alignment], key='n_matches', reverse_sort: bool = True) -> list[Alignment]:
    """
    Sort and filter a collection of alignments based on a specified key. The function sorts
    alignments by the given key and iteratively compares the alignments, keeping those that
    meet certain criteria and removes others.

    Arguments:
        alignments (Iterable[Alignment]): A collection of Alignment objects to be processed.
        key (str): The attribute name used for sorting alignments. Defaults to 'n_matches'.
        reverse_sort (bool): If True, sort alignments in descending order using the specified
            key. Defaults to True.

    Returns:
        list[Alignment]: A filtered and sorted list of alignments, meeting specified criteria.
    """
    kept_alignments = []
    sorted_alignments = sorted(alignments, key=attrgetter(key), reverse=reverse_sort)
    while sorted_alignments:
        kept_alignments.append(sorted_alignments.pop(0))
        sorted_alignments = list(cull(kept_alignments[-1], sorted_alignments))
    return kept_alignments


def parse_cigar(cigar: str) -> Generator[tuple[Literal['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'], int], None, None]:
    """
    Parses the cigar with a regular expression
    :return: A Generator of tuples containing the operation symbol and span
    """
    if cigar:
        for match in _CIGAR_OPERATIONS.finditer(cigar):
            if not match:
                raise AlignmentError(f'Could not parse cigar: {cigar}')
            yield match['operation'], int(match['n'])

