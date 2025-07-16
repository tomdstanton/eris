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
from typing import Iterable, Generator, Literal, Iterator, Any
from itertools import groupby
from operator import attrgetter
from re import compile

from eris.seq import HasLocation, Location, Seq, Feature, Qualifier

# Constants ------------------------------------------------------------------------------------------------------------
_CIGAR_OPERATIONS = compile(r'(?P<n>[0-9]+)(?P<operation>[MIDNSHP=X])')
_QUERY_CONSUMING_OPERATIONS = {"M", "I", "S", "=", "X"}
_TARGET_CONSUMING_OPERATIONS = {"M", "D", "N", "=", "X"}

# Classes --------------------------------------------------------------------------------------------------------------
class AlignmentError(Exception):
    pass


class Alignment(HasLocation):
    """
    Class representing an alignment between a query sequence and target sequence.

    As this class contains a location attribute, it can be used for slicing and extraction of ``Seq``, ``Record`` and
    ``Feature`` objects.

    :param query: Query sequence name
    :param query_start: Query start coordinate (0-based)
    :param query_end: Query end coordinate (0-based)
    :param target: Target sequence name
    :param location: ``Location`` of the alignment on the target sequence
    :param query_length: Query sequence length
    :param target_length: Target sequence length
    :param length: Number bases, including gaps, in the alignment
    :param cigar: CIGAR string
    :param aligned_seqs: Tuple of aligned sequences
    """

    def __init__(self, query: str, query_start: int, query_end: int, target: str, location: Location,
                 query_length: int = 0, target_length: int = 0, length: int = 0, cigar: str = None,
                 n_matches: int = 0, quality: int = 0, tags: dict = None, score: float = 0, E: float = 0,
                 identity: float = 0,  query_coverage: float = 0, target_coverage: float = 0,
                 aligned_seqs: tuple[Seq, Seq] = None):
        super().__init__(location)
        self.query = query  # Query sequence name
        self.query_length = query_length  # Query sequence length
        self.query_start = query_start  # Query start coordinate (0-based)
        self.query_end = query_end  # Query end coordinate (0-based)
        self.target = target  # Target sequence name
        self.target_length = target_length  # Target sequence length
        self.length = length  # Number of residues in the alignment including gaps
        self.cigar = cigar
        self.score = score
        self.E = E
        self.n_matches = n_matches  # Number of matching residues in the alignment
        self.quality = quality  # Mapping quality (0-255 with 255 for missing)
        self.tags = tags or {}  # {tag: value} pairs):
        self.identity = identity
        self.query_coverage = query_coverage
        self.target_coverage = target_coverage
        self.aligned_seqs = aligned_seqs

    def __repr__(self):
        return f'{self.query}:{self.query_start}-{self.query_end} {self.target}:{self.location}'

    def __len__(self):
        return self.length

    @classmethod
    def from_paf(cls, paf_line: str):
        """
        Parse a line in PAF format and return an Alignment object.
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
            if self.location.end + (self.query_length - self.query_end) >= self.target_length:  # + strand only, test -
                self.location.partial_end = True
            if self.location.start < self.query_start:  # + strand only, test -
                self.location.partial_start = True
            return self
        except Exception as e:
            raise AlignmentError(f"Error parsing PAF line: {paf_line}: {e}")

    @classmethod
    def from_sam(cls, sam_line: str):
        """
        Parse a line in SAM format and return an Alignment object.
        """
        raise NotImplementedError
        # TODO: Implement this

    def iter_cigar(self) -> Generator[tuple[Literal['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'], int], None, None]:
        """
        Parses the cigar with a regular expression
        :return: A Generator of tuples containing the operation symbol and span
        """
        if not self.cigar:
            raise AlignmentError(f'No CIGAR computed for {self}')
        for match in _CIGAR_OPERATIONS.finditer(self.cigar):
            if not match:
                raise AlignmentError(f'Could not parse cigar: {self.cigar}')
            yield match['operation'], int(match['n'])

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
            for operation, n in self.iter_cigar():
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
            location=self.location, id_=f"{self.query}|{self.target}|{self.location}", kind=kind,
            qualifiers=[
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
def group_alignments(alignments: Iterable[Alignment], key: str = 'query'
                     ) -> Generator[tuple[Any, Iterator[Alignment]], None, None]:
    """
    Groups alignments by a specified attribute.

    This function takes an iterable of Alignment objects and groups them based on
    a specified attribute (defaulting to 'query'). The alignments are first sorted
    by the specified attribute, and then grouped using itertools.groupby, yielding
    tuples containing the group key and an iterator over the grouped alignments.

    Parameters:
        alignments: Iterable[Alignment]
            A collection of Alignment objects to group.
        key: str, default 'query'
            The attribute of the Alignment objects to use for grouping.

    Yields:
        tuple[Any, Iterator[Alignment]]
            A tuple consisting of the group key and an iterator over the grouped
            Alignment objects.
    """
    yield from groupby(sorted(alignments, key=attrgetter(key)), key=attrgetter(key))


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

