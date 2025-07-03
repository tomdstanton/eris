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
from operator import attrgetter
from typing import Literal, Iterator, Union, Any, Generator, Iterable
from copy import copy
from random import Random
from uuid import uuid4

from eris import Requires, RESOURCES
from eris.alphabet import Alphabet, DNA, _DNA

# Constants ------------------------------------------------------------------------------------------------------------
_TYPE2TAG = {float: 'f', int: 'i', str: 'Z'}  # See: https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#optional-fields

# Classes --------------------------------------------------------------------------------------------------------------
class LocationError(Exception):
    pass


class Location:
    """
    Represents a region on a sequence.

    :param start: Start coordinate (0-based)
    :param end: End coordinate (0-based)
    :param strand: Strand (1 or -1)
    :param partial_start: Boolean indicating if the start is partial
    :param partial_end: Boolean indicating if the end is partial
    :param parent_id: Reference sequence name
    :param joins: List of other locations that are joined to this location
    """
    def __init__(self, start: int, end: int, strand: Literal[1, -1] = 1, partial_start: bool = False,
                 partial_end: bool = False, parent_id: str = None, *joins: 'Location'):
        self.start: int = start
        self.end: int = end
        self.strand: Literal[1, -1] = strand
        self.partial_start: bool = partial_start
        self.partial_end: bool = partial_end
        self.parent_id: str = parent_id
        self.joins: list[Location] = list(joins)
        # TODO: Implement more methods inspired by
        #  https://github.com/sanger-pathogens/Fastaq/blob/master/src/pyfastaq/intervals.py

    def __repr__(self):
        return f"{self.start}-{self.end} {self.strand}"

    def __len__(self):
        return self.end - self.start

    def __iter__(self):
        return iter((self.start, self.end, self.strand))

    def __contains__(self, item: Union['Location', 'HasLocation', int, float, slice]):
        if isinstance(item, slice):
            return self.start <= item.start and self.end >= item.stop
        elif isinstance(item, (Location, HasLocation)):
            if isinstance(item, HasLocation):
                item = item.location
            return self.start <= item.start and self.end >= item.end
        elif isinstance(item, (int, float)):
            return self.start <= item <= self.end
        else:
            raise TypeError(item)

    def __add__(self, other: Union['Location', 'HasLocation']) -> 'Location':
        if isinstance(other, HasLocation):
            other = other.location
        if isinstance(other, Location):
            return Location(min(self.start, other.start), max(self.end, other.end), self.strand)
        else:
            raise TypeError(other)

    def __radd__(self, other: Union['Location', 'HasLocation']) -> 'Location':
        if isinstance(other, HasLocation):
            other = other.location
        if isinstance(other, Location):
            return other.__add__(self)
        else:
            raise TypeError(other)

    def __iadd__(self, other: Union['Location', 'HasLocation']):
        if isinstance(other, HasLocation):
            other = other.location
        if isinstance(other, Location):
            self.start = min(self.start, other.start)
            self.end = max(self.end, other.end)
            return self
        else:
            raise TypeError(other)

    def partial(self) -> bool:
        """
        Returns True if the location is partial, False otherwise
        """
        return self.partial_start or self.partial_end

    def overlap(self, other: Union['Location', 'HasLocation']) -> int:
        """
        Returns the length of the overlap between two locations
        :param other: Location or HasLocation object
        """
        if isinstance(other, HasLocation):
            other = other.location
        if isinstance(other, Location):
            return max(0, min(self.end, other.end) - max(self.start, other.start))
        else:
            raise TypeError(other)

    def extract(self, parent: Union[Union['Seq', 'Record', 'Feature'], dict[str, Union['Seq', 'Record', 'Feature']]]) -> 'Seq':
        """
        Extracts a sequence from a parent object based on the location.

        :param parent: Parent object (``Seq``, ``Record``, ``Feature``) or dictionary of parent objects
        :return: Seq object
        """
        new_seq = []
        if isinstance(parent, dict):
            if not (ref := parent.get(self.parent_id)):
                raise LocationError(f'{self.parent_id=} not in parent dictionary')
            new_seq.append(ref.seq[self] if not isinstance(ref, Seq) else ref[self])
        else:
            if self.parent_id and not isinstance(parent, Seq) and parent.id != self.parent_id:
                raise LocationError(f'{self.parent_id=} does not match {parent.id=}')
            new_seq.append(parent.seq[self] if not isinstance(parent, Seq) else parent[self])
        if self.joins:
            new_seq += (i.extract(parent) for i in self.joins)
        if len(alphabet := {i.alphabet for i in new_seq}) > 1:
            raise LocationError('Cannot extract sequence from mixed alphabet types')
        return Seq(''.join(map(str, new_seq)), alphabet.pop())

    def shift(self, by: int):
        return Location(self.start + by, self.end + by, self.strand)

    @classmethod
    def random(cls, rng: Random = RESOURCES.rng, length: int = None, min_len: int = 1, max_len: int = 10000,
               min_start: int = 0, max_start: int = 1000000):
        if not length:
            length = rng.randint(min_len, max_len)
        start = rng.randint(min_start, max_start - length)
        return cls(start, start + length, rng.choice([1, -1]))

    def _get_points(self, max_head_length: float, shape: Literal['arrow', 'box'] = 'arrow',
                    arrow_shaft_width: float = 0.4, head_width_ratio: float = 1.5, head_length_ratio: float = 0.2):
        """
        Generates a set of points to define either an arrow or box shape along a strand.
        """
        if shape == 'box':
            return [(self.start, 0), (self.end, 0), (self.end, 1), (self.start, 1)]
        elif shape == 'arrow':
            head_length = min(head_length_ratio * len(self), max_head_length)
            head_width = head_width_ratio * arrow_shaft_width  # Head width relative to shaft
            shaft_end = (self.end - head_length) if self.strand == 1 else (self.start + head_length)
            return [
                (self.start if self.strand == 1 else self.end, -arrow_shaft_width / 2),  # Start of shaft
                (shaft_end, -arrow_shaft_width / 2),  # End of shaft, base of head
                (shaft_end, -head_width / 2),  # Head tip, left
                (self.end if self.strand == 1 else self.start, 0),  # Head tip, center
                (shaft_end, head_width / 2),  # Head tip, right
                (shaft_end, arrow_shaft_width / 2),  # End of shaft, base of head (top)
                (self.start if self.strand == 1 else self.end, arrow_shaft_width / 2),  # Start of shaft (top)
            ]
        else:
            raise ValueError(shape)


class HasLocation:
    """
    Base class for objects with a location attribute
    """
    def __init__(self, location: Location):
        self.location = location

    def __len__(self) -> int:
        return self.location.end - self.location.start

    def overlap(self, other: Union[Location, 'HasLocation']) -> int:
        """
        Returns the length of the overlap between two locations
        """
        return self.location.overlap(other)

    def __contains__(self, item: Union[Location, 'HasLocation', int, float]) -> bool:
        """
        Returns True if the object contains the item, False otherwise
        """
        return self.location.__contains__(item)


class SeqError(Exception):
    pass


class TranslationError(SeqError):
    pass


class Seq:
    """
    Class representing a sequence.

    :param seq: Sequence string
    :param alphabet: alphabet (DNA or Amino), attempts to guess if not provided
    """
    def __init__(self, seq: str, alphabet: Literal['DNA', 'Amino'] = None):
        if len(seq := seq.strip().upper()) < 1:
            raise SeqError("Seq must be >=1 character(s)")
        self._seq = seq
        self.alphabet = alphabet or ('DNA' if self._seq in DNA else 'Amino')

    def __repr__(self):
        return f"{self.alphabet}({self._seq})" if len(self._seq) < 13 else f"{self.alphabet}({self._seq[:5]}...{self._seq[-5:]})"

    def __len__(self):
        return len(self._seq)

    def __hash__(self) -> int:
        return hash(self._seq)

    def __eq__(self, other):
        if isinstance(other, Seq):
            return self._seq == other._seq
        return False

    def __str__(self):
        return self._seq

    def __bytes__(self):
        return self._seq.encode()

    def __iter__(self):
        return iter(self._seq)

    def __contains__(self, item: str):
        return item in self._seq

    def __reversed__(self) -> 'Seq':
        return Seq(self._seq[::-1], self.alphabet)

    def __add__(self, other: Union[str, 'Seq']) -> 'Seq':
        if isinstance(other, Seq):
            assert self.alphabet == other.alphabet, SeqError('Both sequences need to be of the same alphabet')
            return Seq(self._seq + str(other), self.alphabet)
        elif isinstance(other, str):
            return Seq(self._seq + other, self.alphabet)
        else:
            raise TypeError(other)

    def __radd__(self, other: Union[str, 'Seq']) -> 'Seq':
        return self.__add__(other)

    def __iadd__(self, other: Union[str, 'Seq']):
        if isinstance(other, Seq):
            assert self.alphabet == other.alphabet, SeqError('Both sequences need to be of the same alphabet')
            self._seq += str(other)
            return self
        elif isinstance(other, str):
            self._seq += other
            return self
        else:
            raise TypeError(other)

    def __getitem__(self, item: Union[slice, int, Location, HasLocation]) -> 'Seq':
        """Quick method of slicing the sequence and reverse complementing if necessary"""
        if isinstance(item, (slice, int)):
            return Seq(self._seq[item], self.alphabet)
        elif isinstance(item, (Location, HasLocation)):
            if isinstance(item, HasLocation):
                item = item.location
            new = Seq(self._seq[item.start:item.end], self.alphabet)
            return new if (item.strand == 1) else new.reverse_complement()
        else:
            raise TypeError(item)

    def __call__(self) -> str:
        return self._seq

    @classmethod
    def random(
            cls, alphabet: Alphabet, rng: Random = RESOURCES.rng, gc: float = 0.5, length: int = None,
            min_len: int = 10, max_len: int = 5000000) -> 'Seq':
        if isinstance(alphabet, _DNA):
            at = (1 - gc) / 2
            gc /= 2
            return cls(''.join(rng.choices(alphabet, weights=[at, gc, at, gc], k=length or rng.randint(min_len, max_len))), 'DNA')
        else:
            return cls(''.join(rng.choice(alphabet) for _ in range(rng.randint(min_len, max_len))), 'Amino')

    def replace(self, *args, **kwargs):
        return Seq(self._seq.replace(*args, **kwargs), self.alphabet)

    def reverse_complement(self) -> 'Seq':
        """Reverse complements the sequence if it is a nucleotide sequence"""
        if self.alphabet == 'DNA':
            return Seq(self._seq.translate(DNA.complement)[::-1], 'DNA')
        return self

    def translate(self, to_stop: bool = True, stop_symbol: str = "*", frame: Literal[0, 1, 2] = 0,
                  gap_character: str = None) -> 'Seq':
        """
        Translates the sequence to amino acid if it is a nucleotide sequence

        :param to_stop: Boolean to stop translation at the first stop codon
        :param stop_symbol: Character to resemble the stop codon
        :param frame: Zero-based frame to begin translation from, must be one of 0, 1 or 2
        :param gap_character: The gap character if performing a gapped-translation
        """
        if self.alphabet == 'DNA':  # Only need to translate DNA
            if len(self._seq) < 3:
                raise TranslationError(f'Cannot translate sequence of length {len(self._seq)}')
            protein = []
            for i in range(frame, len(self), 3):  # Iterate over seq codons (chunks of 3)
                if len(codon := self._seq[i:i + 3]) < 3:
                    break  # We can't translate chunks less than 3 (not codons) so break here
                if gap_character and gap_character in codon:
                    protein.append(gap_character)
                    continue
                protein.append(residue := DNA.codons.get(codon, stop_symbol))  # If lookup fails (get == None), assume stop codon
                if to_stop and residue == stop_symbol:
                    break  # Break if to_stop == True
            if protein == ['*']:
                raise TranslationError('Translation consists of a single stop codon')
            return Seq(''.join(protein), 'Amino')
        return self  # alphabet == 'Amino' so return self

    def GC(self) -> float:
        """
        Returns the GC content of the sequence
        """
        return sum(1 for i in self._seq if i in {'G', 'C'}) / len(self._seq)


class Qualifier:
    """
    Class to represent a qualifier as a key and value

    Attributes:
        key: str
        value: Any
    """
    def __init__(self, key: str, value: Any = None):
        self.key = key
        self.value = value

    def __repr__(self):
        return f"{self.key}={self.value}"

    def __iter__(self) -> Iterator[tuple[str, Any]]:
        return iter((self.key, self.value))

    def __str__(self):
        return f"{self.key}:{_TYPE2TAG[type(self.value)]}:{self.value}"

    def __eq__(self, other):
        if isinstance(other, Qualifier):
            return self.key == other.key and self.value == other.value
        return False


class Record:
    """
    Represents a sequence record.

    Attributes:
        id: str
        seq: The sequence of the record.
        desc: A description of the record.
        qualifiers: A list of Qualifier objects providing additional information about the record.
        features: A list of Feature objects representing features on the sequence.
    """

    def __init__(self, id_: str, seq: Union[Seq, str], desc: str = None, qualifiers: list[Qualifier] = None,
                 features: list['Feature'] = None):
        self.id = id_
        self.seq = seq if isinstance(seq, Seq) else Seq(seq)
        self.desc = desc or ''
        self.qualifiers = qualifiers or []
        self.features = features or []  # TODO: Consider using a heapq for positional lookup

    def __getitem__(self, item: Union[slice, HasLocation, Location]) -> 'Record':
        if isinstance(item, HasLocation):  # Convert HasLocation to Location
            item = item.location
        elif isinstance(item, slice):  # Convert slice to location
            item = Location(item.start or 0, item.stop or len(self), 1)
        elif not isinstance(item, Location):
            raise TypeError(f"Item must be a slice,a Location or object with a location attribute, not {item}")
        new_record = Record(f"{self.id}_{item.start}-{item.end}", self.seq[item], self.desc)
        for feature in self.features:  # Assume features are sorted
            if feature.location in item:
                new_record.features.append(new_feature := feature.shift(-item.start))
                new_feature.location.parent_id = new_record.id
        return new_record

    def __repr__(self) -> str:
        return f"{self.id} {self.seq.__repr__()}"

    def __str__(self) -> str:
        return self.id

    def __len__(self) -> int:
        return len(self.seq)

    def __hash__(self) -> int:
        return hash(self.id)

    def __eq__(self, other) -> bool:
        if isinstance(other, Record):
            return self.id == other.id
        return False

    def __iter__(self) -> Iterator['Feature']:
        return iter(self.features)

    def __format__(self, __format_spec: Literal['fasta', 'fna', 'ffn', 'faa', 'bed', 'gfa'] = '') -> str:
        if __format_spec == '':
            return self.__str__()
        elif __format_spec in {'fasta', 'fna'}:
            return f">{self.id}\n{self.seq}\n"
        elif __format_spec in {'faa', 'ffn'}:
            return ''.join(format(i, __format_spec) for i in self.features)
        elif __format_spec == 'bed':
            return ''.join(f"{self.id}\t{format(i, __format_spec)}" for i in self.features)
        elif __format_spec == 'gfa':
            # See: https://gfa-spec.github.io/GFA-spec/GFA1.html
            return f"S\t{self.id}\t{self.desc}\t{self.seq}\n"
        else:
            raise NotImplementedError(f'Invalid format: {__format_spec}')

    def __add__(self, other: 'Record') -> 'Record':
        if isinstance(other, Record):
            new = Record(f"{self.id}_{other.id}", self.seq + other.seq, self.desc, [], [])
            for feature in self.features:
                new_feat = copy(feature)  # Creates a copy for the new record
                new_feat.location.parent_id = new.id  # Update the ref ID

            for feature in other.features:  # Other feature locations need to be updated
                feature = feature.shift(len(self))  # Creates a new feature and location
                feature.location.parent_id = new.id  # Update the ref ID
                new.features.append(feature)

            new.features.sort(key=attrgetter('location.start'))  # Sort features
            return new
        else:
            raise TypeError(other)

    def __radd__(self, other: 'Record') -> 'Record':
        return other.__add__(self)

    def __iadd__(self, other: 'Record'):
        if isinstance(other, Record):
            self.id = f"{self.id}_{other.id}"  # Update the ref ID
            for feature in self.features:
                feature.location.parent_id = self.id  # Update the ref ID

            for feature in other.features:
                feature = feature.shift(len(self))  # Creates a new feature and location
                feature.location.parent_id = self.id
                self.features.append(feature)  # Update the ref ID

            self.features.sort(key=attrgetter('location.start'))  # Sort features
            self.seq += other.seq  # Now we can update the sequence
            return self
        else:
            raise TypeError(other)

    def __delslice__(self, i, j):
        del self.seq[i:j]  # Remove the part of the sequence
        self.features = []  # Reset the features
        for feature in self.features:
            if feature.location.end <= i:
                continue  # Feature is before the deleted part
            elif feature.location.start >= j:
                self.features.append(feature.shift(-(j - i)))  # Adjust location
            else:
                del feature  # Feature overlaps the deleted part, remove it

    @classmethod
    def random(cls, id_: str = None, alphabet: Alphabet = DNA, rng: Random = RESOURCES.rng, gc: float = 0.5,
               length: int = None, min_len: int = 10, max_len: int = 5000000) -> 'Record':
        """
        Generates a random record for testing purposes.

        :param id_: ID of the record. If not provided, a random UUID will be used.
        :param alphabet: Alphabet of the sequence.
        :param rng: Random number generator.
        :param gc: GC content of the sequence.
        :param length: Length of the sequence. If not provided, a random length will be generated.
        :param min_len: Minimum length of the sequence if length is not specified.
        :param max_len: Maximum length of the sequence if length is not specified.
        :return: A Record instance.
        """
        return cls(id_ or str(uuid4()), Seq.random(alphabet, rng, gc, length, min_len, max_len))

    def shred(self, rng: Random = RESOURCES.rng, n_breaks: int = None, break_points: list[int] = None
              ) -> Generator['Record', None, None]:
        """
        Shreds the record into smaller records at the specified break points.

        :param rng: Random number generator
        :param n_breaks: The number of breaks to make in the record. If not provided, a random number of breaks will be
            made between 1 and half the length of the record.
        :param break_points: A list of break points to use. If not provided, random break points will be generated.
        :return: A generator of smaller records
        """
        if not n_breaks:
            n_breaks = rng.randint(1, len(self) // 2)
        if not break_points:
            break_points = sorted([rng.randint(0, len(self)) for _ in range(n_breaks)])
        previous_end = 0
        for break_point in break_points:
            yield self[previous_end:break_point]
            previous_end = break_point
        yield self[previous_end:]

    def insert(self, other: 'Record', at: int, replace: bool = True) -> 'Record':
        """
        Inserts another record into this record at the specified position.

        :param other: The record to insert.
        :param at: The position to insert the other record at.
        :param replace: Whether to replace the existing sequence at the insertion point with the inserted sequence.
            If False, the inserted sequence will be inserted without removing any existing sequence.
        :return: A new Record instance with the other record inserted.
        """
        if not 0 < at < len(self):
            raise IndexError(f'Cannot insert at {at}, must be between 0 and {len(self)}')
        else:
            return self[:at] + other + self[at if not replace else at + len(other):]

    def translate(self, to_stop: bool = True, stop_symbol: str = "*", frame: Literal[0, 1, 2] = 0,
                  gap_character: str = None) -> 'Record':
        """
        Translates the record to amino acid if it is a nucleotide sequence

        :param to_stop: Boolean to stop translation at the first stop codon
        :param stop_symbol: Character to resemble the stop codon
        :param frame: Zero-based frame to begin translation from, must be one of 0, 1 or 2
        :param gap_character: The gap character if performing a gapped-translation
        """
        if self.seq.alphabet == 'DNA':
            return self
        else:
            return Record(self.id, self.seq.translate(to_stop, stop_symbol, frame, gap_character), self.desc)

    @Requires(requires='matplotlib')
    def plot(
            self, ax: 'plt.Axes' = None, arrow_shaft_width: float = 0.4, head_width_ratio: float = 1.5,
            head_length_ratio: float = 0.2, face_color: str = "lightsalmon",
            line_width: int = 3, edge_color='orangered', label_color: str = 'black', label_offset: float = 0.2,
            record_color: str = "black", show_record: bool = True, label_features: bool = True, label_size: float = 8
    ) -> tuple['Figure', 'plt.Axes']:
        """
        Plots the record and its features.

        :param ax: Matplotlib axes object
        :param arrow_shaft_width: Width of the arrow shaft
        :param head_width_ratio: Ratio of head width to shaft width
        :param head_length_ratio: Ratio of head length to feature length
        :param face_color: Face color of the features
        :param line_width: Line width of the feature edges
        :param edge_color: Edge color of the features
        :param label_color: Color of the feature labels
        :param label_offset: Offset of the feature labels from the features
        :param record_color: Color of the record line
        :param show_record: Boolean to show the record line
        :param label_features: Boolean to label the features
        :param label_size: Size of the feature labels
        :return: Matplotlib figure and axes objects
        """
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection

        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 2))
        else:
            fig = ax.figure  # Get existing figure
        max_head_length = 0.1 * len(self)  # Cap head length at 10% of record length
        if show_record:
            ax.plot([0, len(self)], [0, 0], color=record_color, lw=2, solid_capstyle='butt')

        patches = []
        for f in self.features:
            points = f.location._get_points(max_head_length, 'arrow' if f.kind == 'CDS' else 'box',
                                            arrow_shaft_width, head_width_ratio, head_length_ratio)
            patches.append(
                Polygon(points, closed=True, fill=True, facecolor=face_color, linewidth=line_width, edgecolor=edge_color)
            )
            if label_features:
                label_x = (f.location.start + f.location.end) / 2
                label_y = label_offset  # Position above the feature
                ax.text(label_x, label_y, f.id, ha='center', va='bottom', rotation=45, fontsize=label_size,
                        color=label_color)

        collection = PatchCollection(patches, match_original=True)
        ax.add_collection(collection)
        ax.set_xlim(0, len(self))
        ax.set_ylim(-1, 1.5)  # Adjust ylim to accommodate labels
        ax.set_title(f"{self.id} - {self.desc}")
        ax.set_xlabel("Position")
        ax.set_yticks([])
        return fig, ax


class Feature(HasLocation):
    """
    Represents a feature on a sequence.

    Attributes:
        id: The unique identifier for the feature.
        kind: The type of feature (e.g., gene, CDS, rRNA).
        seq: The sequence of the feature (optional).
        qualifiers: A list of Qualifier objects providing additional information about the feature.
    """

    def __init__(self, location: Location, id_: str = 'unknown', kind: str = 'misc_feature',
                 seq: Union[Seq, str] = None, qualifiers: list[Qualifier] = None):
        super().__init__(location)
        self.id = id_
        self.kind = kind
        self.seq = (seq if isinstance(seq, Seq) else Seq(seq)) if seq else None
        self.qualifiers = qualifiers or []

    def __repr__(self):
        return f"{self.kind}({self.id} {self.location})"

    def __str__(self):
        return self.id

    def __hash__(self) -> int:
        return hash(self.id)

    def __eq__(self, other):
        if isinstance(other, Feature):
            return self.id == other.id
        return False
    
    def __getitem__(self, item: str) -> Union[Any, None]:
        """
        Quick method of getting a qualifier by key using next()
        Warning:
            This will only return the first instance, which in most cases is fine
        """
        if not isinstance(item, str):
            raise TypeError(item)
        return next((v for k, v in self.qualifiers if k == item), None)

    @classmethod
    def random(cls, parent: Record, id_: str = None, rng: Random = RESOURCES.rng, min_len: int = 1,
               max_len: int = 10000, min_start: int = 0, max_start: int = 1000000):
        location = Location.random(rng, min_len=min_len, max_len=min(max_len, len(parent)), min_start=min_start,
                                   max_start=min(max_start, len(parent) - max_len))
        location.parent_id = parent.id
        feature = cls(location, id_ or str(uuid4()), 'CDS')
        feature.extract(parent, store_seq=True)
        return feature

    def translate(
            self, to_stop: bool = True, stop_symbol: str = "*", frame: Literal[0, 1, 2] = 0, gap_character: str = None,
            store_translation: bool = True,
            parent: Union[Union['Seq', 'Record', 'Feature'], dict[str, Union['Seq', 'Record', 'Feature']]] = None,
            store_seq: bool = False
    ) -> Seq:
        """
        Translates the feature to amino acid if it is a nucleotide sequence
        :param to_stop: Boolean to stop translation at the first stop codon
        :param stop_symbol: Character to resemble the stop codon
        :param frame: Zero-based frame to begin translation from
        :param gap_character: The gap character if performing a gapped-translation
        :param store_translation: Boolean to store the translation in the feature
        :param parent: Parent object (``Seq``, ``Record``, ``Feature``) or dictionary of parent objects
        :param store_seq: Boolean to store the translated sequence in the feature
        :return: Translated sequence
        """
        if not (translation := next((v for k, v in self.qualifiers if k == 'translation'), None)):
            if self.seq is None:
                if parent is None:
                    raise AttributeError(f'Cannot translate feature "{self.id}" without a stored sequence or parent to '
                                         f'extract from')
                else:
                    translation = self.extract(parent, store_seq).translate(to_stop, stop_symbol, frame, gap_character)
            else:
                translation = self.seq.translate(to_stop, stop_symbol, frame, gap_character)
            if store_translation:
                self.qualifiers.append(Qualifier('translation', translation))
        return translation

    def __format__(self, __format_spec: Literal['fasta', 'ffn', 'fna', 'faa', 'bed'] = ''):
        if __format_spec == '':
            return self.__str__()
        elif __format_spec in {'fasta', 'ffn', 'fna'}:
            if self.seq is None:
                raise AttributeError(f'No seq extracted for {self=}')
            else:
                return f">{self.id}\n{self.seq}\n"
        elif __format_spec == 'faa':
            return f">{self.id}\n{self.translate()}\n"
        elif __format_spec == 'bed':
            return (f'{self.location.start}\t{self.location.end}\t{self.id}\t'
                    f'{next((v for k, v in self.qualifiers if k == "score"), 0)}\t'
                    f'{"+" if self.location.strand == 1 else "-"}\n')
        else:
            raise NotImplementedError(f'Format "{__format_spec}" not supported')

    def extract(self, parent: Union[Union['Seq', 'Record', 'Feature'], dict[str, Union['Seq', 'Record', 'Feature']]],
                store_seq: bool = False) -> 'Seq':
        seq = self.location.extract(parent)
        if store_seq:
            self.seq = seq
        return seq

    def shift(self, by: int) -> 'Feature':
        return Feature(self.location.shift(by), self.id, self.kind, None, self.qualifiers)


class Cluster:
    """
    Represents a cluster of genes as ``Feature`` objects

    """
    def __init__(self, id_: str = 'unknown', description: str = 'unknown', *members: Feature):
        self.id = id_
        self.description = description
        self._members = {i for i in members}

    def __len__(self):
        return len(self._members)

    def __repr__(self):
        return f"Cluster({self.id} - {self.description} - {len(self._members)} members)"

    def __str__(self):
        return self.id

    def __iter__(self):
        return iter(self._members)

    def __hash__(self) -> int:
        return hash(self.id)

    def __eq__(self, other):
        if isinstance(other, Cluster):
            return self._members == other._members
        return False

    def add(self, other: Any):
        self._members.add(other)

    def __format__(self, __format_spec: Literal['tsv'] = ''):
        if __format_spec == '':
            return self.__str__()
        elif __format_spec == 'tsv':
            return ''.join(f"{self.id}\t{i}\n" for i in self._members)
        else:
            raise NotImplementedError(f'Format "{__format_spec}" not supported')


# Functions ------------------------------------------------------------------------------------------------------------
def merge_locations(
        locations: Iterable[Union[Location, HasLocation]], tolerance: Union[int, float, None] = None
) -> Generator[Location, None, None]:
    """
    Merges a list of locations into a single location, possibly with a tolerance

    :param locations: Iterable of locations or objects wth a location to merge
    :param tolerance: Tolerance for merging locations or None to merge all locations
    :return: A generator of merged locations

    """
    if not locations:
        return None

    locations = sorted((i.location if isinstance(i, HasLocation) else i for i in locations), key=attrgetter('start'))
    if len(locations) == 1:
        yield from locations
        return None

    if tolerance is None:  # Merge locations regardless
        starts, ends, strands = [], [], []
        for i in locations:
            starts.append(i.start); ends.append(i.end); strands.append(i.strand)
        yield Location(min(starts), max(ends), max(set(strands), key=strands.count))

    else:  # Merge locations within tolerance
        current_location = locations[0]  # Start with first location
        strands = [current_location.strand]  # Collect strands to calculate consensus
        # TODO: Think about how to deal with locations with joins
        # TODO: Also think about how to deal with refs
        for start, end, strand in locations[1:]:  # Iterate through the locations and unpack
            if start - tolerance <= current_location.end:  # Overlap, merge the locations
                current_location = Location(
                    current_location.start, max(current_location.end, end), current_location.strand
                )
                strands.append(strand)
            else:  # No overlap, add the current location to the merged list and start a new location
                current_location.strand = max(set(strands), key=strands.count)
                yield current_location  # Yield the current location
                current_location = Location(start, end, strand)  # Start a new location
                strands = [current_location.strand]

        current_location.strand = max(set(strands), key=strands.count)
        yield current_location  # Yield the last location
