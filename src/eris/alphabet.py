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
# Classes --------------------------------------------------------------------------------------------------------------
class AlphabetError(Exception):
    pass


class Alphabet:
    def __init__(self, symbols: str):
        self._symbols = symbols
        self._set: frozenset[str] = frozenset(symbols)
        if len(self._set) < len(self._symbols):
            raise AlphabetError(f'Alphabet symbols "{symbols}" are not unique')

    def __len__(self):
        return len(self._symbols)

    def __repr__(self):
        return self._symbols

    def __hash__(self):
        return hash(self._set)

    def __eq__(self, other):
        if isinstance(other, Alphabet):
            return self._symbols == other._symbols
        return False

    def __contains__(self, item: str):
        return set(item) <= self._set

    def __iter__(self):
        return iter(self._symbols)

    def __getitem__(self, item):
        return self._symbols[item]


class _Amino(Alphabet):
    def __init__(self):
        super().__init__('ACDEFGHIKLMNPQRSTVWY')


class _DNA(Alphabet):
    def __init__(self):
        super().__init__('ACGT')
        self.complement = str.maketrans({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}) # Complement
        self.codons = dict(
            zip([a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'],  # Translation table 11 (prokaryotes)
                'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))


# Constants ------------------------------------------------------------------------------------------------------------
DNA = _DNA()
Amino = _Amino()
