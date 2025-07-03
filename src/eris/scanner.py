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
from pathlib import Path
from typing import Literal, Union, Generator
from concurrent.futures import ThreadPoolExecutor, Executor
from dataclasses import asdict, dataclass
from collections import defaultdict
from operator import attrgetter
from warnings import warn
from itertools import chain
from functools import partial

from eris import Requires, RESOURCES, ErisWarning
from eris.seq import Qualifier, Feature, Record
from eris.io import Genome, SeqFile
from eris.graph import Edge
from eris.utils import group_positions, Config


# Constants ------------------------------------------------------------------------------------------------------------
_DATA_PATH = Path('eris/data')

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


# @dataclass
# class HMMSearchConfig(Config):
#     bias_filter: bool = True
#     null2: bool = True
#     seed: int = 42
#     Z: float = None
#     domZ: float = None
#     F1: float = 0.02
#     F2: float = 1e-3
#     F3: float = 1e-5
#     E: float = 10.0
#     T: float = None
#     domE: float = 10.0
#     domT: float = None
#     incE: float = 0.01
#     incT: float = None
#     incdomE: float = 0.01
#     incdomT: float = None
#     bit_cutoffs: Literal["gathering", "trusted", "noise"] = None


class ISScannerResult:
    def __init__(
            self, genome_id: str
    ):
        self.genome_id = genome_id

    def __iter__(self):
        return iter(self.pieces)

    def __len__(self):
        if not self.pieces:
            return 0
        else:
            return sum(len(i) for i in self.pieces)

    def __repr__(self):
        return f'{self.genome_id}; {len(self.pieces)} piece(s) across {self.n_contigs} contig(s)'

    def __format__(self, __format_spec: Literal['fasta', 'fna', 'ffn', 'faa', 'bed', 'gfa'] = '') -> str:
        if __format_spec == '':
            return self.__str__()
        # elif __format_spec in {'fasta', 'fna', 'ffn', 'faa', 'bed'}:
        #     return ''.join(format(i, __format_spec) for i in self.pieces)
        # elif __format_spec == 'gfa':
        #     return ''.join(format(i, __format_spec) for i in chain(self.pieces, self.edges))
        else:
            raise NotImplementedError(f'Invalid format: {__format_spec}')


class ISScannerWarning(ErisWarning):
    pass


class ISScannerError(Exception):
    pass


# class ISScanner:
#     def __init__(
#             self
#     ):

def main():
    from eris.minimap2 import Minimap2, Minimap2AlignConfig, Minimap2IndexConfig
    from eris.io import SeqFile
    from eris.alignment import cull_all, group_alignments
    queries = list(SeqFile('eris/data/IS.fna'))
    genome = Genome.from_file('test/ERR4920392.gfa')

    aligner = Minimap2(genome, index_config=Minimap2IndexConfig(k=11, w=5))
    aligner = Minimap2(genome)

    alignments = {k: list(cull_all(v)) for k, v in group_alignments(aligner.align(queries), 'target')}

    tirs = []
    for alignment in aligner.align(genome, config=Minimap2AlignConfig(c=True, x='sr', rev_only=True, n=1, m=15)):
        if alignment.query == alignment.target:
            tirs.append(alignment)

