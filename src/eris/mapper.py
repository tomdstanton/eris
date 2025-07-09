# """
# Copyright 2025 Tom Stanton (tomdstanton@gmail.com)
# https://github.com/tomdstanton/eris
#
# This file is part of eris. eris is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version. eris is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details. You should have received a copy of the GNU General Public License along with eris.
# If not, see <https://www.gnu.org/licenses/>.
# """
# from typing import Literal
# from dataclasses import dataclass
#
# from eris import ErisWarning, RESOURCES
# from eris.db import Database
# from eris.io import ReadSet, ReadFile, SeqFile
# from eris.utils import Config
# from eris.alignment import cull_all, group_alignments
# from eris.external import Minimap2, Minimap2AlignConfig, Minimap2IndexConfig
#
#
# # Classes --------------------------------------------------------------------------------------------------------------
# class InsertionSequence:
#     def __init__(self, feature: Feature, contig: Record):
#         self.feature: Feature = feature
#         # self.contig_gc: float = contig.seq.GC()
#         self.contg_copy_number: float = next((v for k, v in contig.qualifiers if k == 'dp'), 0)
#         self.CDS_in_element: set[Feature] = set()
#         self.CDS_flanking_element: set[Feature] = set()
#         self.edges: set[Edge] = set()
#
#     def __repr__(self):
#         return self.feature.__repr__()
#
#     def __str__(self):
#         return self.feature.__str__()
#
#     def __len__(self):
#         return len(self.feature)
#
#     def __format__(self, __format_spec: Literal['tsv'] = '') -> str:
#         if __format_spec == '':
#             return self.__str__()
#         elif __format_spec == 'tsv':
#             results = ['\t'.join(f'{self}', self.feature.location.parent_id, str(self.feature.qualifiers))]
#
#             return '\n'.join(results)
#
#         else:
#             raise NotImplementedError(f'Invalid format: {__format_spec}')
#
#
# class ISMapperResult:
#     def __init__(self, readset_id: str):
#         self.readset_id = readset_id
#         self.insertion_sequences: list[InsertionSequence] = []
#
#     def __iter__(self):
#         return iter(self.insertion_sequences)
#
#     def __len__(self):
#         return len(self.insertion_sequences)
#
#     def __format__(self, __format_spec: Literal['fasta', 'fna', 'ffn', 'faa', 'bed', 'gfa'] = '') -> str:
#         if __format_spec == '':
#             return self.__str__()
#         # elif __format_spec in {'fasta', 'fna', 'ffn', 'faa', 'bed'}:
#         #     return ''.join(format(i, __format_spec) for i in self.pieces)
#         # elif __format_spec == 'gfa':
#         #     return ''.join(format(i, __format_spec) for i in chain(self.pieces, self.edges))
#         else:
#             raise NotImplementedError(f'Invalid format: {__format_spec}')
#
#
# class ISMapperWarning(ErisWarning):
#     pass
#
#
# class ISMapperError(Exception):
#     pass
#
#
# class ISMapper:
#     def __init__(self, align_config: Minimap2AlignConfig = None, index_config: Minimap2IndexConfig = None,
#                  gene_finder_config: GeneFinderConfig = None, pool: Executor = None):
#         self.db = Database()
#         self.align_config = align_config or Minimap2AlignConfig(c=True)
#         self.index_config = index_config or Minimap2IndexConfig()
#
#     def cleanup(self):
#         pass
#
#     def __enter__(self):
#         return self
#
#     def __exit__(self, exc_type, exc_val, exc_tb):
#         self.cleanup()
#
#     def __del__(self):
#         self.cleanup()
#
#     def __call__(self, *readsets: ReadSet):
#         yield from map(self._pipeline, readsets)
#
#     def _pipeline(self, readset: ReadSet) -> Union[ISMapperResult, None]:
#         """
#         The scanner pipeline takes a single bacterial genome, predicts ORFs if necessary, and uses Minimap2 to find IS
#         elements with nucleotide-nucleotide alignment against contigs. Then, ORFs belonging to- or flanked by- the IS
#         elements are identified (including potential graph traversal), and are returned in a result object for later
#         reporting.
#         """
#         if not isinstance(readset, ReadSet):
#             readset = ReadSet.from_file(readset)
#
