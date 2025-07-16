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
from eris.io import Genome
from eris.db import Database
from eris.external import Minimap2
from eris.alignment import group_alignments, cull_all

def _pipeline():
    db = Database()
    pangenome = Genome.from_file('tests/ST17/ST17_bifrost.gfa')
    unitigs_with_is = []
    with Minimap2(targets=pangenome) as aligner:  # Align IS to unitigs
        for unitig, alignments in group_alignments(aligner.align(db), 'target'):
            pangenome[unitig].add_features(*(i.as_feature('mobile_element') for i in cull_all(alignments)))
            unitigs_with_is.append(unitig)  # Add unitig name to list, also indicator that we have alignments

    pangenome_graph = pangenome.as_assembly_graph()


    neighbours = pangenome_graph.get_neighbors(unitig)