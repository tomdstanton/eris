"""
Module for running the pangenome pipeline on bacterial pangenome graphs.
"""
# from eris.io import Genome
# from eris.db import Database
# from eris.external import Minimap2
# from eris.alignment import group_alignments, cull_all
#
# def _pipeline():
#     db = Database()
#     pangenome = Genome.from_file('tests/ST17/ST17_bifrost.gfa')
#     unitigs_with_is = []
#     with Minimap2(targets=pangenome) as aligner:  # Align IS to unitigs
#         for unitig, alignments in group_alignments(aligner.align(db), 'target'):
#             pangenome[unitig].add_features(*(i.as_feature('mobile_element') for i in cull_all(alignments)))
#             unitigs_with_is.append(unitig)  # Add unitig name to list, also indicator that we have alignments
#
#     pangenome_graph = pangenome.as_assembly_graph()
