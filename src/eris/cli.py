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
from argparse import RawDescriptionHelpFormatter, ArgumentParser
# _METADATA = metadata.metadata
# __version__ = _METADATA

from eris.utils import bold, get_logo
from eris.db import Database


# Main CLI Entry Point -------------------------------------------------------------------------------------------------
def main():
    # TODO: See if these arguments can be obtained from the package metadata
    parser = ArgumentParser(
        description=get_logo('A suite for bacterial surface polysaccharides'),
        usage="%(prog)s <module>",
        add_help=False,
        prog="eris",
        formatter_class=RawDescriptionHelpFormatter,
        epilog=f'For more help, visit: {bold("eris.readthedocs.io")}'
    )
    subparsers = parser.add_subparsers(
        title=bold('module'), dest='module', prog='eris', metavar='<module>', help=None, required=True
    )
    typing.add_parser(subparsers)
    db.add_parser(subparsers)
    opts = parser.add_argument_group(bold('Other options'), "")
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', action='help')

    # This doesn't allow the user to use the following CLI modules if the extra packages aren't installed
    # if dependencies.is_available('eris.predict'):
        # predict.add_parser(subparsers)

    args = parser.parse_args()

    if args.module == 'typing':
        typing.main(args)

    elif args.module == 'db':
        db.main(args)


# Functions ------------------------------------------------------------------------------------------------------------
def add_parser(parsers):
    parser = parsers.add_parser(
        'db',
        description=get_logo('eris Databases'),
        add_help=False,
        formatter_class=RawTextHelpFormatter,
        help='Explore, manage or create eris databases',
        usage="eris db <submodule>"
    )
    subparsers = parser.add_subparsers(title=bold('Submodule\n'), dest='submodule', metavar="<submodule>", required=True)
    _add_convert_subparser(subparsers)
    opts = parser.add_argument_group(bold('Other options'), "")
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', action='help')


def _add_convert_subparser(parsers):
    parser = parsers.add_parser(
        'convert',
        description=get_logo('Convert entries from a eris database'),
        # epilog=f'For more help, visit: {bold(_URL)}',
        add_help=False,
        formatter_class=RawTextHelpFormatter,
        help='Convert entries from a eris database',
        usage="eris db convert <db> [formats] [options]"
    )
    add_db_arg(parser)
    opts = parser.add_argument_group(bold('Other options'), "")
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', action='help')


def add_db_arg(parser):
    parser.add_argument(
        'db', metavar='db path/keyword', type=Database,
        help="eris database path or keyword\nUse 'eris db available' to see available databases and keywords"
    )


def add_db_options(parser):
    parser = parser.add_argument_group(bold('Database options'), "")
    parser.add_argument('--locus-regex', type=compile, metavar='',
                      help=f'Python regular-expression to match locus names in db source note')
    parser.add_argument('--type-regex', type=compile, metavar='',
                      help=f'Python regular-expression to match locus types in db source note')
    parser.add_argument('--filter', type=compile, metavar='',
                      help='Python regular-expression to select loci to include in the database')


# Functions ------------------------------------------------------------------------------------------------------------
def add_parser(parsers):
    parser = parsers.add_parser(
        'type',
        description=get_logo('In silico serotyping'),
        add_help=False,
        formatter_class=RawTextHelpFormatter,
        help='In silico serotyping',
        prog='eris type',
        usage="%(prog)s <submodule>",
    )
    subparsers = parser.add_subparsers(title=bold('Submodule\n'), dest='submodule', metavar="<submodule>", required=True)
    _add_assembly_subparser(subparsers)
    _add_read_subparser(subparsers)
    opts = parser.add_argument_group(bold('Other options'), "")
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', action='help')


def _add_assembly_subparser(parsers):
    parser = parsers.add_parser(
        'assemblies',
        description=get_logo('In silico serotyping of assemblies'),
        add_help=False,
        formatter_class=RawTextHelpFormatter,
        help='In silico serotyping of assemblies',
        prog='eris type assemblies',
        usage="%(prog)s <db> <input> ... [options]"
    )
    parser = parser.add_argument_group(bold('Inputs'), "")
    add_db_arg(parser)
    parser.add_argument('input', nargs='+', metavar='input', type=Path, help='Assemblies in FASTA, GFA or Genbank format')

    parser = parser.add_argument_group(bold('Output options'), "\nNote, text outputs accept '-' for stdout")
    parser.add_argument('-o', '--out', metavar='', default=stdout, type=FileType('at'),
                      help='Output file to write/append tabular results to (default: stdout)')
    parser.add_argument('-f', '--fasta', metavar='', nargs='?', default=None, const='.', type=write_to_file_or_directory,
                      help='Turn on fasta output\n'
                           'Accepts a single file or a directory (default: cwd)')
    parser.add_argument('-j', '--json', metavar='', nargs='?', default=None, const='eris_results.json', type=FileType('at'),
                      help='Turn on JSON lines output\n'
                           'Optionally choose file (can be existing) (default: %(const)s)')

    parser = parser.add_argument_group(bold('Scoring options'), "")
    parser.add_argument('--min-cov', type=float, required=False, default=50.0, metavar='',
                      help='Minimum gene %%coverage (blen/q_len*100) to be used for scoring (default: %(default)s)')
    parser.add_argument("--score-metric", metavar='', default=0, type=int, choices=range(4),
                      help="Metric for scoring each locus (default: %(default)s)\n"
                           "  0: AS (alignment score of genes found)\n"
                           "  1: mlen (matching bases of genes found)\n"
                           "  2: blen (aligned bases of genes found)\n"
                           "  3: q_len (query length of genes found)")
    parser.add_argument("--weight-metric", metavar='', default=3, type=int, choices=range(6),
                      help="Weighting for the 1st stage of the scoring algorithm (default: %(default)s)\n"
                           "  0: No weighting\n"
                           "  1: Number of genes found\n"
                           "  2: Number of genes expected\n"
                           "  3: Proportion of genes found\n"
                           "  4: blen (aligned bases of genes found)\n"
                           "  5: q_len (query length of genes found)")
    parser.add_argument('--n-best', type=int, default=2, metavar='',
                      help='Number of best loci from the 1st round of scoring to be\n'
                           'fully aligned to the assembly (default: %(default)s)')

    parser = parser.add_argument_group(bold('Confidence options'), "")
    parser.add_argument("--gene-threshold", type=float, metavar='',
                      help="Species-level locus gene identity threshold (default: database specific)")
    parser.add_argument("--max-other-genes", type=int, metavar='', default=1,
                      help="Typeable if <= other genes (default: %(default)s)")
    parser.add_argument("--percent-expected", type=float, metavar='', default=50,
                      help="Typeable if >= %% expected genes (default: %(default)s)")
    parser.add_argument("--below-threshold", type=bool, default=False, metavar='',
                      help="Typeable if any genes are below threshold (default: %(default)s)")
    add_db_options(parser)
    opts = parser.add_argument_group(bold('Other options'), "")
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', action='help')


def _add_read_subparser(parsers):
    parser = parsers.add_parser(
        'reads',
        description=get_logo('In silico serotyping of reads'),
        add_help=False,
        formatter_class=RawTextHelpFormatter,
        help='In silico serotyping of reads',
        prog='eris type reads',
        usage="%(prog)s <db> <input> ... [options]"
    )
    parser = parser.add_argument_group(bold('Inputs'), "")
    add_db_arg(parser)
    parser.add_argument(
        'input', nargs='+', metavar='reads', type=ReadFile, help='Read sets in FASTQ format'
    )
    add_db_options(parser)
    opts = parser.add_argument_group(bold('Other options'), "")
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', action='help')


# def _other_opts(parser):
#     opts = parser.add_argument_group(bold('Other options'), "")
#     opts.add_argument('--no-header', action='store_true', help='Suppress header line')
#     opts.add_argument('-t', '--threads', type=int, default=system_resources.available_cpus, metavar='',
#                       help="Number of alignment threads or 0 for all available (default: %(default)s)")
#     opts.add_argument('-v', '--version', help='Show version number and exit', metavar='')
#     opts.add_argument('-h', '--help', help='Show this help message and exit', metavar='')


# Main -----------------------------------------------------------------------------------------------------------------
# def main(args: Namespace):
#     pass
    # from eris.typing import AssemblyTyper
    # # with AssemblyTyper('kpsc_k', auxiliary_genes='eris/typing/IS.fna') as typer, open('test.tsv', 'wt') as tsv:
    # #     for result in typer(*Path('eris/typing/test/data').glob("*.fasta")):
    # #         tsv.write(f"{result:tsv}")
    # with AssemblyTyper(args.db) as typer:
    #     for result in typer(*args.input):
    #         pass
