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
from argparse import RawTextHelpFormatter, RawDescriptionHelpFormatter, ArgumentParser
from importlib.metadata import metadata
from sys import stdout

from eris.io import SeqFile, Genome
from eris.utils import bold, get_logo


# Constants ------------------------------------------------------------------------------------------------------------
_PACKAGE_NAME = 'eris'
_METADATA = metadata(_PACKAGE_NAME)
_VERSION = _METADATA.json['version']


# Main CLI Entry Point -------------------------------------------------------------------------------------------------
def main():
    parser = ArgumentParser(
        description=get_logo('Uncovering IS-mediated discord in bacterial genomes'),
        usage="%(prog)s <module>",
        add_help=False,
        prog=_METADATA,
        formatter_class=RawDescriptionHelpFormatter,
        epilog=f'For more help, visit: {bold("eris.readthedocs.io")}'
    )
    subparsers = parser.add_subparsers(
        title=bold('module'), dest='module', prog=_PACKAGE_NAME, metavar='<module>', help=None, required=True
    )

    scan_parser = subparsers.add_parser(
        'scan', description=get_logo('Scan for IS in bacterial genomes'),
        epilog=f'For more help, visit: {bold('eris.readthedocs.io')}',
        add_help=False, formatter_class=RawTextHelpFormatter,
        help='Scan for IS in bacterial genomes', usage="eris scan <genomes> [options]"
    )
    scan_parser.add_argument(
        'genomes', help='Genomes in FASTA, GFA or Genbank format; use - for stdin',
        metavar='<genomes>', nargs='+'
    )
    scan_parser.add_argument_group(bold('Outputs'), "\nNote, text outputs accept '-' for stdout")
    scan_parser.add_argument('--tsv', metavar='', default='-',
                      help='Path to output tabular results (default: stdout)')
    scan_parser.add_argument('--no-tsv-header', action='store_true',
                             help='Suppress header in TSV output')

    parser.add_argument_group(bold('Other options'), '')
    parser.add_argument('-v', '--version', help='Show version number and exit', action='version')
    parser.add_argument('-h', '--help', help='Show this help message and exit', action='help')

    args = parser.parse_args()

    if args.module == 'scan':
        from eris.scanner import Scanner, ScannerResultWriter
        with Scanner() as scanner, ScannerResultWriter(tsv=args.tsv, no_tsv_header=args.no_tsv_header) as writer:
            for result in scanner(*args.genomes):
                writer.write(result)
