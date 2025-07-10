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
from typing import Iterable, Callable, Any, Union, IO
from argparse import RawTextHelpFormatter, RawDescriptionHelpFormatter, ArgumentParser
# from importlib.metadata import metadata
from os import fstat
from io import IOBase
from sys import stdout, stderr, stdin
from pathlib import Path
import time

from eris.io import SeqFile, Genome
from eris.utils import bold, get_logo


# Constants ------------------------------------------------------------------------------------------------------------
# _PACKAGE = Path(__file__).parent.name
_PACKAGE = 'eris'
# _METADATA = metadata(_PACKAGE)
# _VERSION = _METADATA.json['version']


# Classes --------------------------------------------------------------------------------------------------------------
class ResultWriter:
    """
    A class to handle the writing of results to multiple types of files; to be used with the CLI.
    """
    def __init__(self, outputs: tuple[tuple[str, Union[str, Path, IO]]], suffix: str = '_eris_results',
                 no_tsv_header: bool = False, tsv_header: str = '#'):
        """
        :param suffix: Suffix to append to output filenames (default: _eris_results)
        :param no_tsv_header: Suppress header in TSV output (default: False)
        """
        self.files = {}
        self.handles = {}
        self._suffix: str = suffix
        self.tsv_header: str = tsv_header
        self.header_written: bool = no_tsv_header
        for format, argument in outputs:
            if argument is not None:
                if isinstance(argument, str):
                    self.files[format] = Path(argument)
                elif isinstance(argument, Path):
                    self.files[format] = argument
                elif isinstance(argument, IOBase):
                    self.handles[format] = argument
                else:
                    raise TypeError(f"{format=} {argument=} must be a string, Path or IO, not {type(argument)}")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def open(self, *args, **kwargs):
        """Opens the handles that aren't stdout or stderr"""
        for handle in self.handles.values():
            if handle.name not in {'<stdout>', '<stderr>'}:
                handle.open(*args, **kwargs)

    def close(self):
        """Close the handles that aren't stdout or stderr"""
        for handle in self.handles.values():
            if handle.name not in {'<stdout>', '<stderr>'}:
                handle.close()

    def write(self, result: Union['Result', None], open_mode: str = 'wt'):
        """
        Write the typing result to files or file handles
        :param result: A Result instance or None
        """
        if result is None:
            return None

        for fmt, file in self.files.items():  # Write to the file outputs
            with open(file / f'{result.genome_id}{self._suffix}.{fmt}', mode=open_mode) as handle:
                handle.write(format(result, fmt))

        for fmt, handle in self.handles.items():  # Write to the handle outputs
            # Logic for non-repeating TSV header, this limits TSV to a handle but should be the case
            if fmt == 'tsv' and not self.header_written:
                handle.write(self.tsv_header)
                self.header_written = True
            handle.write(format(result, fmt))


class ProgressBar:
    """
    A CLI progress bar that wraps an iterable and calls a function for each item.

    This class acts as an iterator, yielding results from the callable while
    printing a dynamic progress bar to stderr.

    Args:
        iterable (Iterable): An iterable of items to process.
        callable (Callable): A function to call for each item from the iterable.
        desc (str): A description to display before the progress bar.
        bar_length (int): The character length of the progress bar.
        bar_character (str): The character used to fill the progress bar.
    """
    def __init__(self, iterable: Iterable, callable: Callable[[Any], Any],
                 desc: str = "Processing items", bar_length: int = 40, bar_character: str = '█'):
        # Eagerly consume the iterable to get a total count for the progress bar.
        assert len(bar_character) == 1, "Bar character must be a single character"
        self.items = list(iterable)
        self.total = len(self.items)
        self.callable = callable
        self.desc = desc
        self.bar_length = bar_length
        self.bar_character = bar_character
        self._iterator: Iterable = None
        self.start_time: float = None
        self._processed_count: int = 0

    def _format_time(self, seconds: float) -> str:
        """Formats seconds into a HH:MM:SS string."""
        if not isinstance(seconds, (int, float)) or seconds < 0:
            return "00:00:00"
        m, s = divmod(seconds, 60)
        h, m = divmod(m, 60)
        return f"{int(h):02d}:{int(m):02d}:{int(s):02d}"

    def __iter__(self):
        self._iterator = iter(self.items)
        self.start_time = time.time()
        self._processed_count = 0
        self._update_progress()  # Display the initial (0%) bar
        return self

    def __next__(self):
        # The for-loop protocol will handle the StopIteration from next()
        item = next(self._iterator)

        # The core of the wrapper: call the provided callable for one item.
        result = self.callable(item)
        self._processed_count += 1
        self._update_progress()
        return result

    def _update_progress(self):
        """Calculates and prints the progress bar to stderr."""
        if self.total == 0:
            return

        percent_complete = self._processed_count / self.total
        filled_length = int(self.bar_length * percent_complete)
        bar = self.bar_character * filled_length + '-' * (self.bar_length - filled_length)

        elapsed_time = time.time() - self.start_time

        # Calculate Estimated Time of Arrival (ETA)
        if self._processed_count > 0:
            avg_time_per_item = elapsed_time / self._processed_count
            remaining_items = self.total - self._processed_count
            eta = avg_time_per_item * remaining_items
        else:
            eta = float('inf')

        # Format time strings
        elapsed_str = self._format_time(elapsed_time)
        eta_str = self._format_time(eta) if eta != float('inf') else '??:??:??'

        # Use carriage return '\r' to stay on the same line
        progress_line = (f'\r{self.desc}: {int(percent_complete * 100):>3}%|{bar}| '
                         f'{self._processed_count}/{self.total} '
                         f'[{elapsed_str}<{eta_str}]')

        stderr.write(progress_line)

        # When the loop is finished, print a newline to move to the next line
        if self._processed_count == self.total:
            stderr.write('\n')

        stderr.flush()

    def __len__(self):
        return self.total


# Functions ------------------------------------------------------------------------------------------------------------
def scan_parser(subparsers):
    parser = subparsers.add_parser(
        'scan', description=get_logo('Scan for IS in bacterial genomes'),
        epilog=f'For more help, visit: {bold('eris.readthedocs.io')}',
        add_help=False, formatter_class=RawTextHelpFormatter,
        help='Scan for IS in bacterial genomes', usage="%(prog)s <genome> <genome...> [options]"
    )
    inputs = parser.add_argument_group(bold('Inputs'), '')
    inputs.add_argument(
        'genome', nargs='*', default=[stdin],
        help='Genome(s) in FASTA, GFA or Genbank format. If not provided, reads from stdin.'
   )
    outputs = parser.add_argument_group(bold('Outputs'), "")
    outputs.add_argument(
        '--tsv', metavar='', default=stdout, help='Path to output tabular results (default: stdout)'
    )
    outputs.add_argument('--no-tsv-header', action='store_true', help='Suppress header in TSV output')
    opts = parser.add_argument_group(bold('Other options'), '')
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', action='help')


# Main CLI Entry Point -------------------------------------------------------------------------------------------------
def main():
    parser = ArgumentParser(
        description=get_logo('Uncovering IS-mediated discord in bacterial genomes'),
        usage="%(prog)s <command>",
        add_help=False,
        prog=_PACKAGE,
        formatter_class=RawDescriptionHelpFormatter,
        epilog=f'For more help, visit: {bold("eris.readthedocs.io")}'
    )
    subparsers = parser.add_subparsers(
        title=bold('Command'), dest='command', prog=_PACKAGE, metavar='<command>', help=None, required=True
    )
    scan_parser(subparsers)
    opts = parser.add_argument_group(bold('Other options'), '')
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', action='help')

    args = parser.parse_args()

    if args.command == 'scan':
        from eris.scanner import Scanner, TSV_HEADER
        with (
            Scanner() as scanner,
            ResultWriter(
                outputs=(('tsv', args.tsv),),
                tsv_header=TSV_HEADER,
                no_tsv_header=args.no_tsv_header
            ) as writer
        ):
            for result in ProgressBar(args.genome, scanner, "Scanning genomes"):
                writer.write(result)
