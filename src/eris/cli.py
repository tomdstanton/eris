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
from itertools import chain
from typing import Iterable, Callable, Any, Union, IO, Dict
from argparse import RawTextHelpFormatter, RawDescriptionHelpFormatter, ArgumentParser
from os import fstat
from io import IOBase
from sys import stdout, stderr, stdin
from pathlib import Path
import time
# from threading import Lock
from concurrent.futures import Executor, ThreadPoolExecutor

from eris import RESOURCES
from eris.io import SeqFile, Genome
from eris.utils import bold, get_logo, write_to_file_or_directory

# Constants ------------------------------------------------------------------------------------------------------------
_EPILOG = 'For more help, visit: ' + bold(f'{RESOURCES.package}.readthedocs.io')
_SUFFIX = f'_{RESOURCES.package}_results'


# Classes --------------------------------------------------------------------------------------------------------------
class ResultWriter:
    """
    A class to handle the writing of results to multiple types of files; to be used with the CLI.
    """
    def __init__(self, *outputs: tuple[str, Union[str, Path, IO]], suffix: str = _SUFFIX,
                 no_tsv_header: bool = False, tsv_header: str = '#', pool: Executor = None):
        """
        :param suffix: Suffix to append to output filenames
        :param no_tsv_header: Suppress header in TSV output (default: False)
        """
        self._files = {}
        self._handles = {}
        self.suffix: str = suffix
        self.tsv_header: str = tsv_header
        self._header_written: bool = no_tsv_header
        for format, argument in outputs:
            if argument is not None:
                if argument in {'-', 'stdout'}:
                    self._handles[format] = stdout
                    # self._handle_locks[format] = Lock()
                elif isinstance(argument, str):
                    self._files[format] = Path(argument)
                elif isinstance(argument, Path):
                    self._files[format] = argument
                elif isinstance(argument, IOBase):
                    self._handles[format] = argument
                    # self._handle_locks[format] = Lock()
                else:
                    raise TypeError(f"{format=} {argument=} must be a string, Path or IO, not {type(argument)}")

        if not self._files and not self._handles:
            raise ValueError("No outputs specified")

        # self.pool: Executor = pool or ThreadPoolExecutor(
        #     min(len(self._files) + len(self._handles), RESOURCES.available_cpus + 4)
        # )

    def __len__(self):
        return len(self._files) + len(self._handles)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        # Wait for all writes to finish before exiting context
        # self.pool.shutdown(wait=True)

    def __del__(self):
        self.close()
        # Avoid __del__ for resource management; __exit__ is preferred. If it is called, try to clean up.
        # self.pool.shutdown(wait=False)

    def close(self):
        """Close the handles that aren't stdout or stderr"""
        for handle in self._handles.values():
            if handle.name not in {'<stdout>', '<stderr>'}:  # These handles cannot be opened or closed
                handle.close()

    def _write(self, fmt: str, out: Union[str, Path, IOBase], result, open_mode: str) -> int:
        if fmt in self._handles:
            # Acquire lock for the specific handle to ensure thread-safe writing
            # with self._handle_locks[out]:
            #     # Logic for non-repeating TSV header
            #     if fmt == 'tsv' and not self._header_written:
            #         out.write(self.tsv_header)
            #         self._header_written = True
            #     out.write(format(result, fmt))
            # Logic for non-repeating TSV header
            if fmt == 'tsv' and not self._header_written:
                out.write(self.tsv_header)
                self._header_written = True
            return out.write(format(result, fmt))
        else:  # This logic is for self.files
            with open(out / f'{result.genome_id}{self.suffix}.{fmt}', mode=open_mode) as handle:
                return handle.write(format(result, fmt))

    def write(self, result: Union['Result', None], open_mode: str = 'wt'):
        """
        Write the typing result to files or file handles
        :param result: A Result instance or None
        """
        if result is None:
            return None
        # Use the pool to write the result, iterating over the handles and files dictionaries
        # for fmt, out in chain(self._handles.items(), self._files.items()):
        #     self.pool.submit(self._write, result, fmt, out, open_mode)
        for fmt, out in chain(self._handles.items(), self._files.items()):
            self._write(fmt, out, result, open_mode)


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
    name, desc = 'scan', 'Scan for IS in bacterial genomes'
    parser = subparsers.add_parser(
        name, description=get_logo(desc), prog=f'{RESOURCES.package} {name}', epilog=_EPILOG,
        formatter_class=RawTextHelpFormatter, help=desc, usage="%(prog)s <genome> <genome...> [options]", add_help=False
    )
    inputs = parser.add_argument_group(bold('Inputs'), '\nNote, input file(s) may be compressed.')
    inputs.add_argument(
        'genome', nargs='*', default=[stdin], metavar='<genome>',
        help='Genome(s) in FASTA, GFA or Genbank format; reads from stdin by default.\n'
             'Genome(s) in FASTA/GFA format can paired up with GFA/BED\nannotation files with the same prefix.'
   )
    outputs = parser.add_argument_group(
        bold('Outputs'),
        '\nNote, text outputs accept "-" or "stdout" for stdout'
        '\nIf a directory is passed, individual files will be written per input genome'
    )
    outputs.add_argument(  # TSV output can be written to a single file or one file per genome, default=stdout
        '--tsv', metavar='', default=stdout, help='Path to output tabular results (default: stdout)', nargs='?',
        type=write_to_file_or_directory
    )
    outputs.add_argument(  # FFN output can be written to a single file or one file per genome
        '--ffn', metavar='', help=f'Path to output feature DNA sequences in FASTA format\n'
                                  f'Defaults to "./[genome]{_SUFFIX}.ffn" when passed without arguments',
        const='.', nargs='?', type=write_to_file_or_directory, default=None
    )
    outputs.add_argument(  # FAA output can be written to a single file or one file per genome
        '--faa', metavar='', help=f'Path to output feature Amino acid sequences in FASTA format\n'
                                  f'Defaults to "./[genome]{_SUFFIX}.faa" when passed without arguments',
        const='.', nargs='?', type=write_to_file_or_directory, default=None
    )
    outputs.add_argument('--no-tsv-header', action='store_true', help='Suppress header in TSV output')
    opts = parser.add_argument_group(bold('Other options'), '')
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', action='help')


# Main CLI Entry Point -------------------------------------------------------------------------------------------------
def main():
    parser = ArgumentParser(
        description=get_logo('Uncovering IS-mediated discord in bacterial genomes'), epilog=_EPILOG,
        usage="%(prog)s <command>", add_help=False, prog=RESOURCES.package, formatter_class=RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(
        title=bold('Command'), dest='command', metavar='<command>', required=True, help=None,
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
                ('tsv', args.tsv), ('faa', args.faa), ('ffn', args.ffn),
                tsv_header=TSV_HEADER, no_tsv_header=args.no_tsv_header
            ) as writer
        ):
            for result in ProgressBar(args.genome, scanner, "Scanning genomes"):
                writer.write(result)
