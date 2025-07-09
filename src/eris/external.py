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
from tempfile import NamedTemporaryFile
from typing import Iterable, Generator, Union, Literal
from subprocess import Popen, PIPE, DEVNULL
from dataclasses import dataclass, asdict
from warnings import warn

from eris import RESOURCES
from eris.seq import Record, Feature
from eris.io import SeqFile, ReadSet
from eris.utils import Config, find_executable_binaries, is_non_empty_file
from eris.alignment import Alignment


# Classes --------------------------------------------------------------------------------------------------------------
class ExternalProgramError(Exception):
    pass

class ExternalProgram:
    """
    Base class to handle an external program to be executed in subprocesses
    """
    def __init__(self, program: str):
        if (binary := next(find_executable_binaries(program))) is None:
            raise ExternalProgramError(f'Could not find {program}')
        self._program = program
        self._binary = binary
        
    def __repr__(self):
        return f'{self._program}({self._binary})'
        
    def __enter__(self):
        return self

    def _run(self, command, input_: str = None, stdout: bool = True):
        with Popen(f"{self._binary} {command}", shell=True, stderr=PIPE, stdin=PIPE if input_ else DEVNULL,
                   stdout=PIPE if stdout else DEVNULL, universal_newlines=True) as process:
            return process.communicate(input=input_)  


@dataclass
class Minimap2IndexConfig(Config):
    t: int = RESOURCES.available_cpus  # Number of threads [3]
    H: bool = False  # use homopolymer-compressed k-mer (preferable for PacBio)
    k: int = None  # k-mer size (no larger than 28) [15]
    w: int = None  # minimizer window size [10]
    I: int = None  # split index for every ~NUM input bases [8G]


@dataclass
class Minimap2AlignConfig(Config):
    x: Literal[
        'lr:hq', 'splice', 'splice:hq', 'asm5', 'asm10', 'asm20', 'sr', 'map-pb', 'map-hifi', 'map-ont', 'map-iclr',
        'ava-pb', 'ava-ont'
    ] = None
    t: int = RESOURCES.available_cpus  # Number of threads [3]
    f: float = None  # filter out top FLOAT fraction of repetitive minimizers [0.0002]
    g: int = None  # stop chain elongation if there are no minimizers in INT-bp [5000]
    G: int = None  # max intron length (effective with -xsplice; changing -r) [200k]
    F: int = None  # max fragment length (effective with -xsr or in the fragment mode) [800]
    r: int = None  # chaining/alignment bandwidth and long-join bandwidth [500,20000]
    n: int = None  # minimal number of minimizers on a chain [3]
    m: int = None  #   minimal chaining score (matching bases minus log gap penalty) [40]
    X: bool = False  # skip self and dual mappings (for the all-vs-all mode)
    p: int = None  # min secondary-to-primary score ratio [0.8]
    N: int = None  # retain at most INT secondary alignments [5]
    A: int = None  # matching score [2]
    B: int = None  # mismatch penalty (larger value for lower divergence) [4]
    O: int = None  # gap open penalty [4,24]
    E: int = None  # gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [2,1]
    z: int = None  # Z-drop score and inversion Z-drop score [400,200]
    s: int = None  # minimal peak DP alignment score [80]
    u: Literal['f', 'b', 'n'] = None  # how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]
    J: Literal[0, 1] = None  # splice mode. 0: original minimap2 model; 1: miniprot model [1]
    a: bool = False  # output in the SAM format (PAF by default)
    c: bool = False  # output CIGAR in PAF
    cs: bool = False  # output the cs tag; STR is 'short' (if absent) or 'long' [none]
    ds: bool = False  # output the ds tag, which is an extension to cs
    MD: bool = False  # output the MD tag
    eqx: bool = False  # write =/X CIGAR operators
    rev_only: bool = False
    for_only: bool = False

    @classmethod
    def from_reads(cls, reads: ReadSet, **kwargs):
        config = cls(**kwargs)
        if reads.set_type in ['illumina', 'short']:
            config.x = 'sr'
        elif reads.set_type == 'pb':
            config.x = 'map-pb'
        elif reads.set_type == 'ont':
            config.x = 'map-ont'
        else:
            config.x = 'lr:hq'
        return config


class Minimap2Error(ExternalProgramError):
    pass


class Minimap2(ExternalProgram):
    """
    An aligner instance using Minimap2. The class can be instantiated with targets which will create a temporary
    index file, that can be handled safely with context management.

    Example:
        >>> from eris.external import Minimap2
        ... from eris.io import Genome
        ... with Minimap2(targets=Genome.from_file('tests.fasta')) as aligner:  # Build the index from genome 1
        ...     for alignment in aligner(Genome.from_file('test2.gfa')):
        ...         print(alignment)
    """
    def __init__(self, targets: Union[str, Path, Union[Record, Feature], Iterable[Union[Record, Feature]]] = None,
                 align_config: Minimap2AlignConfig = None, index_config: Minimap2IndexConfig = None):
        super().__init__('minimap2')
        self._target_index: Path = None
        self._align_config: Minimap2AlignConfig = align_config or Minimap2AlignConfig()
        self._index_config: Minimap2IndexConfig = index_config or Minimap2IndexConfig()
        if targets:
            self.build_index(targets)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()

    def __del__(self):
        self.cleanup()

    def __call__(self, *args, **kwargs):
        return self.align(*args, **kwargs)

    def cleanup(self):
        """
        Cleans up resources related to the target index if it is set.

        This method checks if there is a target index specified and, if so, closes the
        index to release any held resources. It also resets the target index
        reference to None, ensuring subsequent operations do not attempt to access
        a closed resource.

        Raises:
            No explicit exceptions are raised, but calling operations on a closed
            resource externally could produce runtime errors.
        """
        if self._target_index and not self._target_index.closed:
            self._target_index.close()
            self._target_index = None

    @staticmethod
    def _validate_seqs(*seqs: Union[str, Path, SeqFile, Record, Feature], pipe_arg: str = '-') -> tuple[
        str, [Union[str, None]]]:
        """
        Checks query or target inputs based on the types passed
        Returns:
            A tuple of the argument to be passed to the subprocess and, if necessary, the piped input (or none)
        """
        # If pipe, the input MUST be a mixture of str, Record, Feature, else MUST be Path, SeqFile
        proc_arg, pipe = [], ''
        for i in seqs:
            if isinstance(i, (Path, SeqFile)):
                # Check arg_types to make sure it doesn't contain pipe-able types (Record, Feature, str)
                if pipe:
                    raise ExternalAlignerError('Cannot mix pipe-able inputs with file-like inputs')
                proc_arg.append(str(i))  # Requires SeqFile str representation to be SeqFile.path.__str__
            else:
                if proc_arg:
                    raise ExternalAlignerError('Cannot mix pipe-able inputs with file-like inputs')
                if isinstance(i, str):
                    if not i.startswith('>') or not i.startswith('@'):
                        raise ExternalAlignerError('If input is a string, it must be in sequence format')
                    pipe += i
                else:
                    pipe += format(i, 'fasta')

        return ' '.join(proc_arg) if proc_arg else pipe_arg, pipe or None

    def build_index(
            self,
            targets: Union[str, Path, SeqFile, Record, Feature, Iterable[Union[str, Path, SeqFile, Record, Feature]]],
            config: Minimap2IndexConfig = None
    ) -> Path:
        """
        Builds an index for the targets and writes to a temporary file

        :param targets: The targets to index; can be a fasta-formatted string, path to a fasta file,
                        ``Records`` or ``Features``
        :returns: A Path to the temporary file containing the index
        """
        arg, pipe = self._validate_seqs(*targets)
        if not config:
            config = self._index_config
        if self._target_index is None:
            self._target_index = NamedTemporaryFile()
        _, stderr = self._run(f'{_build_params(config)} -d {self._target_index.name} {arg}', pipe, False)
        if not is_non_empty_file(self._target_index.name):
            self._target_index = None
            raise Minimap2Error(f'Failed to build index, {stderr=}')
        return Path(self._target_index.name)

    def align(
            self,
            queries: Union[str, Path, SeqFile, Record, Feature, Iterable[Union[str, Path, SeqFile, Record, Feature]]],
            targets: Union[str, Path, SeqFile, Record, Feature, Iterable[Union[str, Path, SeqFile, Record, Feature]]] = None,
            config: Minimap2AlignConfig = None
    ) -> Generator[Alignment, None, None]:
        """
        Aligns queries to targets.

        :param queries: Alignment queries; can be a fasta-formatted string, path to a fasta file, ``Records`` or
                        ``Features``

        :param targets: Optional Alignment targets; can be a fasta-formatted string, path to a fasta file, ``Records``
                        or ``Features``; will otherwise use the target index built at instantiation.
                        If ``targets`` is passed and a target index exists, the target index is ignored.
                        If ``query`` is not a path and ``targets`` is not a path, a new target index will need to be
                        built, which will overwrite any existing target index, printing a ``erisAlignmentWarning``.

        :returns: A generator of ``Alignment`` objects
        """
        query_arg, query_pipe = self._validate_seqs(*queries)
        if not config:
            config = self._align_config
        if targets:
            target_arg, target_pipe = self._validate_seqs(*targets)
            if query_pipe and target_pipe:  # Both are stdin, so we need to create the target index
                if self._target_index:  # Existing target index, raise an error rather than overwriting
                    warn('Overwriting existing target index', ExternalAlignerWarning)
                target_arg, target_pipe = self.build_index(targets).name, None
        else:
            if self._target_index:
                target_arg, target_pipe = self._target_index.name, None
            else:
                raise Minimap2Error('Targets must be supplied if no target index has been built')

        stdout, stderr = self._run(f'{_build_params(config)} {target_arg} {query_arg}', query_pipe or target_pipe)

        for line in stdout.splitlines():
            yield Alignment.from_sam(line) if config.a else Alignment.from_paf(line)


# class FragGeneScanRs:
#     """To use as a potential alternative to pyrodigal/prodigal - it is faster and works on reads"""
#     def __init__(self):
#         super().__init__('FragGeneScanRs')


# Functions ------------------------------------------------------------------------------------------------------------
def _build_params(config: Config) -> str:
    """
    Builds CLI parameters from a config object, returning a string.
    Follows the convention that single letter arguments are preceded by a single hyphen and multiple-letter arguments
    are preceded by a double hyphen.
    """
    return ' '.join(f'{"-" if len(k) == 1 else "--"}{k.replace("_", "-")} {v if not isinstance(v, bool) else ""}'.strip()
                    for k, v in asdict(config).items() if v).strip()
