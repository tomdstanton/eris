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
from dataclasses import dataclass, asdict
from os import environ, X_OK, pathsep, access
from pathlib import Path
from sys import stdout
from typing import IO, Generator, Union, Literal
from io import IOBase
from argparse import Namespace
from urllib.request import urlopen, Request
from gzip import (open as gz_open, decompress as gz_decompress)
from bz2 import (open as bz2_open, decompress as bz2_decompress)
from lzma import (open as xz_open, decompress as xz_decompress)

from eris import RESOURCES

# Constants ------------------------------------------------------------------------------------------------------------
_MAGIC_BYTES = {b'\x1f\x8b': 'gz', b'\x42\x5a': 'bz2', b'\xfd7zXZ\x00': 'xz'}
_OPEN = {'gz': gz_open, 'bz2': bz2_open, 'xz': xz_open}
_DECOMPRESS = {'gz': gz_decompress, 'bz2': bz2_decompress, 'xz': xz_decompress}
_MIN_N_BYTES = max(len(i) for i in _MAGIC_BYTES)  # Minimum number of bytes to read in a file to guess the compression
if 'zstandard' in RESOURCES.optional_packages:
    from zstandard import (open as zst_open, decompress as zst_decompress)
    _MAGIC_BYTES[b'\x28\xb5\x2f\xfd'] = 'zst'
    _OPEN['zst'] = zst_open
    _DECOMPRESS['zst'] = zst_decompress

# Classes --------------------------------------------------------------------------------------------------------------
@dataclass  # (kw_only=True) https://medium.com/@aniscampos/python-dataclass-inheritance-finally-686eaf60fbb5
class Config:
    """
    Config parent class that can conveniently set attributes from CLI args
    """

    @classmethod
    def from_args(cls, args: Namespace):
        """
        Sets attributes of the class from a Namespace object (e.g. from argparse)

        Parameters
        ----------
        args : :class:`argparse.Namespace`
            :class:`argparse.Namespace` object containing attributes to set

        Returns
        -------
        cls
            Class instance with attributes set from args

        """
        return cls(**{f: getattr(args, f) for f in asdict(cls) if hasattr(args, f)})

# Functions ------------------------------------------------------------------------------------------------------------
def find_executable_binaries(*programs: str) -> Generator[Path, None, None]:
    """
    Check if programs are installed and executable

    :param programs: List of programs to check
    :return: Generator of Path objects for each program found
    """
    programs, found = set(programs), 0
    while len(programs) > found:
        for path in map(Path, environ["PATH"].split(pathsep)):
            if path.is_dir():
                for binary in path.iterdir():
                    if binary.name in programs and access(binary, X_OK):
                        found += 1
                        yield binary


def is_non_empty_file(file: Union[str, Path], min_size: int = 1) -> bool:
    """
    Checks if a file exists, is a file, and is non-empty (optionally above a minimum size).

    :param file: Path to the file to check.
    :param min_size: Minimum size of the file in bytes.
    :return: True if the file exists, is a file, and is non-empty, False otherwise.
    """
    if not isinstance(file, Path):
        file = Path(file)
    if file.exists() and file.is_file():
        return file.stat().st_size >= min_size
    else:
        return False


def write_to_file_or_directory(path: Union[str, Path, IO], mode: str = 'at') -> Union[Path, IO]:
    """
    Writes to a file or creates a directory based on the provided path.

    If the path is '-' or 'stdout', it returns stdout.
    If the path has a suffix, it's treated as a file and opened for appending.
    If the path has no suffix, it's treated as a directory and created if it doesn't exist.

    :param path: The path to the file or directory.
    :param mode: The mode to open the file in if it's a file.
    :return: A file handle (IO) if a file is specified, or a Path object if a directory is specified.
    """
    if isinstance(path, IOBase):
        return path
    if path in {'-', 'stdout'}:  # If the path is '-', return stdout
        return stdout
    if not isinstance(path, Path):  # Coerce to Path object
        path = Path(path)
    if path.suffix:  # If the path has an extension, it's probably a file
        # NB: We can't use is_file or is_dir because it may not exist yet, `open()` will create or append
        return open(path, mode)  # Open the file
    else:
        path.mkdir(exist_ok=True, parents=True)  # Create the directory if it doesn't exist
    return path


def xopen(
        file: Union[str, Path],
        method: Literal['magic', 'guess', 'uncompressed', 'gz', 'bz2', 'xz', 'zst'] = 'magic',
        **open_args
) -> IO:
    """
    Opens a file with the appropriate open function based on the magic bytes at the beginning of the data.
    Inspired by `xopen <https://pypi.org/project/xopen/>`_

    :param file: Path to file for opening as a string or Path object.
    :param method: Method to determine how to open file; ``guess`` will guess from the file extension, ``magic``
                   will determine from the magic bytes at the beginning of the file, the other options force the file
                   type and should be one of ``uncompressed``, ``gz``, ``bz2``, ``xz`` or ``zst`` (if installed).
    :param open_args: Additional keyword arguments to pass to the open function.
    :return: File handle in text or binary mode depending on the arguments passed to the open function.
    """
    if method not in {'magic', 'guess', 'uncompressed', 'gz', 'bz2', 'xz', 'zst'}:
        raise ValueError(f'Invalid {method=}')

    if method in {'uncompressed', 'gz', 'bz2', 'xz', 'zst'}:
        if method == 'zst' and 'zstandard' not in RESOURCES.optional_packages:
            raise ImportError(f'Package to deal with {method=} not imported, is it installed?')
        return _OPEN.get(method, open)(file, **open_args)

    elif method == 'magic':  # Use the empirical method using magic bytes
        with open(file, 'rb') as f:  # Open the file to read bytes
            first_bytes = f.read(_MIN_N_BYTES)  # Get the bytes necessary to guess the compression type
        method = 'uncompressed'
        for magic, compression in _MAGIC_BYTES.items():
            if first_bytes.startswith(magic):
                method = compression
                break
        return _OPEN.get(method, open)(file, **open_args)

    elif method == 'guess':
        if not isinstance(file, Path):
            file = Path(file)
        return _OPEN.get(file.suffix, open)(file, **open_args)


def decompress(
        buffer: bytes, method: Literal['magic', 'uncompressed', 'gz', 'bz2', 'xz', 'zst'] = 'magic'
) -> bytes:
    """
    Decompresses a buffer of bytes using the appropriate decompression function based on the magic bytes at the
    beginning of the data.

    :param buffer: Buffer of bytes to decompress.
    :param method: Method to determine how to decompress buffer; ``magic`` will determine from the magic bytes at the
                   beginning of the buffer, the other options force the buffer type and should be one of
                   ``uncompressed``, ``gz``, ``bz2``, ``xz`` or ``zst`` (if installed).
    :return: Decompressed buffer of bytes.
    """
    if method not in {'magic', 'uncompressed', 'gz', 'bz2', 'xz', 'zst'}:
        raise ValueError(f'Invalid {method=}')

    if method in {'uncompressed', 'gz', 'bz2', 'xz', 'zst'}:
        if method == 'zst' and 'zstandard' not in RESOURCES.optional_packages:
            raise ImportError(f'Package to deal with {method=} not imported, is it installed?')
        return buffer if method == 'uncompressed' else _DECOMPRESS[method](buffer)

    elif method == 'magic':  # Use the empirical method using magic bytes
        first_bytes = buffer[0:_MIN_N_BYTES]
        method = 'uncompressed'
        for magic, compression in _MAGIC_BYTES.items():
            if first_bytes.startswith(magic):
                method = compression
                break
        return buffer if method == 'uncompressed' else _DECOMPRESS[method](buffer)


def download(url: Union[str, Request], dest: Union[str, Path] = None) -> Union[Path, bytes, None]:
    """
    Downloads a file from a URL to a destination or returns the data as bytes.

    :param url: URL to download from.
    :param dest: Destination path to save the file to. If None, returns the data as bytes.
    :return: Path to the downloaded file or the data as bytes.
    :raises ValueError: If no data is written to the destination file.
    """
    if not isinstance(url, Request):
        url = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    response = urlopen(url)
    if dest:
        with open(dest, mode='wb') as handle:
            handle.write(response.read())
        if (dest := Path(dest)).stat().st_size == 0:
            dest.unlink()
            raise ValueError('No data written')
        return dest
    else:
        return response.read()


def bold(text: str):
    """
    Makes text bold in the terminal.

    :param text: Text to make bold.
    :return: Bold text.
    """
    return f"\033[1m{text}\033[0m"


def get_logo(message: str, width: int = 58) -> str:  # 43 is the width of the logo
    """
    Returns the eris logo with a message centered below it.

    :param message: Message to display below the logo.
    :param width: Width of the logo.
    :return: Formatted logo string.
    """
    return f"\033[1;35m========================|> eris |>========================\n{message.center(width)}\033[0m"
