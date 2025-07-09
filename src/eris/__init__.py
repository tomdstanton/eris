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
from typing import Iterable
from warnings import warn
import functools
from importlib import import_module
from importlib.resources import files

# Classes --------------------------------------------------------------------------------------------------------------
class Resources:
    """
    Holds global resources for eris.

    Attributes:
        available_cpus: Number of available CPUs.
        rng: A random number generator instance.
    """
    def __init__(self):
        self.data = files('eris') / 'data'
        try:
            from os import process_cpu_count as cpu_count
        except ImportError:
            from os import cpu_count
        self.available_cpus = cpu_count()
        from random import Random
        self.rng = Random()


class ErisWarning(Warning):
    """
    eris warning.

    eris should use this warning (or subclasses of it), making it easy to
    silence all our warning messages should you wish to:

    >>> import warnings
    >>> from eris import ErisWarning
    ... warnings.simplefilter('ignore', ErisWarning)

    Consult the warnings module documentation for more details.
    """

    pass


class DependencyWarning(ErisWarning):
    """Custom warning class for missing optional dependencies."""

    pass


class Requires:
    """
    Decorator to gracefully handle missing optional dependencies.

    Args:
        requires (set[str]): The name of the dependency that might be missing.
    """

    def __init__(self, *, requires: Iterable[str]):
        self.requires = set(requires)

    def __call__(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if all(filter(DEPENDENCIES.is_available, self.requires)):
                return func(*args, **kwargs)
            else:
                warn(f"The function '{func.__name__}' is unavailable because one of the following dependencies"
                     f"{', '.join(self.requires)} are missing.", DependencyWarning)
                return None

        return wrapper


class Dependencies:
    """
    Manages optional dependencies for the eris library.
    """

    def __init__(self, *packages: str):
        """
        Initializes the Dependencies object.

        Args:
            packages (Iterable[str]): An iterable of package names to check.
        """
        self.available_packages = set(filter(self._check_module, packages))

    def is_available(self, module_name: str) -> bool:
        """
        Checks if a specific module or submodule is available.

        Will check packages found at initialisation followed by a dynamic check.
        """
        return module_name in self.available_packages or self._check_module(module_name)

    @staticmethod
    def _check_module(module_name: str) -> bool:
        """Checks if a module can be imported.

        Args:
            module_name (str): The name of the module to check.

        Returns:
            bool: True if the module can be imported, False otherwise.
        """
        try:
            import_module(module_name)
            return True
        except ImportError:
            return False

# Constants ------------------------------------------------------------------------------------------------------------
RESOURCES = Resources()
DEPENDENCIES = Dependencies('zstandard', 'matplotlib')
