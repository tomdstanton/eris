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
from typing import Callable
from warnings import warn
from functools import wraps
from importlib import import_module
from importlib.resources import files
from importlib.metadata import metadata
from random import Random
from pathlib import Path
try:
    from os import process_cpu_count as cpu_count
except ImportError:
    from os import cpu_count

# Classes --------------------------------------------------------------------------------------------------------------
class Resources:
    """
    Holds global resources for eris.

    Attributes:
        data: Path
        available_cpus: Number of available CPUs.
        rng: A random number generator instance.
    """
    def __init__(self, *optional_packages: str):
        self.package = Path(__file__).parent.name
        self.data: Path = files(self.package) / 'data'
        self.metadata = metadata(self.package)
        self.available_cpus: int = cpu_count()
        self.rng: Random = Random()
        self.optional_packages: set[str] = set(filter(self._check_module, optional_packages))

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


def require(*packages: str) -> Callable:
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            if missing_deps := [dep for dep in packages if dep not in RESOURCES.optional_packages]:
                warn(
                    f"Function '{func.__name__}' requires the following missing dependencies: "
                    f"{', '.join(missing_deps)}. Skipping execution.",
                    DependencyWarning
                )
                return None
            return func(*args, **kwargs)
        return wrapper
    return decorator


# Constants ------------------------------------------------------------------------------------------------------------
RESOURCES = Resources('pyrodigal', 'matplotlib')
