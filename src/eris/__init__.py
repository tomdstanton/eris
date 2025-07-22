"""
Top-level module, including resource and optional dependency management.
"""
from typing import Callable
from warnings import warn
from functools import wraps
from importlib import import_module
from pathlib import Path


# Classes --------------------------------------------------------------------------------------------------------------
class Resources:
    """
    Holds global resources for this package which are generated on demand.

    Attributes:
        package: Name of the package
        optional_packages: Set of optional packages to check for
    """
    def __init__(self, *optional_packages: str):
        """
        Parameters:
            optional_packages: Optional packages to check for, e.g. 'numpy', 'pandas'
        """
        self.package: str = Path(__file__).parent.name
        self.optional_packages: set[str] = set(filter(self._check_module, optional_packages))
        self._data: 'Traversible' = None  # Generated on demand
        self._metadata: 'PackageMetadata' = None  # Generated on demand
        self._available_cpus: int = None  # Generated on demand
        self._rng: 'Random' = None  # Generated on demand

    @property
    def data(self) -> 'Traversable':
        """Path to the package data"""
        if self._data is None:
            from importlib.resources import files
            self._data = files(self.package) / 'data'
        return self._data

    @property
    def metadata(self) -> 'PackageMetadata':
        """Package metadata"""
        if self._metadata is None:
            from importlib.metadata import metadata
            self._metadata = metadata(self.package)
        return self._metadata

    @property
    def rng(self) -> 'Random':
        """A random number generator instance, can be reused"""
        if self._rng is None:
            from random import Random
            self._rng = Random()
        return self._rng

    @property
    def available_cpus(self):
        """Number of available CPUs"""
        if self._available_cpus is None:
            try:
                from os import process_cpu_count as cpu_count
            except ImportError:
                from os import cpu_count
            self._available_cpus = cpu_count()
        return self._available_cpus

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
    A warning class for this package, making it easy to silence all our warning messages should you wish to.
    Consult the `python.warnings` module documentation for more details.

    Examples:
        >>> import warnings
        >>> from eris import ErisWarning
        ... warnings.simplefilter('ignore', ErisWarning)
    """

    pass


class DependencyWarning(ErisWarning):
    """Custom warning class for missing optional dependencies."""
    pass


def require(*packages: str) -> Callable:
    """
    A decorator to check for required optional packages before executing a function.

    Args:
        *packages: Variable number of package names (strings) that are required.

    Returns:
        A decorator that wraps the function. If any required packages are missing,
        it issues a DependencyWarning and returns None. Otherwise, it executes
        the original function.

    Examples:
        >>> from eris import require
        ... @require('numpy')
        ... def some_numpy_func():
        ... ...
    """
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
