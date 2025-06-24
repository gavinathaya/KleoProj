"""
KleoPy: Custom tools for 216-Kleopatra using the Dipole-Segment model.
===========

Submodules
-----------

constants                     : Physical constants and parameters for 216-Kleopatra
grid_search                   : Numerical methods for finding symmetric periodic orbits
integrator                    : Numerical integrators for solving ODEs
orbital_eq                    : Equations related to the orbital dynamics of 216-Kleopatra
"""
    # kleopy
    #     ├── constants
    #     ├── grid_search
    #     ├── integrator
    #     └── orbital_eq 

import logging

#Package version
__version__ = "0.1.0" 

#Package-level logger
logger = logging.getLogger("kleopy")
logger.addHandler(logging.NullHandler())

# Import key classes/functions for easy access
# from .module1 import MainClass, useful_function
# from .module2 import AnotherClass

submodules = [
    "constants",
    "grid_search",
    "integrator",
    "orbital_eq",
]

__all__ = submodules + [
    "__version__",
    "logger",
]

#Package-level initialization code
def show_version():
    """Print the installed version of kleopy."""
    print(f"kleopy version {__version__}")
