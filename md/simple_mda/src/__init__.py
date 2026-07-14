"""
MD analysis utilities package.

Provides a collection of tools for analyzing molecular dynamics trajectories,
including RMSF calculations, PCA, helix analysis, dihedral angles, and distance
measurements.
"""

import importlib
import pkgutil
from pathlib import Path

__all__ = []

# Dynamically import all public modules and expose their public names
for module_info in pkgutil.iter_modules([str(Path(__file__).parent)]):
    module = importlib.import_module(f".{module_info.name}", package=__name__)
    
    # Grab all public names from the module (those not starting with underscore)
    names = getattr(module, "__all__", [
        n for n in vars(module) if not n.startswith("_")
    ])
    
    # Add module members to package namespace
    for name in names:
        globals()[name] = getattr(module, name)
    __all__.extend(names)