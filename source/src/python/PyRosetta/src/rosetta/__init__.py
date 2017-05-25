"""Pyrosetta.rosetta top-level shim module.

Shim import module exposing the 'pyrosetta.rosetta' extension module as a
top-level 'rosetta' module. This hook imports the 'pyrosetta.rosetta' extension
module then inserts module entries into sys.path for the rosetta namespace.

Not that this technique will not cleanly support 'reload'-ing shim modules,
however as 'pyrosetta.rosetta' is a compiled extension module it does not
support 'reload'-ing in any case. Extension of this mirroring technique to
"reload"-able python extension modules would require implementation of a PEP302
compliant import hook.
"""

import sys
import warnings
warnings.warn(
    "Import of 'rosetta' as a top-level module is deprecated and may be removed in 2018, import via 'pyrosetta.rosetta'.",
    stacklevel=2
)

import pyrosetta.rosetta

pyrosetta_rosetta_modules = ["pyrosetta.rosetta"] + [m for m in sys.modules if m.startswith("pyrosetta.rosetta.")]
remapped_rosetta_modules = { m[len("pyrosetta."):] : sys.modules[m] for m in pyrosetta_rosetta_modules }

sys.modules.update(remapped_rosetta_modules)
