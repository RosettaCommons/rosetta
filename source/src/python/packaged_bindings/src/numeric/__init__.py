from __numeric_all_at_once_ import *

import inspect
from types import MethodType

# Setup xyzVector and xyzMatrix pickle support
def __xyzVector__reduce_ex_(self, protocol_level):
    """xyzVector __reduce__ implementation."""

    return ( self.__class__, (self[0], self[1], self[2]) )

def __xyzMatrix__reduce_ex_(self, protocol_level):
    """xyzMatrix __reduce__ implementation."""

    return ( self.__class__, (), (self.col_x(), self.col_y(), self.col_z()) )

def __xyzMatrix__setstate_(self, cols ):
    col_x, col_y, col_z = cols
    self.col_x(col_x)
    self.col_y(col_y)
    self.col_z(col_z)

numeric_classes = inspect.getmembers(__numeric_all_at_once_, inspect.isclass)

for nname, nclass in numeric_classes:
    if nname.startswith("xyzVector"):
        nclass.__reduce_ex__ = MethodType(__xyzVector__reduce_ex_, None, nclass)
    elif nname.startswith("xyzMatrix"):
        nclass.__setstate__ = MethodType(__xyzMatrix__setstate_, None, nclass)
        nclass.__reduce_ex__ = MethodType(__xyzMatrix__reduce_ex_, None, nclass)

try:
    __import__("rosetta.utility")

    __import__("rosetta.numeric.geometry")
    __import__("rosetta.numeric.geometry.hashing")
except ImportError as e:
    import warnings
    warnings.warn("rosetta.numeric import error: %s" % e, ImportWarning)
    raise

xyzMatrix_double = xyzMatrix_Real
xyzVector_double = xyzVector_Real
