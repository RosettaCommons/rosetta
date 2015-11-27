from __utility_all_at_once_ import *
vector1_double = vector1_Real

from types import MethodType

def _vector_extend(self, other):
    for o in other:
        self.append(o)

for k, v in globals().items():
    if k.startswith("vector1_"):
        v.extend = MethodType(_vector_extend, None, v)

def Vector1(list_in):
    """Creates a Vector1 object, deducing type from the given list."""

    if all([isinstance(x, bool) for x in list_in]):
        t = utility.vector1_bool
    elif all([isinstance(x, int) for x in list_in]):
        t = utility.vector1_int
    elif all([isinstance(x, float) or isinstance(x, int) for x in list_in]):
        t = utility.vector1_double
    elif all([isinstance(x, str) for x in list_in]):
        t = utility.vector1_string
    elif all([isinstance(x, core.id.AtomID) for x in list_in]):
        t = utility.vector1_AtomID
    else:
        raise Exception('Vector1: attemting to create vector of unknow type ' +
                        'or mixed type vector init_list = ' + str(list_in))

    v = t()
    for i in list_in:
        v.append(i)
    return v

def Set(list_in):
    """Creates a Vector1 object, deducing type from the given list."""
    if all([isinstance(x, int) for x in list_in]):
        t = utility.set_int
    elif all([isinstance(x, float) or isinstance(x, int) for x in list_in]):
        t = utility.set_double
    elif all([isinstance(x, str) for x in list_in]):
        t = utility.set_string
    else:
        raise Exception('Set: attemting to create vector of unknow type ' +
                        'or mixed type vector init_list = ' + str(list_in))

    s = t()
    for i in list_in: s.add(i)
    return s


###############################################################################
# Exit callback handler classes.

from .py import PyExitCallback
_python_py_exit_callback = None

class PyRosettaException(Exception):
    def __str__(self):
        return 'PyRosettaException'

class PythonPyExitCallback(PyExitCallback):
    """Exit callback handler, used to shim library-level exit calls."""
    _instance = None

    @classmethod
    def setup(cls):
        cls._instance = cls()
        cls.set_PyExitCallBack(cls._instance)

    def exit_callback(self):
        raise PyRosettaException()

    def __init__(self):
        PyExitCallback.__init__(self)

try:
    __import__("rosetta.utility.excn")
    __import__("rosetta.utility.file")
except ImportError as e:
    import warnings
    warnings.warn("rosetta.utility import error: %s" % e, ImportWarning)
    raise

