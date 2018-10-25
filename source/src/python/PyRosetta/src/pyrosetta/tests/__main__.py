import sys
try:
    import faulthandler
    faulthandler.enable()
except ImportError:
    import warnings
    warnings.warn("Unable to import faulthandler, install for compiled-module traceback support.")
    pass

import pyrosetta.tests

result = pyrosetta.tests.test()
sys.exit(not result.wasSuccessful())
