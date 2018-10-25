import os
from glob import glob

import pickle
import numpy

_root = os.path.join(os.path.dirname(__file__))

if os.path.exists(_root + "/source_structure_residues.npy"):
    source_structure_residues = numpy.load(_root + "/source_structure_residues.npy")
else:
    source_structure_residues = None
