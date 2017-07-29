from __future__ import print_function

import pyrosetta
import pyrosetta.rosetta as rosetta
from pyrosetta import *

pyrosetta.init(extra_options = "-constant_seed -nstruct 5")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

print( pyrosetta.version() )

print( 'out:nstruct: {}'.format( rosetta.basic.options.get_integer_option('out:nstruct') ) )
assert rosetta.basic.options.get_integer_option('out:nstruct') == 5
