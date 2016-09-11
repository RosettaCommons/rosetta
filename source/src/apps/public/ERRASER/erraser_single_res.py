#!/usr/bin/env phenix.python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os.path
import imp

file_path = os.path.split( os.path.abspath(__file__) ) [0]
imp.load_source('erraser_wrapper', file_path + '/erraser_wrapper.py')
imp.load_source('erraser_option', file_path + '/erraser_option.py')
imp.load_source('erraser_util', file_path + '/erraser_util.py')

from erraser_wrapper import erraser_single_res
from erraser_option import erraser_option
from erraser_util import *

option = erraser_option()
option.read_cmdline_erraser_single_res( sys.argv )
option.finalize()
erraser_single_res( option )
