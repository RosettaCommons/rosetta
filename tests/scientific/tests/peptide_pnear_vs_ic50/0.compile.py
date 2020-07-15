#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  peptide_pnear_vs_ic50/0.compile.py
## @brief this script is part of the peptide_pnear_vs_ic50 scientific test.
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).

import benchmark

from benchmark import to_bytes
from benchmark import _S_passed_, _S_failed_, _S_build_failed_, _S_script_failed_

config = benchmark.config()
if {'mpi', 'cxx11thread', 'serialization'}.issubset(config['platform']['extras']) == False :
    benchmark.error(_S_build_failed_, 'This test requires the `mpi`, `cxx11thread` and `serialization` platforms.  Aborting...')

res, output, build_command_line = benchmark.build_rosetta()
with open('build-log.txt', 'wb') as f: f.write( to_bytes(output) )

if res: benchmark.error(_S_build_failed_, f'Building rosetta failed!\n{build_command_line}\n{output}\n')

benchmark.save_variables() # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
