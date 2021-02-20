#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  simple_cycpep_predict/0.compile.py
## @brief this script is part of the simple_cycpep_predict scientific test.
## @author Sergey Lyskov.
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

import benchmark
import os

from benchmark import to_bytes
from benchmark import _S_passed_, _S_failed_, _S_build_failed_, _S_script_failed_

config = benchmark.config()

rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']

#Create a symlink to data if we don't already have it.
# if not os.path.islink('data'):
#     os.symlink(f'{rosetta_dir}/tests/scientific/data', 'data')


if 'mpi' and 'serialization' not in config['platform']['extras']: benchmark.error(_S_build_failed_, 'This test requires the `mpi` and `serialization` platform.  Aborting...')

res, output, build_command_line = benchmark.build_rosetta()
with open('build-log.txt', 'wb') as f: f.write( to_bytes(output) )

if res: benchmark.error(_S_build_failed_, f'Building rosetta failed!\n{build_command_line}\n{output}\n')

benchmark.save_variables() # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
