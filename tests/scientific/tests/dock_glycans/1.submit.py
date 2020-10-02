#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license.
# (c) The Rosetta software is developed by the contributing members of the
# (c) Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to
# (c) University of Washington CoMotion, email: license@uw.edu.

## @file   dock_glycans/1.submit.py
## @brief  This script is part of the dock_glycans scientific test.
## @author Sergey Lyskov
## @author Labonte <JWLabonte@jhu.edu>

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "dock_glycans"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#-mute all
command_line = '''
-database {rosetta_dir}/database

-include_sugars
-lock_rings

-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/{target}.pdb
-in:file:native {rosetta_dir}/tests/scientific/data/{testname}/{target}.pdb
-cst_fa_file {rosetta_dir}/tests/scientific/data/{testname}/{target}.cst

-nstruct {nstruct}
-n_cycles 100

-out:file:scorefile {prefix}/{target}.score

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

nstruct = 2 if debug else 500

targets = '1DIW 1GU3 1MFA 1OF4'.split()
targets = targets[:2] if debug else targets

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:
    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)

    hpc_job_ids.append( hpc_driver.submit_hpc_job(
        name=f'{testname}-{target}',

        executable = f'{rosetta_dir}/source/bin/dock_glycans.{extension}',
        arguments = command_line.format_map(vars()),
        working_dir = prefix,
        jobs_to_queue = min(nstruct, 50),
        log_dir = hpc_logs,
        time=24,
        block=False)
    )


# if not debug:
#     hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
#     time.sleep(64)  # waiting for NFS caching
hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

benchmark.save_variables('debug targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
