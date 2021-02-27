#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  RosettaCM/1.submit.py
## @brief this script is part of RosettaCM scientific test
## @author Sergey Lyskov
## @author Jason Fell

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "RosettaCM"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

command_line = '''
-database {rosetta_dir}/database
-in:file:fasta {rosetta_dir}/tests/scientific/data/{testname}/{target}.fasta
-in:file:native {rosetta_dir}/tests/scientific/data/{testname}/{target}_xray.pdb
-nstruct {nstruct}
-parser:protocol {rosetta_dir}/tests/scientific/data/{testname}/{target}_hybridize.xml
-parser:script_vars ROSETTA={rosetta_dir} TEST={testname}
-out:file:scorefile {prefix}/{target}.score
-out:prefix {prefix}/{target}_
-default_max_cycles 200
-dualspace
-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

nstruct = 1 if debug else 200

targets = '2LVC_A 4EPZ_A 4EZI_A 4F0J_A 4FCZ_A 4FDY_A 4FJ6_A 4FM3_A 4FR9_A 4FS7_A 4GHB_A 4GPV_A 4H1X_A 4H41_A 4IC1_A 4JQ6_A'.split()
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

    hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
        name=f'{testname}-{target}',

        executable = f'{rosetta_dir}/source/bin/rosetta_scripts.{extension}',
        arguments = command_line.format_map(vars()),
        working_dir = prefix,
        jobs_to_queue = 50, # The nstruct is of an appropriate size where we do not need to run replicates.
        log_dir = hpc_logs,
        time=24,
        block=False)
    )

# if not debug:
#     hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
#     time.sleep(64)  # waiting for NFS caching
hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

benchmark.save_variables('debug targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
