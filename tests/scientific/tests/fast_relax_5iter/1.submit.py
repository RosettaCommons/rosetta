#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/1.submit.py
## @brief this script is part of cartesian_relax scientific test
## @author Sergey Lyskov

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "fast_relax_5iter"
debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']

hpc_driver = benchmark.hpc_driver()
extension  = benchmark.calculate_extension()


command_line = '''
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/{target}.pdb
-in:file:native {rosetta_dir}/tests/scientific/data/{testname}/{target}.pdb
-nstruct {nstruct}
-parser:protocol {working_dir}/{testname}.xml
-out:file:scorefile {prefix}/{target}.score

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')


nstruct = 2 if debug else 100

targets = '1POH 1TTZ 2AQD 2DCF 2FKK 2NR7 2OSS 2RK6 3EA6 3ELF 3ESS 3K0M'.split()
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
        executable = f'{rosetta_dir}/source/bin/rosetta_scripts.{extension}',
        arguments = command_line.format_map(vars()),
        working_dir = prefix,
        jobs_to_queue = min(nstruct, 50),
        log_dir = hpc_logs,
        time=24,
        block=False)
    )


if not debug:
    hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
    time.sleep(64)  # waiting for NFS caching


benchmark.save_variables('targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
