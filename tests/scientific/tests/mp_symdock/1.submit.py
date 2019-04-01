#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mp_symdock/1.submit.py
## @brief this script is part of the mp_symdock scientific test
## @author JKLeman

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "mp_symdock"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE
command_line = '''
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/{target}_opm__0001_INPUT.pdb
-in:file:native {rosetta_dir}/tests/scientific/data/{testname}/{target}_opm__0001.pdb
-mp:setup:spanfiles {rosetta_dir}/tests/scientific/data/{testname}/{target}_opm_A.span
-mp:scoring:hbond
-nstruct {nstruct}
-symmetry:symmetry_definition {rosetta_dir}/tests/scientific/data/{testname}/{target}.sym
-symmetry:initialize_rigid_body_dofs
-packing:pack_missing_sidechains 0
-docking:dock_lowres_filter 10 20
-out:file:scorefile {prefix}/{target}.score
-score:weights mpframework_docking_fa_2015.wts
-docking::dock_ppk true
-multiple_processes_writing_to_one_directory true
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

#==> EDIT HERE
nstruct = 2 if debug else 1000

#==> EDIT HERE
targets = '1AFO 1BL8 2MPN 2OAR 2UUH'.split()
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

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/mp_symdock.{extension}',
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

#==> EDIT HERE
benchmark.save_variables('targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
