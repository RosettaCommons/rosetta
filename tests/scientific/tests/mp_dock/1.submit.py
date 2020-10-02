#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  _template_/1.submit.py
## @brief this script is part of mp_dock scientific test
## @author JKLeman

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "mp_dock"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE
command_line = '''
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/{target}_ppk.pdb
-in:file:native {rosetta_dir}/tests/scientific/data/{testname}/{target}.pdb
-mp:setup:spanfiles {rosetta_dir}/tests/scientific/data/{testname}/{target}.span
-out:file:scorefile {prefix}/{target}.score
-out:path:pdb {prefix}
-score:weights mpframework_docking_fa_2015.wts
-mp:scoring:hbond 0
-nstruct {nstruct}
-docking:partners {target[5]}_{target[6]}
-dock_pert 3 8
-restore_talaris_behavior 1
-run:multiple_processes_writing_to_one_directory 1
-packing:pack_missing_sidechains false
-restore_talaris_behavior
-restore_lazaridis_imm_behavior true
-detect_disulf false
'''.replace('\n', ' ').replace('  ', ' ')

#==> EDIT HERE
nstruct = 2 if debug else 1000

#==> EDIT HERE
targets = '1AFO_AB 2K9J_AB 2L35_AB 2NQ2_AB 2QI9_AB 3DIN_CE 3ODU_AB 3RLF_FG 3TUI_AB 4DJH_AB'.split()
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
        executable = f'{rosetta_dir}/source/bin/mp_dock.{extension}',
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

#==> EDIT HERE
benchmark.save_variables('targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
