#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  _template_/1.submit.py
## @brief this script is part of the mp_domain_assembly scientific test
## @author JKLeman

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "mp_domain_assembly"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE
command_line = '''
-database {rosetta_dir}/database
-in:file:fasta {rosetta_dir}/tests/scientific/data/{testname}/{target}_tr_A.fasta
-in:file:frag3 {rosetta_dir}/tests/scientific/data/{testname}/{target}.frags.3.200_v1_3
-in:file:frag9 {rosetta_dir}/tests/scientific/data/{testname}/{target}.frags.9.200_v1_3
-in:file:native {rosetta_dir}/tests/scientific/data/{testname}/{target}_tr_A_opt.pdb
-mp:setup:spanfiles {rosetta_dir}/tests/scientific/data/{testname}/{target}_tr_A.span
-nstruct {nstruct}
-mp:assembly:TM_pose_number {tm_pose}
-mp:assembly:poses {pose}
-out:path:pdb {prefix}
-out:file:scorefile {prefix}/{target}.score
-multiple_processes_writing_to_one_directory true
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

#==> EDIT HERE
nstruct = 2 if debug else 5000

#==> EDIT HERE
targets = '2Z5X 5HK1 5HK2 5HLB 5XSY'.split()
targets = targets[:2] if debug else targets
tm_poses = '2 1 1 1 2'.split()

poses = []
for i in tm_poses:
    if i == '1':
        poses.append( rosetta_dir+'/tests/scientific/data/'+testname+'/XXX_tr_A_opt_tm.pdb '+rosetta_dir+'/tests/scientific/data/'+testname+'/XXX_tr_A_opt_sol.pdb' )
    else:
        poses.append( rosetta_dir+'/tests/scientific/data/'+testname+'/XXX_tr_A_opt_sol.pdb '+rosetta_dir+'/tests/scientific/data/'+testname+'/XXX_tr_A_opt_tm.pdb' )

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for i in range(0, len(targets)):
    target = targets[i]
    pose = poses[i]
    tm_pose = tm_poses[i]

    prefix = f'{working_dir}/output/{target}'
    print (command_line.format_map(vars()).replace( "XXX", target ))
    if not os.path.exists(prefix): os.makedirs(prefix)

    hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
        name=f'{testname}-{target}',

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/mp_domain_assembly.{extension}',
        arguments = command_line.format_map(vars()).replace( "XXX", target ),
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
