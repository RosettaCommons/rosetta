#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  _template_/1.submit.py
## @brief this script is part of <template> scientific test
## @author Sergey Lyskov

import os, sys, time, subprocess
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "antibody_snugdock"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()


#==> EDIT HERE
command_line = '''
-partners HL_G
-check_cdr_chainbreaks false
-spin
-dock_pert 3 8
-ex1
-ex2aro
-extrachi_cutoff 0
-detect_disulf false
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/targets/intial_orientation/{target}_aligned_model_0_0001.pdb
-ensemble1 {prefix}/{target}_h3.list
-ensemble2 {prefix}/{target}_ag.list
-nstruct {nstruct}
-out:file:scorefile {prefix}/{target}.score
-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

if debug: command_line += "-loops:fast true"

#==> EDIT HERE
nstruct = 2 if debug else 500
njobs = 2 if debug else 500

#==> EDIT HERE
targets = '1ahw 1jps 1mlc 1ztx 2aep 2jel'.split()
print(targets)
#targets = targets[:1] if debug else targets

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

    # note the ensemble lists are not with correct paths, replace with sed!!
    subprocess.getoutput( 'sed \'s|{dir_to_replace}|' + f'{rosetta_dir}/tests/scientific/data/{testname}/targets/ab_ensemble|g\' {rosetta_dir}/tests/scientific/data/{testname}/targets/{target}_h3.list > {prefix}/{target}_h3.list' )
    subprocess.getoutput( 'sed \'s|{dir_to_replace}|' + f'{rosetta_dir}/tests/scientific/data/{testname}/targets/ag_ensemble|g\' {rosetta_dir}/tests/scientific/data/{testname}/targets/{target}_ag.list > {prefix}/{target}_ag.list' )

    hpc_job_ids.append( hpc_driver.submit_hpc_job(
        name=f'{testname}-{target}',

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/snugdock.{extension}',
        arguments = command_line.format_map(vars()),
        working_dir = prefix,
        jobs_to_queue = njobs,
        log_dir = hpc_logs,
        time=24,
        block=False)
    )


# if not debug:
#     hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
#     time.sleep(64)  # waiting for NFS caching
hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

#==> EDIT HERE
benchmark.save_variables('debug targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
