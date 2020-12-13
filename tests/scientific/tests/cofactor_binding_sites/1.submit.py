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

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "cofactor_binding_sites"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE
# -nstruct {nstruct}
command_line = '''
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/input/{directory}/{pdb}.pdb
-nstruct 1
-mute protocols.backrub.BackrubMover
-multiple_processes_writing_to_one_directory
-no_color
-ex1
-ex2
-extrachi_cutoff 0
-extra_res_fa {rosetta_dir}/tests/scientific/data/{testname}/input/{directory}/{lig}.params
-resfile {rosetta_dir}/tests/scientific/data/{testname}/input/{directory}/{pdb}.resfile
-coupled_moves::mc_kt 2.4
-coupled_moves::boltzmann_kt 2.4
-coupled_moves::ntrials {ntrials}
-coupled_moves::initial_repack false
-coupled_moves::ligand_mode true
-coupled_moves::ligand_weight 2
-coupled_moves::fix_backbone false
-coupled_moves::bias_sampling true
-coupled_moves::bump_check true
-coupled_moves::backbone_mover backrub
-coupled_moves::exclude_nonclashing_positions true
-out:prefix {prefix}{id4}_
-out:suffix
'''.replace('\n', ' ').replace('  ', ' ')

#==> EDIT HERE
nstruct = 1 if debug else 200
ntrials = 100 if debug else 1000

#==> EDIT HERE
targets = open( f'{rosetta_dir}/tests/scientific/data/{testname}/input/systems.txt' ).readlines()
#targets = targets[:2] if debug else targets


#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []

prefix = f'{working_dir}/output/'
if not os.path.exists(prefix): os.makedirs(prefix)

for target in targets:
    pdb = target.split()[0]
    lig = target.split()[1]
    directory = target.split()[2]

    for n in range( 0, nstruct ):
        id4 = str(n+1).zfill(4)
        hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
            name=f'{testname}-{pdb}-{lig}',
            executable = f'{rosetta_dir}/source/bin/coupled_moves.{extension}',
            arguments = command_line.format_map(vars()),
            working_dir = prefix,
            jobs_to_queue = 1,
            log_dir = hpc_logs,
            time=1,
            block=False)
        )

# if not debug:
#     hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
#     time.sleep(64)  # waiting for NFS caching
hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

#==> EDIT HERE
benchmark.save_variables('targets nstruct working_dir testname rosetta_dir')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
