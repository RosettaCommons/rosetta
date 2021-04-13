#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  ligand_docking/1.submit.py
## @brief this script is part of ligand docking scientific test
## @author Sergey Lyskov
## @author Rocco Moretti
## @author Shannon Smith

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "sb_ligand_docking"
debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']

hpc_driver = benchmark.hpc_driver()
extension  = benchmark.calculate_extension()

sfxns = ['ligand', 'talaris', 'ref2015', 'betanov16']
for sfxn in sfxns:
    command_line = '''
    @ {rosetta_dir}/tests/scientific/data/ligand_docking/flags_{sfxn}
    -database {rosetta_dir}/database
    -in:file:s {rosetta_dir}/tests/scientific/data/ligand_docking/{target}_input.pdb
    -in:file:native {rosetta_dir}/tests/scientific/data/ligand_docking/{target}_native.pdb
    -extra_res_fa {rosetta_dir}/tests/scientific/data/ligand_docking/{target}_ligand.params
    -parser:protocol {rosetta_dir}/tests/scientific/data/ligand_docking/ligand_docking_{sfxn}.xml
    -parser::script_vars startfrom={rosetta_dir}/tests/scientific/data/ligand_docking/{target}_startfrom.pdb
    -out:file:scorefile {working_dir}/output/{sfxn}/{target}/{target}.sc
    -no_color
    -ex1
    -ex2
    -ignore_zero_occupancy false
    -overwrite
    -nstruct {nstruct}
    '''.replace('\n', ' ').replace('  ', ' ')

#    -out:pdb_gz

    nstruct = 1 if debug else 1000
    targets = '3FKL 4QMQ 2VTS 4MCC 3U5K 4BQH 5D7C 3R5N 5AOJ 3K5V 4CHN 5FQV 5D7P 4CJF 5AM0 4A9N 4XT2 4CCU 4GB2 3ZXZ 3PE2 1EYQ 4UWC 5CTU 4GV0 3FLY 3BLL 2HK5 4F1L 4FA2 4J93 3P0Q 4W7T 1ZZL 3TLL 4KTN 2ZC9 4FPJ 4WUA 3T96 3UWK 4I4F 2ZDT 2WXG 4UAL 3WYP 3ZO4 4ZBQ 2BR1 3VRI'.split()
    targets = targets[:2] if debug else targets

    hpc_logs = f'{working_dir}/hpc-logs'
    if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
    hpc_job_ids = []
    for target in targets:
        prefix = f'{working_dir}/output/{sfxn}/{target}'
        if not os.path.exists(prefix): os.makedirs(prefix)

        hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
            name=f'{testname}-{sfxn}-{target}',
            executable = f'{rosetta_dir}/source/bin/rosetta_scripts.{extension}',
            arguments = command_line.format_map(vars()),
            working_dir = prefix,
            jobs_to_queue = 1, # nstruct is low enough that we don't need to do replicates
            log_dir = hpc_logs,
            time=24,
            block=False)
        )

#    if not debug:
#        hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
#        time.sleep(64)  # waiting for NFS caching
hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

benchmark.save_variables('sfxns targets working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
