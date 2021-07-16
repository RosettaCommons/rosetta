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

## @file   glycan_dock/1.submit.py
## @brief  This script is part of the GlycanDock scientific test.
## @author Sergey Lyskov
## @author Morgan Nance (@mlnance; revised for GlycanDock sci test)

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "glycan_dock"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

# GlycanDock was benchmarked using ref2015 ScoreFunction
# and, within GlycanDockProtocol.cc, the following weights were set
# sugar_bb = 0.5
# fa_intra_rep_nonprotein = weight of fa_rep
# so the flags here will only specify the ref2015 sf
# but will rely on GlycanDockProtocol.cc to set weights as needed
command_line = '''
-database {rosetta_dir}/database

-include_sugars
-lock_rings

-score:weights ref2015

-in:file:s {rosetta_dir}/tests/scientific/tests/{testname}/inputs/{target}-input.pdb
-docking:partners A_X
-cst_fa_file {rosetta_dir}/tests/scientific/tests/{testname}/inputs/{target}.cst
-in:file:native {rosetta_dir}/tests/scientific/tests/{testname}/inputs/{target}.pdb

-nstruct {nstruct}
-n_cycles 10

-out:file:scorefile {prefix}/{target}.score
-out:pdb_gz

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

# minimum of nstruct = 5 because using 5-top-scoring for bootstrapping
nstruct = 1000
if debug: nstruct = 5

targets = '1UXX 2RDK 2J1U 5OYE 5ZHO 6R3M'.split()
if debug: targets[:2]

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

        executable = f'{rosetta_dir}/source/bin/GlycanDock.{extension}',
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
