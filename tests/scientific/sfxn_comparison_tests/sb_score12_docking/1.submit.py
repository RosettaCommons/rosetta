#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  docking/1.submit.py
## @brief this script is part of docking scientific test
## @author Sergey Lyskov
## @author Shourya S. Roy Burman

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "sb_score12_docking"
debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']

hpc_driver = benchmark.hpc_driver()
extension  = benchmark.calculate_extension()

def partners(argument):
    switcher = {
        "1AY7": "A_B",
        "1MAH": "A_F",
        "2PCC": "A_B",
        "2SIC": "E_I",
        "1CGI": "E_I",
        "1LFD": "B_A",
        "3EO1": "AB_CF",
        "4IZ7": "A_B",
        "1FQ1": "A_B",
        "2IDO": "A_B",
    }
    return switcher.get(argument, -1)


command_line = '''
-database {rosetta_dir}/database

-in:file:s {rosetta_dir}/tests/scientific/data/docking/targets/{target}.pdb
-unboundrot {rosetta_dir}/tests/scientific/data/docking/targets/{target}.pdb
-in:file:native {rosetta_dir}/tests/scientific/data/docking/natives/{target}.pdb

-nstruct        {nstruct}

-partners       {partners_value}
-dock_pert      3 8
-spin

-docking_low_res_score motif_dock_score
-mh:path:scores_BB_BB {rosetta_dir}/database/additional_protocol_data/motif_dock/xh_16_
-mh:score:use_ss1 false
-mh:score:use_ss2 false
-mh:score:use_aa1 true
-mh:score:use_aa2 true

-score:weights score12
-restore_pre_talaris_2013_behavior true

-detect_disulf true
-rebuild_disulf true
-ignore_zero_occupancy  false
-ex1
-ex2aro

-out:file:scorefile     {prefix}/{target}.score

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')


nstruct = 2 if debug else 5000

targets = '1AY7 1CGI 1FQ1 1LFD 1MAH 2IDO 2PCC 2SIC 3EO1 4IZ7'.split()
targets = targets[:2] if debug else targets

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:
    partners_value = partners(f'{target}')
    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)

    hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
        name=f'docking-{target}',
        executable = f'{rosetta_dir}/source/bin/docking_protocol.{extension}',
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
