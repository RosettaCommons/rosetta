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
## @author Phuong T. Nguyen (tranphuonguns@gmail.com)

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "sb_talaris13_loop_modeling_kic_fragments_12res"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE
gen_command_line = '''
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/loop_modeling_kic_fragments_12res/{target}.pdb
-in:file:native {rosetta_dir}/tests/scientific/data/loop_modeling_kic_fragments_12res/reference/{target}.pdb
-parser:script_vars loop={rosetta_dir}/tests/scientific/data/loop_modeling_kic_fragments_12res/{target}.loop fast=0 loop_res=loop_res_holder
-loops:extended true
-loops:frag_sizes 9 3 1
-loops:frag_files {rosetta_dir}/tests/scientific/data/loop_modeling_kic_fragments_12res/{target}.200.9mers {rosetta_dir}/tests/scientific/data/loop_modeling_kic_fragments_12res/{target}.200.3mers none
-ex1
-ex2
-refine_outer_cycles 5
-score:weights talaris2013
-corrections::restore_talaris_behavior true

-nstruct {nstruct}
-parser:protocol {rosetta_dir}/tests/scientific/sfxn_comparison_tests/{testname}/loop_modeling_kic_fragments_12res.xml
-out:file:scorefile {prefix}/{target}.score

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

#==> EDIT HERE
nstruct = 2 if debug else 1000

#==> EDIT HERE
targets = '2cpl 2ebn 2exo 2pia 2rn2 2sil 2tgi'.split()
targets = targets[:2] if debug else targets

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:
    # Read loop file and extract loop residue numbers
    file = open(f'{rosetta_dir}/tests/scientific/data/loop_modeling_kic_fragments_12res/{target}.loop', 'r')
    data = file.read().split()
    # Conventional loop file is Pose numbering, no need to add chain ID
    loop_resnum = data[1] + '-' + data[2]		
    # Replace loop residue in general commandline above with specific loop residue numbers
    command_line = gen_command_line.replace('loop_res_holder', loop_resnum)
    file.close()
    
    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)

    hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
        name=f'{testname}-{target}',

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/rosetta_scripts.{extension}',
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
benchmark.save_variables('debug targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
