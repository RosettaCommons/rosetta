#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mp_f19_decoy_discrimination/1.submit.py
## @brief this script is part of the franklin19 decoy discrimination scientific test
## @author Sergey Lyskov, Rebecca F. Alford (ralford3@jhu.edu)

import os, sys, time, subprocess
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "mp_f19_decoy_discrimination"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

nstruct = 1
targets = 'vatp brd7 fmr5 rhod'.split(" ")
targets = [ targets[0] ] if debug else targets

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []

for target in targets: 

    # Write a custom path dependent list of decoys
    targetdir = f'{rosetta_dir}/tests/scientific/data/{testname}/{target}'
    decoys = subprocess.check_output(f'ls {targetdir}/decoy.*.pdb', shell=True).decode("utf-8").split('\n')
    decoylist_file = f"{working_dir}/{target}.decoy.list"
    decoylist = []
    with open( decoylist_file, 'w' ) as f: 
        if debug: 
            for d in decoys[0:5]: 
                f.write( d + "\n" )
                decoylist.append( d )
        else: 
            for d in decoys: 
                f.write( d + "\n" )
                decoylist.append( d )

    print(decoylist)
    # Submit an individual job to refine each decoy
    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)
    for d in decoylist: 

        command_line = '''
        -database {rosetta_dir}/database
        -in:file:s {d}
        -in:file:native {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}_native.pdb
        -mp:setup:spanfiles {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.span
        -nstruct {nstruct}
        -parser:protocol {working_dir}/membrane_relax.xml
        -parser:script_vars sfxn_weights=franklin2019
        -out:file:scorefile {prefix}/{target}.score
        -multiple_processes_writing_to_one_directory
        -constrain_relax_to_start_coords true
        -no_color
        '''.replace('\n', ' ').replace('  ', ' ')

  
        hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
            name=f'{testname}-{target}',

            executable = f'{rosetta_dir}/source/bin/rosetta_scripts.{extension}',
            arguments = command_line.format_map(vars()),
            working_dir = prefix,
            jobs_to_queue = min(nstruct, 50),
            log_dir = hpc_logs,
            time=24,
            block=False)
        )

hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

benchmark.save_variables('debug targets nstruct extension rosetta_dir working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
