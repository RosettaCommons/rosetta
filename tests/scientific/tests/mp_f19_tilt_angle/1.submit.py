#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  _template_/1.submit.py
## @brief this script is part of mp_f19_energy_landscape scientific test
## @author Sergey Lyskov

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "mp_f19_energy_landscape"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()
nstruct = 1

command_line1 = '''
-database {rosetta_dir}/database
-overwrite
-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/{target}.pdb
-mp:setup:spanfiles {rosetta_dir}/tests/scientific/data/{testname}/{target}.span
-nstruct 1
-parser:protocol {working_dir}/mp_energy_landscape.xml
-parser:script_vars sfxn_weights=franklin2019
-out:file:scorefile {prefix}/{target}.score
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

#targets = '1a11 1mp6 2nr1 1pje WALP23'.split()
targets = '1a11 1mp6 2nr1 1pje WALP23'.split()
targets = targets[:2] if debug else targets

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:

    command_line2 = []
    start_z = ["-60","-40","-20","0","20","40"]
    end_z = ["-40","-20","0","20","40","60"]

    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)


    if((target=='1a11') or (target=='2nr1')):
        print(target)
        command_line2 = command_line1 + '-mp:lipids:composition DLPC -mp:lipids:temperature 20.0'
    elif((target=='1pje') or (target=='WALP23')):
        print(target)
        command_line2 = command_line1 + '-mp:lipids:composition DOPC -mp:lipids:temperature 30.0'
    elif((target=='1mp6')):
        print(target)
        command_line2 = command_line1 + '-mp:lipids:composition DMPC -mp:lipids:temperature 30.0'

    for j in range(len(start_z)):

        command_line = []
        command_line = command_line2 + ' -parser::script_vars start_z=' + start_z[j] + ' -parser::script_vars end_z=' + end_z[j]

        hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
            name=f'{testname}-{target}_{start_z[j]}_to_{end_z[j]}',

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
