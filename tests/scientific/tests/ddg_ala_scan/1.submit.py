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
import json

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "ddg_ala_scan"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE

command_line = '''
-database {rosetta_dir}/database
-out:file:scorefile {prefix}/out.score
-parser:view -inout:dbms:mode sqlite3 
-inout:dbms:database_name {prefix}/out.db3 -no_optH true -restore_talaris_behavior
-multiple_processes_writing_to_one_directory
-no_color

-in:file:s {target_pdb}
-parser:protocol {working_dir}/{testname}.xml
-parser:script_vars {script_vars}
'''.replace('\n', ' ').replace('  ', ' ')



#==> EDIT HERE
nstruct = 1
with open(working_dir+'/job_dict.json','r') as f:
    targets_dict = json.load(f)
targets = list(targets_dict.keys())
targets_command_line = list(targets_dict.values())

if debug:
    targets = targets[:1] 
    targets_command_line = targets_command_line[:1]

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for n, target in enumerate(targets):
    target = target.replace('/','__')
    prefix = f'{working_dir}/output/{target}'
    print('PREFIX:')
    print(prefix)
    if not os.path.exists(prefix): 
        os.makedirs(prefix)
        print("MADE DIRECTORY: ",prefix)

    target_pdb = targets_command_line[n]["-in:file:s"]
    script_vars = " ".join(targets_command_line[n]["-parser:script_vars"])
    #format_map must be called at least twice, since script_vars contain replacment tags as well
    cmd_line = command_line.format_map(vars())
    cmd_line = cmd_line.format_map(vars())
    print('COMMAND LINE IS:')
    print(cmd_line)

    hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
        name=f'{testname}-{target}',

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/rosetta_scripts.{extension}',
        arguments = cmd_line,
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
