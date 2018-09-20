#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  enzyme_design/1.submit.py
## @brief this script is part of enzyme_design scientific test
## @author Sergey Lyskov
## @author Rocco Moretti

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "enzyme_design"
debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']

hpc_driver = benchmark.hpc_driver()
extension  = benchmark.calculate_extension()


command_line = '''
@ {working_dir}/ENZDES.flags
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/{target}_{target}.pdb
-in:file:native {rosetta_dir}/tests/scientific/data/{testname}/{target}_{target}.pdb
-extra_res_fa {rosetta_dir}/tests/scientific/data/{testname}/{target}.params
-nstruct {nstruct}
-parser:protocol {working_dir}/{testname}.xml
-parser::script_vars PSSM={rosetta_dir}/tests/scientific/data/{testname}/{target}.pssm

-scorefile_format json
-out:file:scorefile {target}.json
-no_color
'''.replace('\n', ' ').replace('  ', ' ')


nstruct = 1 if debug else 1


targets = '1A99 1db1 1fby 1FZQ 1H6H 1hmr 1hsl 1J6Z 1l8b 1LKE 1n4h 1nl5 1nq7 1OPB 1POT 1RBP 1sw1 1TYR 1urg 1USK 1uw1 1wdn 1x7r 1XT8 1XZX 1y2u 1y3n 1y52 1z17 1ZHX 2b3b 2DRI 2e2r 2f5t 2FME 2FQX 2FR3 2GM1 2h6b 2HZQ 2ifb 2ioy 2p0d 2PFY 2Q2Y 2Q89 2qo4 2rct 2RDE 2UYI 3B50'.split()
targets = targets[:1] if debug else targets

metrics = "seqrec pssm_seqrec pssm_delta_seqrec".split()

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

    hpc_job_ids.append( hpc_driver.submit_hpc_job(
        name=f'{testname}-{target}',
        executable = f'{rosetta_dir}/source/bin/rosetta_scripts.{extension}',
        arguments = command_line.format_map(vars()),
        working_dir = prefix,
        jobs_to_queue = 1, # nstruct is low enough that we don't need to do replicates
        log_dir = hpc_logs,
        time=24,
        block=False)
    )

if not debug:
    hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
    time.sleep(64)  # waiting for NFS caching


benchmark.save_variables('targets metrics nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
