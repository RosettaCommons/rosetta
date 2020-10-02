#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/2.analyze.py
## @brief this script is part of cartesian_relax scientific test
## @author Sergey Lyskov

import os, sys, subprocess
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

interface_analysis_flags = '''
-l {prefix}/models.list
-native {rosetta_dir}/tests/scientific/data/{testname}/natives/{target}.pdb
-partners HL_G
-docking_local_refine
-dock_min
-out:file:score_only {prefix}/{target}_rmsd.sc
-multiple_processes_writing_to_one_directory
-nstruct 1
'''.replace('\n', ' ').replace('  ', ' ')

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:
    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)

    # get models for analysis
    subprocess.getoutput( f'ls {prefix}/*.pdb > {prefix}/models.list' )

    hpc_job_ids.append( hpc_driver.submit_hpc_job(
        name=f'{testname}-1b-{target}',

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/docking_protocol.{extension}',
        arguments = interface_analysis_flags.format_map(vars()),
        working_dir = prefix,
        jobs_to_queue = len(targets),
        log_dir = hpc_logs,
        time=24,
        block=False)
    )


# if not debug:
#     hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
#     time.sleep(64)  # waiting for NFS caching
hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

benchmark.save_variables('debug targets nstruct working_dir testname')
