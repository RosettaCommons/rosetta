#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mp_f19_ddG_of_mutation/1.submit.py
## @brief this script is part of mp_f19_ddG_of_mutation scientific test
## @author Sergey Lyskov, Rebecca F. Alford (ralford3@jhu.edu)

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "mp_f19_ddG_of_mutation"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

nstruct = 1
targets = ['OmpLA', 'OmpLA_aro',  'PagP']
targets = [targets[0]] if debug else targets

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:
    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)

    # Submitting PyRosetta job
    hpc_job_ids.append(
        hpc_driver.submit_serial_hpc_job(
            name=f'{testname}-{target}',

            executable = config["python_virtual_environment"]["python"],
            arguments = f'{rosetta_dir}/tests/scientific/tests/{testname}/predict_ddG.py --energy_fxn franklin2019 --implicit_lipids --repack_radius 8.0 --outdir {prefix} --pdb {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.pdb --spanfile {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.span --mutation_list {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.mutations.dat',
            working_dir = prefix,
            jobs_to_queue = min(nstruct, 50),
            log_dir = hpc_logs,
            time = 24,
            block = False,
            shell_wrapper = True,
        )
    )

hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

benchmark.save_variables('debug targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
