#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  simple_cycpep_predict/1.submit.py
## @brief this script is part of simple_cycpep_predict scientific test
## @author Sergey Lyskov.
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#os.system('rm -r data')
#os.system(f'ln -s {rosetta_dir}/tests/scientific/data .')

#==> EDIT HERE
testname    = "glycan_structure_prediction"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

ntasks = 20 if debug else 80
node_time_to_request = 1.5 if debug else 12

#Approximately 3200 CPU-hours for this test

nstruct = 5 if debug else 500

#{rosetta_dir}/tests/scientific/data/{testname}

#==> EDIT HERE
command_line = '''
-ignore_unrecognized_res
-ignore_zero_occupancy false
-load_PDB_components false
-pdb_comments
-out:pdb_gz
-scorefile_format json
-ex1
-ex2
-use_input_sc
-ideal_sugars
-jd2:delete_old_poses
-other_pose_to_scorefile
-input_ab_scheme AHo_Scheme
-output_ab_scheme AHo_Scheme
-include_sugars
-auto_detect_glycan_connections
-alternate_3_letter_codes pdb_sugar
-maintain_links
-write_pdb_link_records
-write_glycan_pdb_codes
-mpi_fraction_outputters .05
-skip_connect_info
-ignore_unrecognized_res
-edensity::score_symm_complex false
-cryst::crystal_refine
-job_definition_file {working_dir}/default_substituted.xml
-nstruct {nstruct}
-database {rosetta_dir}/database
-no_fconfig
-out:path:all decoys
-no_color
'''.replace('\n', ' ').replace('  ', ' ')
#-nblist_autoupdate

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
if not os.path.exists(hpc_logs+"/test"): os.system("touch " + hpc_logs + "/test")
hpc_job_ids = []

# create decoys directory
decoys_dir = f'{working_dir}/decoys'
if not os.path.exists(decoys_dir): os.makedirs(decoys_dir)

hpc_job_ids.append( hpc_driver.submit_mpi_hpc_job(
    name=testname,
    executable = f'{rosetta_dir}/source/bin/rosetta_scripts_jd3.{extension}',
    arguments = command_line.format_map(vars()),
    working_dir = working_dir,
    log_dir = hpc_logs,
    time=node_time_to_request,
    ntasks=ntasks,
    block=False)
)

hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
time.sleep(64)  # waiting for NFS caching

benchmark.save_variables('debug nstruct working_dir rosetta_dir debug testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
