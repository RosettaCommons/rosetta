#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  ligand_docking_ddg/1.submit.py
## @brief this script is part of ligand docking ddg scientific test
## @author Sergey Lyskov
## @author Rocco Moretti
## @author Shannon Smith

import os, sys, subprocess, time
import benchmark
import csv

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "ligand_docking_ddg"
debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']

hpc_driver = benchmark.hpc_driver()
extension  = benchmark.calculate_extension()

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

target_file = f'{rosetta_dir}/tests/scientific/data/{testname}/target_ligand_onesthatweknowwork_test.txt'
# extracting protein-ligand pairs from file
data = open(target_file, 'r')
data_reader = csv.reader(data, delimiter=',')
header = next(data_reader)
target_index = header.index("target")
ligid_index = header.index("ligid")
ligand_index = header.index("ligand")

command_line = '''
@ {rosetta_dir}/tests/scientific/tests/{testname}/{testname}.flags
-database {rosetta_dir}/database
-in:file:s "{rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}_cleaned.pdb {rosetta_dir}/tests/scientific/data/{testname}/{target}/{ligand}_aligned.pdb"
-in:file:native "{rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}_cleaned.pdb {rosetta_dir}/tests/scientific/data/{testname}/{target}/{ligand}_aligned.pdb"
-extra_res_fa {rosetta_dir}/tests/scientific/data/{testname}/{target}/{ligand}.params
-extra_res_mol {rosetta_dir}/tests/scientific/data/{testname}/cofactors/*conf_pH.sdf
-parser:protocol {rosetta_dir}/tests/scientific/tests/{testname}/{testname}.xml
-out:file:scorefile {working_dir}/output/{target}/{target}_{ligand}.sc
-no_color
-ex1
-ex2
-ignore_zero_occupancy false
-overwrite
-nstruct 10
'''.replace('\n', ' ').replace('  ', ' ')

nstruct = 1 if debug else 1

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []

targets = []
ligands = []
for row in data_reader:
    targets.append(row[target_index])
    ligands.append(row[ligand_index])

for target,ligand in zip(targets,ligands):
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


benchmark.save_variables('targets ligands rosetta_dir working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)