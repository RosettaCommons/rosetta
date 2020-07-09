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

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "simple_cycpep_predict"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

run_time = 900 if debug else 9000
nodes_to_request = 1 if debug else 4
node_time_to_request = 1.5 if debug else 4

#==> EDIT HERE
command_line = '''
-database {rosetta_dir}/database
-in:file:native {working_dir}/native.pdb
-ex1
-ex2
-MPI_pnear_lambda 1.25
-MPI_pnear_kbt 0.62
-nstruct 1000000000
-cyclic_peptide:MPI_batchsize_by_level 1
-cyclic_peptide:MPI_auto_2level_distribution true
-out:file:silent out.silent
-multithreading:total_threads 1
-threads_per_worker 1
-cyclic_peptide:compute_rmsd_to_lowest true
-cyclic_peptide:compute_ensemble_sasa_metrics true
-cyclic_peptide:MPI_stop_after_time {run_time}
-out:file:silent out.silent
-cyclic_peptide:sample_cis_pro_frequency 0.3
-cyclic_peptide:MPI_output_fraction 0.001
-cyclic_peptide:sequence_file {working_dir}/seq.txt
-score:symmetric_gly_tables true
-cyclic_peptide:genkic_closure_attempts 250
-cyclic_peptide:genkic_min_solution_count 1
-cyclic_peptide:use_rama_filter true
-cyclic_peptide:rama_cutoff 3.0
-cyclic_peptide:min_genkic_hbonds 2
-cyclic_peptide:min_final_hbonds 2
-mute all
-unmute protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication_MPI_summary protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication_MPI 
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []

prefix = f'{working_dir}/output/'
if not os.path.exists(prefix): os.makedirs(prefix)

hpc_job_ids.append( hpc_driver.submit_mpi_hpc_job(
    name=f'simple_cycpep_predict',
    executable = f'{rosetta_dir}/source/bin/simple_cycpep_predict.{extension}',
    arguments = command_line.format_map(vars()),
    working_dir = prefix,
    log_dir = hpc_logs,
    time=node_time_to_request,
    requested_nodes=nodes_to_request,
    requested_processes_per_node=20,
    block=False)
)

hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
time.sleep(64)  # waiting for NFS caching

benchmark.save_variables('debug nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
