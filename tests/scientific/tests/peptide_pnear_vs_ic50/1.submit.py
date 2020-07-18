#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  peptide_pnear_vs_ic50/1.submit.py
## @brief this script is part of peptide_pnear_vs_ic50 scientific test
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "peptide_pnear_vs_ic50"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

maybe_beta  = '--beta' if 'beta' in config['variants'] else ''

samples = 2000 if debug else 80000
nodes_to_request = 1 if debug else 2
node_time_to_request = 1.0 if debug else 4.0

#==> EDIT HERE
command_line = '''
-database {rosetta_dir}/database
-ex1
-ex2
-MPI_pnear_lambda 1.5
-MPI_pnear_kbt 0.62
-nstruct {samples}
-cyclic_peptide:MPI_batchsize_by_level 1
-cyclic_peptide:MPI_auto_2level_distribution true
-multithreading:total_threads 1
-threads_per_worker 1
-cyclic_peptide:compute_rmsd_to_lowest true
-cyclic_peptide:compute_ensemble_sasa_metrics true
-cyclic_peptide:sample_cis_pro_frequency 0.3
-cyclic_peptide:MPI_output_fraction 0.001
-score:symmetric_gly_tables true
-cyclic_peptide:genkic_closure_attempts 250
-cyclic_peptide:genkic_min_solution_count 1
-cyclic_peptide:use_rama_filter true
-cyclic_peptide:rama_cutoff 3.0
-cyclic_peptide:min_genkic_hbonds 2
-cyclic_peptide:min_final_hbonds 2
{maybe_beta}
-mute all
-unmute protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication_MPI_summary protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication_MPI
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

command_line_additions = '''
-in:file:native {working_dir}/inputs/NDM1i_XXXXX.pdb
-cyclic_peptide:sequence_file {working_dir}/inputs/NDM1i_XXXXX_seq.txt
-out:file:silent out_NDM1i_XXXXX.silent
'''.replace('\n', ' ').replace('  ', ' ')


#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

suffices = ['1A','1B','1C','1D','1E','1F','1G']

for suffix in suffices:
    hpc_logs = f'{working_dir}/hpc-logs_' + suffix
    if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []

prefix = f'{working_dir}/output/'
if not os.path.exists(prefix): os.makedirs(prefix)

for suffix in suffices:
    hpc_job_ids.append( hpc_driver.submit_mpi_hpc_job(
        name=f'peptide_pnear_vs_ic50_NDM1i_' + suffix,
        executable = f'{rosetta_dir}/source/bin/simple_cycpep_predict.{extension}',
        arguments = (command_line + " " + command_line_additions.replace('XXXXX', suffix) ).format_map(vars()),
        working_dir = prefix,
        log_dir = f'{working_dir}/hpc-logs_' + suffix,
        time=node_time_to_request,
        requested_nodes=nodes_to_request,
        requested_processes_per_node=20,
        block=False)
    )

hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
time.sleep(64)  # waiting for NFS caching

benchmark.save_variables('debug samples working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
