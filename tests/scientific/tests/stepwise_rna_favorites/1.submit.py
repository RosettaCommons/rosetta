#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  stepwise_rna_favorites/1.submit.py
## @brief this script is part of stepwise_rna_favorites scientific test
## @author Andy Watkins

import os, sys, time, glob
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "stepwise_rna_favorites"
debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']

hpc_driver = benchmark.hpc_driver()
extension  = benchmark.calculate_extension()

# A target, for us, is a name *and* a set of flags
# Contrast the relax benchmark.
# OK, though, first try is "can we just use wildcards"
command_line = '''
-database {rosetta_dir}/database
-in:file:s {s}
-in:file:native {native}
-in:file:fasta {rosetta_dir}/tests/scientific/data/{testname}/{target}.fasta
-nstruct {nstruct}
-score:weights stepwise/rna/rna_res_level_energy4.wts
-restore_talaris_behavior
-motif_mode
-cycles 400
-save_times
-out:path {prefix}
-out:file:silent {target}.out

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

nstruct = 2 if debug else 5 # 200 split over 40 jobs

targets = 'gagua_pentaloop gcaa_tetraloop gg_mismatch j44a_p4p6 r2_4x4 srl_fixed srl_free_bulgedG srp_domainIV srp_domainIV_fixed tandem_ga_imino tandem_ga_sheared uucg_tetraloop'.split()
#targets = '1POH 1TTZ 2AQD 2DCF 2FKK 2NR7 2OSS 2RK6 3EA6 3ELF 3ESS 3K0M'.split()
targets = targets[:2] if debug else targets

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
    s = " ".join(glob.glob('{rosetta_dir}/tests/scientific/data/{testname}/{target}_HELIX*.pdb'.format_map(vars())))
    native = " ".join(glob.glob('{rosetta_dir}/tests/scientific/data/{testname}/{target}_NATIVE*.pdb'.format_map(vars())))
    jobs = 1 if debug else 40

    #print(command_line)
    #print(command_line.format_map(vars()))
    hpc_job_ids.append( hpc_driver.submit_hpc_job(
        name=f'stepwise_rna_favorites-{target}',
        executable = f'{rosetta_dir}/source/bin/stepwise.{extension}',
        arguments = command_line.format_map(vars()).replace("\n", " "),
        working_dir = prefix,
        jobs_to_queue = jobs,
        log_dir = hpc_logs,
        time=24,
        block=False)
    )


# if not debug:
#     hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
#     time.sleep(64)  # waiting for NFS caching
hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

benchmark.save_variables('targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
