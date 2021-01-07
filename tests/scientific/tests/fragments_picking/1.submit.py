#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  fragments_picking/1.submit.py
## @brief this script is part of fragments_picking scientific test
## @author Sergey Lyskov

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "fragments_picking"

debug       = config['debug']
print("debug",debug)
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE
frag_command_line = '''
-database 			{rosetta_dir}/database
-in::file::vall 		{rosetta_dir}/../tools/fragment_tools/vall.jul19.2011.gz
-frags::frag_sizes		 3 9
-frags::describe_fragments	 {prefix}/{target}.fsc
-out::file::frag_prefix 	 {prefix}/{target}-multiR.3w
-frags::scoring::config 	 {rosetta_dir}/tests/scientific/data/{testname}/scoring-multirama.wghts
-mute 				 core.fragment.picking.VallProvider
-frags::ss_pred 		 {rosetta_dir}/tests/scientific/data/{testname}/{target}.psipred_ss2 psipred {rosetta_dir}/tests/scientific/data/{testname}/{target}.porter.ss2 porter {rosetta_dir}/tests/scientific/data/{testname}/{target}.rdb_ss2 rdb
-in:file::native 		 {rosetta_dir}/tests/scientific/data/{testname}/{target}.pdb
-in::file::checkpoint 		 {rosetta_dir}/tests/scientific/data/{testname}/{target}.checkpoint
-in::file::s 			 {rosetta_dir}/tests/scientific/data/{testname}/{target}.pdb
-frags::denied_pdb 		 {rosetta_dir}/tests/scientific/data/{testname}/{target}.homolog
-frags::n_candidates 	 	 {n_candidates}
-frags::n_frags 		 {n_frags}
-frags::picking::quota_config_file	 {rosetta_dir}/tests/scientific/data/{testname}/quota.def
-mute 				 core.conformation
-mute 				 core.chemical
-out:file:scorefile 		{prefix}/{target}.score

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

#==> EDIT HERE
ab_initio_command_line = '''
-database 			{rosetta_dir}/database
-abinitio::rsd_wt_helix 	0.5
-abinitio::rsd_wt_loop		0.5
-abinitio::rg_reweight		0.5
-abinitio::use_filters		false
-abinitio::relax		false
-abinitio::stage2_patch 	{rosetta_dir}/tests/scientific/data/{testname}/crmsd_patch
-abinitio::stage3a_patch 	{rosetta_dir}/tests/scientific/data/{testname}/crmsd_patch
-abinitio::stage3b_patch 	{rosetta_dir}/tests/scientific/data/{testname}/crmsd_patch
-abinitio::stage4_patch  	{rosetta_dir}/tests/scientific/data/{testname}/crmsd_patch
-silent_gz

-nstruct 			{nstruct}
-increase_cycles 		10

-frag3 				{prefix}/{target}-multiR.3w.200.3mers
-frag9 				{prefix}/{target}-multiR.3w.200.9mers

-in:file::native 		 {rosetta_dir}/tests/scientific/data/{testname}/{target}-ideal.pdb
-in::file::s 			 {rosetta_dir}/tests/scientific/data/{testname}/{target}-ideal.pdb

-out:file:silent 		{prefix}/{target}-abinitio.out
-out:sf				 {prefix}/{target}-abinitio.fsc
-out:file:scorefile		 {prefix}/{target}.score

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')
#-parser:protocol {working_dir}/{testname}.xml //nie wiem czy ta linijka im jest potrzebna, na wszelki ją tu zostawiłam

#==> EDIT HERE
n_frags = 200
n_candidates = 1000
nstruct = 2 if debug else 300

#==> EDIT HERE
targets = '1ptqA 1urnA 1tulA 1bkrA 1r69A 1vlsA 1tenA 1bk2A 5croA 2vikA'.split(" ")
targets = targets[:1] if debug else targets

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

    hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
    name=f'{testname}-{target}-frag',

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/fragment_picker.{extension}',
        arguments = frag_command_line.format_map(vars()),
        working_dir = prefix,
        jobs_to_queue = min(nstruct, 50),
        log_dir = hpc_logs,
        time=24,
        block=False)
    )

hpc_driver.wait_until_complete(hpc_job_ids,silent=True)

hpc_job_ids = []
for target in targets:
    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)

    hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
        name=f'{testname}-{target}-ab-initio',

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/minirosetta.{extension}',
        arguments = ab_initio_command_line.format_map(vars()),
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
benchmark.save_variables('debug targets working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
