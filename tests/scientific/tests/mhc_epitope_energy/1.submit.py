#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mhc_epitope_energy/1.submit.py
## @brief this script is part of mhc_epitope_energy scientific test
## @author Sergey Lyskov
## @author Brahm Yachnin (brahm.yachnin@rutgers.edu)

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "mhc_epitope_energy"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

command_line = '''
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/pdbs/{target}.pdb
-in:file:native {rosetta_dir}/tests/scientific/data/{testname}/pdbs/{target}.pdb
-nstruct {nstruct}
-parser:protocol {working_dir}/{testname}.xml
-parser:script_vars
 aacomp_setup={rosetta_dir}/tests/scientific/data/{testname}/aacompfiles/{target}.comp
 resfile={rosetta_dir}/tests/scientific/data/{testname}/resfiles/{target}_thresh{pssm_thresh}.res
 pdb={rosetta_dir}/tests/scientific/data/{testname}/pdbs/{target}.pdb
-out:file:scorefile {prefix}/{target}.score
-run:preserve_header
-linmem_ig 10
-no_packstat_calculation
-no_optH
-delete_old_poses 1

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

nstruct = 2 if debug else 100
pssm_thresh = 3 if debug else 1 # Note that changing pssm_thresh to anything other than 1 will require changing the cutoffs for all PDBs!

#Leaving out 1k1b_0001, 1yck_0001, 2fc3_0001, and 3mx7_0001, as they behave weirdly wrt netcharge
targets = '1aa2_0001 1dq0_0001 1eer_A_0001 1ewx_0001 1ey4_0001 1f94_0001 1i2t_0001 1jvw_0001 1jyh_0001 1l6p_0001 1ln4_0001 1lwb_0001 1lz1_0001 1mhn_0001 1nkd_0001 1nwa_0001 1pqe_0001 1q5z_0001 1qzn_0001 1r26_0001 1r29_0001 1rc9_0001 1srv_0001 1tp6_0001 1vfq_0001 1zlb_0001 1zlm_0001 1zpw_0001 2a4v_0001 2d7j_0001 2erf_0001 2erw_0001 2f6e_0001 2g5x_0001 2h8e_0001 2icc_0001 2igd_0001 2lis_0001 2oeb_0001 2sak_A_0001 3b02_0001 3cnu_0001 3iu5_0001 3k6y_0001 3lfo_0001 3r2x_C_0001 3rdy_0001 3ri9_0001 3zzp_0001 4eug_0001'.split()
targets = targets[:2] if debug else targets

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:
    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)

    hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
        name=f'{testname}-{target}',

        executable = f'{rosetta_dir}/source/bin/rosetta_scripts.{extension}',
        arguments = command_line.format_map(vars()),
        working_dir = prefix,
        jobs_to_queue = min(nstruct, 50),
        log_dir = hpc_logs,
        time=24,
        block=False)
    )

hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

benchmark.save_variables('debug targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
