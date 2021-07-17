#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  abinitio_RosettaNMR_rdc/1.submit.py
## @brief this script is part of abinitio_RosettaNMR_rdc scientific test
## @author Sergey Lyskov

import os, sys, time, subprocess
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "abinitio_RosettaNMR_rdc"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE
command_line = '''
-database {rosetta_dir}/database
-in:file:native {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.pdb
-in:file:fasta {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.fasta
-in:file:frag3 {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.cs.200.3mers
-in:file:frag9 {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.cs.200.9mers
-out:file:silent {prefix}/{target}.denovo.out
-out:file:silent_struct_type binary
-out:file:scorefile {prefix}/{target}.score
-out:nstruct {nstruct}

-chemical:patch_selectors CENTROID_HA
-abinitio:increase_cycles 1
-abinitio:rg_reweight 0.5
-abinitio:rsd_wt_helix 0.5
-abinitio:rsd_wt_loop 0.5
-abinitio:stage1_patch {rosetta_dir}/tests/scientific/data/{testname}/{target}/rdc.wts_patch
-abinitio:stage2_patch {rosetta_dir}/tests/scientific/data/{testname}/{target}/rdc.wts_patch
-abinitio:stage3a_patch {rosetta_dir}/tests/scientific/data/{testname}/{target}/rdc.wts_patch
-abinitio:stage3b_patch {rosetta_dir}/tests/scientific/data/{testname}/{target}/rdc.wts_patch
-abinitio:stage4_patch {rosetta_dir}/tests/scientific/data/{testname}/{target}/rdc.wts_patch
-relax:fast

-broker:setup {working_dir}/{target}.tbp
-run:protocol broker
-reinitialize_mover_for_each_job
-score:find_neighbors_3dgrid

-evaluation:rmsd NATIVE 2native FULL
-nmr:rdc:input_file {working_dir}/{target}.rdc.inp
-nmr:rdc:multiset_weights 1.0
-nmr:rdc:normalization_type none
-nmr:rdc:correct_sign false

-out:mute core.pose.nmr.rdc.RDCSingle

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

#==> EDIT HERE
nstruct = 2 if debug else 2000

#==> EDIT HERE
targets = '2a7o 2k5u 2klc'.split()
targets = targets[:2] if debug else targets

#==========================
# write input file with proper path
# nothing else works unfortunately
def write_input_file( file2read, file2write, path2add ):
    
    # read file
    lines = (subprocess.getoutput("cat " + file2read)).splitlines()
    
    #file will be written in working directory
    f = open( file2write, "w" )
    
    for l in lines:
         
        # replace PATH with actual path
        l = l.replace( "<PATH>", path2add )
        f.write( l + "\n" )
    f.close()
    
    return

#==========================

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:
    
    # write broker file with current paths
    file2read = f'{rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.tbp'
    file2write = f'{working_dir}/{target}.tbp'
    path2add = f'{rosetta_dir}/tests/scientific/data/{testname}/{target}'
    write_input_file( file2read, file2write, path2add )

    # write rdc input file with current paths
    file2read = f'{rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.rdc.inp'
    file2write = f'{working_dir}/{target}.rdc.inp'
    path2add = f'{rosetta_dir}/tests/scientific/data/{testname}/{target}'
    write_input_file( file2read, file2write, path2add )

    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)

    hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
        name=f'{testname}-{target}',

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/minirosetta.{extension}',
        arguments = command_line.format_map(vars()),
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
benchmark.save_variables('debug targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
