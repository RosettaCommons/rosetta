#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  ligand_scoring_ranking/1.submit.py
## @brief this script is part of ligand scientific test
## @author Sergey Lyskov
## @author Shannon Smith (shannon.t.smith.1@vanderbilt.edu)

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "ligand_scoring_ranking"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE
command_line = '''
@ {rosetta_dir}/tests/scientific/data/{testname}/flags_ligand_scoring_ranking
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/{target}_protein_ligand.pdb
-in:file:native {rosetta_dir}/tests/scientific/data/{testname}/{target}_protein_ligand.pdb
-in:file:extra_res_fa {rosetta_dir}/tests/scientific/data/{testname}/{target}_ligand.params
-nstruct {nstruct}
-parser:protocol {working_dir}/{testname}.xml
-out:file:scorefile {prefix}/{target}.score

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

nstruct = 10

targets = '4m0z 4m0y 3qgy 4qd6 4rfm 4llx 5c28 3uuo 3ui7 5c2h 2v00 3wz8 3pww 3prs 3uri 4cr9 4cra 4x6p 4crc 4ty7 5aba 5a7b 4agn 4agp 4agq 3bgz 3jya 2c3i 4k18 5dwr 3mss 3k5v 3pyy 2v7a 4twp 3wtj 3zdg 3u8k 4qac 3u8n 1a30 2qnq 1g2k 1eby 3o9i 4lzs 3u5j 4wiv 4ogj 3p5o 1ps3 3dx1 3d4z 3dx2 3ejr 3l7b 4eky 3g2n 3syr 3ebp 2w66 2w4x 2wca 2xj7 2vvn 3aru 3arv 3ary 3arq 3arp 4ih5 4ih7 3cj4 4eo8 3gnw 1gpk 1gpn 1h23 1h22 1e66 3f3a 3f3c 4mme 3f3d 3f3e 2wbg 2cbv 2j78 2j7h 2cet 3udh 3rsx 4djv 2vkm 4gid 4jfs 4j28 2wvt 2xii 4pcs 3rr4 1s38 1r5y 3gc5 3ge7 4dli 2zb1 4f9w 3e92 3e93 4owm 3twp 3r88 4gkm 3qqs 3gv9 3gr2 4kz6 4jxs 2r9w 2hb1 1bzc 2qbr 2qbq 2qbp 1q8t 1ydr 1q8u 1ydt 3ag9 3fcq 1z9g 1qf1 5tmn 4tmn 4ddk 4ddh 3ivg 3coz 3coy 3pxf 4eor 2xnb 1pxn 2fvd 4k77 4e5w 4ivb 4ivd 4ivc 4f09 4gfm 4hge 4e6q 4jia 2brb 2br1 3jvr 3jvs 1nvq 3acw 4ea2 2zcr 2zy1 2zcq 1bcu 3bv9 1oyt 2zda 3utu 3u9q 2yfe 3fur 3b1m 2p4y 3uo4 3up2 3e5a 2wtv 3myg 3kgp 1c5z 1o5b 1owh 1sqa 4jsz 3kwa 2weg 3ryj 3dd0 2xdl 3b27 1yc1 3rlr 2yki 1z95 3b68 3b5r 3b65 3g0w 4u4s 1p1q 1syi 1p1n 2al5 3g2z 3g31 4de2 4de3 4de1 1vso 4dld 3gbb 3fv2 3fv1 4mgd 2qe4 1qkt 2pog 2p15 2y5h 1lpg 2xbv 1z6e 1mq6 1nc3 1nc1 1y6r 4f2w 4f3c 1uto 4abg 3gy4 1k1i 1o3f 2yge 2fxs 2iwx 2wer 2vw5 4kzq 4kzu 4j21 4j3l 3kr8 2ymd 2wnc 2xys 2wn9 2x00 3ozt 3ozs 3oe5 3oe4 3nw9 3ao4 3zt2 3zsx 4cig 3zso 3n7a 4ciw 3n86 3n76 2xb8 4bkt 4w9c 4w9l 4w9i 4w9h 3nq9 3ueu 3uev 3uew 3uex 3lka 3ehy 3tsk 3nx7 4gr0 3dxg 3d6q 1w4o 1o0h 1u1b'.split()
## Debugging using first 10 because contains first two clusters. Very quick though.
targets = targets[:10] if debug else targets

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
        name=f'{testname}-{target}',

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/rosetta_scripts.{extension}',
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


"""
# Submitting PyRosetta job
hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
    name=f'{testname}-{PyRosetta-example-job}',

    #==> EDIT HERE, substitute <MyPythonScript.py> and <script-arguments-if-any> with name of your PyRosetta Python script (located inside your test dir) and command line flags for it
    executable = config["python_virtual_environment"]["python"],
    arguments = '<MyPythonScript.py> <script-arguments-if-any>',
    working_dir = prefix,
    jobs_to_queue = min(nstruct, 50),
    log_dir = hpc_logs,
    time=24,
    block=False)
)
"""

#==> EDIT HERE
benchmark.save_variables('debug targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
