#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  antibody_grafting/1.submit.py
## @brief this script is part of antibody_grafting scientific test
## @author Sergey Lyskov

import os, subprocess, sys, time, shutil
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "antibody_grafting"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

# check if blastp is in path, else ..
ncbi_blast_plus = shutil.which("blastp")
if ncbi_blast_plus == None: # overwrite/guess
    ncbi_blast_plus = "/Users/benchmark/prefix/ncbi-blast-2.9.0+/bin/blastp"

#==> EDIT HERE -- change target names with directories
command_line = '''
-antibody:blastp {ncbi_blast_plus}
-fasta {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.truncated.fasta
-in:file:native {rosetta_dir}/tests/scientific/data/{testname}/benchmark_crystals/{target}_trunc.pdb
-database {rosetta_dir}/database/
-antibody:grafting_database {rosetta_dir}/database/additional_protocol_data/antibody/
-antibody:exclude_pdbs {target}
-antibody:n_multi_templates 1
-no_relax
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

#==> EDIT HERE
# exclude 1x9q for now - there is a grafting error and it needs custom regex
targets = '1dlf 1fns 1gig 1jfq 1jpt 1mfa 1mlb 1mqk 1nlb 1oaq 1seq 2adf 2d7t 2e27 2fb4 2fbj 2r8s 2v17 2vxv 2w60 2xwt 2ypv 3e8u 3eo9 3g5y 3giz 3gnm 3go1 3hc4 3hnt 3i9g 3ifl 3liz 3lmj 3m8o 3mxw 3nps 3oz9 3p0y 3t65 3umt 3v0w 4f57 4h0h 4h20 4hpy 4nzu'.split()
targets = targets[:6] if debug else targets

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

# antibody has a special quirk where pdbs need to be unzipped first
# this is usually done by the grafting executable, but when multiple
# processes are run there is a bit of a collision, so let's unzip now
cwd = os.getcwd()
os.chdir(f'{rosetta_dir}/tools/antibody/antibody_database/')
for f in os.listdir("./"):
    # only extract if it is a bz2 and an unzipped version is not present
    if f.endswith('.bz2') and not os.path.isfile(f[:-4]):
        print("Extracting {}".format(f))
        subprocess.call(["bunzip2", "-k", f])
os.chdir(cwd)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:
    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)

    hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
        name=f'{testname}-{target}',

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/antibody.{extension}',
        arguments = command_line.format_map(vars()),
        working_dir = prefix,
        jobs_to_queue = 1,
        log_dir = hpc_logs,
        time=24,
        block=False)
    )


# if not debug:
#     hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
#     time.sleep(64)  # waiting for NFS caching
hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

#==> EDIT HERE
benchmark.save_variables('targets working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
