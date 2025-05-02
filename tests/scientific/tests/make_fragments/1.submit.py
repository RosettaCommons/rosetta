#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  _template_/1.submit.py
## @brief this script is part of <template> scientific test
## @author Sergey Lyskov

import os, sys, time
import gzip
import json
import pathlib
import benchmark
from benchmark.tests import parallel_execute
from concurrent.futures import ProcessPoolExecutor as PPE, wait as PPE_wait
import subprocess
import datetime
import shutil
import copy
from tempfile import TemporaryDirectory as T_D
import urllib.request as request
from contextlib import closing


benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "make_fragments"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = os.path.abspath(config['working_dir'])
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

debug_target_names = ["T0859", "T0860"]
num_cores_to_use = 6

# T0859 takes 2 hrs
# T0950 takes 7 hrs

# WARNING!!!!!!!
# BECAUSE THIS TAKES TOO LONG ON THE TESTING SERVER WE HAD TO REMOVE MANY TARGETS
# TO FULLY TEST THE DATASET UNCOMMENT THE FOLLOWING LINE AND COMMENT THE LINE AFTER
skip_list = []
# skip_list = ["T0963", "T0960", "T0953s1", "T0863", "T0912", "T0904", "T0918", "T0954", "T0902", "T0905", "T0920", "T0945", "T0943", "T0873", "T0966", "T0920", "T0896"]



# WARNING!!
# fdbft = fragment database for test
# but sparksx cant handle paths that long, so we have to make path as short as possible
frag_db_prefix = os.path.abspath(config['prefix'] + '/fdbft')
print("FRAG DB PREFIX", frag_db_prefix)

# >>> IF YOU WANT TO CHANGE THE BLAST DATABASE, ALTER THIS <<<
current_expected_sig = 'make_fragments benchmark v1.0 installed June 11 2019'

try:
    with open(os.path.join(frag_db_prefix, 'signature')) as fh:
        found_sig = fh.read()
except FileNotFoundError:
    found_sig = ''

# SPARKSX requires 2.7
conda_config = copy.deepcopy(config)
conda_config['platform']['python'] = '2.7'

# This is the website for compbiocore https://cbc.brown.edu/ https://github.com/compbiocore
# This is the website for bioconda bioconda: https://bioconda.github.io/
# Unfortunately conda doesn't have the necessary perl libraries we need, so we need these channels
conda_channels_to_add = ['bioconda', 'compbiocore']

conda = benchmark.tests.setup_conda_virtual_environment(
        working_dir,
        conda_config['platform'],
        conda_config)

for channel in conda_channels_to_add:
    benchmark.execute(f'Adding extra extra channel {channel}...',
                      f'{conda.activate} && conda config --env --add channels {channel}' )

conda_extra_packages = 'perl perl-if perl-parallel-forkmanager perl-switch perl-text-balanced perl-json'
benchmark.execute( f'Setting up extra packages {conda_extra_packages}...', f'{conda.activate} && conda install --quiet --yes {conda_extra_packages}' )

# We should be careful this doesn't change.
perl5libdir = os.path.join(conda.root, 'lib/site_perl/5.26.2/')
assert os.path.isdir(perl5libdir)

make_fragments_pl = os.path.join(frag_db_prefix, "make_fragments.pl")

make_frags_environment_commands = (
        f'{conda.activate}'
        f' && export ROSETTA_DATABASE={config["rosetta_dir"]}/database'
        f' && export VALL={config["rosetta_dir"]}/tools/fragment_tools/vall.jul19.2011'
        f' && export FRAGMENT_PICKER={config["rosetta_dir"]}/source/bin/fragment_picker.{extension}'
        f' && export FRAGMENT_PICKER_NUM_CPUS=1'
        f' && export BLAST_NUM_CPUS={num_cores_to_use}'
        f' && export PERL5LIB={perl5libdir}')

# Always copy over new code from install deps/make_frags
if os.path.isdir(frag_db_prefix):
    shutil.copy(f"{config['rosetta_dir']}/tools/fragment_tools/install_dependencies.pl", os.path.join(frag_db_prefix, "install_dependencies.pl"))
    shutil.copy(f"{config['rosetta_dir']}/tools/fragment_tools/make_fragments.pl", os.path.join(frag_db_prefix, "make_fragments.pl"))


if not os.path.isdir(frag_db_prefix): pathlib.Path(frag_db_prefix).mkdir(exist_ok=True, parents=True)

seqres_fn = os.path.join(frag_db_prefix, "pdb_seqres.txt")
if not os.path.isfile(seqres_fn):
    with closing(request.urlopen("https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt")) as r:
            with open(seqres_fn, 'wb') as f:
                        shutil.copyfileobj(r, f)
    # subprocess.call(["curl", "ftp://ftp.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt", "-o", seqres_fn])

print("Checking for frag_db...")
if not os.path.isdir(frag_db_prefix) or found_sig != current_expected_sig:
    print(f"{os.path.isdir(frag_db_prefix)} or {found_sig} != {current_expected_sig} -- rebuilding database")
    pathlib.Path(frag_db_prefix).mkdir(exist_ok=True, parents=True)
    # Make database first time
    shutil.copy(f"{config['rosetta_dir']}/tools/fragment_tools/install_dependencies.pl", os.path.join(frag_db_prefix, "install_dependencies.pl"))
    shutil.copy(f"{config['rosetta_dir']}/tools/fragment_tools/make_fragments.pl", os.path.join(frag_db_prefix, "make_fragments.pl"))

    benchmark.execute('Installing dependencies...',
                      (f'cd {frag_db_prefix}'
                       f' && {conda.activate}'
                       f' && export LOCAL_NR_COPY={config["mounts"]["fragment_picking_database"]}'
                       f' && ./install_dependencies.pl standard localnrcopy'
                       f' && cd {os.getcwd()}'))

    # Run once to double check for dependencies + make sure everything works
    with T_D() as tmpdir:
        with open(os.path.join(tmpdir, "1a19.fasta"), 'w') as fh:
            fh.write(f">1a19.fasta\nMKKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDCLTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGADITIILS\n")
        benchmark.execute(
                'Running test-make_fragments.pl...',
                (f'cd {tmpdir}'
                 f' && {make_frags_environment_commands}'
                 f' && {make_fragments_pl} -n_frags 200 -frag_sizes 3,9 1a19.fasta'))
        assert os.path.isfile(os.path.join(tmpdir, 't001_.200.3mers'))
        assert os.path.isfile(os.path.join(tmpdir, 't001_.200.9mers'))

    with open(os.path.join(frag_db_prefix, "signature"), 'w') as fh:
        fh.write(current_expected_sig)
else:
    print("Found signature matched, not rebuilding db")

with open("make_a_fragfile.sh", 'w') as fh:
    fh.write(f"""#!/usr/bin/env bash
# '''
# usage:
# ./make_a_fragfile.sh {{SEQVENCE}} {{ABSPATH_TO_MAKE_FRAGMENTS.pl}} {{ABSPATH_TO_FINAL_DIRECTORY}}
# '''

one=$1
two=$2
three=$3
shift
shift
shift

{make_frags_environment_commands}

ogdir=$(pwd)
tmpdir=$(mktemp -d)
cd $tmpdir
echo -e ">mysequence\n${{one}}\n" > sequence.fasta
start=$(date +%s)
$two -n_frags 200 -frag_sizes 3,9 sequence.fasta -nohoms
end=$(date +%s)
runtime=$((end-start))
echo "Time to run: ${{runtime}}" >runtime.txt

echo mkdir -p $three
mkdir -p $three
echo  cp -r $(pwd) ${{three}}/
cp -r $(pwd) ${{three}}/
echo  cp t001_.200.3mers $three
cp t001_.200.3mers $three
echo cp t001_.200.9mers $three
cp t001_.200.9mers $three
echo cp sequence.fasta $three
cp sequence.fasta $three
echo cp runtime.txt $three
cp runtime.txt $three
cd $ogdir
""")
subprocess.call(['chmod', '+x', 'make_a_fragfile.sh'])


casp_targets = os.path.join(config["rosetta_dir"], "tests/scientific/data/make_fragments/casp_12-13.json.gz")

hpc_logs = f'{working_dir}/hpc-logs'
pathlib.Path(hpc_logs).mkdir(exist_ok=True, parents=True)

hpc_job_ids = []

with gzip.open(casp_targets, 'rt') as fh:
    json_data = json.loads(fh.read())
print(f"loaded casp_targets from {casp_targets}")

if debug:
    for db_name, db_data in list(json_data.items()):
        for target_name, target_data in list(db_data.items()):
            if target_name not in debug_target_names:
                del json_data[db_name][target_name]
all_names = []
for db_name, db_data in list(json_data.items()):
    for target_name, target_data in list(db_data.items()):
        all_names.append(target_name)
print(f"Is debug?: {debug} names: {' '.join(all_names)}")

print("Writing pdb files and setting workdirs")
for db_name, db_data in json_data.items():
    db_dir = os.path.join(working_dir, db_name)
    shutil.rmtree(db_dir, ignore_errors=True)
    # for x in range(10):
    #     print("WARNING revert this line!")
    for target_name, target_data in db_data.items():
        target_loc = os.path.abspath(os.path.join(db_dir, target_name))
        pathlib.Path(target_loc).mkdir(exist_ok=True, parents=True)
        target_pdb = os.path.join(target_loc, f"{target_name}.pdb")
        target_data["pdb_file"] = target_pdb
        target_data["workdir"] = os.path.abspath(target_loc)
        with open(target_pdb, 'w') as fh:
            fh.write(target_data["pdb_text"])

def fn(driver, max_num_jobs):
    jobs = list(driver.jobs)
    jobs = [job for job in jobs if not driver.complete(job)]
    if len(jobs) < max_num_jobs:
        return False # stops sleeping
    return True # continues sleeping

current_jobs = 0
my_jobs = {}

for db_name, db_data in json_data.items():
    for target_name, target_data in db_data.items():
        workdir = target_data["workdir"]
        # 9mers is -1
        expected_fragfiles = [os.path.abspath(os.path.join(target_data["workdir"], x)) for x in ['t001_.200.3mers', 't001_.200.9mers']]
        print("WARNING This will checkpoint from the last run and should only be used for debugging!!!")
        if all(os.path.isfile(x) for x in expected_fragfiles):
            print("WARNING FOUND ALL", expected_fragfiles, "NOT GOING TO RUN")
            target_data["final_files"] = expected_fragfiles
            continue
        # uncomment for local running
        # my_jobs[f'{testname}-{db_name}-{target_name}'] = f"cd {workdir} && {os.path.join(os.getcwd(), 'make_a_fragfile.sh')} {target_data['sequence']} {make_fragments_pl} {workdir}"
        target_data["run_command"] = f"cd {workdir} && {os.path.join(os.getcwd(), 'make_a_fragfile.sh')} {target_data['sequence']} {make_fragments_pl} {workdir}"

        # hpc_job_ids.append(hpc_driver.submit_serial_hpc_job(
        #     name=f'{testname}-{db_name}-{target_name}',
        #     executable=os.path.join(os.getcwd(), "make_a_fragfile.sh"),
        #     arguments=f"{target_data['sequence']} {make_fragments_pl} {workdir}",
        #     working_dir=workdir,
        #     jobs_to_queue=1,
        #     log_dir=hpc_logs,
        #     time=150,
        #     block=False))
        # hpc_driver.block_until(False, fn, hpc_driver, 8)

        target_data["final_files"] = expected_fragfiles

# hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

parallel_cores = int(config['cpu_count'])
print(f"going to run in parallel {len(my_jobs)} jobs on {parallel_cores} cores")

def subprocess_wrapper(command):
    ret = subprocess.run(command, shell=True)
    if ret.returncode:
        print("FAILED WITH COMMAND:\n{command}\ngot return code {ret}")
        raise RuntimeError("FAILED WITH COMMAND:\n{command}\ngot return code {ret}")


# Keep just in case for debugging
to_delete = []
for db_name, db_data in json_data.items():
    for target_name, target_data in db_data.items():
        if target_name.upper() in skip_list:
            print(f"SKIPPING {target_name} BECAUSE OTHERWISE IT TAKES TOO LONG!")
            to_delete.append((db_name, target_name))
for to_del in to_delete:
    del json_data[to_del[0]][to_del[1]]


with PPE(int(parallel_cores*3)) as ppe:
    to_submit = []
    for db_name, db_data in json_data.items():
        for target_name, target_data in db_data.items():
            if "run_command" not in target_data:
                print(f"WARNING COULDNT FIND 'run_command' IN 'target_data' for {db_name} {target_name}")
            to_submit.append((subprocess_wrapper, target_data['run_command'], target_data["sequence"]))
    # Longest jobs first
    to_submit.sort(key=lambda x: len(x[-1]), reverse=True)
    jobs = [ppe.submit(x[0], x[1]) for x in to_submit]
    PPE_wait(jobs)

    # parallel_execute("make_frags", my_jobs, rosetta_dir, working_dir, int(config['cpu_count']/(num_cores_to_use/2)), time=1000)

with open("1.checkpoint.json", 'w') as fh:
    fh.write(json.dumps(json_data))

#==> EDIT HERE
benchmark.save_variables('debug targets nstruct working_dir testname json_data')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
