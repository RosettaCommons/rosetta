#!/usr/bin/env python
"""
Run cryoAssembly

"""

import argparse
import sys
import json
import typing
from typing import List
import itertools
import copy
import subprocess
import os
import pathlib
import glob

from dask_jobqueue import SLURMCluster
from dask.distributed import Client


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--premade_json", help="The premade json you want to use", default=None)
    parser.add_argument("--submit", action="store_true", help="Should we submit to slurm", default=None)
    parser.add_argument("--name", help="Slurm job name", default="cryo_as_cmdl")
    args = parser.parse_args()
    return vars(args)


def init_config():
    basic_config = {
        "assembly": {
            "connectivity_weight": 1000.0,
            "dc_weight": 1.0,
            "dsm": "",
            "extras": "",
            "final_out": 20,
            "n_iter": 1000,
            "n_out": 10,
            "null_weight": -1150,
            "onebody_weight": 260.0,
            "twobody_weight": 150.0,
            "jobs": 2,
        },
        "general": {
            "dc_file": "",
            "dc_d_add": "2",
            "equal_dc_strength": "true",
            "max_dc_score": "10000",
            "min_dc_score": "2",
            "executable": "",
            "database": "",
            "extras": "",
            "in:file:silent": "",
            "mapfile": "",
            "mute": "true",
            "name": "ssnf",
        },
        "iterover": {"general": {}, "twobody": {}, "assembly": {"twobody_weight": [10, 100]}},
        "scheduler": {
            "cpu_limit": 10000,
            "worker_walltime": 4,
            "cores": 10000,
            "queue": "dimaio",
            "name": "cryo_ass_cmdline",
        },
        "twobody": {
            "constrain_type": "legacy",
            "extras": (
                "-hbond_bb_per_residue_energy"
                " -mh:path:scores_BB_BB /rosetta_location/database/additional_protocol_data_submodule/motif_dock/xh_16"
            ),
            "bb_iter": "0",
            "rb_iter": "20",
            "refine_dens_wt": "10",
            "score_closability": "true",
            "scorefile_out": "sc",
            "slide_domains": "5",
            "legacy_radii": "12",
            "trim_clash": "true",
            "allow_domains_to_minimize": "true",
            "core_split": 15000,
        },
    }
    return basic_config


def get_iterover_params(options: typing.Dict[typing.Any, typing.Any]):
    basic_keys = ["general", "twobody", "assembly"]

    params = []
    for key in basic_keys:
        this_params = []
        keys = list(options["iterover"][key].keys())
        for x in itertools.product(*map(options["iterover"][key].get, keys)):
            this_params.append(dict(zip(keys, x)))
        params.append(this_params)
    final_params = [dict(zip(basic_keys, x)) for x in itertools.product(*params)]
    return final_params


def get_serial_name(name, params, prefix=""):
    """
    {folder_type(assembly, twobody)+(option1=?)+(option2=?)}
    example: assembly:dc_weight = [1,5,10] & twobody:dc_d = [ 20,30,40 ]
    would build:
    twobody+Tdc_d=20, twobody+Tdc_d=30, twobody+Tdc_d=40
    assembly+Tdc_d=20+Adc_weight=1,assembly+Tdc_d=20+Adc_weight=5, assembly+Tdc_d=20+Adc_weight=10
    assembly+Tdc_d=30+Adc_weight=1,assembly+Tdc_d=30+Adc_weight=5, assembly+Tdc_d=30+Adc_weight=10
    assembly+Tdc_d=40+Adc_weight=1,assembly+Tdc_d=40+Adc_weight=5, assembly+Tdc_d=40+Adc_weight=10
    this should be easy to split, and organizie in postprocessing
    index = index after splitting by '+', variable_name = index 0 after selecting a split by '+'
    var_res = index 1 after selecting a split by '+'
    """
    if prefix:
        ret_name = f"{prefix}-{name}"
    else:
        ret_name = name
    for key, val in params.items():
        ret_name = f"{ret_name}+{key}={val}"
    return ret_name


def get_basic_flags(options):
    default_flags = (
        f"{options['general']['executable']}"
        f" -database {options['general']['database']}"
        " -hbond_bb_per_residue_energy true"
        " -ignore_unrecognized_res"
    )
    return default_flags


def add_general_flags(general_params, basic_flags: str):
    for key, val in general_params.items():
        if key in ["executable", "database", "name", "in:file:silent", "score_from_pdbs", "dc_file"]:
            continue
        if key in ["extras"]:
            basic_flags = f"{basic_flags} {val}"
            continue
        basic_flags = f"{basic_flags} -{key} {val}"
    if "dc_file" in general_params:
        if not os.path.isfile(general_params["dc_file"]):
            raise RuntimeError(f"couldn't find dc file: {general_params['dc_file']}")
        basic_flags = f"{basic_flags} -dc_file {general_params['dc_file']}"
    if "*" in general_params["in:file:silent"]:
        basic_flags = (
            f"{basic_flags} -in:file:silent {' '.join([x for x in glob.glob(general_params['in:file:silent'])])}"
        )
    elif general_params["in:file:silent"]:
        basic_flags = f"{basic_flags} -in:file:silent {general_params['in:file:silent']}"
    elif "*" in general_params["score_from_pdbs"]:
        basic_flags = (
            f"{basic_flags} -score_from_pdbs {' '.join([x for x in glob.glob(general_params['score_from_pdbs'])])}"
        )
    elif general_params["score_from_pdbs"]:
        basic_flags = f"{basic_flags} -score_from_pdbs {general_params['score_from_pdbs']}"
    return basic_flags


def _run_twobody_scoring(runcommand: str, expected_file: str, expected_logfile: str):

    if os.path.isfile(expected_file):
        if len(open(expected_file, "rb").read()) < 100:
            os.remove(expected_file)
        else:
            return expected_file

    with open(expected_logfile, "w") as fh:
        ret = subprocess.run(runcommand.split(), stdout=fh, stderr=subprocess.STDOUT)

    if ret.returncode:
        with open(f"{expected_logfile}_brk", "w") as fh:
            fh.write(f"Failure with command: {runcommand}\nlogfile {expected_logfile}")
        print(f"Failure with command: {runcommand}, logfile {expected_logfile}")
        raise RuntimeError(f"Failure with command: {runcommand}")
    if not os.path.isfile(expected_file):
        print(f"Failure to acquire {expected_file} with command: {runcommand}, logfile {expected_logfile}")
        raise RuntimeError(f"Failure to acquire {expected_file} with command: {runcommand}")
    return expected_file


def run_twobody_scoring(
    options: typing.Dict[typing.Any, typing.Any],
    iterover_params,
    base_command: str,
    output_directory: str,
    standardclient: Client,
):
    """
    returns a future of score files! be careful!
    """
    pathlib.Path(output_directory).mkdir(exist_ok=True)
    logpath = os.path.join(output_directory, "logs")

    pathlib.Path(logpath).mkdir(exist_ok=True)
    twobody_opts = copy.deepcopy(options["twobody"])
    twobody_opts.update(iterover_params)
    for key, val in twobody_opts.items():
        if key in ["extras"]:
            base_command = f"{base_command} {val}"
            continue
        if key in ["core_split"]:
            continue
        base_command = f"{base_command} -{key} {val}"

    base_command = f"{base_command} -out:prefix {output_directory}/ -mode score"

    total_split = int(twobody_opts["core_split"])

    output_score_files = []
    output_score_files.append(os.path.join(output_directory, f"{twobody_opts['scorefile_out']}_onebody.sc"))
    for idx in range(1, total_split + 1):
        runscript = f"{base_command} -core_idx {idx} {total_split}"
        if idx == 1:
            with open("twobody_example", "w") as fh:
                fh.write(runscript)
        expected_file = os.path.join(output_directory, f"{twobody_opts['scorefile_out']}_{idx:06}_{total_split:06}.sc")

        logfile = os.path.join(logpath, f"{twobody_opts['scorefile_out']}_{idx:06}_{total_split:06}.log")
        if os.path.isfile(expected_file) and len(open(expected_file, "rb").read()) > 100:
            output_score_files.append(expected_file)
        else:
            if os.path.isfile(expected_file):
                os.remove(expected_file)
            output_score_files.append(standardclient.submit(_run_twobody_scoring, runscript, expected_file, logfile))
    return output_score_files


def _run_assembly(
    runcommand: str, twobody_scorefiles: typing.List[str], expected_file: str, expected_logfile: str = "/dev/null"
):
    runcommand = f"{runcommand} -scorefile_in {' '.join(twobody_scorefiles)}"
    with open(expected_logfile, "w") as fh:
        print(runcommand, file=sys.stderr)
        ret = subprocess.run(runcommand.split(), stdout=fh, stderr=subprocess.STDOUT)

    if ret.returncode:
        raise RuntimeError(f"Failure with logfile: {expected_logfile} Failure with command: {runcommand}")
    if not os.path.isfile(expected_file):
        raise RuntimeError(f"Failure to acquire {expected_file} with command: {runcommand}")
    return expected_file


def _run_compile_assembly(
    runcommand: str,
    twobody_scorefiles: typing.List[str],
    input_files: typing.List[str],
    expected_files: typing.List[str],
    expected_logfile: str = "/dev/null",
):
    runcommand = f"{runcommand} -scorefile_in {' '.join(twobody_scorefiles)}"
    runcommand = f"{runcommand} -result_file_names {' '.join(input_files)}"
    with open(expected_logfile, "w") as fh:
        ret = subprocess.run(runcommand.split(), stdout=fh, stderr=subprocess.STDOUT)

    if ret.returncode:
        with open("failed_compile.sh", "w") as fh:
            fh.write(runcommand)
        raise RuntimeError(f"Failure with command: {runcommand}")
    for expected_file in expected_files:
        if not os.path.isfile(expected_file):
            print(f"Failure to acquire {expected_file} with command: {runcommand}", file=sys.stderr)
    return expected_files


def run_assembly(
    options: typing.Dict[typing.Any, typing.Any],
    iterover_params,
    base_command: str,
    twobody_scorefiles,
    output_directory: str,
    standardclient: Client,
):
    pathlib.Path(output_directory).mkdir(exist_ok=True)
    all_results_folder = os.path.join(output_directory, "all_results")
    pathlib.Path(all_results_folder).mkdir(exist_ok=True)
    top_results_folder = os.path.join(output_directory, "top_results")
    pathlib.Path(top_results_folder).mkdir(exist_ok=True)

    assembly_opts = copy.deepcopy(options["assembly"])
    assembly_opts.update(iterover_params)

    for key, val in assembly_opts.items():
        if key in ["extras"]:
            base_command = f"{base_command} {val}"
            continue
        if key in ["jobs", "final_out"]:
            continue
        base_command = f"{base_command} -{key} {val}"
    base_command = f"{base_command}" " -mode assemble" " -start_random_assembly true"

    total_jobs = int(assembly_opts["jobs"])

    assembly_files = []
    for idx in range(1, total_jobs + 1):
        output_file_flag = os.path.join(all_results_folder, f"assembly_{idx:04}")
        output_file = os.path.join(all_results_folder, f"assembly_{idx:04}.assembly_results")
        logfile = os.path.join(all_results_folder, f"assembly_{idx:04}.log")
        run_command = f"{base_command}" " -dump_as_text" f" -out:file:silent {output_file_flag}"
        if os.path.isfile(output_file):
            assembly_files.append(output_file)
        else:
            assembly_files.append(
                standardclient.submit(_run_assembly, run_command, twobody_scorefiles, output_file, logfile)
            )
    # compile results
    final_file_base_name = os.path.join(top_results_folder, "assembly")
    run_command = f"{base_command}" " -text_to_results" " -core_idx 1 1" f" -out:file:silent {final_file_base_name}"
    top_assembly_files = standardclient.submit(
        _run_compile_assembly,
        run_command,
        twobody_scorefiles,
        assembly_files,
        [f"{final_file_base_name}_{x+1:04}.silent" for x in range(assembly_opts["final_out"])],
        os.path.join(top_results_folder, "compile.log"),
    )
    return top_assembly_files


def run_cryo_assembly(options: typing.Dict[typing.Any, typing.Any], output_directory: str, client: Client) -> List[str]:
    """
    require dictionary stays sorted the same way!
    """
    parameters = get_iterover_params(options)
    jobs = {}
    basic_flags = get_basic_flags(options)

    jobs_to_complete = []
    serial_names = []
    for param_set in parameters:
        general_params = options["general"]
        general_params.update(param_set["general"])
        general_serial_name = get_serial_name("general", param_set["general"])
        if general_serial_name not in jobs:
            jobs[general_serial_name] = add_general_flags(general_params, basic_flags)

        twobody_params = param_set["twobody"]
        twobody_serial_name = get_serial_name("twobody", twobody_params, general_serial_name)
        twobody_output_dir = os.path.join(output_directory, twobody_serial_name)
        if twobody_serial_name not in jobs:
            jobs[twobody_serial_name] = run_twobody_scoring(
                options, twobody_params, jobs[general_serial_name], twobody_output_dir, client
            )

        assembly_params = param_set["assembly"]
        assembly_serial_name = get_serial_name("assembly", assembly_params, twobody_serial_name)
        assembly_output_dir = os.path.join(output_directory, assembly_serial_name)
        if assembly_serial_name not in jobs:
            jobs[assembly_serial_name] = run_assembly(
                options,
                assembly_params,
                jobs[general_serial_name],
                jobs[twobody_serial_name],
                assembly_output_dir,
                client,
            )
            jobs_to_complete.append(jobs[assembly_serial_name])
            serial_names.append(assembly_serial_name)
    found_assembly_files = []
    top_assembly_files = client.gather(jobs_to_complete)
    for (serial_name, expected_files) in zip(serial_names, list(top_assembly_files)):
        for expected_file in expected_files:
            if not os.path.isfile(expected_file):
                print(f"WARNING: {serial_name} COULDN'T FIND FILE {expected_file}")
            else:
                found_assembly_files.append(expected_file)
    return found_assembly_files


def cryo_assembly_from_commandline(options):
    # Start slurm interface
    if not options["premade_json"]:
        if not os.path.isfile("config.json"):
            with open("config.json", "w") as fh:
                fh.write(json.dumps(init_config(), indent=2))
        else:
            raise RuntimeError("config.json exists but not supplied")
        sys.exit(0)
    options = json.loads(open(options["premade_json"]).read())

    with SLURMCluster(
        cores=1,
        processes=1,
        memory=options["scheduler"]["memory"],
        queue=options["scheduler"]["queue"],
        walltime=options["scheduler"]["worker_walltime"],
        name=f"{options['scheduler']['name']}_worker",
        job_name=f"{options['scheduler']['name']}_worker",
        extra=["--no-nanny", "--no-dashboard"],
    ) as cluster:
        cluster.adapt(minimum=0, maximum=options["scheduler"]["cores"], wait_count=400)
        with Client(cluster) as client:
            print("Server info (standard):", client, client.scheduler_info(), file=sys.stderr)
            run_cryo_assembly(options, client)
            print("Finished cryo assembly!")


if __name__ == "__main__":
    args = parseargs()
    if args["submit"]:
        from daimyo.distribution.slurm import ez_run_dask_master

        this_command = ["python", "-u"] + [x for x in sys.argv if x != "--submit"]
        extra = None
        if "depend" in args and args["depend"]:
            extra = [f"#SBATCH --dependency=afterany:{args['depend']}"]
        ez_run_dask_master(this_command, args["name"], hours=20, gb_mem=10, extra=extra)
    else:
        del args["submit"]
        if "depend" in args:
            del args["depend"]
        cryo_assembly_from_commandline(args)
