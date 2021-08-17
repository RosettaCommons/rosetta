#!/usr/bin/env python
"""
A class to run dgdp on a selection of inputs

TODO:
    Maybe move to class?
"""

import os
import sys
import argparse
import json
import sys
import asyncio
from typing import List, Dict, Any, Optional
import uuid
import zlib
import math
import pathlib
import copy
import shutil

try:
    from dask_jobqueue import SLURMCluster
    from dask.distributed import Client, LocalCluster, as_completed
except ModuleNotFoundError:
    print("Warning: Unable to run cryo_dock_and_assemble without dask_jobqueue installed")
    sys.exit()

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from distribution import run_subprocess

from user_input import (
    add_default_rosetta_exes_arguments,
    add_default_jobdist_arguments,
    save_args,
    get_rosetta_exe_and_database_from_args,
)


class DGDPConfig:
    def __init__(self, dgdp_exe: str, rosetta_database: str, outdir: str, basic_flags: Optional[str] = None):
        self.dgdp_exe = dgdp_exe
        self.rosetta_database = rosetta_database
        self.outdir = outdir
        self.basic_flags = basic_flags

    def set_basic_flags(self, input_dict) -> None:
        """
        Notes:
            we ignore mapreso because we set it manually to low because docking showed better
            results then.
        """
        flags_to_ignore = [
            "cores",
            "xml",
            "run",
            "edensity::mapreso",
            None,
            "queue",
            "worker_walltime",
            "worker_memory",
            "executable",
            "database",
            "name",
            "extras",
            "in:file:s",
            "multi_natives",
            "redo",
            "refinement_reso",
            "final_result_name",
            "walltime",
            "memory",
            "dgdp_exe",
            "rosetta_database",
            "rosetta_database_remote",
            "rosetta_bin_location_remote",
            "rosetta_bin_location",
            "mapfiles",
        ]

        basic_flags = f"{self.dgdp_exe} -database {self.rosetta_database}"
        for key, val in input_dict.items():
            if key in flags_to_ignore:
                continue
            basic_flags += f" -{key} {val}"
        # Default flags always on
        basic_flags += " -ignore_unrecognized_res -ignore_zero_occupancy false" f" {input_dict['extras']}"
        self.basic_flags = basic_flags
        print("set basic flags", self.basic_flags)


def parseargs():
    parser = argparse.ArgumentParser()

    basic_job_input_parser = parser.add_argument_group("basic_job_inputs")
    basic_job_input_parser.add_argument("--pdbs", help="Pdbs to dock into density", required=True, nargs="+")
    basic_job_input_parser.add_argument("--mapfiles", nargs="+", help="Electron density mapfiles", required=True)
    basic_job_input_parser.add_argument("--mapreso", help="Electron density mapfile resolution", required=True)
    basic_job_input_parser.add_argument(
        "--multi_native", help="All possible natives (docked already)", required=False, nargs="+"
    )

    algorithm_options_parser = parser.add_argument_group("algorithm_options")
    algorithm_options_parser.add_argument("--clust_radius", help="Radius to cluster in", required=False, type=int)
    algorithm_options_parser.add_argument(
        "--convolute_single_residue",
        help="Should we convolute single residue?",
        required=False,
        action="store_true",
        default=False,
    )
    algorithm_options_parser.add_argument(
        "--laplacian_offset", help="Laplacian offset 0 if none", type=int, required=False, default=None
    )
    algorithm_options_parser.add_argument("--n_filtered", help="Number of results to refine", type=int, required=True)
    algorithm_options_parser.add_argument("--n_output", help="Number of results to output", type=int, required=True)
    algorithm_options_parser.add_argument("--n_to_search", help="Number of points to search", type=int, required=True)
    algorithm_options_parser.add_argument("--point_radius", help="Point radius", type=float, required=True)
    algorithm_options_parser.add_argument(
        "--min", help="Rigid Body minimize? (true/false)", required=False, default=False, action="store_true"
    )
    algorithm_options_parser.add_argument(
        "--min_bb", help="Backbone Body minimize? (true/false)", required=False, default=False, action="store_true"
    )
    algorithm_options_parser.add_argument(
        "--constrain_refinement",
        help="Constrain refinement to x A (int, 0 is don't constrain)",
        type=int,
        default=0,
        required=False,
    )
    algorithm_options_parser.add_argument(
        "--refinement_reso", help="Refinement resolution", type=int, default=None, required=False
    )
    algorithm_options_parser.add_argument(
        "--rot_middle_ca",
        help="In SPHARM docking, will rotate pose on middle CAs (ca closest to COM)",
        action="store_true",
        required=False,
        default=False,
    )
    algorithm_options_parser.add_argument(
        "--extras",
        help='Enter any extra commands you would like to run (blank if none) (use "")',
        type=str,
        required=False,
        default="",
    )
    algorithm_options_parser.add_argument(
        "--max_rot_per_trans",
        help="max rotations to take per translation(empty uses default)",
        type=int,
        required=False,
        default=None,
    )
    algorithm_options_parser.add_argument("--bw", help="SPHARM bandwidth", type=int, required=False, default=None)

    job_distribution_options = parser.add_argument_group("job_distribution_options")
    job_distribution_options.add_argument("--cores", help="Total cores you want to execute on", type=int, required=True)
    job_distribution_options.add_argument(
        "--overwrite",
        help="Overwrite any previous results? (true/false)",
        action="store_true",
        required=False,
        default=None,
    )
    job_distribution_options.add_argument(
        "--redo", help="re run everything", action="store_true", required=False, default=False
    )

    script_options_parser = parser.add_argument_group("script_options")
    script_options_parser.add_argument("-j", "--premade_json", help="The premade json you want to use", default=None)
    script_options_parser.add_argument(
        "--final_result_names", help="Final result name (what to call the silent structs)", required=True, nargs="+"
    )

    parser = add_default_rosetta_exes_arguments(parser)
    parser = add_default_jobdist_arguments(parser)

    args = parser.parse_args()

    exes, rosetta_db = get_rosetta_exe_and_database_from_args(args, ["dgdp"])
    args.dgdp_exe = exes["dgdp"]
    args.rosetta_database = rosetta_db

    save_args()
    return args


def clean_up_args(args: Dict[str, Any]) -> Dict[str, Any]:
    """
    this is to make commandline compatible with the more config focused json which
    is used in assembly
    """
    in_f_s = []
    for pdb_grp in args["pdbs"]:
        curr_s = []
        for pdb_fn in pdb_grp:
            curr_s.append(os.path.abspath(pdb_fn))
        in_f_s.append(curr_s)
    args["in:file:s"] = in_f_s
    del args["pdbs"]

    args["mapfiles"] = [os.path.abspath(x) for x in args["mapfiles"]]

    if "multi_native" in args and "multi_natives" in args:
        raise RuntimeError("Cannot specify 'multi_native' and 'multi_natives' at the same time")
    if "multi_native" in args:
        args["multi_natives"] = [args["multi_native"]]
        del args["multi_native"]
    elif "multi_natives" in args:
        new_natives = []
        for x in args["multi_natives"]:
            new_natives.append([os.path.abspath(y) for y in x])
        args["multi_natives"] = new_natives

    if "mapreso" in args:
        args["edensity::mapreso"] = args["mapreso"]
        del args["mapreso"]
    keys_to_del = []
    for key, val in args.items():
        if val is None and key not in ["submit", "extras", "edensity::mapreso", "refinement_reso"]:
            keys_to_del.append(key)
    for key in keys_to_del:
        del args[key]
    return args


def load_premade_json(premade, skip_name_check=False):
    if len(premade) != 0:
        try:
            input_dict = json.load(open(premade))
            print("Loading your premade input file... ")
            for key, value in input_dict.items():
                print("loaded: ", key, " ", value)
            prev = input_dict["name"]
            if not skip_name_check:
                input_dict["name"] = input("Please enter a name for this job: ")
                if input_dict["name"] == "asdf":
                    print("Secret code entered, naming job ", prev)
                    input_dict["name"] = prev
            json.dump(input_dict, open("config.json", "w"), indent=4, sort_keys=True)
            return input_dict
        except ValueError:
            print("ERROR: premade input file could not be read into the json format")
            print("Check your input file and try again")
            print(f"Could not load: {premade}")
            sys.exit()


def _run_get_points(runcommand: List[str], expected_outputs: List[str], input_files: Dict[str, bytes], task_id: str):
    result_files = run_subprocess(runcommand, expected_outputs, input_files, task_id, save_logs=True)
    return result_files


async def run_get_points(pdb: str, config: DGDPConfig, options, client):
    pdb_no_path = os.path.basename(pdb)
    pdb_no_extension = os.path.splitext(pdb_no_path)[0]
    points_to_search_fname = f"{os.path.splitext(pdb_no_path)[0]}.pts"
    points_to_search_pdb_fname = f"{os.path.splitext(pdb_no_path)[0]}.pts.pdb"

    points_to_search_final_fname = os.path.join(config.outdir, f"{os.path.splitext(pdb_no_path)[0]}.pts")
    points_to_search_pdb_final_fname = os.path.join(config.outdir, f"{os.path.splitext(pdb_no_path)[0]}.pts.pdb")
    runscript = (
        f"{config.basic_flags}"
        " -mode find_pts"
        f" -s {os.path.abspath(pdb)}"
        f" -edensity::mapreso {options['edensity::mapreso']}"
        f" -points_to_search_fname {points_to_search_fname}"
        f" -points_to_search_pdb_fname {points_to_search_pdb_fname}"
        f" -out::file::silent {pdb_no_extension}"
    )
    input_files = {}
    input_files[pdb_no_path] = zlib.compress(open(os.path.abspath(pdb), "rb").read())
    if not os.path.isfile(points_to_search_final_fname) or options["redo"]:
        print("didnt find", points_to_search_final_fname)
        get_points_files = await client.submit(
            _run_get_points,
            runscript.split(),
            [points_to_search_fname, points_to_search_pdb_fname],
            input_files,
            uuid.uuid4().hex,
        )
        with open(points_to_search_final_fname, "wb") as fh:
            fh.write(zlib.decompress(get_points_files[points_to_search_fname]))

        with open(points_to_search_pdb_final_fname, "wb") as fh:
            fh.write(zlib.decompress(get_points_files[points_to_search_pdb_fname]))

        out_logfn = os.path.join(config.outdir, f"log_out/{pdb_no_extension}_get_points.log")
        with open(out_logfn, "wb") as fh:
            fh.write(zlib.decompress(get_points_files["out.log"]))
    return points_to_search_final_fname


def _run_search(runcommand: List[str], expected_outputs: List[str], input_files: Dict[str, bytes], task_id: str):
    print("running search")
    result_files = run_subprocess(runcommand, expected_outputs, input_files, task_id, save_logs=True)
    print("search result done", [key for key in result_files.keys()])
    return result_files, expected_outputs[0]


async def run_search(points_to_search_file: str, pdb: str, config: DGDPConfig, options, client: Client):
    tasks = []
    input_files = {}
    pdb_no_path = os.path.basename(pdb)
    pdb_no_extension = os.path.splitext(pdb_no_path)[0]
    input_files[pdb_no_path] = zlib.compress(open(pdb, "rb").read())
    points_to_search_fname = f"{pdb_no_extension}.pts"
    input_files[points_to_search_fname] = zlib.compress(open(points_to_search_file, "rb").read())
    input_files_future = await client.submit(dict, input_files)
    search_file_names = []
    for core in range(options["cores"]):
        core += 1
        avg_workload = math.floor(int(options["n_to_search"]) / int(options["cores"]))
        search_start = (core - 1) * avg_workload + 1
        search_end = core * avg_workload
        if core == options["cores"]:
            search_end = int(options["n_to_search"])
        pdb_no_extension = pdb.split("/")[-1][:-4]
        output_file = f"{pdb_no_extension}_search_{search_start:07}_{search_end:07}.sscores"
        runscript = (
            f"{config.basic_flags}"
            " -mode search_points"
            f" -core_idx {core} {options['cores']}"
            f" -s {pdb_no_path}"
            f" -points_to_search_fname {points_to_search_fname}"
            f" -point_search_results_fname {output_file}"
            f" -edensity::mapreso {options['edensity::mapreso']}"
            f" -out::file::silent {pdb_no_extension}"
        )
        out_final_fn = os.path.join(config.outdir, output_file)
        if not os.path.isfile(out_final_fn) or options["redo"]:
            tasks.append(
                client.submit(_run_search, runscript.split(), [output_file], input_files_future, uuid.uuid4().hex)
            )
        else:
            search_file_names.append(out_final_fn)
    # There should be 1 return file, and 1 log file
    async for result_fut in as_completed(tasks):
        result_files_expected_output = await result_fut
        result_files, expected_output = result_files_expected_output
        logfile_name = expected_output.replace(".sscores", ".log")

        with open(os.path.join(config.outdir, f"log_out/{logfile_name}"), "wb") as fh:
            fh.write(zlib.decompress(result_files["out.log"]))

        out_final_fn = os.path.join(config.outdir, expected_output)
        with open(out_final_fn, "wb") as fh:
            fh.write(zlib.decompress(result_files[expected_output]))
        search_file_names.append(out_final_fn)
    return search_file_names


def _run_combine_search(
    runcommand: List[str], expected_outputs: List[str], input_files: Dict[str, bytes], task_id: str
):
    result_files = run_subprocess(runcommand, expected_outputs, input_files, task_id, save_logs=True)
    return result_files


async def run_combine_search(
    pdb: str, search_result_file_names: List[str], config: DGDPConfig, options, client: Client
):
    pdb_no_path = os.path.basename(pdb)
    pdb_no_extension = os.path.splitext(pdb_no_path)[0]
    search_result_file_basenames = [os.path.basename(x) for x in search_result_file_names]
    output_file = f"{pdb_no_extension}_rfr.sscores"
    output_final_file = os.path.join(config.outdir, output_file)
    runscript = (
        f"{config.basic_flags}"
        f" -mode combine_search"
        f" -s {pdb_no_path}"
        f" -local_result_files {' '.join(search_result_file_basenames)}"
        f" -edensity::mapreso {options['edensity::mapreso']}"
        f" -out::file::silent {pdb_no_extension}"
        f" -combined_search_results_fname {output_file}"
    )

    input_files = {}
    if not os.path.isfile(output_final_file) or options["redo"]:
        for i, input_file in enumerate(search_result_file_names):
            input_files[search_result_file_basenames[i]] = zlib.compress(open(input_file, "rb").read())

        input_files[pdb_no_path] = zlib.compress(open(pdb, "rb").read())

        fut = client.submit(_run_combine_search, runscript.split(), [output_file], input_files, uuid.uuid4().hex)
        result_files = await fut

        with open(output_final_file, "wb") as fh:
            fh.write(zlib.decompress(result_files[output_file]))
        out_log_file = os.path.join(config.outdir, f"log_out/{pdb_no_extension}_combine_search.log")
        with open(out_log_file, "wb") as fh:
            fh.write(zlib.decompress(result_files["out.log"]))
    return output_final_file


def _run_refinement(runcommand: List[str], expected_outputs: List[str], input_files: Dict[str, bytes], task_id: str):
    result_files = run_subprocess(runcommand, expected_outputs, input_files, task_id, save_logs=True)
    return result_files, expected_outputs[0]


async def run_refinement(pdb: str, combined_search_result_file_name: str, config: DGDPConfig, options, client):
    """
    """
    print("running refinement on pdb:", pdb)
    pdb_no_path = os.path.basename(pdb)
    pdb_no_extension = os.path.splitext(pdb_no_path)[0]
    input_files = {}
    input_files[pdb_no_path] = zlib.compress(open(pdb, "rb").read())
    combined_search_result_file_basename = os.path.basename(combined_search_result_file_name)
    input_files[combined_search_result_file_basename] = zlib.compress(
        open(combined_search_result_file_name, "rb").read()
    )
    input_files_future = await client.submit(dict, input_files)

    # sometimes to total to refine is less thatn n_filtered
    with open(combined_search_result_file_name) as fh:
        total_to_refine = len(fh.readlines()) - 1  # -1 for header

    tasks = []
    refinement_file_names = []
    for core in range(options["cores"]):
        core += 1
        avg_workload = math.floor(total_to_refine / int(options["cores"]))
        refine_start = (core - 1) * avg_workload + 1
        refine_end = core * avg_workload
        if core == options["cores"]:
            refine_end = total_to_refine
        mapreso = 2
        if options["refinement_reso"] is not None:
            mapreso = options["refinement_reso"]
        output_file = f"{pdb_no_extension}_inter_{refine_start:05}_{refine_end:05}_ref.silent"
        runscript = (
            f"{config.basic_flags}"
            " -mode refine"
            f" -core_idx {core} {options['cores']}"
            f" -s {pdb_no_path}"
            f" -local_result_files {combined_search_result_file_basename}"
            f" -edensity::mapreso {mapreso}"
            f" -out::file::silent {output_file}"
        )
        output_final_file = os.path.join(config.outdir, output_file)
        if not os.path.isfile(output_final_file) or options["redo"]:
            tasks.append(
                client.submit(_run_refinement, runscript.split(), [output_file], input_files_future, uuid.uuid4().hex)
            )
        else:
            refinement_file_names.append(output_final_file)

    print("made all refinement tasks for pdb:", pdb)
    # There should be 1 return file, and 1 log file
    async for result_fut in as_completed(tasks):
        result_files_expected_output = await result_fut
        result_files, expected_output = result_files_expected_output

        expected_final_output = os.path.join(config.outdir, expected_output)
        with open(expected_final_output, "wb") as fh:
            fh.write(zlib.decompress(result_files[expected_output]))
        refinement_file_names.append(expected_final_output)

        logfile_name = os.path.join(config.outdir, "log_out", expected_output.replace("_ref.silent", ".log"))
        with open(logfile_name, "wb") as fh:
            fh.write(zlib.decompress(result_files["out.log"]))

    return refinement_file_names


def _run_combine_refinement(
    runcommand: List[str], expected_outputs: List[str], input_files: Dict[str, bytes], task_id: str
):
    result_files = run_subprocess(runcommand, expected_outputs, input_files, task_id, save_logs=True)
    return result_files


async def run_combine_refinement(pdb, refinement_file_names, config, options, standardclient):
    """
    """
    refinement_file_basenames = [os.path.basename(x) for x in refinement_file_names]
    pdb_no_path = os.path.basename(pdb)
    pdb_no_extension = os.path.splitext(pdb_no_path)[0]
    mapreso = 2
    if options["refinement_reso"] is not None:
        mapreso = options["refinement_reso"]

    output_file = f"{pdb_no_extension}_final.silent"
    output_final_file = os.path.join(config.outdir, output_file)

    runscript = (
        f"{config.basic_flags}"
        " -mode combine_refine"
        f" -s {pdb_no_path}"
        f" -out::file::silent {pdb_no_extension}"
        f" -edensity::mapreso {mapreso}"
        f" -in:file:silent {' '.join(refinement_file_basenames)}"
        f" -out:file:silent {output_file}"
    )

    input_files = {}
    if not os.path.isfile(output_final_file) or options["redo"]:

        for file_bn_in, file_in in zip(refinement_file_basenames, refinement_file_names):
            input_files[file_bn_in] = zlib.compress(open(file_in, "rb").read())

        input_files[pdb_no_path] = zlib.compress(open(pdb, "rb").read())
        input_files_future = await standardclient.submit(dict, input_files)
        result_files = await standardclient.submit(
            _run_combine_refinement, runscript.split(), [output_file], input_files_future, uuid.uuid4().hex
        )
        with open(output_final_file, "wb") as fh:
            fh.write(zlib.decompress(result_files[output_file]))

        output_log_fn = os.path.join(config.outdir, f"log_out/{pdb_no_extension}_combine_refinement.log")
        with open(output_log_fn, "wb") as fh:
            fh.write(zlib.decompress(result_files["out.log"]))

    return output_final_file


async def run_dgdp_on_pdb(pdb_file: str, options: Dict[str, Any], config: DGDPConfig, client: Client) -> str:
    print("getting points for", pdb_file)
    points_to_search_file = await run_get_points(pdb_file, config, options, client)
    print("searching points for", pdb_file)
    search_result_file_names = await run_search(points_to_search_file, pdb_file, config, options, client)
    print("combining search for", pdb_file)
    combined_search_result_file_name = await run_combine_search(
        pdb_file, search_result_file_names, config, options, client
    )
    print("running refinement for", pdb_file)
    refinement_file_names = await run_refinement(pdb_file, combined_search_result_file_name, config, options, client)
    print("combining refinement for", pdb_file)
    final_result_silent_file_name = await run_combine_refinement(
        pdb_file, refinement_file_names, config, options, client
    )
    print("done with", pdb_file)
    return final_result_silent_file_name


def _run_cluster_silent(
    runcommand: List[str], expected_outputs: List[str], input_files: Dict[str, bytes], task_id: str
):
    result_files = run_subprocess(runcommand, expected_outputs, input_files, task_id, save_logs=True)
    return result_files


async def run_cluster_silent(
    dgdp_result_filenames: List[str], config: DGDPConfig, options: Dict[str, Any], client: Client
) -> str:
    pathlib.Path(config.outdir).mkdir(exist_ok=True, parents=True)

    # sort of a hack [0] we will see how it works
    mapreso = 2
    if options["refinement_reso"] is not None:
        mapreso = options["refinement_reso"]
    dgdp_result_basenames = [os.path.basename(x) for x in dgdp_result_filenames]
    output_file = f"{options['final_result_name']}_clustered.silent"
    runscript = (
        f"{config.basic_flags}"
        " -mode cluster_silent"
        f" -out::file::silent {output_file}"
        f" -edensity::mapreso {mapreso}"
        f" -in:file:silent {' '.join(dgdp_result_basenames)}"
    )
    output_final_file = os.path.join(config.outdir, output_file)

    if not os.path.isfile(output_final_file):
        input_files = {}
        for input_file_bn, input_file in zip(dgdp_result_basenames, dgdp_result_filenames):
            input_files[input_file_bn] = zlib.compress(open(input_file, "rb").read())

        result_files = await client.submit(
            _run_cluster_silent, runscript.split(), [output_file], input_files, uuid.uuid4().hex
        )

        with open(output_final_file, "wb") as fh:
            fh.write(zlib.decompress(result_files[output_file]))

        logdir = os.path.join(config.outdir, "log_out")
        pathlib.Path(logdir).mkdir(exist_ok=True, parents=True)
        output_log_fn = os.path.join(logdir, f"cluster_{options['final_result_name']}.log")
        with open(output_log_fn, "wb") as fh:
            fh.write(zlib.decompress(result_files["out.log"]))

    return output_final_file


async def dgdp_call_1(options: Dict[str, Any], config: DGDPConfig, client: Client) -> List[str]:
    logdir = os.path.join(config.outdir, "log_out")
    pathlib.Path(logdir).mkdir(exist_ok=True, parents=True)

    dgdp_tasks = []
    for pdb in options["in:file:s"]:
        dgdp_tasks.append(asyncio.create_task(run_dgdp_on_pdb(pdb, options, config, client)))
    all_final_result_files = await asyncio.gather(*dgdp_tasks)
    print("running cluster silent")
    return list(all_final_result_files)


def check_docking_inputs(docking_dict: Dict[str, Any]) -> None:
    pdbs = docking_dict["in:file:s"]
    if len(pdbs) == 0:
        raise RuntimeError('you must set the "pdbs" value in the docking dictionary')

    final_result_names = docking_dict["final_result_names"]
    if len(final_result_names) != len(pdbs):
        raise RuntimeError('you must set the "final_result_names" value to be the same length as the "pdbs" value')

    if "multi_natives" in docking_dict:
        multi_natives = docking_dict["multi_natives"]
        if multi_natives:
            if len(multi_natives) != len(pdbs):
                raise RuntimeError(
                    'you must set the "multi_natives" value to bethe same length as the "pdbs" value' " (can be empty)"
                )
        else:
            del docking_dict["multi_natives"]


async def setup_docking_jobs_2(docking_dict: Dict[str, Any], config: DGDPConfig, i: int, client: Client):
    pdb_group = docking_dict["in:file:s"][i]
    multi_natives = []
    if "multi_natives" in docking_dict:
        multi_natives = docking_dict["multi_natives"][i]
    final_result_name = docking_dict["final_result_names"][i]

    current_options = copy.deepcopy(docking_dict)
    del current_options["in:file:s"]
    if "multi_natives" in current_options:
        del current_options["multi_natives"]
    del current_options["final_result_names"]
    del current_options["mapfiles"]

    per_map_tasks = []
    for mapfile in docking_dict["mapfiles"]:
        current_map_options = copy.deepcopy(current_options)
        current_map_options["in:file:s"] = pdb_group
        if multi_natives:
            current_map_options["multi_native"] = " ".join(multi_natives)
        current_map_options["final_result_name"] = final_result_name
        current_map_options["mapfile"] = mapfile
        mapfile_noext = os.path.splitext(os.path.basename(mapfile))[0]
        current_output_dir = os.path.join(config.outdir, f"{final_result_name}-{mapfile_noext}")

        c_conf = copy.deepcopy(config)
        c_conf.outdir = current_output_dir
        c_conf.set_basic_flags(current_map_options)

        per_map_tasks.append(asyncio.create_task(dgdp_call_1(current_map_options, c_conf, client)))
    per_map_tasks = await asyncio.gather(*per_map_tasks)
    all_result_files = [y for x in per_map_tasks for y in x]

    final_cluster_options = copy.deepcopy(current_options)

    final_cluster_options["in:file:s"] = pdb_group
    if multi_natives:
        final_cluster_options["multi_native"] = " ".join(multi_natives)
    final_cluster_options["final_result_name"] = final_result_name
    # we don't rescore here this is just so rosetta doesn't choke
    final_cluster_options["mapfile"] = docking_dict["mapfiles"][0]

    final_output_dir = os.path.join(config.outdir, "docking_outputs")
    c_conf = copy.deepcopy(config)
    c_conf.outdir = final_output_dir
    c_conf.set_basic_flags(final_cluster_options)

    final_clustered_results_fn = await run_cluster_silent(all_result_files, c_conf, final_cluster_options, client)
    return final_clustered_results_fn


async def setup_docking_jobs_1(docking_dict: Dict[str, Any], config: DGDPConfig, client: Client):
    check_docking_inputs(docking_dict)
    jobs = []
    print(docking_dict["in:file:s"])
    for i, pdb_group in enumerate(docking_dict["in:file:s"]):
        output_pdb_directory = os.path.join(config.outdir, "input_pdbs")
        pathlib.Path(output_pdb_directory).mkdir(exist_ok=True, parents=True)
        for j, pdb_fn in enumerate(pdb_group):
            basename = os.path.basename(pdb_fn)
            shutil.copyfile(pdb_fn, os.path.join(output_pdb_directory, f"{i:02}_{j:02}_{basename}"))
        jobs.append(asyncio.create_task(setup_docking_jobs_2(docking_dict, config, i, client)))
    jobs = await asyncio.gather(*jobs)
    return jobs


async def dgdp_from_commandline(args: argparse.Namespace):
    """
    """
    # Start slurm interface
    async with SLURMCluster(
        cores=1,
        memory=f"{args.memory}GB",
        queue=args.queue,
        job_name=args.name,
        processes=1,
        job_cpu=1,
        walltime=args.walltime,
        extra=["--no-nanny", "--no-bokeh"],
        asynchronous=True,
    ) as cluster:
        async with Client(cluster, asynchronous=True) as client:
            print("Server info:", client, client.scheduler_info(), file=sys.stderr)

            cwd = os.getcwd()
            config = DGDPConfig(args.dgdp_exe, args.rosetta_database, cwd)
            args.pdbs = [args.pdbs]
            options = clean_up_args(vars(args))

            await setup_docking_jobs_1(options, config, client)


if __name__ == "__main__":
    args = parseargs()
    if args.submit:
        from daimyo.distribution.slurm import ez_run_dask_master

        this_command = ["python", "-u"] + [x for x in sys.argv if x != "--submit"]
        extra = None
        ez_run_dask_master(this_command, args.name, hours=20, gb_mem=10, extra=extra)
    else:
        del args.submit
        loop = asyncio.get_event_loop()
        loop.run_until_complete(dgdp_from_commandline(args))
        loop.close()
