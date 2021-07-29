#!/usr/bin/env python
import argparse
import json
import asyncio
from typing import Dict, Any, List, Union
import copy
import os
import sys
import pathlib


from dask_jobqueue import SLURMCluster
from dask.distributed import Client, LocalCluster

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from dgdp import setup_docking_jobs_1, DGDPConfig, clean_up_args
from cryo_assembly import run_cryo_assembly

from user_input import (
    add_default_rosetta_exes_arguments,
    add_default_jobdist_arguments,
    save_args,
    get_rosetta_exe_and_database_from_args,
)


def example_config():
    return {
        "assembly": {
            "dsm": "",
            "connectivity_weight": 1000.0,
            "extras": "",
            "final_out": 20,
            "n_iter": 200,
            "n_out": 10,
            "null_weight": -1150,
            "onebody_weight": 260.0,
            "twobody_weight": 150.0,
            "dc_weight": 0,
            "overlap_weight": 100,
            "vdw_weight": 1.0,
            "rg_weight": 0,
            "cen_env_smooth_weight": 2.0,
            "cen_pair_smooth_weight": 1.0,
            "hbond_sr_bb_weight": 1,
            "hbond_lr_bb_weight": 1,
            "rama_weight": 0.2,
            "omega_weight": 0.2,
            "cbeta_smooth_weight": 1.0,
            "cenpack_smooth_weight": 1.0,
            "motif_dock_weight": 0,
            "interchain_vdw_weight": 0,
            "neighbor_count_weight": 0,
            "jobs": 50,
        },
        "general": {
            "extras": "",
            "mapfile": "/home/danpf/collaborations/best_collaborators/P2_J418_map_resampled.mrc",
            "mute": "true",
            "name": "vps",
            "rank_range": "1 100",
        },
        "iterover": {"general": {}, "twobody": {}, "assembly": {}},
        "scheduler": {"cores": 400, "name": "vps", "queue": "dimaio", "memory": "12G", "worker_walltime": "6:00:00"},
        "twobody": {
            "scorefile_out": "sc",
            "scorefile": "sc",
            "constrain_type": "legacy",
            "extras": "-hbond_bb_per_residue_energy",
            "bb_iter": "0",
            "rb_iter": "20",
            "refine_dens_wt": "10",
            "core_split": 500,
            "score_closability": "true",
            "slide_domains": "5",
            "legacy_radii": "12",
            "trim_clash": "true",
            "allow_domains_to_minimize": "true",
        },
        "docking": {
            "pdbs": [],  # groups of pdbs which can be docked together
            "mapfiles": [],
            "mapreso": 4.0,
            "multi_natives": [],  # multi_natives must be same length as pdbs but can be empty
            "clust_radius": 10.0,
            "convolute_single_residue": True,
            "laplacian_offset": 0,
            "n_to_search": 20000,
            "n_filtered": 1000,
            "n_output": 200,
            "point_radius": 2.0,
            "min": True,
            "min_bb": False,
            "constrain_refinement": 10.0,
            "refinement_reso": 4.0,
            "rot_middle_ca": False,
            "extras": "",
            "max_rot_per_trans": 11,
            "bw": None,
            "final_result_names": [],  # final_result_name pdbs same length as "pdbs"
            "cores": 4,
        },
    }


def parseargs() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--premade_json", help="The premade json you want to use", required=True)
    parser.add_argument("-o", "--output_dir", help="The directory to output to", required=True)
    parser.add_argument("-do", "--dock_only", help="dock only -- don't assemble", default=False, action="store_true")
    parser.add_argument("-rl", "--run_locally", help="run locally dont use slurm", default=False, action="store_true")
    parser = add_default_rosetta_exes_arguments(parser)
    parser = add_default_jobdist_arguments(parser)
    args = parser.parse_args()
    with open(args.premade_json) as fh:
        args.premade_json = json.load(fh)

    exes, rosetta_db = get_rosetta_exe_and_database_from_args(args, ["dgdp", "cryo_assembly", "extract_pdbs"])
    args.dgdp_exe = exes["dgdp"]
    args.cryo_assembly_exe = exes["cryo_assembly"]
    args.extract_pdbs_exe = exes["extract_pdbs"]
    args.rosetta_database = rosetta_db
    save_args()
    return args


def check_docking_inputs(docking_dict: Dict[str, Any]) -> None:
    pdbs = docking_dict["pdbs"]
    if len(pdbs) == 0:
        raise RuntimeError('you must set the "pdbs" value in the docking dictionary')

    final_result_names = docking_dict["final_result_names"]
    if len(final_result_names) != len(pdbs):
        print(len(final_result_names))
        print(len(pdbs))

        raise RuntimeError('you must set the "final_result_names" value to be the same length as the "pdbs" value')

    multi_natives = docking_dict["multi_natives"]
    if len(multi_natives) != len(pdbs):
        raise RuntimeError(
            'you must set the "multi_natives" value to bethe same length as the "pdbs" value' " (can be empty)"
        )


async def run_docking(
    docking_dict: Dict[str, Any], dgdp_exe: str, rosetta_database: str, output_dir: str, client: Client
) -> List[str]:
    dock_dict = copy.deepcopy(docking_dict)
    check_docking_inputs(dock_dict)

    config = DGDPConfig(dgdp_exe, rosetta_database, output_dir)
    options = clean_up_args(dock_dict)
    config.set_basic_flags(options)

    results = await setup_docking_jobs_1(options, config, client)
    return list(results)


async def run_assembly(
    docking_result_fns: List[str],
    premade_json: Dict[str, Any],
    cryo_assembly_exe: str,
    rosetta_database: str,
    output_dir: str,
    client: Client,
) -> List[str]:
    raise
    if 1 + 1 == 2:
        premade_json["general"]["executable"] = "/home/danpf/Rosetta/cryo_assembly/source/bin/cryo_assembly"
        premade_json["general"]["database"] = "/home/danpf/Rosetta/cryo_assembly/database"
    else:
        print("hacking json for local running")
        print("hacking json for local running")
        print("hacking json for local running")
        print("hacking json for local running")
        premade_json["general"]["executable"] = cryo_assembly_exe
        premade_json["general"]["database"] = rosetta_database

    premade_json["general"]["in:file:silent"] = " ".join(docking_result_fns)
    cryo_results = run_cryo_assembly(premade_json, output_dir, client)
    return cryo_results


async def main(args: argparse.Namespace, client: Client, sync_client: Client) -> None:
    pathlib.Path(args.output_dir).mkdir(exist_ok=True, parents=True)
    with open(os.path.join(args.output_dir, "input.json"), "w") as fh:
        fh.write(json.dumps(args.premade_json))

    docking_job_fns = await run_docking(
        args.premade_json["docking"], args.dgdp_exe, args.rosetta_database, args.output_dir, client
    )
    if args.dock_only:
        print("docking only... finishing up!")
        return
    docking_multi_natives = args.premade_json["docking"].get("multi_natives", [])
    assembly_multi_natives = args.premade_json["cryo_assembly"]["general"].get("multi_natives", [])
    if len(docking_multi_natives) != 0 and len(assembly_multi_natives) == 0:
        args.premade_json["cryo_assembly"]["general"]["multi_native"] = " ".join(
            [y for x in docking_multi_natives for y in x]
        )

    cryo_a_result_filenames = await run_assembly(
        docking_job_fns,
        args.premade_json["cryo_assembly"],
        args.dgdp_exe,
        args.rosetta_database,
        args.output_dir,
        sync_client,
    )
    print("results in", cryo_a_result_filenames)


def get_sync_cluster(args: argparse.Namespace) -> Union[SLURMCluster, LocalCluster]:
    if args.run_locally:
        cluster = LocalCluster(n_workers=1, threads_per_worker=1, processes=True, asynchronous=True)
    else:
        cluster = SLURMCluster(
            cores=1,
            memory=f"{args.memory}GB",
            queue=args.queue,
            job_name=args.name,
            processes=1,
            job_cpu=1,
            walltime=args.walltime,
            extra=["--no-nanny", "--no-bokeh"],
            asynchronous=True,
        )
        cluster.adapt(minimum=0, maximum=999, wait_count=400)
    return cluster


def get_work_cluster(args: argparse.Namespace) -> Union[SLURMCluster, LocalCluster]:
    if args.run_locally:
        cluster = LocalCluster(n_workers=1, threads_per_worker=1, processes=True, asynchronous=True)
    else:
        cluster = SLURMCluster(
            cores=1,
            memory=f"{args.memory}GB",
            queue=args.queue,
            job_name=args.name,
            processes=1,
            job_cpu=1,
            walltime=args.walltime,
            extra=["--no-nanny", "--no-bokeh"],
            asynchronous=True,
        )
        cluster.adapt(minimum=0, maximum=999, wait_count=400)
    return cluster


async def commandline_main(args: argparse.Namespace) -> None:
    async with get_sync_cluster(args) as sync_cluster:
        async with Client(sync_cluster, asynchronous=True) as sync_client:
            async with get_work_cluster(args) as cluster:
                async with Client(cluster, asynchronous=True) as client:
                    print("Server info:", client, client.scheduler_info(), file=sys.stderr)
                    await main(args, client, sync_client)


if __name__ == "__main__":
    args = parseargs()
    loop = asyncio.get_event_loop()
    loop.run_until_complete(commandline_main(args))
    loop.close()
