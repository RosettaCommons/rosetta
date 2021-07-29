"""
A collection of ways to get input from users via an interactive prompt.
"""

from typing import List, Tuple, Dict, Optional
import os
import sys
import argparse


def save_args() -> None:
    x = 1
    previous = None
    while True:
        out_fn = f"saved_cmdline_{x:03}.cmd"
        if os.path.isfile(out_fn):
            x += 1
            previous = out_fn
            continue
        out_args_string = " ".join(sys.argv)
        previous_text = None
        if previous:
            with open(previous) as fh:
                previous_text = fh.read()
        if previous_text and (previous_text.strip() == out_args_string.strip()):
            break
        with open(out_fn, "w") as fh:
            fh.write(" ".join(sys.argv))
        break


def add_default_jobdist_arguments(parser) -> argparse.ArgumentParser:
    parser.add_argument("-na", "--name", help="name for this job", required=True)
    parser.add_argument("--submit", help="Submit to scheduler", action="store_true", default=False)
    parser.add_argument("--queue", help="what queue to submit workers to", required=False, type=str, default="medium")

    parser.add_argument(
        "--walltime", help="what walltime for workers to have", required=False, type=str, default="12:00:00"
    )

    parser.add_argument("--memory", help="what memory (GB) for workers to have", required=False, type=int, default=5)
    return parser


def add_default_jobdist_high_low_arguments(parser):
    parser.add_argument("-na", "--name", help="name for this job", required=True)
    parser.add_argument("--submit", help="Submit to scheduler", action="store_true", default=False)
    parser.add_argument(
        "--high-memory-queue",
        help="""what queue to submit high memory worker to (for fragment picking)""",
        required=False,
        type=str,
        default="medium",
    )

    parser.add_argument(
        "--low-memory-queue",
        help="what queue to submit low memory worker to",
        required=False,
        type=str,
        default="medium",
    )

    parser.add_argument(
        "--high-memory-walltime",
        help="what walltime for high memory workers to have",
        required=False,
        type=str,
        default="12:00:00",
    )

    parser.add_argument(
        "--low-memory-walltime",
        help="what walltime for low memory workers to have",
        required=False,
        type=str,
        default="12:00:00",
    )

    parser.add_argument(
        "--high-memory-memory",
        help="what memory (GB) for high memory workers to have",
        required=False,
        type=int,
        default=5,
    )

    parser.add_argument(
        "--low-memory-memory",
        help="what memory (GB) for low memory workers to have",
        required=False,
        type=int,
        default=3,
    )
    return parser


def check_rosetta_bin_directory(bin_directory: str, exes_to_check: List[str]) -> str:
    found = [False for _ in exes_to_check]
    files_and_folders = os.listdir(bin_directory)
    for x in files_and_folders:
        for i, y in enumerate(exes_to_check):
            if found[i]:
                continue
            if x.startswith(y):
                found[i] = True
        if all(found):
            break

    if not all(found):
        unfound_exes = []
        for i, x in enumerate(found):
            if not x:
                unfound_exes.append(exes_to_check[i])
        raise ValueError(f"Unable to find '{' '.join(unfound_exes)}' in {bin_directory}.")
    return bin_directory


def check_rosetta_bin_directory_defaults(bin_directory: str):
    return check_rosetta_bin_directory(bin_directory, ["rosetta_scripts", "grower_prep", "extract_pdbs"])


def check_rosetta_database_location(database_directory: str):
    if not os.path.isdir(database_directory):
        raise ValueError(f"Unable to locate directory {database_directory}.")
    return database_directory


def add_default_rosetta_exes_arguments(parser):
    rosetta_exes_args = parser.add_mutually_exclusive_group(required=True)
    rosetta_exes_args.add_argument(
        "--rosetta_bin_location", help="directory of the rosetta executables", type=check_rosetta_bin_directory_defaults
    )
    rosetta_exes_args.add_argument(
        "--rosetta_bin_location_remote",
        help="directory of the rosetta executables remote so not error checked",
        type=str,
    )

    rosetta_database_args = parser.add_mutually_exclusive_group(required=True)
    rosetta_database_args.add_argument("--rosetta_database", help="location of the rosetta database", type=str)
    rosetta_database_args.add_argument(
        "--rosetta_database_remote", help="location of the rosetta database, remote so not error checked", type=str
    )
    return parser


def get_rosetta_exe_and_database_from_args(
    args: argparse.Namespace, exes: Optional[List[str]] = None
) -> Tuple[Dict[str, str], str]:
    if exes is None:
        exes = ["rosetta_scripts", "grower_prep", "extract_pdbs"]
    ret = {}
    if args.rosetta_bin_location:
        for exe in exes:
            ret[exe] = os.path.join(args.rosetta_bin_location, exe)
    else:
        for exe in exes:
            ret[exe] = os.path.join(args.rosetta_bin_location_remote, exe)

    if args.rosetta_database:
        rosetta_database = args.rosetta_database
    else:
        rosetta_database = args.rosetta_database_remote
    return ret, rosetta_database
