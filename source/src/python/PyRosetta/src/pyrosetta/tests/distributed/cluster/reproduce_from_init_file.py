# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import argparse
import ast
import inspect
import os
import pyrosetta.distributed.io as io
import sys
import tempfile
import textwrap
import types

from functools import partial
from pyrosetta import init_from_file
from pyrosetta.distributed.cluster import get_scores_dict, reproduce

sys.path.insert(0, os.path.dirname(__file__))
try:
    from test_reproducibility import TestReproducibilityMulti
except ImportError as ex:
    raise ImportError(ex)
test_suite = globals().get("TestReproducibilityMulti")


def get_protocols(scores_dict):
    """Get PyRosettaCluster protocols from the source code."""
    protocol_names = scores_dict["metadata"]["protocols"]
    test_case = test_suite.test_reproducibility_from_reproduce
    source_code = textwrap.dedent(inspect.getsource(test_case))
    source_code_lines = source_code.splitlines()
    module = types.ModuleType("_protocols")
    for node in ast.walk(ast.parse(source_code)):
        if isinstance(node, ast.FunctionDef) and node.name in protocol_names:
            exec(
                textwrap.dedent(os.linesep.join(source_code_lines[node.lineno - 1: node.end_lineno])),
                module.__dict__,
            )
    _get_protocol = partial(getattr, module)
    protocols = list(map(_get_protocol, protocol_names))

    return protocols


def reproduce_init_from_file_test(input_file, scorefile_name, input_init_file, sequence):
    """Reproduce decoy from '.pdb.bz2' file with a '.init' file."""
    skip_corrections = False  # Do not skip corrections since not using results for another reproduction
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Initialize PyRosetta on the host node before instantiating `input_pose`
        io.init_from_file(
            input_init_file,
            output_dir=os.path.join(tmp_dir, "reproduce_init_from_file_test"),
            dry_run=False,
            skip_corrections=skip_corrections,
            relative_paths=False,
            max_decompressed_bytes=100_000,
            database=None,
            verbose=True,
            set_logging_handler="logging",
            notebook=None,
            silent=False,
        )
        # Instantiate original input pose
        input_pose = io.to_pose(io.pose_from_sequence(sequence))
        # Get protocols
        scores_dict = get_scores_dict(input_file)
        protocols = get_protocols(scores_dict)
        # Reproduce
        output_init_file = os.path.join(tmp_dir, "pyrosetta.init")
        compressed = False
        reproduce(
            input_file=input_file,
            scorefile=None,
            decoy_name=None,
            protocols=protocols,
            input_packed_pose=input_pose,
            clients=None,
            instance_kwargs={
                "sha1": None,
                "scorefile_name": scorefile_name,
                "output_init_file": output_init_file, # Test `IO._dump_init_file` with custom path
                "compressed": compressed,
            },
            skip_corrections=skip_corrections,
            init_from_file_kwargs={},
        )
        # Test output '.init' file compression
        if compressed:
            output_init_file += ".bz2"
        assert os.path.isfile(output_init_file)


def reproduce_test(input_file, scorefile_name, input_init_file, skip_corrections, init_from_file_skip_corrections):
    """Reproduce decoy from a '.init' file."""
    print(
        "Running test case:",
        f"`reproduce(skip_corrections={skip_corrections}, "
        + f"init_from_file_kwargs=dict(skip_corrections={init_from_file_skip_corrections}))`",
    )
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Get protocols
        scores_dict = get_scores_dict(input_file)
        protocols = get_protocols(scores_dict)
        # Reproduce
        output_init_file = os.path.join(tmp_dir, "pyrosetta.init")
        compressed = True
        init_from_file_kwargs = dict(
            output_dir=os.path.join(tmp_dir, "reproduce_test"),
            dry_run=False,
            skip_corrections=init_from_file_skip_corrections,
            relative_paths=False,
            max_decompressed_bytes=100_000,
            database=None,
            verbose=True,
            set_logging_handler="logging",
            notebook=None,
            silent=False,
        )
        reproduce(
            input_file=input_init_file,
            scorefile=None,
            decoy_name=None,
            protocols=protocols,
            input_packed_pose=None,
            clients=None,
            instance_kwargs={
                "sha1": None,
                "scorefile_name": scorefile_name,
                "output_init_file": output_init_file, # Test `IO._dump_init_file` with custom path
                "compressed": compressed,
            },
            skip_corrections=skip_corrections,
            init_from_file_kwargs=init_from_file_kwargs,
        )
        # Test `get_scores_dict()` with '.init' file syntax after PyRosetta initialization
        _scores_dict = get_scores_dict(input_init_file)
        assert scores_dict == _scores_dict, f"Scores dictionaries differ: {scores_dict} != {_scores_dict}"
        # Test output '.init' file compression
        if compressed:
            output_init_file += ".bz2"
        assert os.path.isfile(output_init_file)


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--scorefile_name', type=str)
    parser.add_argument('--input_init_file', type=str)
    parser.add_argument('--sequence', type=str)
    parser.add_argument('--test_case', type=int)
    args = parser.parse_args()
    if args.test_case == 0:
        reproduce_init_from_file_test(args.input_file, args.scorefile_name, args.input_init_file, args.sequence)
    elif args.test_case == 1:
        reproduce_test(args.input_file, args.scorefile_name, args.input_init_file, False, False)
    elif args.test_case == 2:
        reproduce_test(args.input_file, args.scorefile_name, args.input_init_file, True, True)
    elif args.test_case == 3:
        reproduce_test(args.input_file, args.scorefile_name, args.input_init_file, True, False)
    elif args.test_case == 4:
        reproduce_test(args.input_file, args.scorefile_name, args.input_init_file, False, True)
