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

from pyrosetta import init_from_file
from pyrosetta.distributed.cluster import get_scores_dict, reproduce

sys.path.insert(0, os.path.dirname(__file__))
try:
    from test_reproducibility import TestReproducibilityMulti
except ImportError as ex:
    raise ImportError(ex)
test_suite = globals().get("TestReproducibilityMulti")


def get_protocols(*protocol_names):
    """Get original user-provided PyRosetta protocols from source code."""
    test_case = test_suite.test_reproducibility_from_reproduce
    source_code = textwrap.dedent(inspect.getsource(test_case))
    source_code_lines = source_code.splitlines()
    for node in ast.walk(ast.parse(source_code)):
        if isinstance(node, ast.FunctionDef) and node.name in protocol_names:
            exec(textwrap.dedent(os.linesep.join(source_code_lines[node.lineno - 1: node.end_lineno])))
    _locals = locals()
    protocols = list(map(_locals.get, protocol_names))

    return protocols


def main(input_file, scorefile_name, input_init_file, sequence):
    """Reproduce decoy from .pdb.bz2 file with a '.init' file."""
    skip_corrections = False # Do not skip corrections since not using results for another reproduction
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Initialize PyRosetta on the host node before instantiating `input_pose`
        init_from_file(
            input_init_file,
            output_dir=os.path.join(tmp_dir, "pyrosetta_init_input_files"),
            dry_run=False,
            skip_corrections=skip_corrections,
            relative_paths=True,
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
        protocol_names = scores_dict["metadata"]["protocols"]
        protocols = get_protocols(*protocol_names)
        # Reproduce
        reproduce(
            input_file=input_file,
            scorefile=None,
            decoy_name=None,
            protocols=protocols,
            input_packed_pose=input_pose,
            client=None,
            instance_kwargs={
                "sha1": None,
                "scorefile_name": scorefile_name,
                "output_init_file": os.path.join(tmp_dir, "pyrosetta.init"), # Test `dump_init_file` with custom path
            },
            input_init_file=None, # Skip `init_from_file` since 'input_packed_pose' is provided
            skip_corrections=skip_corrections,
        )


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--scorefile_name', type=str)
    parser.add_argument('--input_init_file', type=str)
    parser.add_argument('--sequence', type=str)
    args = parser.parse_args()
    main(args.input_file, args.scorefile_name, args.input_init_file, args.sequence)
