"""
PyRosettaCluster reproducibility unit test suite using the `unittest` framework.
"""
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

import glob
import itertools
import os
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import random
import sys
import tempfile
import unittest

try:
    import toolz
except ImportError:
    print(
        "Importing 'pyrosetta.tests.distributed.cluster.test_reproducibility_task_updates' requires the "
        + "third-party package 'toolz' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/toolz/\n",
        flush=True,
    )
    raise

from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    reproduce,
)
from pyrosetta.distributed.cluster.io import secure_read_pickle
from pyrosetta.tests.distributed.cluster.unittest_utils import score_function_is_available


class TestReproducibilityTaskUpdates(unittest.TestCase):
    """Test case for decoy reproducibility with user-defined task dictionary updates."""

    def setUp(self):
        random.seed(111)
        self.tmpdir = tempfile.TemporaryDirectory()
        self.workdir = self.tmpdir.name

    def tearDown(self):
        self.tmpdir.cleanup()
        sys.stdout.flush()

    def reproduce_task_updates(self, norm_task_options=False, with_init_file=False, verbose=False):
        """
        Test for PyRosettaCluster decoy reproducibility with updated task dictionaries
        per user-provided PyRosetta protocol.
        """

        pyrosetta.distributed.init(
            options=f"-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )

        # Protocol number to Rosetta command-line options map
        protocol_options = [
            "-score:weights ref2015",
            "-beta_nov16 1",
            "-beta_nov16_cart 1",
        ]
        if score_function_is_available("beta_jan25"):
            protocol_options.append("-beta_jan25 1")
        else:
            protocol_options.append("-score:weights ref2015")
        protocol_options = {str(k): v for k, v in enumerate(protocol_options, start=0)} # Make JSON-serializable
        for k, v in toolz.dicttoolz.keymap(int, protocol_options).items():
            self.assertFalse(
                (v == "-beta_jan25 1") and (protocol_options.get(str(k + 1), None) == "-score:weights ref2015"),
                msg="Cannot transport `PackedPose` objects (with Cys residues) from 'beta_jan25' to 'ref2015' scorefunction environments."
            )
        if verbose:
            print(f"Protocol options: {protocol_options}", flush=True)

        def create_tasks(verbose=verbose):
            custom_options = protocol_options["0"]
            constant_options = "-out:level 200 -multithreading:total_threads 1" # Constant flags for each protocol
            # For stability of concatenation of the values for the 'options' and 'extra_options' keys, if updating Rosetta command-line flags
            # between PyRosetta protocols, either:
            # (1) exclude the 'options' key (using the default '-ex1 -ex2aro' flags of pyrosetta.init) and only set a value for the 'extra_options' key
            # (2) or set the value of the 'options' key to an empty string ('') and only set a value for the 'extra_options' key
            # (3) or exclude the 'extra_options' key, and only set a value for the 'options' key
            options = ""
            extra_options = f"{custom_options} {constant_options}"
            if verbose:
                print(f"Input task 'options' value: '{options}'", flush=True)
                print(f"Input task 'extra_options' value: '{extra_options}'", flush=True)
            yield {
                "options": options,
                "extra_options": extra_options,
                "set_logging_handler": "logging",
                "protocol_options": protocol_options,
                "verbose": verbose,
                "test_key": None,
            }

        def my_protocol(packed_pose, **kwargs):
            import pyrosetta
            import random

            from pyrosetta.distributed.cluster.init_files import PackedPoseHasher

            # Initialize seed
            seed = kwargs["PyRosettaCluster_seed"]
            random.seed(seed)
            # Unpack scores
            cache = packed_pose.pose.cache
            protocol_scorefxn_names = cache.get("protocol_scorefxn_names", None) or {}
            protocol_total_scores = cache.get("protocol_total_scores", None) or {}
            protocol_n_res = cache.get("protocol_n_res", None) or {}
            protocol_sequence = cache.get("protocol_sequence", None) or {}
            protocol_pose_hash = cache.get("protocol_pose_hash", None) or {}
            protocol_pose_hash_cache = cache.get("protocol_pose_hash_cache", None) or {}
            protocol_pose_hash_cache_comments = cache.get("protocol_pose_hash_cache_comments", None) or {}
            protocol_seeds = cache.get("protocol_seeds", None) or {}
            protocol_random_states = cache.get("protocol_random_states", None) or {}
            # Add current values
            scorefxn = pyrosetta.get_score_function()
            protocol_number = kwargs["PyRosettaCluster_protocol_number"]
            protocol_scorefxn_names[protocol_number] = scorefxn.get_name()
            protocol_total_scores[protocol_number] = scorefxn(packed_pose.pose)
            protocol_n_res[protocol_number] = packed_pose.pose.size()
            protocol_sequence[protocol_number] = packed_pose.pose.sequence()
            protocol_pose_hash[protocol_number] = hash(PackedPoseHasher(packed_pose, include_cache=False, include_comments=False).digest())
            protocol_pose_hash_cache[protocol_number] = hash(PackedPoseHasher(packed_pose, include_cache=True, include_comments=False).digest())
            protocol_pose_hash_cache_comments[protocol_number] = hash(PackedPoseHasher(packed_pose, include_cache=True, include_comments=True).digest())
            protocol_seeds[protocol_number] = pyrosetta.rosetta.numeric.random.rg().get_seed()
            protocol_random_states[protocol_number] = hash(str(random.getstate()))
            # Update scores
            packed_pose = packed_pose.update_scores(
                protocol_scorefxn_names=protocol_scorefxn_names,
                protocol_total_scores=protocol_total_scores,
                protocol_n_res=protocol_n_res,
                protocol_sequence=protocol_sequence,
                protocol_pose_hash=protocol_pose_hash,
                protocol_pose_hash_cache=protocol_pose_hash_cache,
                protocol_pose_hash_cache_comments=protocol_pose_hash_cache_comments,
                protocol_seeds=protocol_seeds,
                protocol_random_states=protocol_random_states,
            )
            # Maybe update task dictionary for next PyRosetta protocol, which gets validated and kept
            if protocol_number + 1 < len(kwargs["protocol_options"]):
                kwargs["options"] = ""
                kwargs["extra_options"] = kwargs["protocol_options"][str(protocol_number + 1)]

            # Test updating another existing key, which gets validated and kept
            assert "test_key" in kwargs, (
                f"Existing key not found in keyword arguments dictionary: 'test_key'"
            )
            kwargs["test_key"] = random.choice(
                ["foo", "bar", "baz", "qux", "quux", frozenset([1, 2, 3]), complex(8, 9), True, False]
            )
            # Test adding a new key to current task dictionary, which gets validated and kept
            key_template = "foobar_{0}"
            for i in range(protocol_number - 1, -1, -1):
                assert key_template.format(i) in kwargs, (
                    f"Previously added key not found in keyword arguments dictionary: '{key_template.format(i)}'"
                )
            kwargs[key_template.format(protocol_number)] = "baz"
            # Test adding a reserved key to current task dictionary, which gets automatically removed
            reserved_key = "PyRosettaCluster_foobar"
            assert reserved_key not in kwargs, (
                f"Reserved key found in keyword arguments dictionary: '{reserved_key}'"
            )
            kwargs[reserved_key] = "baz"
            # Test removing reserved keys from current task dictionary, which get automatically added back
            kwargs.pop("PyRosettaCluster_task")
            kwargs.pop("PyRosettaCluster_client_repr")
            kwargs.pop("PyRosettaCluster_datetime_start")
            # Test updating reserved keys in current task dictionary, which get automatically reverted
            kwargs["PyRosettaCluster_protocol_name"] = "xyzzy"
            kwargs["PyRosettaCluster_protocol_number"] = random.getrandbits(14 * 8).to_bytes(14, "little")
            # Maybe print
            if kwargs["verbose"]:
                print(f"PyRosetta protocol number {protocol_number} Pose.cache:", packed_pose.pose.cache, flush=True)

            return packed_pose, kwargs

        pyrosetta.secure_unpickle.add_secure_package("pandas")
        pyrosetta.secure_unpickle.add_secure_package("pyarrow")
        sequence = "ADEFGHIKLMNPQRSTVWY"
        input_pose = io.to_pose(io.pose_from_sequence(sequence))
        scorefile_name = "test_scores.json"
        compressed = True
        protocols = [my_protocol] * len(protocol_options)
        output_path = os.path.join(self.workdir, f"outputs_norm_task_options_{int(norm_task_options)}")
        if verbose:
            print(f"Running original simulation for `norm_task_options={norm_task_options}` and `with_init_file={with_init_file}`.", flush=True)
        PyRosettaCluster(
            tasks=create_tasks,
            input_packed_pose=input_pose,
            output_path=output_path,
            scratch_dir=self.workdir,
            simulation_records_in_scorefile=True,
            scorefile_name=scorefile_name,
            compressed=compressed,
            norm_task_options=norm_task_options,
            sha1=None,
            project_name="PyRosettaCluster",
            simulation_name="update_tasks",
            output_decoy_types=[".pdb", ".init"],
            output_scorefile_types=[".json", ".xz"],
            output_init_file=None,
            min_workers=1,
            max_workers=2,
            cooldown_time=2.5,
        ).distribute(protocols=protocols)

        scorefile_path = os.path.join(output_path, os.path.splitext(scorefile_name)[0] + ".xz")
        df = secure_read_pickle(scorefile_path, compression="infer")
        self.assertEqual(df.index.size, 1)
        original_record = df.iloc[0]

        if verbose:
            print("Original record:", original_record["instance"], sep="\n", flush=True)
        if verbose:
            logging_file = original_record["metadata"]["logging_file"]
            _logging_files = glob.glob(os.path.join(os.path.dirname(logging_file), "*"))
            for _logging_file in _logging_files:
                print("Output from logging file:", _logging_file, flush=True)
                with open(_logging_file, "r") as f:
                    print(f.read(), flush=True)

        # Reproduce decoy
        if verbose:
            print(f"Original simulation complete for `norm_task_options={norm_task_options}` and `with_init_file={with_init_file}`.", flush=True)
            print(f"Running reproduction simulation for `norm_task_options={norm_task_options}` and `with_init_file={with_init_file}`.", flush=True)
        if with_init_file: # Reproduce from a PyRosetta initialization file
            input_file = original_record["metadata"]["output_file"].split(os.extsep)[0] + (".init.bz2" if compressed else ".init")
            scorefile = None
            decoy_name = None
            input_packed_pose = None
        else: # Reproduce from a pickled `pandas.DataFrame` scorefile
            input_file = None
            scorefile = scorefile_path
            decoy_name = original_record["metadata"]["decoy_name"]
            input_packed_pose = input_pose
        output_path = os.path.join(self.workdir, f"outputs_reproduce_norm_task_options_{int(norm_task_options)}_{with_init_file}")
        reproduce_scorefile_name = f"reproduce_test_scores_{with_init_file}.json"
        skip_corrections = False
        reproduce_kwargs = dict(
            input_file=input_file,
            scorefile=scorefile,
            decoy_name=decoy_name,
            input_packed_pose=input_packed_pose,
            protocols=protocols,
            client=None,
            clients=None,
            instance_kwargs={
                "output_path": output_path,
                "scratch_dir": self.workdir,
                "simulation_records_in_scorefile": True,
                "scorefile_name": reproduce_scorefile_name,
                "compressed": compressed,
                "norm_task_options": norm_task_options,
                "sha1": None,
                "output_decoy_types": [".pdb", ".init"],
                "output_scorefile_types": [".json", ".xz"],
                "output_init_file": None, # Skip `dump_init_file`
                "min_workers": 1,
                "max_workers": 2,
                "cooldown_time": 2.5,
            },
            skip_corrections=skip_corrections,
            init_from_file_kwargs=dict(
                dry_run=None,
                output_dir=os.path.join(self.workdir, f"reproduce_pyrosetta_init_files_{with_init_file}"),
                skip_corrections=skip_corrections,
                relative_paths=None,
                max_decompressed_bytes=1_000_000,
                restore_rg_state=None,
                database=None,
                verbose=verbose,
                set_logging_handler=None,
                notebook=None,
                silent=None,
            ),
        )
        if with_init_file:
            with self.assertWarns(UserWarning) as cm:
                # PyRosetta is already initialized on the head node process but we reproduce from a PyRosetta initialization file
                reproduce(**reproduce_kwargs)
            self.assertTrue(
                str(cm.warning).startswith(
                    "Skipping ScoreFunction corrections for the PyRosettaCluster task"
                    if skip_corrections
                    else "Preserving ScoreFunction corrections for the PyRosettaCluster task"
                )
            )
        else:
            reproduce(**reproduce_kwargs)

        reproduce_scorefile_path = os.path.join(output_path, os.path.splitext(reproduce_scorefile_name)[0] + ".xz")
        df_reproduce = secure_read_pickle(reproduce_scorefile_path, compression="infer")
        self.assertEqual(df_reproduce.index.size, 1)
        reproduce_record = df_reproduce.iloc[0]

        # Assert identical protocol results across original versus reproduction
        for protocol_number in range(len(protocols)):
            self.assertEqual(
                original_record["scores"]["protocol_scorefxn_names"][protocol_number],
                reproduce_record["scores"]["protocol_scorefxn_names"][protocol_number],
                msg=f"Protocol number {protocol_number} score function names differ.",
            )
            self.assertEqual(
                original_record["scores"]["protocol_total_scores"][protocol_number],
                reproduce_record["scores"]["protocol_total_scores"][protocol_number],
                msg=f"Protocol number {protocol_number} total scores differ.",
            )
            self.assertEqual(
                original_record["scores"]["protocol_n_res"][protocol_number],
                reproduce_record["scores"]["protocol_n_res"][protocol_number],
                msg=f"Protocol number {protocol_number} number of residues differ.",
            )
            self.assertEqual(
                original_record["scores"]["protocol_sequence"][protocol_number],
                reproduce_record["scores"]["protocol_sequence"][protocol_number],
                msg=f"Protocol number {protocol_number} sequences differ.",
            )
            self.assertEqual(
                original_record["scores"]["protocol_pose_hash"][protocol_number],
                reproduce_record["scores"]["protocol_pose_hash"][protocol_number],
                msg=f"Protocol number {protocol_number} `Pose` hashes differ.",
            )
            self.assertEqual(
                original_record["scores"]["protocol_pose_hash_cache"][protocol_number],
                reproduce_record["scores"]["protocol_pose_hash_cache"][protocol_number],
                msg=f"Protocol number {protocol_number} `Pose` hashes (with `Pose.cache`) differ.",
            )
            self.assertEqual(
                original_record["scores"]["protocol_pose_hash_cache_comments"][protocol_number],
                reproduce_record["scores"]["protocol_pose_hash_cache_comments"][protocol_number],
                msg=f"Protocol number {protocol_number} `Pose` hashes (with `Pose.cache` and `Pose` comments) differ.",
            )
            self.assertEqual(
                original_record["scores"]["protocol_seeds"][protocol_number],
                reproduce_record["scores"]["protocol_seeds"][protocol_number],
                msg=f"Protocol number {protocol_number} PyRosetta RNG seeds differ.",
            )
            self.assertEqual(
                original_record["scores"]["protocol_random_states"][protocol_number],
                reproduce_record["scores"]["protocol_random_states"][protocol_number],
                msg=f"Protocol number {protocol_number} `random` module states differ.",
            )
        # Assert unique scorefunction results within original and reproduction
        for record in (original_record, reproduce_record):
            scorefxn_names = record["scores"]["protocol_scorefxn_names"]
            total_scores = record["scores"]["protocol_total_scores"]
            scorefxn_name_to_total_score_dict = {scorefxn_names[k]: total_scores[k] for k in scorefxn_names}
            for (name_i, total_score_i), (name_j, total_score_j) in itertools.combinations(scorefxn_name_to_total_score_dict.items(), 2):
                self.assertNotEqual(
                    total_score_i,
                    total_score_j,
                    msg=f"Scorefunctions '{name_i}' and '{name_j}' resulted in identical total score: {total_score_i}",
                )
        # Assert identical sequence through original and reproduction simulations
        for record in (original_record, reproduce_record):
            self.assertEqual(
                io.pose_from_file(record["metadata"]["output_file"]).pose.sequence(),
                sequence,
                msg="Pose sequence diverged."
            )
        if verbose:
            print(f"Successfully validated reproduction simulation for `norm_task_options={norm_task_options}` and `with_init_file={with_init_file}`.", flush=True)

    def test_reproduce_task_updates(self):
        return self.reproduce_task_updates(norm_task_options=False, with_init_file=False)

    def test_reproduce_task_updates_norm_task_options(self):
        return self.reproduce_task_updates(norm_task_options=True, with_init_file=False)

    def test_reproduce_task_updates_with_init_file(self):
        return self.reproduce_task_updates(norm_task_options=False, with_init_file=True)

    def test_reproduce_task_updates_norm_task_options_with_init_file(self):
        return self.reproduce_task_updates(norm_task_options=True, with_init_file=True)
