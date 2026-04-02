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
import random
import sys
import tempfile
import unittest

try:
    import toolz
except ImportError:
    print(
        "Importing 'pyrosetta.tests.distributed.cluster.test_reproducibility_remodel_task_updates' requires the "
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


class TestReproducibilityRemodelTaskUpdates(unittest.TestCase):
    """Test case for decoy reproducibility with user-defined task dictionary updates with RosettaRemodel."""

    _dummy_seq: str = "G"

    def setUp(self):
        random.seed(111)
        self.tmpdir = tempfile.TemporaryDirectory()
        self.workdir = self.tmpdir.name
        self.cwd = os.getcwd()
        self.remodel_1_pdb_file = os.path.join(self.cwd, "1.pdb")
        self.assertFalse(os.path.isfile(self.remodel_1_pdb_file))
        os.chdir(self.workdir) # Prevent RosettaRemodel dumping '1.pdb' in current working directory

    def tearDown(self):
        if os.path.isfile(self.remodel_1_pdb_file):
            try:
                os.remove(self.remodel_1_pdb_file)
            except Exception:
                print(f"Warning: RosettaRemodel output file still exists: {self.remodel_1_pdb_file}", flush=True)
        os.chdir(self.cwd)
        self.tmpdir.cleanup()
        sys.stdout.flush()

    @staticmethod
    def format_stable_float(x):
        dig = sys.float_info.dig # 15
        return format(x, f".{dig}g")

    @staticmethod
    def get_remodel_options():
        return {
            "-remodel:num_trajectory": "1",
            "-remodel:quick_and_dirty": "1",
            "-remodel:dr_cycles": str(random.randint(1, 5)),
            "-remodel:use_clusters": str(random.choice([0, 1])),
            "-remodel:core_cutoff": str(random.randint(10, 20)),
            "-remodel:boundary_cutoff": str(random.randint(5, 15)),
            "-remodel:generic_aa": random.choice(list("ADEFGHIKLMNPQRSTVWY")),
            "-ex1": str(random.choice([0, 1])),
            "-ex2": str(random.choice([0, 1])),
            "-ex3": str(random.choice([0, 1])),
            "-ex4": str(random.choice([0, 1])),
            "-ex1:level": str(random.choice([0, 1])),
            "-ex2:level": str(random.choice([0, 1])),
            "-ex3:level": str(random.choice([0, 1])),
            "-ex4:level": str(random.choice([0, 1])),
            "-extrachi_cutoff": str(random.choice([0, 18])),
            "-use_input_sc": str(random.choice([0, 1])),
            # For unit tests specifying `norm_task_options=True`, here we serialize double precision Python `float` values for roundtrip stability
            # through Rosetta C++ double precision `core::Real` values, both of which have at most 15 digits of reliable numerical precision.
            # This string formatting ensures that the `pyrosetta.get_init_options` function outputs precisely the same weights that are input here:
            "-remodel:hb_srbb": TestReproducibilityRemodelTaskUpdates.format_stable_float(random.random()),
            "-remodel:vdw": TestReproducibilityRemodelTaskUpdates.format_stable_float(random.random()),
            "-remodel:rama": TestReproducibilityRemodelTaskUpdates.format_stable_float(random.random()),
            "-remodel:cbeta": TestReproducibilityRemodelTaskUpdates.format_stable_float(random.random()),
            "-remodel:cenpack": TestReproducibilityRemodelTaskUpdates.format_stable_float(random.random()),
            "-remodel:rg": TestReproducibilityRemodelTaskUpdates.format_stable_float(random.random()),
        }

    @staticmethod
    def get_random_options():
        return TestReproducibilityRemodelTaskUpdates.get_remodel_options()

    def reproduce_remodel_task_updates(self, norm_task_options=False, with_init_file=False, verbose=False):
        """
        Test for PyRosettaCluster decoy reproducibility with updated task dictionaries
        using RosettaRemodel per user-provided PyRosetta protocol.
        """
        pyrosetta.distributed.init(
            options=f"-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        _n_protocols = 3

        def create_tasks(verbose=verbose):
            custom_options = TestReproducibilityRemodelTaskUpdates.get_random_options()
            constant_options = {
                "-multithreading:total_threads": "1",
                "-out:level": "200",
                "-out:levels": " ".join(
                    [
                        "core.fragment.picking_old.vall.vall_io:100",
                        "protocols.forge.components.VarLengthBuild:100",
                        "core.fragment.picking_old.vall.VallLibrarian:100",
                        "core.fragment.picking_old.vall.eval.IdentityEval:100",
                        "protocols.forge.remodel.RemodelMover:100",
                        "protocols.forge.remodel.RemodelDesignMover:100",
                        "protocols.forge.remodel.RemodelWorkingSet:100",
                        "protocols.forge.remodel.RemodelData:100",
                        "protocols.forge.build.BuildManager:100",
                        "protocols.cluster:100",
                    ],
                )
            }
            # For stability of concatenation of the values for the 'options' and 'extra_options' keys, if updating Rosetta command-line flags
            # between PyRosetta protocols, either:
            # (1) exclude the 'options' key (using the default '-ex1 -ex2aro' flags of pyrosetta.init) and only set a value for the 'extra_options' key
            # (2) or set the value of the 'options' key to an empty string ('') and only set a value for the 'extra_options' key
            # (3) or exclude the 'extra_options' key, and only set a value for the 'options' key
            options = ""
            extra_options = {**custom_options, **constant_options}
            if verbose:
                print(f"Input task 'options' value: '{options}'", flush=True)
                print("Input task 'extra_options' value:", extra_options, flush=True)
            yield {
                "options": options,
                "extra_options": extra_options,
                "set_logging_handler": "logging",
                "verbose": verbose,
                "n_protocols": _n_protocols,
                "constant_options": constant_options,
            }

        def write_remodel_blueprint(packed_pose, workdir, n_res_cterm_ext=4, n_term_threshold=2, verbose=False):
            line_template = "{res} {name1} {ss}{resfile_cmd}"
            pose = packed_pose.pose
            lines = []
            for res in range(1, pose.size() + 1):
                name1 = pose.residue(res).name1()
                ss = random.choice(list("HEL")) if res > max((pose.size() - n_term_threshold), 0) else "."
                resfile_cmd = " NATAA" if ss != "." else ""
                line = line_template.format(res=res, name1=name1, ss=ss, resfile_cmd=resfile_cmd)
                lines.append(line)
            for _ in range(n_res_cterm_ext):
                res = 0
                name1 = "X"
                ss = random.choice(list("HEL"))
                resfile_cmd = " PIKAA " + random.choice(list("ADEFGHIKLMNPQRSTVWY"))
                line = line_template.format(res=res, name1=name1, ss=ss, resfile_cmd=resfile_cmd)
                lines.append(line)
            blueprint_file = os.path.join(workdir, "tmp.bp")
            with open(blueprint_file, "w") as f:
                f.write("\n".join(lines))
            if verbose:
                with open(blueprint_file, "r") as f:
                    print("Wrote Remodel blueprint file:", f.read(), sep="\n", flush=True)

            return blueprint_file

        def my_remodel_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io
            import random

            from pyrosetta.distributed.cluster.init_files import PackedPoseHasher
            from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

            verbose = kwargs["verbose"]
            tmp_path = kwargs["PyRosettaCluster_tmp_path"]
            protocol_number = kwargs["PyRosettaCluster_protocol_number"]
            seed = kwargs["PyRosettaCluster_seed"]
            random.seed(seed)
            # Maybe print
            if verbose:
                if "options" in kwargs:
                    print(f"PyRosetta protocol number {protocol_number} value of 'options' key:", kwargs["options"], flush=True)
                if "extra_options" in kwargs:
                    print(f"PyRosetta protocol number {protocol_number} value of 'extra_options' key:", kwargs["extra_options"], flush=True)
            # Setup PackedPose
            if packed_pose is None:
                # Instantiate first `PackedPose` object
                assert protocol_number == 0, f"Protocol number must be `0`. Received: {protocol_number}"
                packed_pose = io.pose_from_sequence("AA")
            elif "current_pose_pdbstring" in kwargs and "current_pose_pdbstring" in kwargs:
                # Re-instantiate current `PackedPose` object from PDB format serialization
                if verbose:
                    print("Re-instantiating current `PackedPose` object from PDB format serialization.")
                assert protocol_number > 0, f"Protocol number must be greater than `0`. Received: {protocol_number}"
                assert packed_pose.pose.sequence() == TestReproducibilityRemodelTaskUpdates._dummy_seq, (
                    f"Received `Pose` object does not have dummy sequence: {packed_pose.pose.sequence()}"
                )
                del packed_pose # Delete dummy `PackedPose` object (optional)
                packed_pose = io.pose_from_pdbstring(kwargs.pop("current_pose_pdbstring")).update_scores(kwargs.pop("current_pose_cache"))
            # Make Remodel blueprint file
            n_res_cterm_ext = random.randint(3, 5)
            n_term_threshold = random.randint(2, 3)
            blueprint_file = write_remodel_blueprint(
                packed_pose,
                tmp_path,
                n_res_cterm_ext=n_res_cterm_ext,
                n_term_threshold=n_term_threshold,
                verbose=verbose,
            )
            # Run RemodelMover
            pose = packed_pose.pose
            xml_obj = XmlObjects.create_from_string(
                f"""<MOVERS><RemodelMover name="remodel" blueprint="{blueprint_file}"/></MOVERS>"""
            ).get_mover("remodel")
            if verbose:
                print(f"Running Remodel in PyRosetta protocol number {protocol_number}.", flush=True)
            xml_obj.apply(pose)
            packed_pose = io.to_packed(pose)
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
            # Maybe print
            if verbose:
                print(f"PyRosetta protocol number {protocol_number} Pose.cache:", packed_pose.pose.cache, flush=True)
            # Maybe update task dictionary for next PyRosetta protocol, which gets validated and kept
            if protocol_number + 1 < kwargs["n_protocols"]:
                if verbose:
                    print(f"PyRosetta protocol number {protocol_number} random RNG state: {hash(str(random.getstate()))}", flush=True)
                kwargs["options"] = ""
                kwargs["extra_options"] = {
                    **TestReproducibilityRemodelTaskUpdates.get_random_options(),
                    **kwargs["constant_options"],
                }
                _reserved: set[str] = pyrosetta.Pose().cache._reserved
                _is_not_reserved = lambda k: k not in _reserved
                kwargs["current_pose_cache"] = toolz.dicttoolz.keyfilter(_is_not_reserved, packed_pose.pose.cache) # Save current `Pose.cache` dictionary
                kwargs["current_pose_pdbstring"] = io.to_pdbstring(packed_pose) # Save current `PackedPose` object
                packed_pose = io.pose_from_sequence(TestReproducibilityRemodelTaskUpdates._dummy_seq) # Dummy `PackedPose` object

            return packed_pose, kwargs

        pyrosetta.secure_unpickle.add_secure_package("pandas")
        pyrosetta.secure_unpickle.add_secure_package("pyarrow")
        scorefile_name = "test_scores.json"
        input_packed_pose = None
        compressed = True
        protocols = [my_remodel_protocol] * _n_protocols
        output_path = os.path.join(self.workdir, f"outputs_norm_task_options_{int(norm_task_options)}")
        if verbose:
            print(f"Running original simulation for `norm_task_options={norm_task_options}` and `with_init_file={with_init_file}`.", flush=True)
        PyRosettaCluster(
            tasks=create_tasks,
            input_packed_pose=input_packed_pose,
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
            input_packed_pose = input_packed_pose
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
            # Initialization from file does not run since PyRosetta is already initialized
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
        # Assert unique scorefunction results within original/reproduce
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
        if verbose:
            print(f"Successfully validated reproduction simulation for `norm_task_options={norm_task_options}` and `with_init_file={with_init_file}`.", flush=True)

    def test_reproduce_remodel_task_updates(self):
        return self.reproduce_remodel_task_updates(norm_task_options=False, with_init_file=False)

    def test_reproduce_remodel_task_updates_norm_task_options(self):
        return self.reproduce_remodel_task_updates(norm_task_options=True, with_init_file=False)

    def test_reproduce_remodel_task_updates_with_init_file(self):
        return self.reproduce_remodel_task_updates(norm_task_options=False, with_init_file=True)

    def test_reproduce_remodel_task_updates_norm_task_options_with_init_file(self):
        return self.reproduce_remodel_task_updates(norm_task_options=True, with_init_file=True)
