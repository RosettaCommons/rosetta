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

import functools
import os
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import sys
import tempfile
import time
import unittest

from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    reproduce,
)
from pyrosetta.distributed.cluster.io import secure_read_pickle
from pyrosetta.tests.distributed.cluster.setup_inputs import get_test_pdb_file


class TestReproducibilityPoseDataFrame(unittest.TestCase):
    """Test decoy reproducibility from pickled `pandas.DataFrame` scorefiles."""
    @classmethod
    def setUpClass(cls):
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 300 -ignore_unrecognized_res 1 -load_PDB_components 0",
            set_logging_handler="logging",
        )
        cls.workdir = tempfile.TemporaryDirectory()
        cls.input_pdb_file = get_test_pdb_file(cls.workdir.name)

    @classmethod
    def tearDownClass(cls):
        cls.workdir.cleanup()

    def tearDown(self):
        sys.stdout.flush()

    @staticmethod
    def timeit(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            import inspect
            t0 = time.time()
            if inspect.isgeneratorfunction(func):
                results = list(func(*args, **kwargs))
            else:
                results = func(*args, **kwargs)
            dt = time.time() - t0
            print(f"{__class__.__name__} finished protocol '{func.__name__}' in {dt:0.3f} seconds.", flush=True)
            return results
        return wrapper

    @staticmethod
    def get_numpy_random_seed(**kwargs):
        protocol_number = kwargs["PyRosettaCluster_protocol_number"]
        seeds = kwargs["PyRosettaCluster_seeds"]
        seed = abs(int(seeds[protocol_number][1]))
        while seed > 2**32 - 1:
            seed //= 10
        return seed

    def assert_atom_coordinates(self, pose1, pose2):
        self.assertEqual(pose1.size(), pose2.size())
        for res in range(1, pose1.size() + 1):
            res1 = pose1.residue(res)
            res2 = pose2.residue(res)
            self.assertEqual(res1.name(), res2.name())
            self.assertEqual(res1.natoms(), res2.natoms())
            for atom in range(1, res1.natoms() + 1):
                self.assertEqual(res1.atom_name(atom), res2.atom_name(atom))
                for axis in "xyz":
                    self.assertEqual(
                        float(getattr(res1.atom(atom).xyz(), axis)),
                        float(getattr(res2.atom(atom).xyz(), axis)),
                    )

    def create_tasks(self):
        yield {
            "options": "-ex1 0 -ex1aro 0 -ex1aro_exposed 0 -ex2 0 -ex2aro 0 -ex2aro_exposed 0 -ex3 0 -ex4 0 -lazy_ig 1",
            "extra_options": "-out:level 300 -multithreading:total_threads 1 -ignore_unrecognized_res 1 -load_PDB_components 0",
            "set_logging_handler": "logging",
        }

    @staticmethod
    @timeit.__func__
    def my_first_protocol(packed_pose, **kwargs):
        """In my_first_protocol, the desired decoy_id to keep is 0."""
        import numpy
        import pyrosetta
        import pyrosetta.distributed.io as io
        from pyrosetta.rosetta.protocols.minimization_packing import (
            PackRotamersMover,
        )

        seed = TestReproducibilityPoseDataFrame.get_numpy_random_seed(**kwargs)
        numpy.random.seed(seed)

        pose = io.to_pose(packed_pose)
        pack_rotamers = PackRotamersMover(
            scorefxn=pyrosetta.create_score_function(numpy.random.choice(["ref2015.wts", "ref2015_cst.wts", "score12.wts"])),
            task=pyrosetta.standard_packer_task(pose),
            nloop=numpy.random.randint(1, 21),
        )
        pack_rotamers.apply(pose)
        dummy_pose = io.to_pose(io.pose_from_sequence("W" * 5))

        return pose.clone(), dummy_pose.clone()

    @staticmethod
    @timeit.__func__
    def my_second_protocol(packed_pose, **kwargs):
        """In my_second_protocol, the desired decoy_id to keep is 1."""
        import numpy
        import pyrosetta
        import pyrosetta.distributed.io as io
        from pyrosetta.rosetta.protocols.simple_moves import VirtualRootMover
        from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

        seed = TestReproducibilityPoseDataFrame.get_numpy_random_seed(**kwargs)
        numpy.random.seed(seed)

        pose = io.to_pose(packed_pose)
        scorefxn = pyrosetta.create_score_function(numpy.random.choice(["ref2015.wts", "ref2015_cst.wts", "score12.wts"]))
        scorefxn(pose)

        virtual_root = VirtualRootMover()
        virtual_root.set_removable(True)
        virtual_root.set_remove(False)
        virtual_root.apply(pose)

        xml = XmlObjects.create_from_string(
            f"""
            <ROSETTASCRIPTS>
            <SCOREFXNS>
                <ScoreFunction name="default_cst" weights="{numpy.random.choice(["ref2015.wts", "ref2015_cst.wts"])}">
                    <Reweight scoretype="atom_pair_constraint" weight="{numpy.random.uniform(5.0, 20.0)}" />
                </ScoreFunction>
            </SCOREFXNS>
            <RESIDUE_SELECTORS>
                <Index name="n" resnums="1"/>
                <Index name="c" resnums="{pose.size()}"/>
            </RESIDUE_SELECTORS>
            <MOVERS>
                <AddConstraints name="add_csts" >
                    <DistanceConstraintGenerator name="dist_cst"
                        residue_selector1="n"
                        residue_selector2="c"
                        function="HARMONIC {numpy.random.uniform(2.0, 3.5)} {numpy.random.uniform(0.1, 1.0)}" />
                </AddConstraints>
                <RemoveConstraints name="rm_csts" constraint_generators="dist_cst" />
                <MinMover name="min"
                    scorefxn="default_cst"
                    chi="1"
                    bb="1"
                    type="dfpmin_armijo_nonmonotone"
                    tolerance="{numpy.random.uniform(0.0001, 0.1)}"
                    max_iter="{numpy.random.randint(100, 501)}" >
                    <MoveMap name="mm"
                        bb="{numpy.random.choice([0, 1])}"
                        chi="{numpy.random.choice([0, 1])}"
                        jump="{numpy.random.choice([0, 1])}"/>
                </MinMover>
            </MOVERS>
            <PROTOCOLS>
                <Add mover="add_csts"/>
                <Add mover="min"/>
                <Add mover="rm_csts"/>
            </PROTOCOLS>
            </ROSETTASCRIPTS>
            """
        ).get_mover("ParsedProtocol")
        xml.apply(pose)

        virtual_root.set_remove(True)
        virtual_root.apply(pose)

        scorefxn(pose)
        dummy_pose = io.to_pose(io.pose_from_sequence("W" * 10))

        return dummy_pose.clone(), pose, kwargs

    @staticmethod
    @timeit.__func__
    def my_third_protocol(packed_pose, **kwargs):
        """In my_third_protocol, the desired decoy_id to keep is 2."""
        import numpy
        import pyrosetta
        import pyrosetta.distributed.io as io
        from pyrosetta.rosetta.protocols.minimization_packing import (
            PackRotamersMover,
        )

        seed = TestReproducibilityPoseDataFrame.get_numpy_random_seed(**kwargs)
        numpy.random.seed(seed)

        pose = io.to_pose(packed_pose)
        pack_rotamers = PackRotamersMover(
            scorefxn=pyrosetta.create_score_function(
                numpy.random.choice(["ref2015.wts", "ref2015_cst.wts", "score12.wts"])
            ),
            task=pyrosetta.standard_packer_task(pose),
            nloop=numpy.random.randint(1, 21),
        )
        pack_rotamers.apply(pose)

        dummy_pose = io.to_pose(io.pose_from_sequence("W" * 15))
        for p in [dummy_pose.clone(), dummy_pose.clone(), pose, dummy_pose.clone()]:
            pyrosetta.rosetta.core.pose.setPoseExtraScore(
                p, "SEQUENCE", p.sequence()
            )
            pyrosetta.rosetta.core.pose.setPoseExtraScore(
                p, "VALUE", numpy.random.randint(1e9),
            )
            yield p
        yield kwargs

    def test_reproducibility_from_reproduce(self):
        """
        Test for PyRosettaCluster decoy reproducibility from instance kwargs
        from '.pose' file and `pandas.DataFrame` scorefile.
        """
        original_output_path = os.path.join(self.workdir.name, "original_outputs")
        self.assertTrue(os.path.isfile(self.input_pdb_file))
        input_pose = io.to_pose(io.pose_from_file(self.input_pdb_file))
        tasks = list(self.create_tasks())
        decoy_dir_name = "test_decoys"
        scorefile_name = "test_scores.json"
        protocols = [
            TestReproducibilityPoseDataFrame.my_first_protocol,
            TestReproducibilityPoseDataFrame.my_second_protocol,
            TestReproducibilityPoseDataFrame.my_third_protocol,
        ]

        instance_kwargs = dict(
            tasks=tasks,
            input_packed_pose=input_pose,
            seeds=None,
            decoy_ids=None,
            client=None,
            scheduler=None,
            scratch_dir=self.workdir.name,
            cores=None,
            processes=None,
            memory=None,
            min_workers=2,
            max_workers=8,
            nstruct=1,
            dashboard_address=None,
            compressed=False,
            compression=True,
            logging_level="DEBUG",
            scorefile_name=scorefile_name,
            project_name="PyRosettaCluster_Tests",
            simulation_name=None,
            environment=None,
            output_path=original_output_path,
            simulation_records_in_scorefile=True,
            decoy_dir_name=decoy_dir_name,
            logs_dir_name="logs",
            ignore_errors=False,
            timeout=0.1,
            max_delay_time=1.0,
            sha1=None,
            dry_run=False,
            save_all=False,
            system_info=None,
            pyrosetta_build=None,
            output_decoy_types=[".pdb", ".pkl_pose", ".b64_pose"],
            output_scorefile_types=[".json", ".gz"],
        )

        if "pandas" not in pyrosetta.secure_unpickle.get_secure_packages():
            with self.assertRaises(AssertionError): # output_scorefile_types=[".gz", ...] requires 'pandas' as a secure package
                _ = PyRosettaCluster(**instance_kwargs)

        pyrosetta.secure_unpickle.add_secure_package("pandas")
        prc = PyRosettaCluster(**instance_kwargs)
        prc.distribute(
            protocols=protocols,
            clients_indices=None,
            resources=None,
        )

        scorefile_path = os.path.join(original_output_path, os.path.splitext(scorefile_name)[0] + ".gz")
        df = secure_read_pickle(scorefile_path, compression="infer")
        selected_decoy_ids = [0, 1, 2]
        self.assertEqual(len(selected_decoy_ids), len(protocols))
        original_record = None
        for index in df.index:
            decoy_ids = df.at[index, "instance"]["decoy_ids"]
            self.assertEqual(len(decoy_ids), len(protocols))
            if all(decoy_ids[i] == selected_decoy_ids[i] for i in range(len(protocols))):
                original_record = df.loc[index]
                self.assertEqual(len(original_record["instance"]["seeds"]), len(protocols))
                break
        self.assertIsNotNone(original_record, msg=f"Could not locate decoy with selected 'decoy_ids'.")

        # Reproduce decoy from .b64_pose or .pkl_pose file
        reproduce_output_path = os.path.join(self.workdir.name, "reproduce_outputs")
        reproduce_scorefile_name = "reproduce_test_scores.json"
        reproduce_input_file_extension = ".b64_pose" # ".pkl_pose"
        self.assertIn(reproduce_input_file_extension, (".b64_pose", ".pkl_pose"))
        reproduce_input_file = os.path.splitext(original_record["metadata"]["output_file"])[0] + reproduce_input_file_extension

        reproduce(
            input_file=reproduce_input_file,
            scorefile=None,
            decoy_name=None,
            protocols=protocols,
            input_packed_pose=input_pose,
            client=None,
            instance_kwargs={
                "output_path": reproduce_output_path,
                "sha1": None,
                "scorefile_name": reproduce_scorefile_name,
                "output_decoy_types": [".pdb", ".b64_pose"],
                "output_scorefile_types": [".json", ".xz"],
            },
        )

        reproduce_scorefile_path = os.path.join(
            reproduce_output_path, os.path.splitext(reproduce_scorefile_name)[0] + ".xz"
        )
        df_reproduce = secure_read_pickle(reproduce_scorefile_path, compression="infer")
        self.assertEqual(df_reproduce.index.size, 1)
        reproduce_record = df_reproduce.iloc[0]
        self.assertListEqual(original_record["instance"]["seeds"], reproduce_record["instance"]["seeds"])
        reproduce_output_file = os.path.splitext(reproduce_record["metadata"]["output_file"])[0] + ".b64_pose"

        self.assertEqual(
            original_record["scores"]["SEQUENCE"],
            reproduce_record["scores"]["SEQUENCE"],
        )
        self.assertEqual(
            original_record["scores"]["VALUE"],
            reproduce_record["scores"]["VALUE"],
        )
        self.assertEqual(
            original_record["scores"]["total_score"],
            reproduce_record["scores"]["total_score"],
        )
        self.assertListEqual(
            original_record["instance"]["seeds"],
            reproduce_record["instance"]["seeds"],
        )
        self.assertListEqual(
            original_record["instance"]["decoy_ids"],
            reproduce_record["instance"]["decoy_ids"],
        )

        if reproduce_input_file_extension == ".b64_pose":
            original_pose = io.pose_from_base64(reproduce_input_file).pose
        elif reproduce_input_file_extension == ".pkl_pose":
            original_pose = io.pose_from_pickle(reproduce_input_file).pose
        reproduce_pose = io.pose_from_base64(reproduce_output_file).pose
        self.assert_atom_coordinates(original_pose, reproduce_pose)

        # Reproduce decoy from `pandas.DataFrame` scorefile and decoy_name
        reproduce2_output_path = os.path.join(self.workdir.name, "reproduce2_outputs")
        reproduce2_scorefile_name = "reproduce2_test_scores.json"
        reproduce(
            input_file=None,
            scorefile=os.path.join(
                reproduce_record["instance"]["output_path"],
                os.path.splitext(reproduce_record["instance"]["scorefile_name"])[0] + ".xz",
            ),
            decoy_name=reproduce_record["metadata"]["decoy_name"],
            input_packed_pose=input_pose,
            protocols=protocols,
            client=None,
            instance_kwargs={
                "output_path": reproduce2_output_path,
                "sha1": None,
                "scorefile_name": reproduce2_scorefile_name,
                "output_decoy_types": [".b64_pose"],
                "output_scorefile_types": [".tar"],
            },
        )

        reproduce2_scorefile_path = os.path.join(
            reproduce2_output_path, os.path.splitext(reproduce2_scorefile_name)[0] + ".tar"
        )
        df_reproduce2 = secure_read_pickle(reproduce2_scorefile_path, compression="infer")
        self.assertEqual(df_reproduce2.index.size, 1)
        reproduce2_record = df_reproduce2.iloc[0]
        self.assertListEqual(reproduce_record["instance"]["decoy_ids"], reproduce2_record["instance"]["decoy_ids"])
        self.assertListEqual(reproduce_record["instance"]["seeds"], reproduce2_record["instance"]["seeds"])
        reproduce2_output_file = os.path.splitext(reproduce2_record["metadata"]["output_file"])[0] + ".b64_pose"

        self.assertEqual(
            reproduce_record["scores"]["SEQUENCE"],
            reproduce2_record["scores"]["SEQUENCE"],
        )
        self.assertEqual(
            original_record["scores"]["VALUE"],
            reproduce_record["scores"]["VALUE"],
        )
        self.assertEqual(
            reproduce_record["scores"]["total_score"],
            reproduce2_record["scores"]["total_score"],
        )
        self.assertListEqual(
            reproduce_record["instance"]["seeds"],
            reproduce2_record["instance"]["seeds"],
        )
        self.assertListEqual(
            reproduce_record["instance"]["decoy_ids"],
            reproduce2_record["instance"]["decoy_ids"],
        )

        reproduce2_pose = io.pose_from_base64(reproduce2_output_file).pose
        self.assert_atom_coordinates(reproduce_pose, reproduce2_pose)
        self.assert_atom_coordinates(original_pose, reproduce2_pose)
