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
import glob
import json
import os
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest
import uuid

from contextlib import contextmanager
from pathlib import Path

from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    export_init_file,
    recreate_environment,
    requires_packed_pose,
    reserve_scores,
    reproduce,
)
from pyrosetta.distributed.cluster.config import get_environment_var
from pyrosetta.distributed.cluster.io import secure_read_pickle
from pyrosetta.tests.distributed.cluster.setup_inputs import get_test_params_file, get_test_pdb_file


class TestReproducibility(unittest.TestCase):
    def test_reproducibility_packer_nstruct(self, filter_results=False):
        """Test for PyRosettaCluster decoy reproducibility with an nstruct of 2."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 300",
            set_logging_handler="logging",
        )

        def my_pyrosetta_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io
            from pyrosetta.rosetta.protocols.minimization_packing import (
                PackRotamersMover,
            )

            pose = io.to_pose(packed_pose)
            pack_rotamers = PackRotamersMover(
                scorefxn=pyrosetta.create_score_function("ref2015.wts"),
                task=pyrosetta.standard_packer_task(pose),
                nloop=25,
            )
            pack_rotamers.apply(pose)
            pyrosetta.rosetta.core.pose.setPoseExtraScore(
                pose, "SEQUENCE", pose.sequence()
            )
            self.assertEqual(type(pose), type(pyrosetta.Pose()))
            return io.to_packed(pose), None, pyrosetta.Pose(), io.to_packed(pyrosetta.Pose())

        with tempfile.TemporaryDirectory() as workdir:
            nstruct = 2
            output_path = os.path.join(workdir, "outputs")
            tasks = [
                {
                    "options": "-ex1",
                    "extra_options": "-out:level 300 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                },
            ]

            cluster = PyRosettaCluster(
                tasks=tasks,
                input_packed_pose=io.pose_from_sequence("TESTING"),
                seeds=1234567,
                decoy_ids=None,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=nstruct,
                dashboard_address=None,
                compressed=True,
                logging_level="INFO",
                scorefile_name="scores.json",
                project_name="PyRosettaCluster_Tests",
                simulation_name=None,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=True,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=True,
                timeout=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                max_delay_time=0.0 if filter_results else 1.0,
                filter_results=filter_results,
                norm_task_options=None,
                output_init_file=None,
            )

            cluster.distribute(protocols=[my_pyrosetta_protocol])

            with open(os.path.join(output_path, "scores.json"), "r") as f:
                data = [json.loads(line) for line in f]
            self.assertEqual(len(data), nstruct)
            self.assertEqual(
                data[0]["scores"]["SEQUENCE"], data[1]["scores"]["SEQUENCE"]
            )
            self.assertEqual(
                data[0]["scores"]["total_score"], data[1]["scores"]["total_score"]
            )
            self.assertListEqual(
                data[0]["instance"]["seeds"], data[1]["instance"]["seeds"]
            )

    def test_reproducibility_minimizer_nstruct(self, filter_results=False):
        """Test for PyRosettaCluster decoy reproducibility with an nstruct of 2."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 300",
            set_logging_handler="logging",
        )

        with tempfile.TemporaryDirectory() as workdir:

            nstruct = 2
            output_path = os.path.abspath(os.path.join(workdir, "outputs"))
            tasks = [
                {
                    "options": "-ex1",
                    "extra_options": "-out:level 300 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                },
            ]

            cluster = PyRosettaCluster(
                tasks=tasks,
                input_packed_pose=io.pose_from_sequence("ACDEFGHIKLMNPQRSTVWY"),
                seeds=-897897897,
                decoy_ids=None,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=nstruct,
                dashboard_address=None,
                compressed=True,
                logging_level="INFO",
                scorefile_name="scores.json",
                project_name="PyRosettaCluster_Tests",
                simulation_name=None,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=True,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=True,
                timeout=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                max_delay_time=0.0 if filter_results else 1.0,
                filter_results=filter_results,
                norm_task_options=None,
                output_init_file=None,
            )

            def my_pyrosetta_protocol(packed_pose, **kwargs):
                import pyrosetta
                import pyrosetta.distributed.io as io
                from pyrosetta.rosetta.core.scoring import ScoreType
                from pyrosetta.rosetta.protocols.minimization_packing import MinMover
                from pyrosetta.rosetta.protocols.simple_moves import VirtualRootMover
                from pyrosetta.rosetta.protocols.constraint_generator import (
                    AddConstraints,
                    RemoveConstraints,
                    TerminiConstraintGenerator,
                )

                pose = io.to_pose(packed_pose)
                scorefxn = pyrosetta.create_score_function("ref2015_cst.wts")
                scorefxn.set_weight(ScoreType.atom_pair_constraint, 20.0)

                scorefxn(pose)
                virtual_root = VirtualRootMover()
                virtual_root.set_removable(True)
                virtual_root.set_remove(False)
                virtual_root.apply(pose)

                cst_gen = TerminiConstraintGenerator()
                cst_gen.set_max_distance(3.0)
                cst_gen.set_min_distance(2.0)
                cst_gen.set_sd(0.1)
                cst_gen.set_weight(20.0)

                add_constraints = AddConstraints()
                add_constraints.add_generator(cst_gen)
                add_constraints.apply(pose)

                mm = pyrosetta.MoveMap()
                mm.set_bb(True)
                mm.set_chi(True)

                min_mover = MinMover()
                min_mover.movemap(mm)
                min_mover.score_function(scorefxn)
                min_mover.min_type("dfpmin_armijo_nonmonotone")
                min_mover.tolerance(0.0001)
                min_mover.max_iter(10000)
                min_mover.apply(pose)

                remove_constraints = RemoveConstraints()
                remove_constraints.add_generator(cst_gen)
                remove_constraints.apply(pose)

                virtual_root.set_remove(True)
                virtual_root.apply(pose)

                scorefxn(pose)

                return pose, None, pyrosetta.Pose(), io.to_packed(pyrosetta.Pose())

            cluster.distribute(my_pyrosetta_protocol,)

            with open(os.path.join(output_path, "scores.json"), "r") as f:
                data = [json.loads(line) for line in f]
            self.assertEqual(len(data), nstruct)
            self.assertAlmostEqual(
                data[0]["scores"]["total_score"],
                data[1]["scores"]["total_score"],
                places=14,
                msg="Minimizer did not converge using identical seeds!",
            )
            self.assertListEqual(
                data[0]["instance"]["seeds"], data[1]["instance"]["seeds"]
            )

    def test_reproducibility_packer_separate(self, filter_results=False):
        """
        Test for PyRosettaCluster decoy reproducibility from two
        separate PyRosettaCluster instantiations.
        """
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 300",
            set_logging_handler="logging",
        )

        def create_tasks():
            yield {
                "options": "-ex1",
                "extra_options": "-out:level 300 -multithreading:total_threads 1",
                "set_logging_handler": "logging",
            }

        def my_custom_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io
            from pyrosetta.rosetta.protocols.minimization_packing import (
                PackRotamersMover,
            )

            pose = io.to_pose(packed_pose)
            pack_rotamers = PackRotamersMover(
                scorefxn=pyrosetta.create_score_function("ref2015.wts"),
                task=pyrosetta.standard_packer_task(pose),
                nloop=25,
            )
            pack_rotamers.apply(pose)
            pyrosetta.rosetta.core.pose.setPoseExtraScore(
                pose, "SEQUENCE", pose.sequence()
            )

            return pose, None, pyrosetta.Pose(), io.to_packed(pyrosetta.Pose())

        with tempfile.TemporaryDirectory() as workdir:

            input_pose = io.to_pose(io.pose_from_sequence("TESTING"))
            output_path = os.path.abspath(os.path.join(workdir, "outputs"))
            seed = 987654321
            nstruct = 1

            cluster = PyRosettaCluster(
                tasks=create_tasks,
                input_packed_pose=input_pose.clone(),
                seeds=seed,
                decoy_ids=None,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=nstruct,
                dashboard_address=None,
                compressed=True,
                logging_level="INFO",
                scorefile_name="scores.json",
                project_name="PyRosettaCluster_Tests",
                simulation_name=None,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=True,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=True,
                timeout=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                max_delay_time=0.0 if filter_results else 1.0,
                filter_results=filter_results,
                norm_task_options=None,
                output_init_file=None,
            )

            def sample_func(packed_pose, **kwargs):
                return my_custom_protocol(packed_pose, **kwargs)

            cluster.distribute(sample_func,)

            independent_cluster_instance = PyRosettaCluster(
                tasks=create_tasks,
                input_packed_pose=input_pose.clone(),
                seeds=seed,
                decoy_ids=None,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=nstruct,
                dashboard_address=None,
                compressed=True,
                logging_level="INFO",
                scorefile_name="scores.json",
                project_name="PyRosettaCluster_Tests",
                simulation_name=None,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=True,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=True,
                timeout=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                max_delay_time=0.0 if filter_results else 1.0,
                filter_results=filter_results,
                norm_task_options=None,
                output_init_file=None,
            )

            independent_cluster_instance.distribute(sample_func)

            with open(os.path.join(output_path, "scores.json"), "r") as f:
                data = [json.loads(line) for line in f]
            self.assertEqual(len(data), nstruct * 2)
            self.assertEqual(
                data[0]["scores"]["SEQUENCE"], data[1]["scores"]["SEQUENCE"]
            )
            self.assertEqual(
                data[0]["scores"]["total_score"], data[1]["scores"]["total_score"]
            )
            self.assertListEqual(
                data[0]["instance"]["seeds"], data[1]["instance"]["seeds"]
            )

    def test_reproducibility_packer_nstruct_filter_results(self, filter_results=True):
        return self.test_reproducibility_packer_nstruct(filter_results=filter_results)

    def test_reproducibility_minimizer_nstruct_filter_results(self, filter_results=True):
        return self.test_reproducibility_minimizer_nstruct(filter_results=filter_results)

    def test_reproducibility_packer_separate_filter_results(self, filter_results=True):
        return self.test_reproducibility_packer_separate(filter_results=filter_results)


class TestReproducibilityMulti(unittest.TestCase):
    def test_reproducibility_packer_nstruct_multi(self, filter_results=False):
        """
        Test for PyRosettaCluster decoy reproducibility with an nstruct of 2
         with multiple protocols.
        """
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 300",
            set_logging_handler="logging",
        )

        def create_tasks():
            yield {
                "options": "-ex1",
                "extra_options": "-out:level 300 -multithreading:total_threads 1",
                "set_logging_handler": "logging",
            }

        def my_first_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io
            from pyrosetta.rosetta.protocols.minimization_packing import (
                PackRotamersMover,
            )

            pose = io.to_pose(packed_pose)
            pyrosetta.rosetta.core.pose.setPoseExtraScore(
                pose, "test_setPoseExtraScore", 123
            )
            pack_rotamers = PackRotamersMover(
                scorefxn=pyrosetta.create_score_function("ref2015.wts"),
                task=pyrosetta.standard_packer_task(pose),
                nloop=10,
            )
            pack_rotamers.apply(pose)

            return pose, None, pyrosetta.Pose(), io.to_packed(pyrosetta.Pose())

        @reserve_scores
        @requires_packed_pose
        def my_second_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io
            from pyrosetta.rosetta.protocols.minimization_packing import (
                PackRotamersMover,
            )

            if packed_pose.pose.empty() or packed_pose is None:
                raise ValueError(
                    "The user-provided PyRosetta protocol is decorated with `@requires_packed_pose`."
                )
            assert packed_pose.pose.size() >= 1

            self.assertEqual(
                dict(packed_pose.scores),
                {**dict(packed_pose.scores), **{"test_setPoseExtraScore": 123}},
            )
            packed_pose.scores.clear()
            self.assertDictEqual({}, packed_pose.scores)

            pose = io.to_pose(packed_pose)
            pack_rotamers = PackRotamersMover(
                scorefxn=pyrosetta.create_score_function("ref2015.wts"),
                task=pyrosetta.standard_packer_task(pose),
                nloop=10,
            )
            pack_rotamers.apply(pose)

            return pose, None, pyrosetta.Pose(), io.to_packed(pyrosetta.Pose())

        def my_third_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io
            from pyrosetta.rosetta.protocols.minimization_packing import (
                PackRotamersMover,
            )

            if packed_pose.pose.empty():
                assert filter_results == False
                return None

            self.assertEqual(
                dict(packed_pose.scores),
                {**dict(packed_pose.scores), **{"test_setPoseExtraScore": 123}},
            )
            pose = io.to_pose(packed_pose)
            pack_rotamers = PackRotamersMover(
                scorefxn=pyrosetta.create_score_function("ref2015.wts"),
                task=pyrosetta.standard_packer_task(pose),
                nloop=10,
            )
            pack_rotamers.apply(pose)
            pyrosetta.rosetta.core.pose.setPoseExtraScore(
                pose, "SEQUENCE", pose.sequence()
            )
            return pose, None, pyrosetta.Pose(), io.to_packed(pyrosetta.Pose())

        with tempfile.TemporaryDirectory() as workdir:

            output_path = os.path.abspath(os.path.join(workdir, "outputs"))
            nstruct = 2
            input_pose = io.to_pose(io.pose_from_sequence("ACDEFGHIKLMNPQRSTVWY"))
            tasks = list(create_tasks())
            seeds = [987654, -1234567, 87654321]
            protocols = [my_first_protocol, my_second_protocol, my_third_protocol]
            self.assertEqual(len(seeds), len(protocols))

            _ = list(
                PyRosettaCluster(
                    tasks=tasks,
                    input_packed_pose=input_pose,
                    seeds=seeds,
                    decoy_ids=None,
                    client=None,
                    scheduler=None,
                    scratch_dir=workdir,
                    cores=None,
                    processes=None,
                    memory=None,
                    min_workers=1,
                    max_workers=1,
                    nstruct=nstruct,
                    dashboard_address=None,
                    compressed=True,
                    logging_level="INFO",
                    scorefile_name=None,
                    project_name="PyRosettaCluster_Tests",
                    simulation_name=None,
                    environment=None,
                    output_path=output_path,
                    simulation_records_in_scorefile=True,
                    decoy_dir_name="decoys",
                    logs_dir_name="logs",
                    ignore_errors=True,
                    timeout=1.0,
                    sha1=None,
                    dry_run=False,
                    save_all=False,
                    system_info=None,
                    pyrosetta_build=None,
                    max_delay_time=0.0 if filter_results else 1.0,
                    filter_results=filter_results,
                    norm_task_options=None,
                    output_init_file=None,
                ).generate(*protocols)
            )

            with open(os.path.join(output_path, "scores.json"), "r") as f:
                data = [json.loads(line) for line in f]
            self.assertEqual(len(data), nstruct)
            self.assertEqual(
                data[0]["scores"]["SEQUENCE"], data[1]["scores"]["SEQUENCE"]
            )
            self.assertEqual(
                data[0]["scores"]["total_score"], data[1]["scores"]["total_score"]
            )
            self.assertListEqual(
                data[0]["instance"]["seeds"], data[1]["instance"]["seeds"]
            )

    def test_reproducibility_packer_nstruct_multi_decoy_ids(self, filter_results=False):
        """Test for PyRosettaCluster decoy reproducibility with an nstruct of 2
        with multiple protocols and a fixed decoy_ids list."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 300",
            set_logging_handler="logging",
        )

        def create_tasks():
            yield {
                "options": "-ex1",
                "extra_options": "-out:level 300 -multithreading:total_threads 1",
                "set_logging_handler": "logging",
            }

        def my_first_protocol(packed_pose, **kwargs):
            """In my_first_protocol, the desired decoy_id to keep is 0."""
            import pyrosetta
            import pyrosetta.distributed.io as io
            from pyrosetta.rosetta.protocols.minimization_packing import (
                PackRotamersMover,
            )

            pose = io.to_pose(packed_pose)
            pack_rotamers = PackRotamersMover(
                scorefxn=pyrosetta.create_score_function("ref2015.wts"),
                task=pyrosetta.standard_packer_task(pose),
                nloop=10,
            )
            pack_rotamers.apply(pose)
            dummy_pose = io.to_pose(io.pose_from_sequence("W" * 6))

            return (
                pose.clone(),
                dummy_pose.clone(),
                dummy_pose.clone(),
                None,
                pyrosetta.Pose(),
                io.to_packed(pyrosetta.Pose()),
            )

        @requires_packed_pose
        def my_second_protocol(packed_pose, **kwargs):
            """In my_second_protocol, the desired decoy_id to keep is 1."""
            import pyrosetta
            import pyrosetta.distributed.io as io
            from pyrosetta.rosetta.protocols.simple_moves import VirtualRootMover
            from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

            if packed_pose.empty() or packed_pose is None:
                raise ValueError(
                    "The user-provided PyRosetta protocol is decorated with `@requires_packed_pose`."
                )
            assert packed_pose.pose.size() >= 1

            pose = io.to_pose(packed_pose)
            scorefxn = pyrosetta.create_score_function("ref2015_cst.wts")
            scorefxn(pose)

            virtual_root = VirtualRootMover()
            virtual_root.set_removable(True)
            virtual_root.set_remove(False)
            virtual_root.apply(pose)

            xml = XmlObjects.create_from_string(
                """
                <ROSETTASCRIPTS>
                <SCOREFXNS>
                    <ScoreFunction name="default_cst" weights="ref2015.wts">
                        <Reweight scoretype="atom_pair_constraint" weight="20.0" />
                    </ScoreFunction>
                </SCOREFXNS>
                <RESIDUE_SELECTORS>
                    <Index name="n" resnums="1"/>
                    <Index name="c" resnums="20"/>
                </RESIDUE_SELECTORS>
                <MOVERS>
                    <AddConstraints name="add_csts" >
                        <DistanceConstraintGenerator name="dist_cst"
                            residue_selector1="n"
                            residue_selector2="c"
                            function="HARMONIC 2.5 0.1" />
                    </AddConstraints>
                    <RemoveConstraints name="rm_csts" constraint_generators="dist_cst" />
                    <MinMover name="min"
                        scorefxn="default_cst"
                        chi="1"
                        bb="1"
                        type="dfpmin_armijo_nonmonotone"
                        tolerance="0.0001"
                        max_iter="500" >
                        <MoveMap name="mm" bb="1" chi="1" jump="1" />
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
            dummy_pose = io.to_pose(io.pose_from_sequence("W" * 6))

            return (
                dummy_pose.clone(),
                pose,
                dummy_pose.clone(),
                None,
                pyrosetta.Pose(),
                io.to_packed(pyrosetta.Pose()),
            )

        def my_third_protocol(packed_pose, **kwargs):
            """In my_third_protocol, the desired decoy_id to keep is 2."""
            import pyrosetta
            import pyrosetta.distributed.io as io
            from pyrosetta.rosetta.protocols.minimization_packing import (
                PackRotamersMover,
            )

            if packed_pose.empty():
                assert packed_pose.pose.empty()
                assert filter_results == False
                return None

            pose = io.to_pose(packed_pose)
            pack_rotamers = PackRotamersMover(
                scorefxn=pyrosetta.create_score_function("ref2015.wts"),
                task=pyrosetta.standard_packer_task(pose),
                nloop=10,
            )
            pack_rotamers.apply(pose)

            dummy_pose = io.to_pose(io.pose_from_sequence("W" * 6))
            for p in [dummy_pose.clone(), dummy_pose.clone(), pose]:
                pyrosetta.rosetta.core.pose.setPoseExtraScore(
                    p, "SEQUENCE", p.sequence()
                )
                yield p

        with tempfile.TemporaryDirectory() as workdir:

            output_path = os.path.abspath(os.path.join(workdir, "outputs"))
            nstruct = 2
            input_pose = io.to_pose(io.pose_from_sequence("ACDEFGHIKLMNPQRSTVWY"))
            tasks = list(create_tasks())
            seeds = [-111111111, -222222222, -333333333]
            decoy_ids = [0, 1, 2]
            protocols = [my_first_protocol, my_second_protocol, my_third_protocol]
            self.assertEqual(len(seeds), len(protocols))
            self.assertEqual(len(decoy_ids), len(protocols))

            PyRosettaCluster(
                tasks=tasks,
                input_packed_pose=input_pose,
                seeds=seeds,
                decoy_ids=decoy_ids,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=nstruct,
                dashboard_address=None,
                compressed=True,
                logging_level="INFO",
                scorefile_name=None,
                project_name="PyRosettaCluster_Tests",
                simulation_name=None,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=True,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=True,
                timeout=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                max_delay_time=0.0 if filter_results else 1.0,
                filter_results=filter_results,
                norm_task_options=None,
                output_init_file=None,
            ).distribute(protocols=protocols)

            with open(os.path.join(output_path, "scores.json"), "r") as f:
                data = [json.loads(line) for line in f]
            self.assertEqual(len(data), nstruct)
            self.assertEqual(
                data[0]["scores"]["SEQUENCE"], data[1]["scores"]["SEQUENCE"]
            )
            self.assertEqual(
                data[0]["scores"]["total_score"], data[1]["scores"]["total_score"]
            )
            self.assertListEqual(
                data[0]["instance"]["seeds"], data[1]["instance"]["seeds"]
            )
            self.assertListEqual(
                data[0]["instance"]["decoy_ids"], data[1]["instance"]["decoy_ids"]
            )

    def test_reproducibility_from_reproduce(self, filter_results=False, verbose=False):
        """Test for PyRosettaCluster decoy reproducibility from instance kwargs."""
        params_dir = tempfile.TemporaryDirectory(prefix="tmp_params_")
        params_file = get_test_params_file(params_dir.name)
        pyrosetta.distributed.init(
            options=f"-run:constant_seed 1 -multithreading:total_threads 1 -extra_res_fa {params_file}",
            extra_options="-out:level 300",
            set_logging_handler="logging",
        )

        def create_tasks():
            yield {
                "options": f"-ex1 -extra_res_fa {params_file}",
                "extra_options": "-out:level 300 -multithreading:total_threads 1",
                "set_logging_handler": "logging",
            }

        def my_first_protocol(packed_pose, **kwargs):
            """In my_first_protocol, the desired decoy_id to keep is 0."""
            import pyrosetta
            import pyrosetta.distributed.io as io
            from pyrosetta.rosetta.protocols.minimization_packing import (
                PackRotamersMover,
            )

            pose = io.to_pose(packed_pose)
            pack_rotamers = PackRotamersMover(
                scorefxn=pyrosetta.create_score_function("ref2015.wts"),
                task=pyrosetta.standard_packer_task(pose),
                nloop=10,
            )
            pack_rotamers.apply(pose)
            dummy_pose = io.to_pose(io.pose_from_sequence("W" * 6))

            return (
                pose.clone(),
                dummy_pose.clone(),
                None,
                pyrosetta.Pose(),
                io.to_packed(pyrosetta.Pose()),
                dummy_pose.clone(),
            )

        def my_second_protocol(packed_pose, **kwargs):
            """In my_second_protocol, the desired decoy_id to keep is 1."""
            import pyrosetta
            import pyrosetta.distributed.io as io
            from pyrosetta.rosetta.protocols.simple_moves import VirtualRootMover
            from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

            if packed_pose.pose.empty():
                assert filter_results == False
                return None

            pose = io.to_pose(packed_pose)
            scorefxn = pyrosetta.create_score_function("ref2015_cst.wts")
            scorefxn(pose)

            virtual_root = VirtualRootMover()
            virtual_root.set_removable(True)
            virtual_root.set_remove(False)
            virtual_root.apply(pose)

            xml = XmlObjects.create_from_string(
                """
                <ROSETTASCRIPTS>
                <SCOREFXNS>
                    <ScoreFunction name="default_cst" weights="ref2015.wts">
                        <Reweight scoretype="atom_pair_constraint" weight="20.0" />
                    </ScoreFunction>
                </SCOREFXNS>
                <RESIDUE_SELECTORS>
                    <Index name="n" resnums="1"/>
                    <Index name="c" resnums="20"/>
                </RESIDUE_SELECTORS>
                <MOVERS>
                    <AddConstraints name="add_csts" >
                        <DistanceConstraintGenerator name="dist_cst"
                            residue_selector1="n"
                            residue_selector2="c"
                            function="HARMONIC 2.5 0.1" />
                    </AddConstraints>
                    <RemoveConstraints name="rm_csts" constraint_generators="dist_cst" />
                    <MinMover name="min"
                        scorefxn="default_cst"
                        chi="1"
                        bb="1"
                        type="dfpmin_armijo_nonmonotone"
                        tolerance="0.0001"
                        max_iter="500" >
                        <MoveMap name="mm" bb="1" chi="1" jump="1" />
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
            dummy_pose = io.to_pose(io.pose_from_sequence("W" * 6))

            return (
                dummy_pose.clone(),
                pose,
                dummy_pose.clone(),
                None,
                pyrosetta.Pose(),
                io.to_packed(pyrosetta.Pose()),
            )

        @requires_packed_pose
        def my_third_protocol(packed_pose, **kwargs):
            """In my_third_protocol, the desired decoy_id to keep is 2."""
            import pyrosetta
            import pyrosetta.distributed.io as io
            from pyrosetta.rosetta.protocols.minimization_packing import (
                PackRotamersMover,
            )

            if packed_pose.pose.empty() or packed_pose is None:
                raise ValueError(
                    "The user-provided PyRosetta protocol is decorated with `@requires_packed_pose`."
                )
            assert packed_pose.pose.size() >= 1

            pose = io.to_pose(packed_pose)
            pack_rotamers = PackRotamersMover(
                scorefxn=pyrosetta.create_score_function("ref2015.wts"),
                task=pyrosetta.standard_packer_task(pose),
                nloop=10,
            )
            pack_rotamers.apply(pose)

            dummy_pose = io.to_pose(io.pose_from_sequence("W" * 6))
            for p in [dummy_pose.clone(), dummy_pose.clone(), pose]:
                pyrosetta.rosetta.core.pose.setPoseExtraScore(
                    p, "SEQUENCE", p.sequence()
                )
                yield p

        def run_subprocess(cmd, module_dir=None):
            print("Running command:", cmd)
            if module_dir:
                env = os.environ.copy()
                env["PYTHONPATH"] = f"{module_dir}{os.pathsep}{os.environ.get('PYTHONPATH', '')}"
            else:
                env = None
            p = subprocess.run(shlex.split(cmd), env=env, check=True, shell=False, stderr=subprocess.PIPE)
            print("Return code: {0}".format(p.returncode))

            return p.returncode

        with tempfile.TemporaryDirectory() as workdir:
            sequence = "ACDEFGHIKLMNPQRSTVWY"
            output_path = os.path.join(workdir, "outputs")
            output_init_file = os.path.join(output_path, "pyrosetta.init")
            input_pose = io.to_pose(io.pose_from_sequence(sequence))
            tasks = list(create_tasks())
            seeds = [-77777777, 888888888, -999999999]
            decoy_ids = [0, 1, 2]
            decoy_dir_name = "test_decoys"
            scorefile_name = "test_scores.json"
            protocols = [my_first_protocol, my_second_protocol, my_third_protocol]
            author = "Username"
            email = "test@example"
            license = "LICENSE.PyRosetta.md"
            self.assertEqual(len(seeds), len(protocols))
            self.assertEqual(len(decoy_ids), len(protocols))

            PyRosettaCluster(
                tasks=tasks,
                input_packed_pose=input_pose,
                seeds=seeds,
                decoy_ids=decoy_ids,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=1,
                dashboard_address=None,
                compressed=True,
                logging_level="INFO",
                scorefile_name=scorefile_name,
                project_name="PyRosettaCluster_Tests",
                simulation_name=None,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=True,
                decoy_dir_name=decoy_dir_name,
                logs_dir_name="logs",
                ignore_errors=False,
                timeout=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                author=author,
                email=email,
                license=license,
                max_delay_time=0.0 if filter_results else 1.0,
                filter_results=filter_results,
                norm_task_options=None,
                output_init_file=output_init_file,
                output_decoy_types=[".pdb", ".b64_pose"],
                output_scorefile_types=[".json"],
            ).distribute(*protocols)

            scorefile_path = os.path.join(output_path, scorefile_name)
            with open(scorefile_path, "r") as f:
                data = [json.loads(line) for line in f]
            self.assertEqual(len(data), 1)
            original_record = data[0]

            if verbose:
                logging_file = original_record["metadata"]["logging_file"]
                _logging_files = glob.glob(os.path.join(os.path.dirname(logging_file), "*"))
                for _logging_file in _logging_files:
                    print("Output from logging file:", _logging_file)
                    with open(_logging_file, "r") as f:
                        print(f.read())

            # Reproduce decoy from .pdb.bz2 file
            run_reproduce_from_pdb_bz2_file = True
            if run_reproduce_from_pdb_bz2_file:
                reproduce_scorefile_name = "reproduce_test_scores.json"
                reproduce(
                    input_file=original_record["metadata"]["output_file"],
                    scorefile=None,
                    decoy_name=None,
                    protocols=protocols,
                    input_packed_pose=input_pose,
                    client=None,
                    instance_kwargs={
                        "sha1": None,
                        "scorefile_name": reproduce_scorefile_name,
                        "output_init_file": None, # Skip `dump_init_file`
                    },
                    skip_corrections=True, # Skip corrections to reuse task result for reproduction
                    init_from_file_kwargs=None,
                )

                reproduce_scorefile_path = os.path.join(
                    output_path, reproduce_scorefile_name
                )
                with open(reproduce_scorefile_path, "r") as f:
                    reproduce_data = [json.loads(line) for line in f]
                self.assertEqual(len(reproduce_data), 1)
                reproduce_record = reproduce_data[0]

                self.assertEqual(
                    original_record["scores"]["SEQUENCE"],
                    reproduce_record["scores"]["SEQUENCE"],
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
                for key in ("author", "email", "license"):
                    self.assertEqual(
                        original_record["instance"][key],
                        reproduce_record["instance"][key],
                    )

            # Reproduce decoy from scorefile and decoy_name
            run_reproduce_from_scorefile = True
            if run_reproduce_from_scorefile:
                reproduce2_scorefile_name = "reproduce2_test_scores.json"
                reproduce(
                    input_file=None,
                    scorefile=os.path.join(
                        reproduce_record["instance"]["output_path"],
                        reproduce_record["instance"]["scorefile_name"],
                    ),
                    decoy_name=reproduce_record["metadata"]["decoy_name"],
                    input_packed_pose=input_pose,
                    protocols=None, # Test detecting protocols in current scope
                    client=None,
                    instance_kwargs={
                        "sha1": None,
                        "scorefile_name": reproduce2_scorefile_name,
                        "output_init_file": "", # Skip `dump_init_file`
                    },
                    skip_corrections=True, # Skip corrections to reuse task result for reproduction
                    init_from_file_kwargs=None,
                )

                reproduce2_scorefile_path = os.path.join(
                    output_path, reproduce2_scorefile_name
                )
                with open(reproduce2_scorefile_path, "r") as f:
                    reproduce2_data = [json.loads(line) for line in f]
                self.assertEqual(len(reproduce2_data), 1)
                reproduce2_record = reproduce2_data[0]

                _records = [original_record]
                if run_reproduce_from_pdb_bz2_file:
                    _records.append(reproduce_record)
                for _record in _records:
                    self.assertEqual(
                        _record["scores"]["SEQUENCE"],
                        reproduce2_record["scores"]["SEQUENCE"],
                    )
                    self.assertEqual(
                        _record["scores"]["total_score"],
                        reproduce2_record["scores"]["total_score"],
                    )
                    self.assertListEqual(
                        _record["instance"]["seeds"],
                        reproduce2_record["instance"]["seeds"],
                    )
                    self.assertListEqual(
                        _record["instance"]["decoy_ids"],
                        reproduce2_record["instance"]["decoy_ids"],
                    )
                    for key in ("author", "email", "license"):
                        self.assertEqual(
                            _record["instance"][key],
                            reproduce2_record["instance"][key],
                        )

            # Test raised exceptions:
            with self.assertRaises(TypeError):
                # Test invalid 'skip_corrections' parameter
                reproduce(skip_corrections=None)

            self.assertIn("init_file", original_record["metadata"])
            with self.assertRaises(ValueError):
                # Test PyRosetta intialization file without an output decoy
                reproduce(input_file=original_record["metadata"]["init_file"])

            # Reproduce decoy with a .init file
            original_init_file = original_record["metadata"]["init_file"]  # Without output decoy
            original_compressed = original_record["instance"]["compressed"]
            self.assertTrue(os.path.isfile(original_init_file))
            input_file = original_record["metadata"]["output_file"]
            self.assertTrue(os.path.isfile(input_file))
            # Export
            _init_dict_before_export = io.read_init_file(original_init_file)
            with self.assertRaises(FileExistsError):
                # Cannot overwrite original .init or .init.bz2 file
                if original_compressed:
                    output_init_file = os.path.splitext(original_init_file)[0]
                else:
                    output_init_file = original_init_file
                export_init_file(
                    input_file,
                    output_init_file=output_init_file,
                    compressed=original_compressed,
                )
            export_init_file(
                input_file,
                output_init_file=None,  # Default
                compressed=original_compressed,
            )
            default_output_init_file = os.path.splitext(input_file)[0]
            if original_compressed:
                default_output_init_file = os.path.splitext(default_output_init_file)[0]
            default_output_init_file += ".init"
            if original_compressed:
                default_output_init_file += ".bz2"
            self.assertTrue(os.path.isfile(default_output_init_file), msg=f"File not found: {default_output_init_file}")
            _init_dict_after_export = io.read_init_file(default_output_init_file)
            self.assertNotEqual(_init_dict_before_export, _init_dict_after_export)
            self.assertNotEqual(
                _init_dict_before_export["metadata"]["sha256"],
                _init_dict_after_export["metadata"]["sha256"],
            )
            self.assertNotEqual(
                _init_dict_before_export["metadata"]["signature"],
                _init_dict_after_export["metadata"]["signature"],
            )

            with self.assertRaises(TypeError):
                # Test PyRosetta initialization file with the output decoy, but wrong input decoy
                reproduce(
                    input_file=default_output_init_file,
                    input_packed_pose=io.pose_from_sequence("INVALID")
                )

            test_script = os.path.join(os.path.dirname(__file__), "reproduce_from_init_file.py")
            module = os.path.splitext(os.path.basename(test_script))[0]
            for test_case in range(5):
                reproduce3_scorefile_name = f"reproduce_test_scores_init_file_{test_case}.json"
                cmd = "{0} -m {1} --input_file '{2}' --scorefile_name '{3}' --input_init_file '{4}' --sequence '{5}' --test_case {6}".format(
                    sys.executable,
                    module,
                    input_file,
                    reproduce3_scorefile_name,
                    default_output_init_file,
                    sequence,
                    test_case,
                )
                returncode = run_subprocess(cmd, module_dir=os.path.dirname(test_script))
                self.assertEqual(returncode, 0, msg=f"Test script failed: {test_script}")

                reproduce3_scorefile_path = os.path.join(
                    output_path, reproduce3_scorefile_name
                )
                with open(reproduce3_scorefile_path, "r") as f:
                    reproduce3_data = [json.loads(line) for line in f]
                self.assertEqual(len(reproduce3_data), 1)
                reproduce3_record = reproduce3_data[0]

                _records = [original_record]
                if run_reproduce_from_pdb_bz2_file:
                    _records.append(reproduce_record)
                if run_reproduce_from_scorefile:
                    _records.append(reproduce2_record)
                for _record in _records:
                    self.assertEqual(
                        _record["scores"]["SEQUENCE"],
                        reproduce3_record["scores"]["SEQUENCE"],
                    )
                    self.assertEqual(
                        _record["scores"]["total_score"],
                        reproduce3_record["scores"]["total_score"],
                    )
                    self.assertListEqual(
                        _record["instance"]["seeds"],
                        reproduce3_record["instance"]["seeds"],
                    )
                    self.assertListEqual(
                        _record["instance"]["decoy_ids"],
                        reproduce3_record["instance"]["decoy_ids"],
                    )
                    for key in ("author", "email", "license"):
                        self.assertEqual(
                            _record["instance"][key],
                            reproduce3_record["instance"][key],
                        )
        params_dir.cleanup()

    def test_reproducibility_packer_nstruct_multi_filter_results(self):
        return self.test_reproducibility_packer_nstruct_multi(filter_results=True)

    def test_reproducibility_packer_nstruct_multi_decoy_ids_filter_results(self):
        return self.test_reproducibility_packer_nstruct_multi_decoy_ids(filter_results=True)

    def test_reproducibility_from_reproduce_filter_results(self):
        return self.test_reproducibility_from_reproduce(filter_results=True)


class TestReproducibilityPoseDataFrame(unittest.TestCase):
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
            print(f"{__class__.__name__} finished protocol '{func.__name__}' in {dt:0.3f} seconds.")
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


class TestEnvironmentReproducibility(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 300 -ignore_unrecognized_res 1 -load_PDB_components 0",
            set_logging_handler="logging",
        )
        cls.workdir = tempfile.TemporaryDirectory()
        cls.run_tag = uuid.uuid4().hex[:12]

    @classmethod
    def tearDownClass(cls):
        try:
            cls.workdir.cleanup()
        except Exception as ex:
            print(f"Warning: failed to cleanup temporary directory: {ex}... Continuing.")
        os.environ.pop(get_environment_var(), None)

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

    @staticmethod
    def run_subprocess(cmd, module_dir=None, cwd=None):
        print("Running command:", cmd)
        if module_dir:
            env = os.environ.copy()
            pythonpath = os.environ.get("PYTHONPATH")
            env["PYTHONPATH"] = f"{module_dir}{os.pathsep + pythonpath if pythonpath else ''}"
        else:
            env = None
        try:
            # Use live output streaming for GitHub Actions visibility
            process = subprocess.Popen(
                shlex.split(cmd),
                cwd=cwd,
                env=env,
                shell=False,
                stdout=sys.stdout,
                stderr=sys.stderr,
                text=True,
            )
            returncode = process.wait()
            if returncode != 0:
                raise subprocess.CalledProcessError(returncode, cmd)
        except subprocess.CalledProcessError as ex:
            print(f"Subprocess command failed (return code: {ex.returncode}): {cmd}", flush=True)
            raise
        except Exception as ex:
            print(f"Unexpected error in subprocess: {ex}", flush=True)
            raise
        else:
            print(f"Return code: {returncode}", flush=True)

        return returncode

    def recreate_environment_test(self, environment_manager="conda"):
        """Test for PyRosettaCluster decoy reproducibility in a recreated virtual environment."""
        self.assertIn(environment_manager, ("conda", "mamba", "uv", "pixi"))

        test_script = os.path.join(os.path.dirname(__file__), "recreate_environment_test_runs.py")

        # Create new environment
        original_env_name = f"{environment_manager}_env_{self.run_tag}"
        original_env_dir = os.path.join(self.workdir.name, original_env_name)
        setup_env_script = os.path.join(os.path.dirname(__file__), "setup_envs.py")
        module = os.path.splitext(os.path.basename(setup_env_script))[0]
        cmd = "{0} -m {1} --env_manager '{2}' --env_dir '{3}'".format(
            sys.executable,
            module,
            environment_manager,
            original_env_dir,
        )
        returncode = TestEnvironmentReproducibility.run_subprocess(
            cmd,
            module_dir=os.path.dirname(setup_env_script),
            cwd=None,
        )
        self.assertEqual(returncode, 0, msg=f"Subprocess command failed: {cmd}")

        # Run original simulation inside new environment
        original_output_path = os.path.join(original_env_dir, f"{environment_manager}_original_outputs")
        original_scorefile_name = "test_scores.json"
        if environment_manager == "pixi":
            cmd = "pixi run python {0} --env_manager '{1}' --output_path '{2}' --scorefile_name '{3}'".format(
                test_script,
                environment_manager,
                original_output_path,
                original_scorefile_name,
            )
        elif environment_manager == "uv":
            cmd = "uv run -p {0} python {1} --env_manager '{2}' --output_path '{3}' --scorefile_name '{4}'".format(
                original_env_dir,
                test_script,
                environment_manager,
                original_output_path,
                original_scorefile_name,
            )
        elif environment_manager in ("conda", "mamba"):
            cmd = "conda run -p {0} python {1} --env_manager '{2}' --output_path '{3}' --scorefile_name '{4}'".format(
                original_env_dir,
                test_script,
                environment_manager,
                original_output_path,
                original_scorefile_name,
            )
        returncode = TestEnvironmentReproducibility.run_subprocess(
            cmd,
            module_dir=None,
            # For pixi, activate the original pixi environment context
            # For conda/mamba/uv, run from environment directory for consistency with pixi workflow
            cwd=original_env_dir,
        )
        self.assertEqual(returncode, 0, msg=f"Subprocess command failed: {cmd}")

        # Recreate new environment from output scorefile
        original_scorefile_path = os.path.join(original_output_path, original_scorefile_name)
        self.assertTrue(os.path.isfile(original_scorefile_path), msg=f"Missing original output scorefile: {original_scorefile_path}")
        with open(original_scorefile_path, "r") as f:
            original_data = [json.loads(line) for line in f]
        self.assertEqual(len(original_data), 1)
        original_record = original_data[0]
        self.assertIn("environment_manager", original_record["metadata"])
        self.assertEqual(original_record["metadata"]["environment_manager"], environment_manager)
        self.assertIn("decoy_name", original_record["metadata"])
        original_decoy_name = original_record["metadata"]["decoy_name"]
        # Set environment manager
        os.environ[get_environment_var()] = environment_manager
        # Recreate environment
        reproduce_env_name = f"{original_env_name}_reproduce"
        recreate_environment(
            environment_name=reproduce_env_name,
            input_file=None,
            scorefile=original_scorefile_path,
            decoy_name=original_decoy_name,
            timeout=999,
            base_dir=self.workdir.name,
        )
        reproduce_env_dir = os.path.join(self.workdir.name, reproduce_env_name)
        self.assertTrue(
            os.path.isdir(reproduce_env_dir),
            f"Reproduced '{environment_manager}' environment directory was not created: '{reproduce_env_dir}'",
        )
        if environment_manager == "uv":
            # The recreated uv environment uses the PyPI 'pyrosetta-installer' package, which does not allow specifying PyRosetta version.
            # Therefore, installing the correct PyRosetta version in the recreated uv environment depends fortuitously on a prompt
            # uv environment recreation after the original uv environment creation.
            print("Running PyRosetta installer in recreated uv environment...")
            # Run PyRosetta installer with mirror fallback
            install_script = textwrap.dedent("""
                import pyrosetta_installer
                try:
                    pyrosetta_installer.install_pyrosetta(
                        distributed=False,
                        serialization=True,
                        skip_if_installed=True,
                        mirror=0
                    )
                except Exception as e:
                    print(f"Recreated PyRosetta installation with 'mirror=0' failed: {e}. Retrying with 'mirror=1'.")
                    pyrosetta_installer.install_pyrosetta(
                        distributed=False,
                        serialization=True,
                        skip_if_installed=True,
                        mirror=1
                    )
            """)
            subprocess.run(
                ["uv", "run", "-p", str(reproduce_env_dir), "python", "-c", install_script],
                check=True,
            )

        # Run reproduction simulation inside recreated environment
        reproduce_output_path = os.path.join(reproduce_env_dir, f"{environment_manager}_reproduce_outputs")
        reproduce_scorefile_name = "test_scores.json"
        if environment_manager == "pixi":
            cmd = (
                f"pixi run python {test_script} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{reproduce_output_path}' "
                f"--scorefile_name '{reproduce_scorefile_name}' "
                f"--original_scorefile '{original_scorefile_path}' "
                f"--original_decoy_name '{original_decoy_name}' "
                "--reproduce"
            )
        elif environment_manager == "uv":
            cmd = (
                f"uv run -p {reproduce_env_dir} python {test_script} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{reproduce_output_path}' "
                f"--scorefile_name '{reproduce_scorefile_name}' "
                f"--original_scorefile '{original_scorefile_path}' "
                f"--original_decoy_name '{original_decoy_name}' "
                "--reproduce"
            )
        elif environment_manager in ("conda", "mamba"):
            cmd = (
                f"conda run -p {reproduce_env_dir} python {test_script} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{reproduce_output_path}' "
                f"--scorefile_name '{reproduce_scorefile_name}' "
                f"--original_scorefile '{original_scorefile_path}' "
                f"--original_decoy_name '{original_decoy_name}' "
                "--reproduce"
            )
        returncode = TestEnvironmentReproducibility.run_subprocess(
            cmd,
            module_dir=None,
            # For pixi, activate the recreated pixi environment context
            # For conda/mamba/uv, run from recreated environment directory for consistency with pixi workflow
            cwd=reproduce_env_dir,
        )
        self.assertEqual(returncode, 0, msg=f"Subprocess command failed: {cmd}")

        # Validate reproduced decoy is identical to original decoy
        reproduce_scorefile_path = os.path.join(reproduce_output_path, reproduce_scorefile_name)
        self.assertTrue(os.path.isfile(reproduce_scorefile_path), msg=f"Missing reproduced output scorefile: {reproduce_scorefile_path}")
        with open(reproduce_scorefile_path, "r") as f:
            reproduce_data = [json.loads(line) for line in f]
        self.assertEqual(len(reproduce_data), 1)
        reproduce_record = reproduce_data[0]
        self.assertIn("environment_manager", reproduce_record["metadata"])
        self.assertEqual(reproduce_record["metadata"]["environment_manager"], environment_manager)

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
        self.assertNotEqual(
            original_record["metadata"]["author"],
            reproduce_record["metadata"]["author"],
        )
        self.assertNotEqual(
            original_record["metadata"]["decoy_name"],
            reproduce_record["metadata"]["decoy_name"],
        )
        original_pose = io.pose_from_file(original_record["metadata"]["output_file"])
        reproduce_pose = io.pose_from_file(reproduce_record["metadata"]["output_file"])
        self.assert_atom_coordinates(original_pose, reproduce_pose)

    @unittest.skipIf(shutil.which("conda") is None, "The executable 'conda' is not available.")
    def test_recreate_environment_conda(self):
        return self.recreate_environment_test(environment_manager="conda")

    @unittest.skipIf(shutil.which("mamba") is None, "The executable 'mamba' is not available.")
    def test_recreate_environment_mamba(self):
        return self.recreate_environment_test(environment_manager="mamba")

    @unittest.skipIf(shutil.which("uv") is None, "The executable 'uv' is not available.")
    def test_recreate_environment_uv(self):
        return self.recreate_environment_test(environment_manager="uv")

    @unittest.skipIf(shutil.which("pixi") is None, "The executable 'pixi' is not available.")
    def test_recreate_environment_pixi(self):
        return self.recreate_environment_test(environment_manager="pixi")


# if __name__ == "__main__":
#     unittest.main(verbosity=2)
