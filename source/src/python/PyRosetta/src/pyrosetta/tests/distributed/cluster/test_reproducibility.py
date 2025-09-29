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


try:
    import pandas
except ImportError:
    print(
        "Importing 'pyrosetta.tests.distributed.cluster.test_reproducibility' requires the "
        + "third-party package 'pandas' as a dependency!\n"
        + "Please install this packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/pandas/\n"
    )
    raise

import functools
import json
import os
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import tempfile
import time
import unittest

from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    requires_packed_pose,
    reserve_scores,
    reproduce,
)
from pyrosetta.distributed.cluster.io import secure_read_pickle


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

    def test_reproducibility_from_reproduce(self, filter_results=False):
        """Test for PyRosettaCluster decoy reproducibility from instance kwargs."""
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

        with tempfile.TemporaryDirectory() as workdir:

            output_path = os.path.join(workdir, "outputs")
            input_pose = io.to_pose(io.pose_from_sequence("ACDEFGHIKLMNPQRSTVWY"))
            tasks = list(create_tasks())
            seeds = [-77777777, 888888888, -999999999]
            decoy_ids = [0, 1, 2]
            decoy_dir_name = "test_decoys"
            scorefile_name = "test_scores.json"
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
                ignore_errors=True,
                timeout=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                max_delay_time=0.0 if filter_results else 1.0,
                filter_results=filter_results,
            ).distribute(*protocols)

            scorefile_path = os.path.join(output_path, scorefile_name)
            with open(scorefile_path, "r") as f:
                data = [json.loads(line) for line in f]
            self.assertEqual(len(data), 1)
            original_record = data[0]

            # Reproduce decoy from .pdb.bz2 file
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
                },
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

            # Reproduce decoy from scorefile and decoy_name
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
                },
            )

            reproduce2_scorefile_path = os.path.join(
                output_path, reproduce2_scorefile_name
            )
            with open(reproduce2_scorefile_path, "r") as f:
                reproduce2_data = [json.loads(line) for line in f]
            self.assertEqual(len(reproduce2_data), 1)
            reproduce2_record = reproduce2_data[0]

            self.assertEqual(
                reproduce_record["scores"]["SEQUENCE"],
                reproduce2_record["scores"]["SEQUENCE"],
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
        cls.input_pdb_file = os.path.join(os.path.dirname(__file__), "data", "1crn.pdb")

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


# if __name__ == "__main__":
#     unittest.main(verbosity=2)
