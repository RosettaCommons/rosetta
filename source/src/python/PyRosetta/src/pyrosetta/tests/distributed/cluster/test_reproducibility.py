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

import json
import os
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import sys
import tempfile
import unittest

from pyrosetta.distributed.cluster import PyRosettaCluster


class TestReproducibility(unittest.TestCase):
    def tearDown(self):
        sys.stdout.flush()

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
                ignore_errors=False,
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
                ignore_errors=False,
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
                ignore_errors=False,
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
                ignore_errors=False,
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
