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
import json
import os
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import shlex
import signal
import subprocess
import sys
import tempfile
import time
import unittest

from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    export_init_file,
    requires_packed_pose,
    reserve_scores,
    reproduce,
)
from pyrosetta.tests.distributed.cluster.setup_inputs import get_test_params_file


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

            self.assertIsInstance(packed_pose.scores, dict)
            self.assertEqual(packed_pose.scores, {})
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.all_keys)
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.extra)
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.extra.real)
            self.assertEqual(packed_pose.pose.cache.extra.real["test_setPoseExtraScore"], 123.0)
            self.assertEqual(packed_pose.pose.cache.extra["test_setPoseExtraScore"], 123.0)
            self.assertEqual(packed_pose.pose.cache["test_setPoseExtraScore"], 123.0)
            packed_pose.scores.clear()
            self.assertDictEqual({}, packed_pose.scores)
            pose = io.to_pose(packed_pose)
            pose.cache.clear()
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

            self.assertIsInstance(packed_pose.scores, dict)
            self.assertEqual(packed_pose.scores, {})
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.all_keys)
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.metrics)
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.metrics.real)
            self.assertEqual(packed_pose.pose.cache.metrics.real["test_setPoseExtraScore"], 123.0)
            self.assertEqual(packed_pose.pose.cache.metrics["test_setPoseExtraScore"], 123.0)
            self.assertEqual(packed_pose.pose.cache["test_setPoseExtraScore"], 123.0)
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
        """
        Test for PyRosettaCluster decoy reproducibility with an nstruct of 2
        with multiple protocols and a fixed `decoy_ids` list.
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

        def run_subprocess(cmd, module_dir=None, timeout=900):
            print("Running command:", cmd, flush=True)

            if module_dir:
                env = os.environ.copy()
                env["PYTHONPATH"] = f"{module_dir}{os.pathsep}{os.environ.get('PYTHONPATH', '')}"
            else:
                env = None

            process = subprocess.Popen(
                shlex.split(cmd),
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                start_new_session=True,
            )

            try:
                try:
                    stdout, stderr = process.communicate(timeout=timeout)
                except subprocess.TimeoutExpired:
                    print(f"Subprocess timeout! Terminating process group: {process.pid}", flush=True)
                    os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                    stdout, stderr = process.communicate()
                    raise RuntimeError(f"Subprocess timeout while running command: {cmd}")

                print(stdout, end="", flush=True)
                if stderr:
                    print(stderr, end="", flush=True)

                if process.returncode != 0:
                    raise subprocess.CalledProcessError(
                        process.returncode, cmd, output=stdout, stderr=stderr
                    )

                return process.returncode

            finally:
                if process.poll() is None:
                    try:
                        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
                        time.sleep(1)
                        os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                    except Exception:
                        pass

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
                    print("Output from logging file:", _logging_file, flush=True)
                    with open(_logging_file, "r") as f:
                        print(f.read(), flush=True)

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
