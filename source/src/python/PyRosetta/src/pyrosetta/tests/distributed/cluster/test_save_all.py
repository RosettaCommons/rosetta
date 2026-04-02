"""
PyRosettaCluster smoke tests using the `unittest` framework.
"""
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

import bz2
import json
import os
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import sys
import tempfile
import unittest
import uuid

from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.utility.initialization import PyRosettaInitFileReader
from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    export_init_file,
)
from pyrosetta.distributed.cluster.init_files import InitFileSigner
from pyrosetta.distributed.cluster.io import (
    METADATA_INPUT_DECOY_KEY,
    METADATA_OUTPUT_DECOY_KEY,
)


class SaveAllTest(unittest.TestCase):
    def tearDown(self):
        sys.stdout.flush()

    def test_save_all(self):
        """Smoke test for PyRosettaCluster usage with the save_all attribute enabled."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )

        _total_tasks = 4
        _total_protocols = 4

        def create_tasks():
            for i in range(_total_tasks):
                yield {
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                    "seq": "LYELL" * (i + 1),
                }

        def my_pyrosetta_protocol(input_packed_pose, **kwargs):
            return input_packed_pose

        with tempfile.TemporaryDirectory() as workdir:
            output_path = os.path.join(workdir, "outputs")
            scorefile_name = "test_save_all.json"
            nstruct = 1
            author = "Username"
            email = "test@example"
            license = "LICENSE.PyRosetta.md"
            project_name = None
            simulation_name = "SaveAllTest"
            input_packed_pose = io.pose_from_sequence("TESTING")
            output_decoy_types = [".pdb", ".init"]
            output_scorefile_types = [".json"]
            compressed = True
            cluster = PyRosettaCluster(
                tasks=create_tasks,
                input_packed_pose=input_packed_pose,
                seeds=None,
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
                compressed=compressed,
                logging_level="CRITICAL",
                scorefile_name=scorefile_name,
                project_name=project_name,
                simulation_name=simulation_name,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=True,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=True,
                timeout=1.0,
                sha1=None,
                dry_run=False,
                save_all=True,
                system_info=None,
                pyrosetta_build=None,
                author=author,
                email=email,
                license=license,
                filter_results=True,
                output_decoy_types=output_decoy_types,
                output_scorefile_types=output_scorefile_types,
                norm_task_options=None,
            ) # Test 'output_init_file' default
            protocol_args = [my_pyrosetta_protocol] * _total_protocols
            cluster.distribute(*protocol_args)

            with open(os.path.join(output_path, scorefile_name), "r") as f:
                data = [json.loads(line) for line in f]
            self.assertEqual(len(data), _total_tasks * _total_protocols)

            _project_name = "PyRosettaCluster" if project_name is None else project_name
            _simulation_name = "PyRosettaCluster" if simulation_name is None else simulation_name
            init_file = os.path.join(output_path, f"{_project_name}_{_simulation_name}_pyrosetta.init")
            if compressed:
                init_file += ".bz2"
            self.assertTrue(os.path.isfile(init_file))
            for record in data:
                self.assertTrue(os.path.samefile(init_file, record["metadata"]["init_file"]))
            if compressed:
                init_data = io.read_init_file(init_file)
            else:
                init_data = PyRosettaInitFileReader.read_json(init_file)

            # Verify output PyRosetta initialization file
            self.assertEqual(init_data["author"], author)
            self.assertEqual(init_data["email"], email)
            self.assertEqual(init_data["license"], license)
            self.assertEqual(init_data["pyrosetta_build"], pyrosetta._build_signature())
            self.assertListEqual(init_data["options"]["run:constant_seed"], ["true"])

            # Verify metadata in output PyRosetta initialization file
            metadata = init_data["metadata"]
            self.assertIsInstance(metadata, dict)
            self.assertEqual(metadata["comment"], "Generated by PyRosettaCluster")
            self.assertEqual(metadata["version"], pyrosetta.distributed.cluster.__version__)
            if isinstance(input_packed_pose, (Pose, PackedPose)):
                self.assertIn(METADATA_INPUT_DECOY_KEY, metadata)
            else:
                self.assertNotIn(METADATA_INPUT_DECOY_KEY, metadata)
            self.assertNotIn(METADATA_OUTPUT_DECOY_KEY, metadata)

            self.assertIn("sha256", metadata)
            sha256 = metadata.pop("sha256", None)
            self.assertIn("signature", metadata)
            signature = metadata.pop("signature", None)
            signer = InitFileSigner(
                input_packed_pose=input_packed_pose,
                output_packed_pose=None,
                metadata=metadata,
            )
            self.assertTrue(signer.verify_sha256(sha256))
            self.assertTrue(signer.verify_signature(signature))

            # Export PyRosetta initialization file with output decoy
            record = data[0] # Use first output decoy
            output_file = record["metadata"]["output_file"]

            if ".init" in output_decoy_types:
                if compressed:
                    decoy_init_file = os.path.splitext(os.path.splitext(output_file)[0])[0] + ".init.bz2"
                    self.assertTrue(os.path.isfile(decoy_init_file))
                    with open(decoy_init_file, "rb") as fbz2:
                        decoy_init_data = PyRosettaInitFileReader.from_json(bz2.decompress(fbz2.read()).decode())
                    _decoy_init_data = io.read_init_file(decoy_init_file)
                    self.assertDictEqual(decoy_init_data, _decoy_init_data)
                else:
                    decoy_init_file = os.path.splitext(output_file)[0] + ".init"
                    self.assertTrue(os.path.isfile(decoy_init_file))
                    decoy_init_data = PyRosettaInitFileReader.read_json(decoy_init_file)
                    _decoy_init_data = io.read_init_file(decoy_init_file)
                    self.assertDictEqual(decoy_init_data, _decoy_init_data)

                # Verify decoy output PyRosetta initialization file
                self.assertEqual(decoy_init_data["author"], author)
                self.assertEqual(decoy_init_data["email"], email)
                self.assertEqual(decoy_init_data["license"], license)
                self.assertEqual(decoy_init_data["pyrosetta_build"], pyrosetta._build_signature())
                self.assertListEqual(decoy_init_data["options"]["run:constant_seed"], ["true"])

                # Verify metadata in decoy output PyRosetta initialization file
                decoy_metadata = decoy_init_data["metadata"]
                self.assertIsInstance(decoy_metadata, dict)
                self.assertEqual(decoy_metadata["comment"], "Generated by PyRosettaCluster")
                self.assertEqual(decoy_metadata["version"], pyrosetta.distributed.cluster.__version__)
                if isinstance(input_packed_pose, (Pose, PackedPose)):
                    self.assertIn(METADATA_INPUT_DECOY_KEY, decoy_metadata)
                else:
                    self.assertNotIn(METADATA_INPUT_DECOY_KEY, decoy_metadata)
                self.assertIn(METADATA_OUTPUT_DECOY_KEY, decoy_metadata)
                poses_output_idx = decoy_metadata[METADATA_OUTPUT_DECOY_KEY]
                self.assertEqual(poses_output_idx, 0)

                self.assertIn("sha256", decoy_metadata)
                decoy_sha256 = decoy_metadata.pop("sha256")
                self.assertIn("signature", decoy_metadata)
                decoy_signature = decoy_metadata.pop("signature")
                decoy_signer_1 = InitFileSigner(
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=io.to_packed(io.to_pose(decoy_init_data["poses"][poses_output_idx])),
                    metadata=decoy_metadata,
                )
                self.assertTrue(decoy_signer_1.verify_sha256(decoy_sha256))
                self.assertTrue(decoy_signer_1.verify_signature(decoy_signature))
                decoy_signer_2 = InitFileSigner(
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=io.pose_from_init_file(decoy_init_file),
                    metadata=decoy_metadata,
                )
                self.assertTrue(decoy_signer_2.verify_sha256(decoy_sha256))
                self.assertTrue(decoy_signer_2.verify_signature(decoy_signature))
                decoy_signer_3 = InitFileSigner(
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=io.poses_from_init_file(decoy_init_file)[poses_output_idx],
                    metadata=decoy_metadata,
                )
                self.assertTrue(decoy_signer_3.verify_sha256(decoy_sha256))
                self.assertTrue(decoy_signer_3.verify_signature(decoy_signature))

            custom_output_init_file = os.path.join(workdir, "custom_output.init")
            export_init_file(
                output_file,
                output_init_file=custom_output_init_file,  # Test custom output file path
                compressed=False,  # Test not compressed
            )
            self.assertTrue(os.path.isfile(custom_output_init_file))
            export_init_file(
                output_file,
                output_init_file=custom_output_init_file,  # Test custom output file path
                compressed=True,  # Test compressed
            )
            custom_output_init_bz2_file = custom_output_init_file + ".bz2"
            self.assertTrue(os.path.isfile(custom_output_init_bz2_file))
            if ".init" not in output_decoy_types:
                export_init_file(
                    output_file,
                    output_init_file=None,  # Test default output file path
                    compressed=compressed,
                )
                exported_init_file = (
                    output_file[: -len(f"{output_decoy_types[0]}.bz2")]
                    if compressed
                    else output_file[: -len(output_decoy_types[0])]
                ) + (".init.bz2" if compressed else ".init")
            else:
                exported_init_file = (
                    output_file[: -len(f"{output_decoy_types[0]}.bz2")]
                    if compressed
                    else output_file[: -len(output_decoy_types[0])]
                ) + "_my_custom.init"
                export_init_file(
                    output_file,
                    output_init_file=exported_init_file,
                    compressed=compressed,
                )
                if compressed:
                    exported_init_file += ".bz2"
            self.assertTrue(os.path.isfile(exported_init_file))
            exported_init_data = io.read_init_file(exported_init_file)

            # Verify exported PyRosetta initialization file
            self.assertEqual(exported_init_data["author"], author)
            self.assertEqual(exported_init_data["email"], email)
            self.assertEqual(exported_init_data["license"], license)
            self.assertEqual(exported_init_data["pyrosetta_build"], pyrosetta._build_signature())
            self.assertListEqual(exported_init_data["options"]["run:constant_seed"], ["true"])

            # Verify metadata in exported PyRosetta initialization file
            exported_metadata = exported_init_data["metadata"]
            self.assertIsInstance(exported_metadata, dict)
            self.assertEqual(exported_metadata["comment"], "Generated by PyRosettaCluster")
            self.assertEqual(exported_metadata["version"], pyrosetta.distributed.cluster.__version__)
            if isinstance(input_packed_pose, (Pose, PackedPose)):
                self.assertIn(METADATA_INPUT_DECOY_KEY, exported_metadata)
            else:
                self.assertNotIn(METADATA_INPUT_DECOY_KEY, exported_metadata)
            self.assertIn(METADATA_OUTPUT_DECOY_KEY, exported_metadata)
            self.assertIn("sha256", exported_metadata)
            exported_sha256 = exported_metadata.pop("sha256", None)
            self.assertIn("signature", exported_metadata)
            exported_signature = exported_metadata.pop("signature", None)
            exported_signer = InitFileSigner(
                input_packed_pose=input_packed_pose,
                output_packed_pose=io.pose_from_init_file(exported_init_file), # Output decoy is first object in 'poses' key value
                metadata=exported_metadata,
            )
            self.assertTrue(exported_signer.verify_sha256(exported_sha256))
            self.assertTrue(exported_signer.verify_signature(exported_signature))

            _decoy_names = set()
            for record in data:
                self.assertEqual(init_data["pyrosetta_build"], record["instance"]["pyrosetta_build"])
                self.assertEqual(exported_init_data["pyrosetta_build"], record["instance"]["pyrosetta_build"])
                for key in ("author", "email", "license"):
                    self.assertIn(key, record["instance"])
                    self.assertNotIn(key, record["metadata"])
                    self.assertNotIn(key, record["scores"])
                self.assertEqual(record["instance"]["author"], author)
                self.assertEqual(record["instance"]["email"], email)
                self.assertEqual(record["instance"]["license"], license)
                self.assertNotIn("init_file", record["instance"])
                self.assertIn("init_file", record["metadata"])
                self.assertTrue(os.path.samefile(record["metadata"]["init_file"], init_file))
                self.assertDictEqual(record["scores"], {})
                _decoy_name = record["metadata"]["decoy_name"]
                self.assertNotIn(_decoy_name, _decoy_names)
                _decoy_names.add(_decoy_name)
                self.assertEqual(record["instance"]["nstruct"], nstruct)

    def test_save_all_dry_run(self):
        """
        Smoke test for PyRosettaCluster usage with the save_all attribute and
        dry_run attribute enabled.
        """
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        _total_tasks = 2
        _total_protocols = 5
        _nstruct = 2

        def create_tasks():
            for i in range(_total_tasks):
                yield {
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                    "seq": "LYELL" * (i + 1),
                }

        def my_pyrosetta_protocol(input_packed_pose, **kwargs):
            return input_packed_pose

        with tempfile.TemporaryDirectory() as workdir:
            output_path = os.path.join(workdir, "outputs")
            scorefile_name = "test_save_all.json"
            decoy_dir_name = "test_decoys"
            logs_dir_name = "test_logs"
            init_file = os.path.join(output_path, "test.init")
            PyRosettaCluster(
                tasks=create_tasks,
                input_packed_pose=io.pose_from_sequence("TESTING"),
                seeds=None,
                decoy_ids=None,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=_nstruct,
                dashboard_address=None,
                compressed=True,
                logging_level="CRITICAL",
                scorefile_name=scorefile_name,
                project_name="PyRosettaCluster_SaveAllTest",
                simulation_name=uuid.uuid4().hex,
                environment=None,
                output_path=output_path,
                decoy_dir_name=decoy_dir_name,
                logs_dir_name=logs_dir_name,
                simulation_records_in_scorefile=True,
                ignore_errors=True,
                timeout=1.0,
                max_delay_time=3.0,
                sha1=None,
                dry_run=True,
                save_all=True,
                system_info=None,
                pyrosetta_build=None,
                filter_results=True,
                norm_task_options=None,
                output_init_file=init_file,
            ).distribute(protocols=[my_pyrosetta_protocol] * _total_protocols)

            self.assertFalse(os.path.exists(os.path.join(output_path, scorefile_name)))
            self.assertFalse(os.path.exists(init_file))
            self.assertTrue(os.path.exists(os.path.join(output_path, decoy_dir_name)))
            self.assertTrue(os.path.exists(os.path.join(output_path, logs_dir_name)))
            self.assertEqual(
                os.listdir(os.path.join(output_path, decoy_dir_name)),
                [],
            )
            self.assertNotEqual(
                os.listdir(os.path.join(output_path, logs_dir_name)),
                [],
            )
