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

import numpy
import os
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import sys
import tempfile
import unittest

from pyrosetta.utility.initialization import (
    PyRosettaInitFileReader,
    PyRosettaInitFileSerializer,
    PyRosettaInitFileWriter,
)
from pyrosetta.distributed.cluster.init_files import InitFileSigner
from pyrosetta.distributed.cluster.io import (
    METADATA_INPUT_DECOY_KEY,
    METADATA_OUTPUT_DECOY_KEY,
    get_poses_from_init_file,
    sign_init_file_metadata_and_poses,
    verify_init_file,
)


class TestInitFileSigner(unittest.TestCase):
    def tearDown(self):
        sys.stdout.flush()

    def test_init_file_signer(self):
        with tempfile.TemporaryDirectory() as workdir:
            output_init_file = os.path.join(workdir, "pyrosetta.init")
            input_packed_pose = io.pose_from_sequence("START")
            output_packed_pose = io.pose_from_sequence("FINAL")
            output_pose = output_packed_pose.pose
            pyrosetta.rosetta.core.pose.add_comment(
                output_pose,
                "REMARK PyRosettaCluster:",
                "{}",
            )
            output_packed_pose = io.to_packed(output_pose)
            poses = [output_packed_pose, input_packed_pose, input_packed_pose.clone()]
            metadata = {"comment": "foo", METADATA_INPUT_DECOY_KEY: "bar", METADATA_OUTPUT_DECOY_KEY: 10, "sha256": "test", "signature": 123}
            with self.assertRaises(ValueError):  # Fails 'sha256' and 'signature' verification
                verify_init_file(
                    output_init_file,
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=output_packed_pose,
                    metadata=metadata,
                )

            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                verbose=False,
            )
            with self.assertRaises(TypeError):  # Fails because metadata `METADATA_INPUT_DECOY_KEY` value is a `str` object
                _ = get_poses_from_init_file(output_init_file, verify=False)
            _poses = io.poses_from_init_file(output_init_file)
            self.assertEqual(len(_poses), 3)
            self.assertEqual(_poses[0].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertEqual(_poses[1].pose.sequence(), input_packed_pose.pose.sequence())

            metadata = {"comment": "bar", METADATA_INPUT_DECOY_KEY: -5, METADATA_OUTPUT_DECOY_KEY: 10, "sha256": None, "signature": None}
            verify_init_file(
                output_init_file,
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
                metadata=metadata,
            )  # Bypass verification with 'sha256' and 'signature' as `NoneType`
            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                overwrite=True,
                verbose=False,
            )
            with self.assertRaises(IndexError):  # Fails since metadata `METADATA_INPUT_DECOY_KEY` and `METADATA_OUTPUT_DECOY_KEY` values are out of range
                _ = get_poses_from_init_file(output_init_file, verify=False)

            metadata = {"comment": "baz", METADATA_INPUT_DECOY_KEY: 1, METADATA_OUTPUT_DECOY_KEY: 1, "sha256": 256}
            with self.assertRaises(ValueError):  # Fails SHA256 verification
                verify_init_file(
                    output_init_file,
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=output_packed_pose,
                    metadata=metadata,
                )
            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                overwrite=True,
                verbose=False,
            )
            _poses = get_poses_from_init_file(output_init_file, verify=False)  # Skip verification test
            self.assertEqual(len(_poses), 2)
            self.assertEqual(_poses[0].pose.sequence(), input_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[0].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertEqual(_poses[1].pose.sequence(), input_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[1].pose.sequence(), output_packed_pose.pose.sequence())

            metadata = {"comment": "baz", METADATA_INPUT_DECOY_KEY: 0, METADATA_OUTPUT_DECOY_KEY: 0}
            verify_init_file(
                output_init_file,
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
                metadata=metadata,
            )  # Bypass verification with 'sha256' and 'signature' missing
            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                overwrite=True,
                verbose=False,
            )
            _poses = get_poses_from_init_file(output_init_file, verify=False)
            self.assertEqual(len(_poses), 2)
            self.assertEqual(_poses[0].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[0].pose.sequence(), input_packed_pose.pose.sequence())
            self.assertEqual(_poses[1].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[1].pose.sequence(), input_packed_pose.pose.sequence())

            metadata = {"comment": "baz", METADATA_INPUT_DECOY_KEY: 1, METADATA_OUTPUT_DECOY_KEY: 0, "signature": "test"}
            with self.assertRaises(ValueError):  # Fails signature verification
                verify_init_file(
                    output_init_file,
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=output_packed_pose,
                    metadata=metadata,
                )
            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                overwrite=True,
                verbose=False,
            )
            with self.assertRaises(ValueError):
                _ = get_poses_from_init_file(output_init_file, verify=True)  # Fails verification
            _poses = get_poses_from_init_file(output_init_file, verify=False)  # Skip verification
            self.assertEqual(len(_poses), 2)

            self.assertEqual(_poses[0].pose.sequence(), input_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[0].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertEqual(_poses[1].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[1].pose.sequence(), input_packed_pose.pose.sequence())

            metadata = {
                "comment": "Generated by PyRosettaCluster",
                METADATA_INPUT_DECOY_KEY: 1,
                METADATA_OUTPUT_DECOY_KEY: 0,
                "version": pyrosetta.distributed.cluster.__version__,
            }
            verify_init_file(
                output_init_file,
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
                metadata=metadata,
            )  # Passes signature verification
            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                overwrite=True,
                verbose=False,
            )
            _poses = get_poses_from_init_file(output_init_file, verify=True)  # Passes verification

            metadata_signed, poses_signed = sign_init_file_metadata_and_poses(
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
            )
            init_dict = PyRosettaInitFileReader.read_json(output_init_file)
            init_dict["metadata"] = metadata_signed
            init_dict["poses"] = [io.to_base64(p) for p in poses_signed]
            init_dict["md5"] = PyRosettaInitFileSerializer.get_md5(init_dict)
            signed_init_file = os.path.join(workdir, "custom_pyrosetta.init")
            PyRosettaInitFileWriter.write_json(init_dict, signed_init_file)
            verify_init_file(
                output_init_file,
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
                metadata=dict(
                    filter(lambda kv: kv[0] not in ("sha256", "signature"), metadata_signed.items())
                ),
            )  # Passes verification

            for _metadata in [
                "string",
                b"Bytes",
                numpy.pi,
                dict(enumerate(range(20))),
                metadata,
                metadata_signed,
            ]:
                if isinstance(_metadata, bytes):
                    with self.assertRaises(TypeError):  # Object of type bytes is not JSON serializable
                        signer = InitFileSigner(
                            input_packed_pose=input_packed_pose,
                            output_packed_pose=output_packed_pose,
                            metadata=_metadata,
                        )
                else:
                    signer = InitFileSigner(
                        input_packed_pose=input_packed_pose,
                        output_packed_pose=output_packed_pose,
                        metadata=_metadata,
                    )
                signatures_dict = signer.sign()
                self.assertIsInstance(signatures_dict, dict)
                signer.verify_sha256(signatures_dict.get("sha256"))
                signer.verify_signature(signatures_dict.get("signature"))
                signer.verify(*signatures_dict.values())
