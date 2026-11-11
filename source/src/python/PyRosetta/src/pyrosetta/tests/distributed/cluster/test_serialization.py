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
import pyrosetta.distributed.io as io
import sys
import unittest

from pyrosetta.distributed.cluster import (
    Serialization,
    update_scores,
)


class SerializationTest(unittest.TestCase):
    """Test case for the use of serialization in PyRosettaCluster."""

    def tearDown(self):
        sys.stdout.flush()

    def test_serialization(self):
        """Smoke test for PyRosettaCluster PackedPose serialization round-trip."""

        for _compression in ("xz", "zlib", "bz2", True, False, None):
            scores = {"test_str": "foo", "test_int": 123, "test_float": numpy.pi}
            for _test_case in range(3):
                input_packed_pose = io.pose_from_sequence("A" * 100)
                if _test_case == 0:
                    # Test scores not cached in Pose
                    input_packed_pose.scores = scores
                elif _test_case == 1:
                    # Test scores not cached in Pose, then update scores
                    input_packed_pose.scores = scores
                    input_packed_pose = update_scores(input_packed_pose)
                elif _test_case == 2:
                    # Test scores cached in Pose
                    input_packed_pose = input_packed_pose.update_scores(scores)

                serializer = Serialization(compression=_compression)
                compressed_packed_pose = serializer.compress_packed_pose(input_packed_pose)
                output_packed_pose = serializer.decompress_packed_pose(compressed_packed_pose)

                _error_msg = f"Failed on test case {_test_case} with compression {_compression}"

                if _compression in (False, None):
                    # Test PackedPose size
                    if _test_case == 0:
                        self.assertGreater(
                            sys.getsizeof(compressed_packed_pose.pickled_pose),
                            sys.getsizeof(input_packed_pose.pickled_pose),
                            msg=_error_msg,
                        ) # `PackedPose.scores` get cached in `Pose.cache`
                    elif _test_case in (1, 2):
                        self.assertEqual(
                            sys.getsizeof(compressed_packed_pose.pickled_pose),
                            sys.getsizeof(input_packed_pose.pickled_pose),
                            msg=_error_msg,
                        ) # `PackedPose.scores` are already cached in `Pose.cache`
                    self.assertEqual(
                        sys.getsizeof(compressed_packed_pose.pickled_pose),
                        sys.getsizeof(output_packed_pose.pickled_pose),
                        msg=_error_msg,
                    ) # `PackedPose.scores` are already cached in `Pose.cache`
                    # Test PackedPose identity
                    self.assertEqual(id(compressed_packed_pose), id(output_packed_pose), msg=_error_msg)
                else:
                    # Test PackedPose size
                    self.assertLess(
                        sys.getsizeof(compressed_packed_pose),
                        sys.getsizeof(input_packed_pose.pickled_pose),
                        msg=_error_msg,
                    )
                    self.assertLess(
                        sys.getsizeof(compressed_packed_pose),
                        sys.getsizeof(output_packed_pose.pickled_pose),
                        msg=_error_msg,
                    )
                    # Test PackedPose identity
                    self.assertNotEqual(id(compressed_packed_pose), id(output_packed_pose), msg=_error_msg)
                # Test PackedPose identity
                self.assertNotEqual(id(input_packed_pose), id(compressed_packed_pose), msg=_error_msg)
                self.assertNotEqual(id(input_packed_pose), id(output_packed_pose), msg=_error_msg)
                # Test PackedPose scores
                if _test_case == 0: # Compression has no effect when user doesn't cache scores in Pose object
                    self.assertNotEqual(input_packed_pose.pose.cache, scores, msg=_error_msg)
                    self.assertEqual(input_packed_pose.scores, scores, msg=_error_msg)
                    self.assertEqual(input_packed_pose.pose.cache, {}, msg=_error_msg)
                    self.assertEqual(output_packed_pose.scores, {}, msg=_error_msg)
                    self.assertEqual(output_packed_pose.pose.cache, scores, msg=_error_msg)
                    self.assertEqual(input_packed_pose.scores, output_packed_pose.pose.cache, msg=_error_msg)
                    self.assertEqual(output_packed_pose.pose.cache, scores, msg=_error_msg)
                    self.assertSetEqual(
                        set(input_packed_pose.scores.keys()),
                        set(output_packed_pose.pose.cache.all_keys),
                        msg=_error_msg,
                    )
                    for scoretype in input_packed_pose.scores.keys():
                        input_value = input_packed_pose.scores[scoretype]
                        output_value = output_packed_pose.pose.cache[scoretype]
                        self.assertEqual(input_value, output_value, msg=_error_msg)
                    self.assertNotEqual(
                        input_packed_pose.pickled_pose,
                        output_packed_pose.pickled_pose,
                        msg=_error_msg,
                    )
                elif _test_case in (1, 2): # Compression has no effect when user caches scores in Pose object (preferred syntax)
                    self.assertNotEqual(input_packed_pose.scores, scores, msg=_error_msg)
                    self.assertEqual(input_packed_pose.scores, {}, msg=_error_msg)
                    self.assertEqual(input_packed_pose.pose.cache, scores, msg=_error_msg)
                    self.assertEqual(output_packed_pose.scores, {}, msg=_error_msg)
                    self.assertEqual(input_packed_pose.scores, output_packed_pose.scores, msg=_error_msg)
                    self.assertEqual(output_packed_pose.pose.cache, scores, msg=_error_msg)
                    self.assertEqual(input_packed_pose.pose.cache, output_packed_pose.pose.cache, msg=_error_msg)
                    self.assertSetEqual(
                        set(input_packed_pose.pose.cache.all_keys),
                        set(output_packed_pose.pose.cache.all_keys),
                        msg=_error_msg,
                    )
                    for scoretype in input_packed_pose.pose.cache.all_keys:
                        input_value = input_packed_pose.pose.cache[scoretype]
                        output_value = output_packed_pose.pose.cache[scoretype]
                        self.assertEqual(input_value, output_value, msg=_error_msg)
                    self.assertEqual(
                        input_packed_pose.pickled_pose,
                        output_packed_pose.pickled_pose,
                        msg=_error_msg,
                    )
