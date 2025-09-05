# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import argparse
import os
import pyrosetta.distributed.io as io
import sys

from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.rosetta.basic import was_init_called

sys.path.insert(0, os.path.dirname(__file__))
try:
    from dump_init_file import POSE_SEQUENCE, POSE_SCORE
except ImportError as ex:
    raise ImportError(ex)


def main(tmp_dir):
    os.chdir(tmp_dir)
    init_file = os.path.join(tmp_dir, "my.init")
    assert not was_init_called(), "PyRosetta is already initialized before calling `io.poses_from_init_file`."

    assert io.poses_from_init_file(None) is None, "NoneType I/O failed."
    packed_poses = io.poses_from_init_file(init_file)

    assert was_init_called(), "PyRosetta is not initialized after calling `io.poses_from_init_file`."
    assert len(packed_poses) == 1, f"Number of PackedPose objects failed: {len(packed_poses)}"
    packed_pose = packed_poses[0]
    assert isinstance(packed_pose, PackedPose), f"PackedPose object type failed: {packed_pose}"
    pose = packed_pose.pose
    sequence = pose.sequence()
    assert sequence == POSE_SEQUENCE, f"Pose sequence failed: {sequence}"
    score = pose.cache["foo"]
    assert isinstance(score, type(POSE_SCORE)), f"Pose score type failed: {type(score)}"
    assert score == POSE_SCORE, f"Pose score identity failed: {score}"

    packed_poses_copy = io.poses_from_init_file(init_file)

    pose_copy = packed_poses_copy[0].pose
    sequence_copy = pose_copy.sequence()
    assert sequence_copy == sequence, f"Pose sequence failed: {sequence_copy}"
    assert score == POSE_SCORE, f"Pose score identity failed: {score}"
    assert id(pose_copy) != id(pose), "Pose memory addresses failed."


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--tmp_dir', type=str)
    args = parser.parse_args()
    main(args.tmp_dir)
