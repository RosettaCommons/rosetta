"""IO routines operating on PackedPose representations."""
import functools

import pyrosetta.io
import pyrosetta.rosetta.utility as utility
import pyrosetta.rosetta.core.pose as pose
import pyrosetta.rosetta.core.import_pose as import_pose

from pyrosetta.distributed import requires_init
from pyrosetta.distributed.packed_pose import (
    pack_result, to_pose, to_packed, to_dict, register_container_traversal
)

__all__ = [
    "pose_from_file",
    "pose_from_sequence",
    "poses_from_silent",
    "pose_from_pdbstring",
    "to_silent",
    "to_pdbstring",
    "to_output_record"
]

# Input functions

pose_from_file = requires_init(pack_result(
    pyrosetta.io.pose_from_file))
pose_from_sequence = requires_init(pack_result(
    pyrosetta.io.pose_from_sequence))
poses_from_silent = requires_init(pack_result(
    pyrosetta.io.poses_from_silent))


@functools.wraps(import_pose.pose_from_pdbstring)
@requires_init
@pack_result
def pose_from_pdbstring(*args, **kwargs):
    result = pose.Pose()
    import_pose.pose_from_pdbstring(result, *args, **kwargs)
    return result


@functools.singledispatch
@requires_init
def to_silent(inp, output_filename):
    """Takes a Pose, PackedPose, or a list and outputs them as a silent file.
    This currently only outputs to a binary silent file.

    Inputs:
    poses: Pose object, PackedPose object or list of either.
    output_filename: The desired name of the output silent file.

    Example:
    pyrosetta.distributed.io.to_silent(poses, "mydesigns.silent")

    The decoy name in your silent file is take from pose.pdb_info().name()
    To set a different decoy name, change it in your pose before calling this function.
    To change the name, you must have a Pose object, not a PackedPose
    Example:
    pose = pyrosetta.distributed.packed_pose.to_pose(packed_pose)
    pose.pdb_info().name("my_tag")

    @srgerb
    """
    if isinstance(inp, (list, tuple, set)):
        pose = [to_pose(p) for p in inp]
    else:
        pose = to_pose(inp)
    pyrosetta.io.poses_to_silent(pose, output_filename)


@functools.singledispatch
@requires_init
def to_pdbstring(inp):
    """Convert to pdb-formatted string with score and energy data.
    """
    from pyrosetta.rosetta.core.io.pdb import dump_pdb
    from pyrosetta.rosetta.core.io import StructFileRepOptions
    from pyrosetta.rosetta.std import ostringstream

    sfro = StructFileRepOptions()
    sfro.set_output_pose_cache_data(True)
    sfro.set_output_pose_energies_table(True)

    oss = ostringstream()
    dump_pdb(to_pose(inp), oss, sfro)

    return oss.bytes().decode()


register_container_traversal(to_pdbstring, lambda d: to_pdbstring(to_packed(d)))


@functools.singledispatch
def to_output_record(inp):
    """Convert to an "archive" output record with scores, pdb string and packed string.
    """

    output_record = dict(to_dict(inp))
    output_record["pdb_pose"] = to_pdbstring(inp)
    output_record["pickled_version"] = utility.Version.version()

    return output_record


register_container_traversal(to_output_record, lambda d: to_output_record(to_packed(d)))
