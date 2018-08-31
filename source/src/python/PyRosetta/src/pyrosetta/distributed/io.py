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
    "pose_from_pdbstring",
    "poses_from_silent"
]

# Input functions

pose_from_file = requires_init(pack_result(
    pyrosetta.io.pose_from_file))
pose_from_sequence = requires_init(pack_result(
    pyrosetta.io.pose_from_sequence))


@requires_init
@pack_result
def poses_from_dir(dirname, loglevel=logging.WARN):
    """Construct and return PackedPoses for every PDB file in a directory.

    Args:
        dirname (str): name of the directory that holds the PDB files to use.

    Yields:
       pyrosetta.distributed.packed_pose.PackedPose: PackedPose instance for each PDB file in the directory.
    """
    from os import listdir, path
    import pyrosetta.distributed.utility.log
    with pyrosetta.distributed.utility.log.LoggingContext(logging.getLogger("rosetta"), level=loglevel):
        yield from (pose_from_file(path.join(dirname, fn)) for fn in listdir(dirname) if fn.endswith(".pdb"))


@functools.wraps(import_pose.pose_from_pdbstring)
@requires_init
@pack_result
def pose_from_pdbstring(*args, **kwargs):
    result = pose.Pose()
    import_pose.pose_from_pdbstring(result, *args, **kwargs)
    return result


poses_from_silent = requires_init(pack_result(pyrosetta.io.poses_from_silent))


@functools.singledispatch
@requires_init
def to_pdbstring(inp):
    """Convert to pdb-formatted string with score and energy data."""
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
    """Convert to an "archive" output record with scores, pdb string and packed string."""

    output_record = dict(to_dict(inp))
    output_record["pdb_pose"] = to_pdbstring(inp)
    output_record["pickled_version"] = utility.Version.version()

    return output_record


register_container_traversal(to_output_record, lambda d: to_output_record(to_packed(d)))
