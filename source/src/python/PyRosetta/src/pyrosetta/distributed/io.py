"""IO routines operating on PackedPose representations."""
import bz2
import functools
import os

import pyrosetta.io
import pyrosetta.rosetta.utility as utility
import pyrosetta.rosetta.core.pose as pose
import pyrosetta.rosetta.core.import_pose as import_pose
import pyrosetta.rosetta.core.scoring as scoring
import pyrosetta.distributed.tasks.score as score

from pyrosetta.distributed import requires_init
from pyrosetta.distributed.packed_pose import (
    pack_result, pose_result, to_pose, to_packed, to_dict, to_base64, to_pickle, register_container_traversal
)

__all__ = [
    "pack_result",
    "pose_result",
    "pose_from_file",
    "pose_from_sequence",
    "pose_from_pdbstring",
    "pose_from_pdb_bz2",
    "pose_from_base64",
    "pose_from_pickle",
    "poses_from_files",
    "poses_from_sequences",
    "poses_from_multimodel_pdb",
    "poses_from_silent",
    "to_base64",
    "to_pickle",
    "to_silent",
    "to_pdbstring",
    "to_output_record"
    "dump_file",
    "dump_pdb",
    "dump_scored_pdb",
    "dump_multimodel_pdb",
    "dump_cif",
    "dump_mmtf",
    "dump_pdb_bz2",
    "dump_base64",
    "dump_pickle",
]

# Input functions

@requires_init
@pack_result
def pose_from_file(*args, **kwargs):
    """
    Uses the input filename from `*args` or `**kwargs` and returns a `PackedPose` object from it,
    deserializing:
        - bz2-encoded files ending with file extensions: (".pdb.bz2", ".bz2")
        - base64-encoded files ending with file extensions: (".base64", ".b64", ".B64", ".pose")
        - pickle-encoded files ending with file extensions: (".pickle", ".pickled_pose")
    Otherwise, implements `io.to_packed(pyrosetta.io.pose_from_file(*args, **kwargs))`.

    @klimaj
    """
    try:
        filename = kwargs.get("filename", False) or next(iter(filter(os.path.isfile, args)))
    except:
        raise FileNotFoundError(
            f"Could not find filename in arguments '{args}' or keyword arguments '{kwargs}'."
        )

    if filename.endswith((".pdb.bz2", ".bz2")):
        pose = pose_from_pdb_bz2(filename)
    elif filename.endswith((".base64", ".b64", ".B64", ".pose")):
        pose = pose_from_base64(filename)
    elif filename.endswith((".pickle", ".pickled_pose")):
        pose = pose_from_pickle(filename)
    else:
        pose = pyrosetta.io.pose_from_file(*args, **kwargs)

    return pose

pose_from_sequence = requires_init(pack_result(
    pyrosetta.io.pose_from_sequence))
poses_from_files = requires_init(pack_result(
    pyrosetta.io.poses_from_files))
poses_from_sequences = requires_init(pack_result(
    pyrosetta.io.poses_from_sequences))
poses_from_silent = requires_init(pack_result(
    pyrosetta.io.poses_from_silent))
poses_from_multimodel_pdb = requires_init(pack_result(
    pyrosetta.io.poses_from_multimodel_pdb))


@functools.wraps(import_pose.pose_from_pdbstring)
@requires_init
@pack_result
def pose_from_pdbstring(*args, **kwargs):
    result = pose.Pose()
    import_pose.pose_from_pdbstring(result, *args, **kwargs)
    return result


@functools.singledispatch
@requires_init
def pose_from_pdb_bz2(filename):
    """Load a `PackedPose` object from a bz2-encoded PDB file.

    @klimaj
    """
    with open(filename, "rb") as f:
        pdbstring = bz2.decompress(f.read()).decode()
    return pose_from_pdbstring(pdbstring)

@pose_from_pdb_bz2.register(type(None))
def pose_from_none(none):
    return None


@functools.singledispatch
@requires_init
@pack_result
@pose_result
def pose_from_base64(filename):
    """
    Load a `PackedPose` object from a base64-encoded file.

    To load a `PackedPose` object from an input base64-encoded string,
    use `io.to_packed(string)` or `io.to_packed(io.to_pose(string))`.

    @klimaj
    """
    with open(filename, "r") as f:
        return f.read()

@pose_from_base64.register(type(None))
def pose_from_none(none):
    return None


@functools.singledispatch
@requires_init
@pack_result
@pose_result
def pose_from_pickle(filename):
    """
    Load a `PackedPose` object from a pickle-encoded binary file.

    To load a `PackedPose` object from an input pickle-encoded bytestring,
    use `io.to_packed(bytestring)` or `io.to_packed(io.to_pose(bytestring))`.

    @klimaj
    """
    with open(filename, "rb") as f:
        return f.read()

@pose_from_pickle.register(type(None))
def pose_from_none(none):
    return None


# Output functions

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
    from pyrosetta.rosetta.core.io import StructFileRepOptions
    from pyrosetta.rosetta.std import ostringstream

    sfro = StructFileRepOptions()
    sfro.set_output_pose_cache_data(True)
    sfro.set_output_pose_energies_table(True)

    oss = ostringstream()
    pyrosetta.io.dump_pdb(to_pose(inp), oss, sfro)

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


@requires_init
def dump_file(inp, output_filename):
    """Dump a file from a `PackedPose` or `Pose` object and output filename.

    @klimaj
    """
    return pyrosetta.io.dump_file(to_pose(inp), output_filename)

@requires_init
def dump_pdb(inp, output_filename):
    """Dump a PDB file from a `PackedPose` or `Pose` object and output filename.

    @klimaj
    """
    return pyrosetta.io.dump_pdb(to_pose(inp), output_filename)

@requires_init
def dump_scored_pdb(inp, output_filename, scorefxn):
    """
    Dump a scored PDB file from a `PackedPose` or `Pose` object, output filename
    and score function. The score function may be a `str` object representing the
    weights passed to `ScorePoseTask`, a `ScorePoseTask` instance, or a
    `ScoreFunction` instance.

    @klimaj
    """
    if isinstance(scorefxn, str):
        scorefxn = score.ScorePoseTask(weights=scorefxn, patch=None)
    if isinstance(scorefxn, score.ScorePoseTask):
        return pyrosetta.io.dump_pdb(scorefxn(inp).pose, output_filename)
    elif isinstance(scorefxn, scoring.ScoreFunction):
        return pyrosetta.io.dump_scored_pdb(to_pose(inp), output_filename, scorefxn)
    else:
        raise ValueError(f"Unsupported argument parameter type for 'scorefxn': {type(scorefxn)}")

@requires_init
def dump_multimodel_pdb(inp, output_filename):
    """
    Dump a multimodel PDB file from a `list`, `tuple`, or `set` of `PackedPose`
    or `Pose` objects, and an output filename.

    @klimaj
    """
    if isinstance(inp, (list, tuple, set)):
        poses = [to_pose(p) for p in inp]
    else:
        poses = to_pose(inp)
    return pyrosetta.io.dump_multimodel_pdb(poses, output_filename)

@requires_init
def dump_cif(inp, output_filename):
    """Dump a CIF file from a `PackedPose` or `Pose` object and output filename.

    @klimaj
    """
    return pyrosetta.io.dump_cif(to_pose(inp), output_filename)

@requires_init
def dump_mmtf(inp, output_filename):
    """Dump a MMTF file from a `PackedPose` or `Pose` object and output filename.

    @klimaj
    """
    return pyrosetta.io.dump_mmtf(to_pose(inp), output_filename)

@requires_init
def dump_pdb_bz2(inp, output_filename):
    """Dump a bz2-encoded PDB file from a `PackedPose` or `Pose` object and output filename.

    @klimaj
    """
    with open(output_filename, "wb") as f:
        f.write(bz2.compress(str.encode(to_pdbstring(inp))))
    return True

@requires_init
def dump_base64(inp, output_filename):
    """Dump a base64-encoded file from a `PackedPose` or `Pose` object and output filename.

    @klimaj
    """
    with open(output_filename, "w") as f:
        f.write(to_base64(inp))
    return True

@requires_init
def dump_pickle(inp, output_filename):
    """Dump a pickle-encoded file from a `PackedPose` or `Pose` object and output filename.

    @klimaj
    """
    with open(output_filename, "wb") as f:
        f.write(to_pickle(inp))
    return True
