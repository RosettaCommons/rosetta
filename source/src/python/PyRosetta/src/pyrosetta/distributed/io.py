"""IO routines operating on PackedPose representations."""
import bz2
import functools
import gzip
import os
import sys
import warnings

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

try:
    import lzma as xz
except ImportError:
    pass


__all__ = [
    "pack_result",
    "pose_result",
    "pose_from_file",
    "pose_from_sequence",
    "pose_from_pdbstring",
    "pose_from_pdb",
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
    "to_output_record",
    "dump_file",
    "dump_pdb",
    "dump_scored_pdb",
    "dump_multimodel_pdb",
    "dump_cif",
    "dump_mmtf",
    "dump_base64",
    "dump_pickle",
    "create_score_function",
    "get_fa_scorefxn",
    "get_score_function",
]

# Input functions

@requires_init
@pack_result
def pose_from_file(*args, **kwargs):
    """
    Uses the input filename from `*args` or `**kwargs` and returns a `PackedPose` object from it,
    deserializing:
        - bz2-encoded files ending with file extensions: (".pdb.bz2", ".bz2")
        - gzip-encoded files ending with file extensions: (".pdb.gz", ".gz")
        - xz-encoded files ending with file extensions: (".pdb.xz", ".xz")
        - base64-encoded files ending with file extension: ".b64_pose"
            *Warning*: This function uses the pickle module to deserialize ".b64_pose" files.
            Use of the pickle module is not secure. Only unpickle data you trust.
        - pickle-encoded files ending with file extension: ".pkl_pose"
            *Warning*: This function uses the pickle module to deserialize ".pkl_pose" files.
            Use of the pickle module is not secure. Only unpickle data you trust.
    Otherwise, implements `io.to_packed(pyrosetta.io.pose_from_file(*args, **kwargs))`.

    @klimaj
    """
    try:
        filename = kwargs.get("filename", False) \
            or next(iter(filter(lambda arg: isinstance(arg, str) and os.path.isfile, args)))
    except:
        raise FileNotFoundError(
            f"Could not find filename in arguments '{args}' or keyword arguments '{kwargs}'."
        )

    if filename.endswith(".b64_pose"):
        pack_or_pose = pose_from_base64(filename)  # returns `PackedPose` object
    elif filename.endswith(".pkl_pose"):
        pack_or_pose = pose_from_pickle(filename)  # returns `PackedPose` object
    else:
        pack_or_pose = pyrosetta.io.pose_from_file(*args, **kwargs)  # returns `Pose` object

    return pack_or_pose


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


@requires_init
@pack_result
@functools.wraps(pyrosetta.io.pose_from_pdbstring)
def pose_from_pdbstring(*args, **kwargs):
    return pyrosetta.io.pose_from_pdbstring(*args, **kwargs)


@functools.singledispatch
@requires_init
@pack_result
@functools.wraps(pose_from_file)
def pose_from_pdb(*args, **kwargs):
    return pose_from_file(*args, **kwargs)

@pose_from_pdb.register(type(None))
def _pose_from_none(none):
    return None


@functools.singledispatch
def pose_from_base64(filename):
    """
    Load a `PackedPose` object from a base64-encoded file.

    *Warning*: This function uses the pickle module to deserialize the input file.
    Use of the pickle module is not secure. Only unpickle data you trust.

    To load a `PackedPose` object from an input base64-encoded string,
    use `io.to_packed(string)` or `io.to_packed(io.to_pose(string))`.

    @klimaj
    """
    raise FileNotFoundError(
        f"The input filename must be an instance of `str`. Recieved: {type(filename)}"
    )

@pose_from_base64.register(type(None))
def _pose_from_none(none):
    return None

@pose_from_base64.register(str)
@requires_init
@pack_result
@pose_result
def _pose_from_str(filename):
    with open(filename, "r") as f:
        return f.read()


@functools.singledispatch
def pose_from_pickle(filename):
    """
    Load a `PackedPose` object from a pickle-encoded binary file.

    *Warning*: This function uses the pickle module to deserialize the input file.
    Use of the pickle module is not secure. Only unpickle data you trust.

    To load a `PackedPose` object from an input pickle-encoded bytestring,
    use `io.to_packed(bytestring)` or `io.to_packed(io.to_pose(bytestring))`.

    @klimaj
    """
    raise FileNotFoundError(
        f"The input filename must be an instance of `str`. Recieved: {type(filename)}"
    )

@pose_from_pickle.register(type(None))
def _pose_from_none(none):
    return None

@pose_from_pickle.register(str)
@requires_init
@pack_result
@pose_result
def _pose_from_str(filename):
    with open(filename, "rb") as f:
        return f.read()


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
    return pyrosetta.io.to_pdbstring(to_pose(inp))


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
    if not output_filename.endswith((".cif", ".mmcif")):
        warnings.warn(
            "Output filename does not end with '.cif' or '.mmcif', "
            + "which `pyrosetta.distributed.io.pose_from_file` expects."
        )
    return pyrosetta.io.dump_cif(to_pose(inp), output_filename)


@requires_init
def dump_mmtf(inp, output_filename):
    """Dump a MMTF file from a `PackedPose` or `Pose` object and output filename.

    @klimaj
    """
    if not output_filename.endswith(".mmtf"):
        warnings.warn(
            "Output filename does not end with '.mmtf', "
            + "which `pyrosetta.distributed.io.pose_from_file` expects."
        )
    return pyrosetta.io.dump_mmtf(to_pose(inp), output_filename)


@requires_init
def dump_base64(inp, output_filename):
    """Dump a base64-encoded file from a `PackedPose` or `Pose` object and output filename.

    @klimaj
    """
    if not output_filename.endswith(".b64_pose"):
        warnings.warn(
            "Output filename does not end with '.b64_pose', "
            + "which `pyrosetta.distributed.io.pose_from_file` expects."
        )
    with open(output_filename, "w") as f:
        f.write(to_base64(inp))
    return True


@requires_init
def dump_pickle(inp, output_filename):
    """Dump a pickle-encoded file from a `PackedPose` or `Pose` object and output filename.

    @klimaj
    """
    if not output_filename.endswith(".pkl_pose"):
        warnings.warn(
            "Output filename does not end with '.pkl_pose', "
            + "which `pyrosetta.distributed.io.pose_from_file` expects."
        )
    with open(output_filename, "wb") as f:
        f.write(to_pickle(inp))
    return True


def create_score_function(*args, **kwargs):
    """
    Returns a `ScorePoseTask` instance.
    Input `*args` and `**kwargs` are passed to `ScorePoseTask`.

    @klimaj
    """
    return score.ScorePoseTask(*args, **kwargs)


def get_fa_scorefxn():
    """
    Returns a `ScorePoseTask` instance with `weights`
    from `pyrosetta.rosetta.protocols.loops.get_fa_scorefxn()`.

    @klimaj
    """
    weights = pyrosetta.io.get_fa_scorefxn().get_name()
    return score.ScorePoseTask(weights=weights, patch=None)


def get_score_function():
    """
    Returns a `ScorePoseTask` instance with `weights` and `patch`
    from `pyrosetta.rosetta.core.scoring.get_score_function()`.

    @klimaj
    """
    return score.ScorePoseTask(weights=None, patch=None)
