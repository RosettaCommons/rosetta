# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

import bz2
import functools
import gzip
import pyrosetta.rosetta as rosetta
import sys
import warnings

from pyrosetta.rosetta.core.import_pose import pose_from_file, pose_from_pdbstring
from pyrosetta.rosetta.core.io.mmcif import dump_cif
from pyrosetta.rosetta.core.io.mmtf import dump_mmtf
from pyrosetta.rosetta.core.pose import make_pose_from_sequence, Pose
from pyrosetta.io.silent_file_map import SilentFileMap

try:
    import lzma as xz
except ImportError:
    pass


# for backward-compatibility
# This needs to be here because of all the people using the Workshops and
# because otherwise, it wrecks a lot of people's scripts.  ~Labonte
def pose_from_pdb(filename):
    return pose_from_file(filename)


def pose_from_sequence(seq, res_type="fa_standard", auto_termini=True):
    """
    Returns a Pose object generated from a single-letter sequence of amino acid
    residues in <seq> using the <res_type> ResidueType and creates N- and C-
    termini if <auto_termini> is set to True.

    Unlike make_pose_from_sequence(), this method generates a default PDBInfo
    and sets all torsion angles to 180 degrees.

    Example:
        pose = pose_from_sequence("THANKSEVAN")
    See also:
        Pose
        make_pose_from_sequence()
        pose_from_file()
        pose_from_rcsb()
    """
    pose = Pose()
    make_pose_from_sequence(pose, seq, res_type, auto_termini)

    for i in range(0, pose.total_residue()):
        res = pose.residue(i + 1)
        if not res.is_protein() or res.is_peptoid() or res.is_carbohydrate():
            continue

        pose.set_phi(i + 1, 180)
        pose.set_psi(i + 1, 180)
        pose.set_omega(i + 1, 180)

    # Empty PDBInfo (rosetta.core.pose.PDBInfo()) is not correct here;
    # we have to reserve space for atoms....
    pose.pdb_info(rosetta.core.pose.PDBInfo(pose))
    pose.pdb_info().name(seq[:8])

    return pose


def poses_from_files(objs, *args, **kwargs):
    """
    Returns an iterator object which is composed of Pose objects from input files.

    Example:
    import glob
    poses = pyrosetta.io.poses_from_files(glob.glob("path/to/pdbs/*.pdb"))
    """
    for obj in objs:
        yield pose_from_file(obj, *args, **kwargs)


def poses_from_sequences(objs, *args, **kwargs):
    """
    Returns an iterator object which is composed of Pose objects with input sequences.

    Example:
    poses = pyrosetta.io.poses_from_sequences(["TEST", "TESTING [ATP]", "ACDEFGHIKLMNPQRSTVWY"])
    """
    for obj in objs:
        yield pose_from_sequence(obj, *args, **kwargs)


def poses_from_silent(silent_filename):
    """Returns an Iterator object which is composed of Pose objects from a silent file.

    @atom-moyer
    """
    sfd = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd.read_file(silent_filename)
    for tag in sfd.tags():
        ss = sfd.get_structure(tag)
        pose = Pose()
        ss.fill_pose(pose)
        pose.pdb_info().name(tag)
        yield pose


def poses_from_multimodel_pdb(filename, *args, **kwargs):
    """Returns an Iterator object which is composed of Pose objects from an input
    multimodel PDB file. Poses are generated in the order in which they appear in
    the input multimodel PDB file. Any models without lines are skipped with a warning.
    Additional input `*args` and `**kwargs` are passed to `pose_from_pdbstring`.

    Example:
    poses = pyrosetta.io.poses_from_multimodel_pdb(filename)
    for pose in poses:
        name = pose.pdb_info().name()

    @klimaj
    """
    with open(filename, "r") as f:
        found_model = False
        for line in f:
            if line.startswith("MODEL"):
                found_model = True
                model_lines = []
                model_number = line.split()[-1]
            elif line.startswith("ENDMDL"):
                if not found_model:
                    raise ValueError(f"Multimodel PDB file is not formatted correctly: {filename}")
                if model_lines:
                    pose = Pose()
                    pdbstring = "".join(model_lines)
                    pose_from_pdbstring(pose, pdbstring, *args, **kwargs)
                    pose.pdb_info().name(f"{filename} MODEL {model_number}")
                    yield pose
                else: # Otherwise raises: "ERROR: Assertion `! lines.empty()` failed."
                    warnings.warn(f"Skipping empty model {model_number} in multimodel PDB file: {filename}")
                found_model = False
            elif found_model:
                model_lines.append(line)


def poses_to_silent(poses, output_filename):
    """Takes a Pose or list of poses and outputs them as a binary silent file.
    This method requires a Pose object.
    If you are using a PackedPose, use pyrosetta.distributed.io.to_silent()

    Inputs:
    poses: Pose or list of poses. This function automatically detects which one.
    output_filename: The desired name of the output silent file.

    Example:
    poses_to_silent(poses, "mydesigns.silent")

    The decoy name written to your silent file is take from pose.pdb_info().name()
    To set a different decoy name, change it in your pose before calling this function.
    Example:
    pose.pdb_info().name("my_tag.pdb")

    @srgerb
    """

    silentOptions = rosetta.core.io.silent.SilentFileOptions()
    silentOptions.set_binary_output(True)
    silentFile = rosetta.core.io.silent.SilentFileData(silentOptions)
    silentStruct = rosetta.core.io.silent.BinarySilentStruct(silentOptions)

    def output_silent(pose):
        decoy_tag = pose.pdb_info().name()
        silentStruct.fill_struct(pose, tag=decoy_tag)
        silentFile.write_silent_struct(silentStruct, filename=output_filename)

    if isinstance(poses, (list, tuple, set)):
        for pose in poses:
            output_silent(pose=pose)
    else:
        output_silent(pose=poses)


def to_pdbstring(pose):
    """Convert to pdb-formatted string with score and energy data.
    """
    from pyrosetta.rosetta.core.io import StructFileRepOptions
    from pyrosetta.rosetta.core.io.pdb import dump_pdb as _dump_pdb
    from pyrosetta.rosetta.std import ostringstream

    sfro = StructFileRepOptions()
    sfro.set_output_pose_cache_data(True)
    sfro.set_output_pose_energies_table(True)

    oss = ostringstream()
    _dump_pdb(pose, oss, sfro)

    return oss.bytes().decode()


@functools.wraps(Pose.dump_file)
def dump_file(pose, output_filename):
    return pose.dump_file(output_filename)


@functools.wraps(Pose.dump_scored_pdb)
def dump_scored_pdb(pose, output_filename, scorefxn):
    return pose.dump_scored_pdb(output_filename, scorefxn)


def dump_pdb(pose, output_filename):
    """
    Dump a PDB file from a `Pose` object and output filename.
    If the output filename ends with ".pdb.bz2" or ".bz2", then dump a bz2-encoded PDB file.
    If the output filename ends with ".pdb.gz" or ".gz", then dump a gzip-encoded PDB file.
    If the output filename ends with ".pdb.xz" or ".xz", then dump a xz-encoded PDB file.

    @klimaj
    """
    from pyrosetta.rosetta.core.io.pdb import dump_pdb as _dump_pdb

    if not output_filename.endswith((".pdb", ".pdb.bz2", ".bz2", ".pdb.gz", ".gz", ".pdb.xz", ".xz")):
        warnings.warn(
            "Output filename does not end with '.pdb', '.pdb.bz2', '.bz2', '.pdb.gz', '.gz', '.pdb.xz', or '.xz', "
            + "which `pyrosetta.io.pose_from_file` expects."
        )

    if output_filename.endswith((".pdb.bz2", ".bz2")):
        with open(output_filename, "wb") as f:
            f.write(bz2.compress(str.encode(to_pdbstring(pose))))
    elif output_filename.endswith((".pdb.gz", ".gz")):
        with gzip.open(output_filename, mode="wt", compresslevel=9) as gz:
            gz.write(to_pdbstring(pose))
    elif output_filename.endswith((".pdb.xz", ".xz")):
        if "lzma" not in sys.modules:
            raise ImportError(
                (
                    "Using 'xz' for compression requires installing the 'xz' package into your python environment. "
                    + "For installation instructions, visit:\n"
                    + "https://anaconda.org/anaconda/xz\n"
                )
            )
        with open(output_filename, "wb") as f:
            f.write(xz.compress(str.encode(to_pdbstring(pose))))
    else:
        return _dump_pdb(pose, output_filename)

    return True


def dump_multimodel_pdb(poses, output_filename):
    """
    Dump a multimodel PDB file from a `list`, `tuple`, or `set` of `Pose` objects,
    and an output filename.

    Inputs:
    poses: `Pose` or iterable of `Pose` objects. This function automatically detects which one.
    output_filename: The desired name of the output PDB file.

    Example:
    dump_multimodel_pdb(poses, "designs.pdb")

    @klimaj
    """
    if not output_filename.endswith(".pdb"):
        warnings.warn("Output filename does not end with '.pdb'.")
    v1 = rosetta.utility.vector1_std_shared_ptr_const_core_pose_Pose_t()
    if isinstance(poses, (list, tuple, set)):
        for p in poses:
            v1.append(p)
    else:
        v1.append(poses)

    return rosetta.core.io.pdb.dump_multimodel_pdb(v1, output_filename)
