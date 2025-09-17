# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import argparse
import glob
import json
import os
import pyrosetta

POSE_SEQUENCE = "DATA"
POSE_SCORE = 1 + 0.5j


def truncate_str(value, _index=50):
    return value[:_index] + "..."

def truncate_init_options(init_options):
    trunc_init_options = {}
    for key, values in init_options.items():
        trunc_values = []
        for value in values:
            if isinstance(value, dict):
                data = {}
                for k1, v1 in value.items():
                    if isinstance(v1, dict):
                        data[k1] = {k2: truncate_str(v2) for k2, v2 in v1.items()}
                    else:
                        data[k1] = truncate_str(v1)
                trunc_values.append(data)
            else:
                trunc_values.append(value)
        trunc_init_options[key] = trunc_values

    return trunc_init_options

def main(tmp_dir):
    os.chdir(tmp_dir)
    pdb_files = sorted(glob.glob(os.path.join(tmp_dir, "tmp_*.pdb")) + glob.glob(os.path.join(tmp_dir, "tmp_*.pdb.gz")))
    assert len(pdb_files) > 0, "PDB files do not exist."
    list_file = os.path.join(tmp_dir, "my_file.list")
    assert os.path.isfile(list_file), "List file does not exist."
    extra_res_fa_files = sorted(glob.glob(os.path.join(tmp_dir, "*.params")))
    assert len(extra_res_fa_files) > 0, "Extra residue files do not exist."
    patch_files = sorted(glob.glob(os.path.join(tmp_dir, "*.txt")))
    assert len(patch_files) > 0, "Patch files do not exist."
    bcl_dir = os.path.join(tmp_dir, "bcl_rosetta")
    os.makedirs(bcl_dir, exist_ok=False)
    if not pyrosetta.rosetta.basic.was_init_called():
        pyrosetta.init(
            options=(
                "-run:constant_seed 1 "
                "-run:jran 1234567 "
                "-override_database_params 1 "
                "-load_PDB_components 0 "
                "-out:levels core.init:0 basic.random.init_random_generator:0 "
                "-parser:script_vars foo=bar baz=12.34"
            ),
            extra_options="-s {0} -l {1} -extra_res_fa {2} -extra_res_fa {3} -extra_patch_fa {4} -bcl {5}".format(
                " ".join(pdb_files),
                list_file,
                " ".join(extra_res_fa_files[:-1]),
                extra_res_fa_files[-1],
                " ".join(patch_files),
                bcl_dir,
            ),
            set_logging_handler="logging",
            notebook=None,
            silent=True,
        )
    else:
        raise RuntimeError("PyRosetta is already initialized.")

    init_options = pyrosetta.get_init_options(compressed=False, as_dict=False)
    assert isinstance(init_options, str)
    assert "-in:path:database" in init_options
    print("Uncompressed PyRosetta initialization options as `str`:", init_options, sep=os.linesep)

    init_options = pyrosetta.get_init_options(compressed=False, as_dict=True)
    assert isinstance(init_options, dict)
    assert "in:path:database" in init_options.keys()
    print("Uncompressed PyRosetta initialization options as `dict`:", init_options, sep=os.linesep)

    init_options = pyrosetta.get_init_options(compressed=True, as_dict=True)
    assert isinstance(init_options, dict)
    assert "in:path:database" in init_options.keys()
    print("Compressed PyRosetta initialization options as `dict`:", truncate_init_options(init_options), sep=os.linesep)

    try:
        pyrosetta.get_init_options(compressed=True, as_dict=False)
        ex = None
    except NotImplementedError as e:
        ex = e
    finally:
        if ex is None:
            raise RuntimeError(f"Did not catch `NotImplementedError`.")
        else:
            print(f"Successfully caught `{type(ex).__name__}: {ex}`")

    init_file = os.path.join(tmp_dir, "my.init")
    metadata = [
        "Contains IGU, GNP, CYX, CED, R1A, and T3P Rosetta residue topology files and the 3prime5prime_methyl_phosphate patch file",
        "Also contains an input pose for my project."
        "Version 2.0",
    ]
    pyrosetta.dump_init_file(
        init_file,
        poses=None,
        author="Username",
        email="test@example",
        license="LICENSE.PyRosetta.md",
        metadata=metadata,
        overwrite=False,
        dry_run=True,
        verbose=True,
    )
    assert not os.path.isfile(init_file), "'.init' file was dumped with `dry_run=True`."

    pose = pyrosetta.pose_from_sequence(POSE_SEQUENCE)
    pose.cache["foo"] = POSE_SCORE
    pyrosetta.dump_init_file(
        init_file,
        poses=pose,
        author="Username",
        email="test@example",
        license="LICENSE.PyRosetta.md",
        metadata=metadata,
        overwrite=False,
        dry_run=False,
        verbose=True,
    )

    pose = pyrosetta.Pose()
    base_res_set = pose.conformation().modifiable_residue_type_set_for_conf().base_residue_types()
    name3_set = set(base_res_set.pop().name3() for _ in range(base_res_set.capacity()))
    with open(os.path.join(tmp_dir, "res_types.json"), "w") as f:
        json.dump(list(name3_set), f)

    _print_init_file_contents = False
    if _print_init_file_contents:
        _delim = "*" * 100
        with open(init_file, "r") as f:
            print(_delim, f.read(), _delim, sep=os.linesep)


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--tmp_dir', type=str)
    args = parser.parse_args()
    main(args.tmp_dir)
