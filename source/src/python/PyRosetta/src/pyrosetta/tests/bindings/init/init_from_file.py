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
import sys

sys.path.insert(0, os.path.dirname(__file__))
try:
    from write_test_files import DATABASE_FILES
except ImportError as ex:
    raise ImportError(ex)


def get_mode(filename):
    return "rb" if filename.endswith(".gz") else "r"

def filter_pdb_rotamers_line(string):
    return os.linesep.join(filter(lambda x: not x.startswith("PDB_ROTAMERS"), string.splitlines()))

def get_pdb_rotamers_file(string, parent_dir=None):
    for line in reversed(string.splitlines()):
        if line.startswith("PDB_ROTAMERS"):
            entry = line.split("PDB_ROTAMERS ")[-1].strip()
            if parent_dir is not None and os.path.isfile(os.path.join(parent_dir, entry)):
                file = os.path.join(parent_dir, entry)
                break
            elif os.path.isfile(entry):
                file = entry
                break
            else:
                raise FileNotFoundError(entry)
    else:
        file = None
    return file

def main(tmp_dir):
    os.chdir(tmp_dir)
    init_file = os.path.join(tmp_dir, "my.init")

    init_options = pyrosetta.get_init_options_from_file(
        init_file,
        dry_run=True,
        output_dir=None,
        relative_paths=False,
        max_decompressed_bytes=None,
        database=None,
        as_dict=True,
    )
    assert isinstance(init_options, dict)
    assert "in:path:database" in init_options
    assert not pyrosetta.rosetta.basic.was_init_called(), "PyRosetta was initialized with `pyrosetta.get_init_options_from_file`"
    print("Dry run PyRosetta initialization from '.init' file options with absolute paths as `dict`:", init_options, sep=os.linesep)

    init_options = pyrosetta.get_init_options_from_file(
        init_file,
        dry_run=True,
        output_dir=None,
        relative_paths=True,
        max_decompressed_bytes=1_000_000_000,
        database=None,
        as_dict=False,
    )
    assert isinstance(init_options, str)
    assert "-in:path:database" in init_options
    assert not pyrosetta.rosetta.basic.was_init_called(), "PyRosetta was initialized with `pyrosetta.get_init_options_from_file`"
    print("Dry run PyRosetta initialization from '.init' file options with relative paths as `str`:", init_options, sep=os.linesep)

    tmp_init_dir = os.path.join(tmp_dir, "tmp1_pyrosetta_init_files")
    init_options = pyrosetta.get_init_options_from_file(
        init_file,
        dry_run=False,
        output_dir=tmp_init_dir,
        relative_paths=True,
        max_decompressed_bytes=10_000_000,
        database=None,
        as_dict=True,
    )
    assert isinstance(init_options, dict)
    assert "in:path:database" in init_options
    assert os.listdir(tmp_init_dir) != []
    assert not pyrosetta.rosetta.basic.was_init_called(), "PyRosetta was initialized with `pyrosetta.get_init_options_from_file`"
    print("PyRosetta initialization from '.init' file options with relative paths as `dict`:", init_options, sep=os.linesep)

    tmp_init_dir = os.path.join(tmp_dir, "tmp2_pyrosetta_init_files")
    init_options = pyrosetta.get_init_options_from_file(
        init_file,
        dry_run=False,
        output_dir=tmp_init_dir,
        relative_paths=None,
        max_decompressed_bytes=1_000_000,
        database=None,
        as_dict=False,
    )
    assert isinstance(init_options, str)
    assert "-in:path:database" in init_options
    assert os.listdir(tmp_init_dir) != []
    assert not pyrosetta.rosetta.basic.was_init_called(), "PyRosetta was initialized with `pyrosetta.get_init_options_from_file`"
    print("PyRosetta initialization from '.init' file options with absolute paths as `str`:", init_options, sep=os.linesep)


    init_dir = os.path.join(tmp_dir, "pyrosetta_init_files")
    try:
        pyrosetta.init_from_file(
            init_file,
            output_dir=init_dir,
            skip_corrections=True,
            relative_paths=True,
            max_decompressed_bytes=1_000, # Testing 1 KB decompressed buffer size limit
            dry_run=True,
            database=None,
            set_logging_handler=None,
            notebook=None,
            silent=False,
        )
        raise RuntimeError("`pyrosetta.init_from_file` did not raise `BufferError` with `max_decompressed_bytes=1_000`.")
    except BufferError as ex:
        print(f"Successfully caught `BufferError` with `max_decompressed_bytes=1_000`: {type(ex).__name__}: {ex}")
    finally:
        assert not pyrosetta.rosetta.basic.was_init_called(), "PyRosetta was initialized with `dry_run=True`"

    pyrosetta.init_from_file(
        init_file,
        output_dir=init_dir,
        skip_corrections=True,
        relative_paths=True,
        max_decompressed_bytes=800_000,
        dry_run=True,
        database=None,
        set_logging_handler=None,
        notebook=None,
        silent=False,
    )
    assert not pyrosetta.rosetta.basic.was_init_called(), "PyRosetta was initialized with `dry_run=True`"

    pyrosetta.init_from_file(
        init_file,
        output_dir=init_dir,
        skip_corrections=False,
        relative_paths=False,
        max_decompressed_bytes=800_000,
        dry_run=True,
        database=os.path.relpath(pyrosetta._rosetta_database_from_env()),
        set_logging_handler=None,
        notebook=None,
        silent=False,
    )
    assert not pyrosetta.rosetta.basic.was_init_called(), "PyRosetta was initialized with `dry_run=True`"

    pyrosetta.init_from_file(
        init_file,
        output_dir=init_dir,
        skip_corrections=None,
        relative_paths=True,
        max_decompressed_bytes=800_000,
        dry_run=False,
        database=None,
        set_logging_handler=None,
        notebook=None,
        silent=False,
    )

    pose = pyrosetta.Pose()
    base_res_set = pose.conformation().modifiable_residue_type_set_for_conf().base_residue_types()
    name3_set = set(base_res_set.pop().name3() for _ in range(base_res_set.capacity()))
    with open(os.path.join(tmp_dir, "res_types.json"), "r") as f:
        original_name3_set = set(json.load(f))
    assert name3_set == original_name3_set, "Residue type sets are not identical."

    def print_identical_files(file1, file2, relpath=tmp_dir):
        print(
            "Files are identical:",
            f"'./{os.path.relpath(file1, start=relpath)}'",
            "==",
            f"'./{os.path.relpath(file2, start=relpath)}'",
        )

    file_counter = 0
    input_files = list(filter(os.path.isfile, glob.glob(os.path.join(init_dir, "**", "*"), recursive=True)))
    for input_file in input_files:
        if os.path.basename(input_file).startswith("tmp_"):
            original_file = os.path.join(tmp_dir, os.path.basename(input_file))
            assert os.path.isfile(original_file), f"Original file does not exist: {original_file}"
            mode = get_mode(input_file)
            with open(original_file, mode) as f1, open(input_file, mode) as f2:
                assert f1.read() == f2.read(), f"{original_file} != {input_file}"
                file_counter += 1
                print_identical_files(original_file, input_file)
        elif os.path.basename(input_file).endswith(".list"):
            original_file = os.path.join(tmp_dir, os.path.basename(input_file))
            assert os.path.isfile(original_file), f"Original file does not exist: {original_file}"
            with open(original_file, "r") as f1, open(input_file, "r") as f2:
                l1, l2 = f1.read().splitlines(), f2.read().splitlines()
                list_file_counter = 0
                assert len(l1) == len(l2), "List files do not contain an identical number of lines."
                for i in range(len(l1)):
                    original_listed_file, input_listed_file = l1[i], l2[i]
                    mode1, mode2 = get_mode(original_listed_file), get_mode(input_listed_file)
                    assert mode1 == mode2, f"File types are not identical: {original_listed_file} {input_listed_file}"
                    with open(original_listed_file, mode1) as f3, open(input_listed_file, mode2) as f4:
                        assert f3.read() == f4.read(), f"{original_listed_file} != {input_listed_file}"
                        print_identical_files(original_listed_file, input_listed_file)
                        file_counter += 1
                        list_file_counter += 1
                assert list_file_counter == len(l1) == len(l2), "List files are not identical."
                file_counter += 1
                print_identical_files(original_file, input_file)
        else:
            for database_file in filter(lambda f: f.endswith((".params", ".txt")), DATABASE_FILES):
                if os.path.basename(database_file) == os.path.basename(input_file):
                    original_file = os.path.join(tmp_dir, os.path.basename(database_file))
                    assert os.path.isfile(original_file), f"Original file does not exist: {original_file}"
                    with open(original_file, "r") as f1, open(input_file, "r") as f2:
                        s1, s2 = f1.read(), f2.read()
                        assert filter_pdb_rotamers_line(s1) == filter_pdb_rotamers_line(s2), f"{database_file} != {input_file}"
                        print_identical_files(original_file, input_file)
                        file_counter += 1
                        original_pdb_rotamers_file = get_pdb_rotamers_file(s1, parent_dir=tmp_dir)
                        input_pdb_rotamers_file = get_pdb_rotamers_file(s2, parent_dir=None)
                        if all(x is not None for x in (original_pdb_rotamers_file, input_pdb_rotamers_file)):
                            with open(original_pdb_rotamers_file, "r") as f3, open(input_pdb_rotamers_file, "r") as f4:
                                assert f3.read() == f4.read(), f"{original_pdb_rotamers_file} != {input_pdb_rotamers_file}"
                                print_identical_files(original_pdb_rotamers_file, input_pdb_rotamers_file)
                                file_counter += 1
    print(f"Successfully tested that {file_counter}/{len(input_files)} cached files are identical to their originals.")

    init_dir_original = os.path.join(tmp_dir, "pyrosetta_init_files")
    init_options_original = pyrosetta.get_init_options_from_file(
        init_file,
        dry_run=True,
        output_dir=init_dir_original,
        relative_paths=True,
        max_decompressed_bytes=800_000,
        database=None,
        as_dict=True,
    )
    init_options_reproduce = pyrosetta.get_init_options(compressed=False, as_dict=True)
    # from pprint import pprint
    # pprint(init_options_original)
    # pprint(init_options_reproduce)
    for option_name, original_values in init_options_original.items():
        reproduce_values = init_options_reproduce[option_name]
        assert len(reproduce_values) == len(original_values), f"{reproduce_values} != {original_values}"
        for i in range(len(original_values)):
            original_value, reproduce_value = original_values[i], reproduce_values[i]
            if all(map(os.path.isdir, (original_value, reproduce_value))):
                original_value_rel, reproduce_value_rel = map(os.path.relpath, (original_value, reproduce_value))
                assert reproduce_value_rel == original_value_rel, f"{reproduce_value_rel} != {original_value_rel}"
            else:
                assert reproduce_value == original_value, f"{reproduce_value} != {original_value}"
    print(f"Successfully tested that original '.init' file options and reproduced PyRosetta initialization options are identical.")


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--tmp_dir', type=str)
    args = parser.parse_args()
    main(args.tmp_dir)
