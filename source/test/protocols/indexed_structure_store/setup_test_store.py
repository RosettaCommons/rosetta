from os import path
import shutil

import numpy
import numpy.lib.recfunctions

import tables

import rosetta; rosetta.init()

import interface_fragment_matching
from interface_fragment_matching.structure_database import StructureDatabase

import fragment_fitting
from fragment_fitting import FragmentDatabase, FragmentSpecification

test_pose = rosetta.pose_from_pdb("test_structure.pdb")
test_pose_residues = StructureDatabase.extract_residue_entries_from_pose(test_pose)

test_fragment_spec = FragmentSpecification(5, ("N", "CA", "C"))
_, test_fragments = test_fragment_spec.fragments_from_source_residues(test_pose_residues)

# append_fields does not propery handle structured arrays with nested structures, convert to coordinate array
test_fragments = numpy.lib.recfunctions.append_fields(
        test_fragment_spec.fragment_data_to_coordinate(test_fragments),
        "threshold_distance",
        numpy.repeat([.1], len(test_fragments)),
        dtypes=float, usemask=False)

if path.exists("test_store"):
    shutil.rmtree("test_store")

with FragmentDatabase(tables.openFile("test_store.h5", "w")) as db:
    db.setup()
    db.add_fragments("full_test_fragments", test_fragments, test_fragment_spec)
    db.write_binary_fragments("full_test_fragments", "test_store")

    db.add_fragments("partial_test_fragments", test_fragments[:-1], test_fragment_spec)
    db.write_binary_fragments("partial_test_fragments", "test_store")

