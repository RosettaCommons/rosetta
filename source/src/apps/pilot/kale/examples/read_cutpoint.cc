// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/kinematic_closure/types.hh>

// Utility headers
#include <devel/init.hh>
#include <utility/exit.hh>
#include <boost/foreach.hpp>
#include <iostream>

using namespace std;
using namespace core;
using core::conformation::Conformation;
using core::id::AtomID;
using core::import_pose::pose_from_file;
using core::kinematics::FoldTree;
using core::pose::Pose;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::loops::fold_tree_from_loops;
using protocols::kinematic_closure::CoordinateList;
using protocols::kinematic_closure::ParameterList;

AtomID id_from_index(Size index, Size first_residue) {
	return AtomID(
			((index - 1) % 3) + 1,
			((index - 1) / 3) + first_residue - 1);
}

int main(int argc, char** argv) {
	devel::init(argc, argv);

	//Pose pose; pose_from_file(pose, "structures/linear/6.1ubq.pdb", core::import_pose::PDB_file);
	//Loop pivots(3, 4, 4);
	Pose pose; pose_from_file(pose, "structures/kic/1srp.pdb", core::import_pose::PDB_file);
	Loop pivots(308, 319, 319);
	Loops pivots_wrapper; pivots_wrapper.add_loop(pivots);

	Size first_residue = pivots.start();
	Size last_residue = pivots.stop();
	Size num_residues = pivots.length();
	Size num_atoms = 3 * (num_residues + 2);
	Size num_bond_lengths = num_atoms - 1;
	Size num_bond_angles = num_atoms - 2;
	Size num_torsion_angles = num_atoms - 3;

	CoordinateList atom_xyzs(num_atoms);
	ParameterList bond_lengths(num_atoms, 0);
	ParameterList bond_angles(num_atoms, 0);
	ParameterList torsion_angles(num_atoms, 0);

	FoldTree tree; fold_tree_from_loops(pose, pivots_wrapper, tree, true);
	pose.fold_tree(tree);
	pose.fold_tree().show(cout);

	Conformation & conformation = pose.conformation();
	AtomID ids[4];

	ids[0] = AtomID(1, 319);
	ids[1] = AtomID(2, 319);
	ids[2] = AtomID(3, 319);
	ids[3] = AtomID(1, 320);

	//conformation.torsion_angle(ids[0], ids[1], ids[2], ids[3]);
	//conformation.set_torsion_angle(ids[0], ids[1], ids[2], ids[3], 0);

	conformation.bond_angle(ids[1], ids[2], ids[3]);
	conformation.set_bond_angle(ids[1], ids[2], ids[3], 0);

	return 0;

	// Read the cartesian coordinates.
	/*

	Size index = 1;

	for (Size i = first_residue - 1; i <= last_residue + 1; i++) {
		core::conformation::Residue const &residue = pose.residue(i);

		for (Size j = 1; j <= 3; j++) {
			atom_xyzs[index].resize(3);
			atom_xyzs[index][1] = residue.xyz(j).x();
			atom_xyzs[index][2] = residue.xyz(j).y();
			atom_xyzs[index][3] = residue.xyz(j).z();
			index++;
		}
	}

	// Read the internal coordinates.

	for (Size index = 1; index <= num_atoms; index++) {
		ids[0] = id_from_index(index + 0, first_residue);
		ids[1] = id_from_index(index + 1, first_residue);
		ids[2] = id_from_index(index + 2, first_residue);
		ids[3] = id_from_index(index + 3, first_residue);

		bond_lengths[index] = (index > num_bond_lengths) ?
			0 : conformation.bond_length(ids[0], ids[1]);

		bond_angles[(index % num_atoms) + 1] = (index > num_bond_angles) ?
			0 : conformation.bond_angle(ids[0], ids[1], ids[2]);

		torsion_angles[(index % num_atoms) + 1] = (index > num_torsion_angles) ?
			0 : conformation.torsion_angle(ids[0], ids[1], ids[2], ids[3]);
	}

	// Set the internal coordinates.

	for (Size index = 1; index <= num_bond_lengths; index++) {
		ids[0] = id_from_index(index + 0, first_residue);
		ids[1] = id_from_index(index + 1, first_residue);

		conformation.set_bond_length(
				ids[0], ids[1], bond_lengths[index]);
	}

	for (Size index = 1; index <= num_bond_angles; index++) {
		ids[0] = id_from_index(index + 0, first_residue);
		ids[1] = id_from_index(index + 1, first_residue);
		ids[2] = id_from_index(index + 2, first_residue);

		conformation.set_bond_angle(
				ids[0], ids[1], ids[2], bond_angles[index + 1]);
	}

	for (Size index = 1; index <= num_torsion_angles; index++) {
		ids[0] = id_from_index(index + 0, first_residue);
		ids[1] = id_from_index(index + 1, first_residue);
		ids[2] = id_from_index(index + 2, first_residue);
		ids[3] = id_from_index(index + 3, first_residue);

		conformation.set_torsion_angle(
					ids[0], ids[1], ids[2], ids[3], torsion_angles[index + 1] + 1);
	}

	// Dump the final structure.

	pose.dump_pdb("1ubq.out.pdb");
	*/
}

