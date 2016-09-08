// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/id/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/import_pose/import_pose.hh>

using namespace std;
using namespace core;

int main(int argc, char** argv) {
	devel::init(argc, argv);

	pose::Pose pose;
	//import_pose::pose_from_file(pose, "structures/contrived/minimal.gly.pdb", core::import_pose::PDB_file);
	import_pose::pose_from_file(pose, "structures/contrived/one-residue.pdb", core::import_pose::PDB_file);
	//import_pose::pose_from_file(pose, "structures/contrived/linear.pdb", core::import_pose::PDB_file);
	pose.dump_pdb("not_messed_up.pdb");

	id::NamedAtomID n_name ("N", 1);
	id::NamedAtomID ca_name ("CA", 1);
	id::NamedAtomID c_name ("C", 1);

	id::AtomID n_id = pose::named_atom_id_to_atom_id(n_name, pose);
	id::AtomID ca_id = pose::named_atom_id_to_atom_id(ca_name, pose);
	id::AtomID c_id = pose::named_atom_id_to_atom_id(c_name, pose);

	cout << "AtomID" << endl;
	cout << n_id.atomno() << endl;
	cout << ca_id.atomno() << endl;
	cout << c_id.atomno() << endl;
	cout << endl;

	cout << "Bond Lengths" << endl;
	cout << pose.conformation().bond_length(n_id, ca_id) << endl;
	cout << pose.conformation().bond_length(ca_id, c_id) << endl;
	cout << endl;

	cout << "Bond Angles" << endl;
	cout << 180 * pose.conformation().bond_angle(n_id, ca_id, c_id) / 3.141592 << endl;
	cout << endl;

	pose.conformation().set_bond_length(n_id, ca_id, 4);
	pose.conformation().set_bond_length(ca_id, c_id, 2);
	pose.conformation().set_bond_angle(n_id, ca_id, c_id, 3.141592);

	pose.dump_pdb("messed_up.pdb");

	/*
	for (Size residue = 1; residue <= pose.size(); residue++) {
		string atoms[] = {"N", "CA", "C"};

		for (Size atom = 0; atom < 3; atom++) {

			id::DOF_ID torsion_angle (atom_id, id::PHI);
			id::DOF_ID bond_angle (atom_id, id::THETA);
			id::DOF_ID bond_length (atom_id, id::D);

			cout << "Torsion Angle (" << residue << "/" << atoms[atom] << "): ";
			if (pose.has_dof(torsion_angle)) cout << pose.dof(torsion_angle) << endl;
			else cout << "N/A" << endl;

			cout << "Bond Angle (" << residue << "/" << atoms[atom] << "): ";
			if (pose.has_dof(bond_angle)) cout << pose.dof(bond_angle) << endl;
			else cout << "N/A" << endl;

			cout << "Bond Length (" << residue << "/" << atoms[atom] << "): ";
			if (pose.has_dof(bond_length)) cout << pose.dof(bond_length) << endl;
			else cout << "N/A" << endl;

			cout << endl;
		}
	}
	*/
}
