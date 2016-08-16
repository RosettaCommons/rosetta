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
#include <core/id/AtomID.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>
#include <numeric/numeric.functions.hh>


using namespace std;
using numeric::constants::r::pi;
using numeric::mod;
using numeric::modulo;

double Mod(double x, double y) { // {{{1
	// Floating-point modulo
	// The result (the remainder) has same sign as the divisor.
	// Similar to matlab's mod(); Not similar to fmod() -   Mod(-3,4)= 1   fmod(-3,4)= -3

	if (0. == y)
		return x;

	double m= x - y * floor(x/y);

	// handle boundary cases resulted from floating-point cut off:

	if (y > 0) {         // modulo range: [0..y)
		if (m >= y)        // Mod(-1e-16, 360.): m= 360.
			return 0;
		if (m < 0) {
			if (y+m == y)
				return 0;      // just in case...
			else
				return y+m;    // Mod(106.81415022205296 , _doubleWO_PI ): m= -1.421e-14
		}
	}
	else {               // modulo range: (y..0]
		if (m<=y)          // Mod(1e-16, -360.): m= -360.
			return 0;
		if (m>0 ) {
			if (y+m == y)
				return 0  ;    // just in case...
			else
				return y+m;    // Mod(-106.81415022205296, -_doubleWO_PI): m= 1.421e-14
		}
	}

	return m;
}

inline double WrapPosNegPI(double fAng) { // {{{1
	// wrap [rad] angle to [-PI..PI)
	return Mod(fAng + pi, 2 * pi) - pi;
}

inline double WrapTwoPI(double fAng) { // {{{1
	// wrap [rad] angle to [0..TWO_PI)
	return Mod(fAng, 2 * pi);
}

inline double WrapPosNeg180(double fAng) { // {{{1
	// wrap [deg] angle to [-180..180)
	return Mod(fAng + 180., 360.) - 180.;
}

inline double Wrap360(double fAng) { // {{{1
	// wrap [deg] angle to [0..360)
	return Mod(fAng ,360.);
}
// }}}1

using core::Size;
using core::Real;
using core::id::AtomID;
using core::pose::Pose;
using core::conformation::Conformation;
using core::PointPosition;


int main(int argc, char** argv) {

	double angles[] = {-420, -360, -240, -120, 0, 120, 240, 360, 420};

	for (int i = 0; i < 9; i++) cout << angles[i] << ", ";
	cout << endl;
	for (int i = 0; i < 9; i++) cout << mod<double>(angles[i], 360.) << ", ";
	cout << endl;
	for (int i = 0; i < 9; i++) cout << modulo<double>(angles[i], 360.) << ", ";
	cout << endl;
	for (int i = 0; i < 9; i++) cout << Wrap360(angles[i]) << ", ";
	cout << endl;

	cout << endl;

	for (int i = 0; i < 9; i++) cout << angles[i] << ", ";
	cout << endl;
	for (int i = 0; i < 9; i++) cout << mod<double>(angles[i] + 180, 360.) - 180 << ", ";
	cout << endl;
	for (int i = 0; i < 9; i++) cout << modulo<double>(angles[i] + 180, 360.) - 180 << ", ";
	cout << endl;
	for (int i = 0; i < 9; i++) cout << WrapPosNeg180(angles[i]) << ", ";
	cout << endl;
}


int not_main(int argc, char** argv) {

	devel::init(argc, argv);

	core::pose::Pose pose_a;
	core::pose::Pose pose_b;

	// 5 Residue
	//core::import_pose::pose_from_file(pose_a,
	//	"/home/kale/rosetta/main/apps/monte_carlo/jobs/tests/5-residue/trajectory/013.pdb");

	// 6 Residue
	//core::import_pose::pose_from_file(pose_a,
	//	"/home/kale/rosetta/main/apps/monte_carlo/jobs/tests/6-residue/trajectory/003.pdb");

	// 7 Residue
	core::import_pose::pose_from_file(pose_a,
		"/home/kale/rosetta/main/apps/monte_carlo/jobs/tests/7-residue/trajectory/005.pdb");

	Conformation &conformation_a = pose_a.conformation();

	///// End of Boilerplate ////////////////////////////////////////////////////

	AtomID ids[4];

	Size num_atoms = 3 * conformation_a.size();
	Size num_bond_lengths = num_atoms - 1;
	Size num_bond_angles = num_atoms - 2;
	Size num_torsion_angles = num_atoms - 3;

	for (Size index = 1; index <= num_atoms; index++) {
		ids[0] = AtomID(1 + (index - 1) % 3, 1 + (index - 1) / 3);
		PointPosition position = conformation_a.xyz(ids[0]);

		cout << "\t\tatom_xyzs[" << index << "] ";
		if (index < 10) cout << " "; cout << "= PointPosition(";
		cout << position.x() << ", ";
		cout << position.y() << ", ";
		cout << position.z() << ");" << endl;
	}
	cout << endl;

	/*
	for (Size index = 1; index <= num_bond_lengths; index++) {
		ids[0] = AtomID(1 + (index - 1) % 3, 1 + (index - 1) / 3);
		ids[1] = AtomID(1 + (index + 0) % 3, 1 + (index + 0) / 3);

		cout << "\t\tlanding_manager->bond_lengths[" << index << "] ";
		if (index < 10) cout << " "; cout << "= ";
		cout << conformation_b.bond_length(ids[0], ids[1]) << ";" << endl;
	}
	cout << endl;

	for (Size index = 1; index <= num_bond_angles; index++) {
		ids[0] = AtomID(1 + (index - 1) % 3, 1 + (index - 1) / 3);
		ids[1] = AtomID(1 + (index + 0) % 3, 1 + (index + 0) / 3);
		ids[2] = AtomID(1 + (index + 1) % 3, 1 + (index + 1) / 3);

		cout << "\t\tlanding_manager->bond_angles[" << index << "] ";
		if (index < 10) cout << " "; cout << "= ";
		cout << WrapTwoPI(conformation_b.bond_angle(ids[0], ids[1], ids[2])) << ";" << endl;
	}
	cout << endl;

	for (Size index = 1; index <= num_torsion_angles; index++) {
		ids[0] = AtomID(1 + (index - 1) % 3, 1 + (index - 1) / 3);
		ids[1] = AtomID(1 + (index + 0) % 3, 1 + (index + 0) / 3);
		ids[2] = AtomID(1 + (index + 1) % 3, 1 + (index + 1) / 3);
		ids[3] = AtomID(1 + (index + 2) % 3, 1 + (index + 2) / 3);

		cout << "\t\tlanding_manager->torsion_angles[" << index << "] ";
		if (index < 10) cout << " "; cout << "= ";
		cout << WrapTwoPI(conformation_b.torsion_angle(ids[0], ids[1], ids[2], ids[3])) << ";" << endl;
	}
	cout << endl;
	*/
}
