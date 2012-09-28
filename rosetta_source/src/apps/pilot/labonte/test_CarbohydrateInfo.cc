// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test_CarbohydrateInfo.cc
/// @brief   Pilot application source code for testing CarbohydrateInfo.
/// @author  labonte


// Project headers
#include <devel/init.hh>
//#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
//#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

// C++ headers
#include <iostream>


int main(int argc, char *argv[])
{
	using namespace std;
	using namespace core;
	using namespace pose;
	using namespace import_pose;
	using namespace conformation;

	// Initialize core.
	devel::init(argc, argv);

	// Declare variables.
	Pose pose;

	// Import test carbohydrate pose.
	pose_from_pdb(pose, "/home/labonte/Workspace/Carbohydrates/heparin-6-mer.pdb");

	for (Size i = 1; i <= pose.total_residue(); ++i) {
		Residue res = pose.residue(i);
		cout << *(res.carbohydrate_info()) << endl << endl;
	}
}
