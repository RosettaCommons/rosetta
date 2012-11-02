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


using namespace std;
using namespace core;
using namespace pose;
using namespace import_pose;
using namespace conformation;


void
test_sugar(Pose & sugar)
{
	cout << "Sequences:" << endl;
	cout << sugar.sequence() << endl;
	cout << sugar.chain_sequence(1) << endl;

	cout << "Residue Info:" << endl;
	for (Size i = 1; i <= sugar.total_residue(); ++i) {
		Residue res = sugar.residue(i);
		cout << res << endl << endl;
	}
}


int
main(int argc, char *argv[])
{
	// Initialize core.
	devel::init(argc, argv);

	// Declare variables.
	Pose maltotriose, isomaltose;

	cout << "------------------------------------------------------------" << endl;
	cout << "Importing maltotriose:" << endl;

	pose_from_pdb(maltotriose, "/home/labonte/Workspace/Carbohydrates/maltotriose.pdb");

	test_sugar(maltotriose);

	cout << "------------------------------------------------------------" << endl;
	cout << "Importing isomaltose:" << endl;

	pose_from_pdb(isomaltose, "/home/labonte/Workspace/Carbohydrates/isomaltose.pdb");

	test_sugar(isomaltose);
}
