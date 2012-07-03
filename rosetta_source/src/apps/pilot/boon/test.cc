// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file   <filename_of_pilot_app>.cc
/// @brief  This is simply a generic pilot app for testing changes.
/// @author <your_name>

// includes
#include <iostream>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/moves/TrialMover.hh>

int main(int argc, char *argv[])
{
	using namespace core;
	using namespace import_pose;
	using namespace pose;

	// initialize core
	devel::init(argc, argv);

	// declare variables
	Pose test_pose;

	// import a test pose
	pose_from_pdb(test_pose, "/home/boon/data/test.pdb");

	std::cout << "Hello, Rosetta World!" << std::endl;
	std::cout << "I just imported my first pose into Rosetta." << std::endl;
	std::cout << "It has " << test_pose.total_residue() << " total residues." << std::endl;

	//std::cout << test_pose.energies() << std::endl;

	protocols::moves::TrialMover trialmover;
	std::cout << trialmover << std::endl;


}
