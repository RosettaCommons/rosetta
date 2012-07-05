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
#include <string>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
/*
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
*/
#include <protocols/simple_moves/ReturnSidechainMover.hh>

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
/*	protocols::simple_moves::SmallMoverOP smallmover = new protocols::simple_moves::SmallMover;
	protocols::simple_moves::ShearMoverOP shearmover = new protocols::simple_moves::ShearMover;
	protocols::moves::SequenceMover mover;
	core::Real num = 1.0;
	mover.add_mover(smallmover, num);
	mover.add_mover(shearmover, num);
	std::cout << mover << std::endl;
	core::Size index = 0;
	std::cout << mover.get_mover(index) << std::endl;
	std::cout << "try" << mover.get_mover(index+1) << std::endl;*/

	/*protocols::simple_moves::SwitchResidueTypeSetMover switchmover;
	std::cout << switchmover << std::endl;
	protocols::simple_moves::SwitchResidueTypeSetMover switch2 = protocols::simple_moves::SwitchResidueTypeSetMover("centroid");
	std::cout << switch2 << std::endl;*/

	protocols::simple_moves::ReturnSidechainMover returnmover = protocols::simple_moves::ReturnSidechainMover(test_pose) ;
	std::cout << returnmover << std::endl;
	protocols::simple_moves::ReturnSidechainMover returnmover2 = protocols::simple_moves::ReturnSidechainMover(test_pose,10,20) ;
	std::cout << returnmover2 << std::endl;

}
