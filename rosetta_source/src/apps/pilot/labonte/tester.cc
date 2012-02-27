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

/// @file   tester.cc
/// @brief  This is simply a generic pilot app for testing changes.
/// @author Labonte

// includes
#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/vector1.hh>
#include <core/import_pose/import_pose.hh>

int
main(int argc, char *argv[])
{
	using namespace core;
	using namespace import_pose;
	using namespace utility;
	using namespace pose;
	using namespace kinematics;

	// initialize core
	devel::init(argc, argv);

	// declare variables
	Pose test_pose;

	// import a test pose
	pose_from_pdb(test_pose,
			"/home/labonte/Code/Fusion_Project/host_insert_fusion.pdb");


	// test new set_chi_true_range() method
	MoveMap mm = MoveMap();
	mm.set_chi_true_range(5, 10);

	mm.show(20);


	// test modified split_by_chain() method
	Pose new_pose;

	vector1<pose::PoseOP> pose_collection;

	pose_collection = test_pose.split_by_chain();

	new_pose = *pose_collection[1];
	new_pose.dump_pdb("/home/labonte/Code/Fusion_Project/test.pdb", "");
}




