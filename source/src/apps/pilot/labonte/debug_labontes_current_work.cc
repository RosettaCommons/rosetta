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

/// @file   debug_labontes_current_work.cc
/// @brief  This is simply a generic pilot app for testing changes.
/// @author Labonte

// includes
#include <iostream>
//#include <algorithm>

#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
//#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/excn/Exceptions.hh>
//#include <utility/vector1.hh>
//#include <utility/vector0.hh>

#include <protocols/simple_moves/BackboneMover.hh>

int main(int argc, char *argv[])
{
    try {
		using namespace std;
		using namespace core;
		using namespace import_pose;
		using namespace pose;
		using namespace utility;
		using namespace kinematics;
		using namespace protocols;

		// initialize core
		devel::init(argc, argv);

		// declare variables
		Pose pose;

		// import a test pose
		//pose_from_pdb(pose, "/home/labonte/Workspace/test_input/test.pdb");

		pose_from_pdb(pose, "/home/labonte/Workspace/Carbohydrates/test.pdb");

		cout << pose << endl << endl;

		Size n_res = pose.total_residue();
		for (core::uint i = 1; i <= n_res; ++i) {
			cout << pose.residue(i) << endl << endl;
		}

		pose.dump_pdb("/home/labonte/Workspace/test_output/modified_sugar.pdb", "");
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}
