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
/// @brief  This is simply a generic pilot app for testing things.
/// @author Labonte



// Package headers
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
//#include <core/pose/annotated_sequence.hh>
//#include <core/pose/PDBInfo.hh>
//#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
//#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
//#include <core/id/TorsionID.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/task/TaskFactory.hh>

#include <protocols/simple_moves/MinMover.hh>
//#include <protocols/simple_moves/BackboneMover.hh>
//#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/util.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <iostream>
//#include <algorithm>

// Construct random-number generator.


int main(int argc, char *argv[])
{
    try {
		using namespace std;
		using namespace utility;
		using namespace core;
		using namespace kinematics;
		using namespace scoring;
		using namespace import_pose;
		using namespace pose;
		using namespace protocols;
		using namespace simple_moves;

		// initialize core
		devel::init(argc, argv);

		// declare variables
		Pose pose, ref;

		// Make a test pose.
		//make_pose_from_sequence( pose, "AAAAAAAAAA", "fa_standard" );
		pose_from_pdb( pose, "/home/labonte/Workspace/Carbohydrates/MBP-G4_ref.pdb" );

		vector1< int > movable_jumps( 1, 1 );
		docking::setup_foldtree( pose, "A_B", movable_jumps );

		cout << pose << endl << endl;

		ScoreFunctionOP sf( get_score_function() );

		cout << "Initial Score: " << (*sf)( pose ) << endl;

		MoveMapOP mm( MoveMapOP( new MoveMap() ) );
		mm->set_jump( 1, true );

		MinMover jump_minimizer;
		jump_minimizer.movemap( mm );
		jump_minimizer.score_function( sf );

		rigid::RigidBodyPerturbMover perturber( 1, 2.0, 0.5 );

		jump_minimizer.apply( pose );

		cout << "Minimization Score Before Any Moves: " << (*sf)( pose ) << endl;

		for ( core::uint i = 1; i <= 100; ++i ) {
			perturber.apply( pose );
			jump_minimizer.apply( pose );

			cout << "Minimization Score After Rigid Move " << i << ": " << (*sf)( pose ) << endl;
		}

		//pose.dump_pdb( "" );
   } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}
