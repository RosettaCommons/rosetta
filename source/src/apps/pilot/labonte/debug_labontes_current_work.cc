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

#include <devel/init.hh>
#include <basic/Tracer.hh>

#include <core/types.hh>
//#include <core/chemical/RingConformerSet.hh>
//#include <core/chemical/carbohydrates/database_io.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/TorsionID.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>

//#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
//#include <utility/vector0.hh>

//#include <protocols/simple_moves/BackboneMover.hh>
//#include <protocols/simple_moves/PackRotamersMover.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Construct random-number generator.
static numeric::random::RandomGenerator RG( 21 );  // the 6th triangular number


int main(int argc, char *argv[])
{
    try {
		using namespace std;
		using namespace core;
		using namespace protocols::loops;
		using namespace protocols::loops::loop_closure::ccd;
		using namespace chemical;
		using namespace import_pose;
		using namespace pose;
		using namespace utility;
		using namespace kinematics;
		//using namespace protocols::simple_moves;
		//using namespace core::scoring;
		//using namespace core::pack::task;
		using namespace core::conformation;

		// initialize core
		devel::init(argc, argv);

		// declare variables
		Pose pose;

		// Make a test pose.

		// pose_from_pdb( pose, "/Volumes/MacintoshSSD/work/weitzner/loop_features_test/heavyAlignClean/1g7lB.clean.pdb" );
		// Loop const loop( pose.pdb_info()->pdb2pose( 'H', 107 ), pose.pdb_info()->pdb2pose( 'H', 138 ),
		//	pose.pdb_info()->pdb2pose( 'H', 111 ) );

		make_pose_from_sequence( pose, "AAAAAAAAAA", "fa_standard" );
		Loop const loop( 2, 9 ,5 );

		MoveMapOP mm = new MoveMap;
		mm->set_bb( true );
		for ( core::uint i = 1; i < pose.total_residue(); ++i ) {
			if (pose.total_residue() == 10 ) { pose.set_omega( i, 180.0 ); } // We only need to do this for the pose from seq
			mm->set( id::TorsionID( i, id::BB, 3 ), false );
		}

		set_single_loop_fold_tree( pose, loop );
		add_single_cutpoint_variant( pose, loop );

		cout << pose << endl << endl;

		for ( core::uint i = loop.start(); i <= loop.stop(); ++i ) {
			pose.set_phi( i, RG.uniform() * 360 );
			pose.set_psi( i, RG.uniform() * 360 );
		}

		pose.dump_pdb( "/Users/weitzner/ccd_test_start.pdb" );

		CCDLoopClosureMover mover( loop, mm );
		mover.check_rama_scores( false );
		mover.max_per_move_torsion_delta_per_residue( 180., 180., 180. );
		mover.max_total_torsion_delta_per_residue( 180., 180., 180. );
		mover.max_cycles( 500 );
		mover.apply( pose );

		pose.dump_pdb( "/Users/weitzner/ccd_test_end.pdb" );
   } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}
