// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/doug/score_min_rtmin_test.cc
/// @brief Simple test of scoring and minimization
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/MinMover.hh>

// core headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("score_min_rtmin_test");

int
main( int argc, char * argv [] )
{
	try {

		devel::init( argc, argv );

		TR << "SCORE FUNCTION" << std::endl;

		// create score function
		//core::scoring::ScoreFunctionOP score_fxn( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::MM_STD_WTS ) );
		core::scoring::ScoreFunctionOP score_fxn( core::scoring::get_score_function() );

		// make a poly-leucine pose
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "LLLLLLLLLLLLLLLLLLLL", "fa_standard" );

		// set dihedrals
		for ( core::Size i(1); i <= 20; ++i ) {
			pose.set_phi( i, -60.0 );
			pose.set_psi( i, -40.0 );
			pose.set_omega( i, 180.0 );
			pose.set_chi( 1, i, 180.0 );
			pose.set_chi( 2, i, 180.0 );
		}

		// copy pose
		core::pose::Pose pose_min_nbt( pose ), pose_min_nbf( pose ), pose_rtmin( pose );

		// dump pose to file
		std::string before_filename( "pose_before.pdb" );
		std::string after_min_nbt_filename( "pose_min_nbt.pdb" );
		std::string after_min_nbf_filename( "pose_min_nbf.pdb" );
		std::string after_rtmin_filename( "pose_rtmin.pdb" );

		pose.dump_pdb( before_filename );

		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		movemap->set_chi( true );

		// minimize nbt
		protocols::simple_moves::MinMoverOP min_mover_nbt(new protocols::simple_moves::MinMover( movemap, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, true ) );
		TR << "Score before min NBT: " << (*score_fxn)( pose_min_nbt ) << std::endl;
		TR << "Minimize NBT" << std::endl;
		min_mover_nbt->apply( pose_min_nbt );
		TR << "Score after min NBT: " << (*score_fxn)( pose_min_nbt ) << std::endl;
		pose_min_nbt.dump_pdb( after_min_nbt_filename );

		// minimize nbf
		protocols::simple_moves::MinMoverOP min_mover_nbf(new protocols::simple_moves::MinMover( movemap, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, false ) );
		TR << "Score before min NBF: " << (*score_fxn)( pose_min_nbf ) << std::endl;
		TR << "Minimize NBF" << std::endl;
		min_mover_nbf->apply( pose_min_nbf );
		TR << "Score after min NBF: " << (*score_fxn)( pose_min_nbf ) << std::endl;
		pose_min_nbf.dump_pdb( after_min_nbf_filename );

		// rt min
		core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
		core::pack::task::operation::RestrictToRepackingOP rtr( new  core::pack::task::operation::RestrictToRepacking() );
		tf->push_back( rtr );

		protocols::simple_moves::RotamerTrialsMinMoverOP rtmm( new protocols::simple_moves::RotamerTrialsMinMover( score_fxn, tf ) );

		TR << "Score before RTMin: " << (*score_fxn)( pose_rtmin ) << std::endl;
		for ( core::Size i(1); i <= 100; ++i ) {
			TR << "RTMin num: " << i << std::endl;
			rtmm->apply( pose_rtmin );
			TR << "Score after RTMin num: " << i << " " << (*score_fxn)( pose_rtmin ) << std::endl;
		}

		pose_rtmin.dump_pdb( after_rtmin_filename );


		TR << "************************************d**o**n**e***********************************" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
