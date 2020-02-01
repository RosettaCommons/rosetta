// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/relax/AcceptToBestMover.cxxtest.hh
/// @brief  test for AcceptToBestMover
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pdb1rpb.hh>

// Unit header
#include <protocols/relax/AcceptToBestMover.hh>

// project headers
#include <core/types.hh>

#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <utility/vector1.hh>

class AcceptToBestMoverTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	bool
	scores_are_close( core::Real s1, core::Real s2 ){
		return ( s1 - s2 < 0.01 ) && ( s2 - s1 < 0.01 );
	}

	void test_accept_to_best() {
#ifdef SERIALIZATION
		using namespace protocols::relax;

		//First, unit test the unit test
		TS_ASSERT( scores_are_close( 2.001, 2.002 ) );
		TS_ASSERT( scores_are_close( 2.003, 2.002 ) );
		TS_ASSERT( ! scores_are_close( 2.103, 2.002 ) );
		TS_ASSERT( ! scores_are_close( 2.003, 2.102 ) );

		protocols::relax::AcceptToBestMover atb_mover;
		auto sfxn = core::scoring::get_score_function();
		atb_mover.set_sfxn( sfxn );

		protocols::simple_moves::MutateResidue mutate( 10, "ARG" );
		protocols::minimization_packing::MinMover minmover;
		minmover.score_function( sfxn );

		//Do this 3 times to make sure the mover properly resets
		for ( core::Size ii = 1; ii <= 3; ++ii ) {
			core::pose::Pose pose(pdb1rpb_pose());

			core::Real const first_score = (*sfxn)(pose);
			atb_mover.apply( pose );
			//Make sure the atb_mover reconizes this as a new pose and does not overwrite with one from the previous loop:
			TS_ASSERT( scores_are_close( (*sfxn)(pose), first_score ) );

			//OPERATION 1
			minmover.apply( pose );
			core::Real const score_after_min = (*sfxn)(pose);
			atb_mover.apply( pose );
			core::Real const second_score = (*sfxn)(pose);

			if ( score_after_min < first_score ) {
				TS_ASSERT( scores_are_close( second_score, score_after_min ) );
			} else {
				//check for reversion to first score
				TS_ASSERT( scores_are_close( second_score, first_score ) );
			}

			//OPERATION 2
			for ( core::Size jj = 1; jj <= 10; ++jj ) {
				mutate.apply( pose );
				core::Real const score_after_mutate = (*sfxn)(pose);
				atb_mover.apply( pose );
				core::Real const third_score = (*sfxn)(pose);

				if ( score_after_mutate < second_score ) {
					TS_ASSERT( scores_are_close( third_score, score_after_mutate ) );
				} else {
					TS_ASSERT( scores_are_close( third_score, second_score ) );
				}
			}
		}


#else
		TS_ASSERT( true );
#endif
	}

	void test_accept_to_best_no_minmover() {
#ifdef SERIALIZATION
		using namespace protocols::relax;

		//First, unit test the unit test
		TS_ASSERT( scores_are_close( 2.001, 2.002 ) );
		TS_ASSERT( scores_are_close( 2.003, 2.002 ) );
		TS_ASSERT( ! scores_are_close( 2.103, 2.002 ) );
		TS_ASSERT( ! scores_are_close( 2.003, 2.102 ) );

		protocols::relax::AcceptToBestMover atb_mover;
		auto sfxn = core::scoring::get_score_function();
		atb_mover.set_sfxn( sfxn );

		protocols::simple_moves::MutateResidue mutate( 10, "ARG" );
		protocols::minimization_packing::MinMover minmover;
		minmover.score_function( sfxn );

		//Do this 3 times to make sure the mover properly resets
		for ( core::Size ii = 1; ii <= 3; ++ii ) {
			core::pose::Pose pose(pdb1rpb_pose());

			core::Real const first_score = (*sfxn)(pose);
			atb_mover.apply( pose );
			//Make sure the atb_mover reconizes this as a new pose and does not overwrite with one from the previous loop:
			TS_ASSERT( scores_are_close( (*sfxn)(pose), first_score ) );

			core::Real current_score = first_score;
			for ( core::Size jj = 1; jj <= 10; ++jj ) {
				mutate.apply( pose );
				core::Real const score_after_mutate = (*sfxn)(pose);
				atb_mover.apply( pose );
				core::Real const second_score = (*sfxn)(pose);

				if ( score_after_mutate < current_score ) {
					TS_ASSERT( scores_are_close( second_score, score_after_mutate ) );
				} else {
					TS_ASSERT( scores_are_close( second_score, current_score ) );
				}

				current_score = second_score;
			}
		}


#else
		TS_ASSERT( true );
#endif
	}
};
