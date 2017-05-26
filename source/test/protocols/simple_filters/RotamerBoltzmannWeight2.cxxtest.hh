// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/RotamerBoltzmannWeight2.cxxtest.hh
/// @brief  test for RotamerBoltzmannWeight2 filter
/// @author Tom Linsky ( tlinsky at uw dot edu )

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>

// Protocol Headers
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_filters/RotamerBoltzmannWeight.hh>
#include <protocols/simple_filters/RotamerBoltzmannWeight2.hh>
#include <protocols/toolbox/EnergyLandscapeEvaluator.hh>

// Core Headers
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <basic/Tracer.hh>

// Boost Headers
#include <boost/assign.hpp>

// C++ Headers


static basic::Tracer TR("protocols.simple_filters.RotamerBoltzmannWeight2.cxxtest.hh");

// --------------- Test Class --------------- //

class RotamerBoltzmannWeight2Tests : public CxxTest::TestSuite {

public:

	void setUp() {
		protocols_init();
	}

	void tearDown() {
	}

	void test_calculator_ids() {
		using namespace protocols::simple_filters;
		RotamerBoltzmannWeight2 rbw;
		RotamerBoltzmannWeight2 rbw2;
		TS_ASSERT_DIFFERS( rbw.calculator_id(), rbw2.calculator_id() );
		RotamerBoltzmannWeight2 rbw3 = rbw2;
		TS_ASSERT_DIFFERS( rbw3.calculator_id(), rbw2.calculator_id() );
	}

	// Note: If the default score function, or FastRelax behavior should change, the
	// hard coded constants (expected_rotamer_scores, rbw2_pnear) are written out to
	// tracers when you run this test. In the case of the expected rotamer scores,
	// look for the "RotamerBoltzCalculator: [Initial input pose/This rotamer] has score XXX"
	// lines; for the rbw2_pnear value, it is output right before the TS_ASSERT_EQUALS
	// statement.
	void test_rbw2_vs_rbw() {
		using namespace protocols::simple_filters;
		using namespace protocols::toolbox;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		// relax input once cycle
		protocols::relax::FastRelax relax( core::scoring::get_score_function(), 2 );
		relax.apply( trpcage );

		trpcage.dump_pdb( "trpcage.pdb" );

		RotamerBoltzmannWeight rbw;
		rbw.rb_jump( 1 );
		rbw.unbound( 0 );
		rbw.skip_ala_scan( false );
		rbw.target_residues( "7" );
		rbw.scorefxn( core::scoring::get_score_function() );

		core::Real const rbw1_score = rbw.compute( trpcage );

		utility::vector1< core::Real > const expected_rotamer_scores = boost::assign::list_of
			(-34.1334) (-38.9484) (-39.1986);

		// initialize scorermspoints object with initial pose score
		ScoreRmsPoints score_rms( ScoreRmsPoint( -39.1986, 0.0 ) );

		// add test scores
		for ( utility::vector1< core::Real >::const_iterator rs=expected_rotamer_scores.begin(); rs!=expected_rotamer_scores.end(); ++rs ) {
			score_rms.push_back( ScoreRmsPoint( *rs, 0.0 ) );
		}

		RotamerBoltzmannWeightEvaluator rbw_evaluator( 0.8, false );
		core::Real const evaluator_score = rbw_evaluator.compute( score_rms );
		TS_ASSERT_DELTA( evaluator_score, -rbw1_score, 1e-3 );
		TR << "RBW1 Score: " << rbw1_score << std::endl;
		TR << "Evaluator Score: " << evaluator_score << std::endl;

		// hpark, May17 2017: rbw2_score is not exact match to rbw1_score and difference may depend on rotamer library...
		RotamerBoltzmannWeight2 rbw2;
		rbw2.set_scorefxn( core::scoring::get_score_function() );
		rbw2.set_residue_selector( core::select::residue_selector::ResidueSelectorCOP( new core::select::residue_selector::ResidueIndexSelector( "7" ) ) );

		core::Real const rbw2_score = rbw2.report_sm( trpcage );
		TR << "RBW2 Score: " << rbw2_score << std::endl;
		rbw2.report( TR, trpcage );
		TS_ASSERT_DELTA( rbw2_score, rbw1_score, 2.0e-2 );

		// switch energy landscape evaluator to test PNear
		EnergyLandscapeEvaluatorOP pnear( new MulliganPNearEvaluator( 0.8, 0.5 ) );
		rbw2.set_energy_landscape_evaluator( pnear );
		core::Real const rbw2_pnear = rbw2.report_sm( trpcage );
		TR << "RBW2 PNear: " << rbw2_pnear << std::endl;
		TS_ASSERT_DELTA( rbw2_pnear, -0.9328, 1e-3 );
	}

};
