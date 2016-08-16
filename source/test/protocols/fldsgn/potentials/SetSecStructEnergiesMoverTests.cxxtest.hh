// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/fldsgn/potentials/SetSecStructEnergiesMoverTests.cxxtest.hh
/// @brief  Test suite for SetSecStructEnergies mover
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/fldsgn/potentials/SetSecStructEnergies.hh>

// Core Headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("SetSecStructEnergiesMoverTests");


class SetSecStructEnergiesMoverTests : public CxxTest::TestSuite {
	//Define Variables
	typedef protocols::fldsgn::potentials::SetSecStructEnergies SetSecStructEnergies;
	typedef core::scoring::ScoreFunctionFactory ScoreFunctionFactory;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;

public:

	void setUp()
	{
		core_init();
	}

	void tearDown()
	{
	}

	void
	compare_fldsgn_cen( ScoreFunctionCOP orig, ScoreFunctionCOP sfxn ) const
	{
		TS_ASSERT_EQUALS( (*orig)[ core::scoring::rg ], (*sfxn)[ core::scoring::rg ] );
		TS_ASSERT_EQUALS( (*orig)[ core::scoring::vdw ], (*sfxn)[ core::scoring::vdw ] );
		TS_ASSERT_DELTA( (*orig)[ core::scoring::natbias_ss ], 0.0, 1e-6 );
		TS_ASSERT_DELTA( (*sfxn)[ core::scoring::natbias_ss ], 1.0, 1e-6 );
		TS_ASSERT_DELTA( (*orig)[ core::scoring::natbias_hh ], 0.0, 1e-6 );
		TS_ASSERT_DELTA( (*sfxn)[ core::scoring::natbias_hh ], 1.0, 1e-6 );
		TS_ASSERT_DELTA( (*orig)[ core::scoring::natbias_hs ], 0.0, 1e-6 );
		TS_ASSERT_DELTA( (*sfxn)[ core::scoring::natbias_hs ], 1.0, 1e-6 );
		TS_ASSERT_DELTA( (*orig)[ core::scoring::ss_pair ], 1.0, 1e-6 );
		TS_ASSERT_DELTA( (*sfxn)[ core::scoring::ss_pair ], 0.0, 1e-6 );
		TS_ASSERT_DELTA( (*orig)[ core::scoring::hs_pair ], 1.0, 1e-6 );
		TS_ASSERT_DELTA( (*sfxn)[ core::scoring::hs_pair ], 0.0, 1e-6 );
		TS_ASSERT_DELTA( (*orig)[ core::scoring::rsigma ], 1.0, 1e-6 );
		TS_ASSERT_DELTA( (*sfxn)[ core::scoring::rsigma ], 0.0, 1e-6 );
		TS_ASSERT_DELTA( (*orig)[ core::scoring::natbias_stwist ], 0.0, 1e-6 );
		TS_ASSERT_DELTA( (*sfxn)[ core::scoring::natbias_stwist ], 1.0, 1e-6 );
	}

	void test_modify_score_function()
	{
		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is( pose, "protocols/fldsgn/filters/test_sheet.pdb" );

		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "fldsgn_cen" );
		ScoreFunctionCOP orig = sfxn->clone();

		TR << *sfxn << std::endl;

		SetSecStructEnergies set_ssene;
		set_ssene.set_scorefunction_ptr( sfxn );
		set_ssene.apply( pose );
		compare_fldsgn_cen( orig, sfxn );

		// Should be able to modify twice in a row
		set_ssene.apply( pose );
		compare_fldsgn_cen( orig, sfxn );

		// Should be able to modify new scorefunction
		ScoreFunctionOP sfxn2 = ScoreFunctionFactory::create_score_function( "fldsgn_cen" );
		set_ssene.set_scorefunction_ptr( sfxn2 );
		set_ssene.apply( pose );
		compare_fldsgn_cen( orig, sfxn2 );
	}

};



