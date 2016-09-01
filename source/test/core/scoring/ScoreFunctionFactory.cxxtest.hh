// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreFunctionFactory.cxxtest.hh
/// @brief  unit testing for core/scoring/ScoreFunctionFactory.cc (just one tiny bit for now)
/// @author Rocco Moretti (rmoretti@u.washington.edu) (ScoreFunctionUtility)
/// @author Steven Lewis (smlewi@gmail.com) (check for beta_nov15)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>

// Package headers
#include <basic/database/open.hh>

#include <basic/options/option.hh>

//Auto Headers

//#include <basic/Tracer.hh>

//using basic::T;
//using basic::Error;
//using basic::Warning;

//static basic::Tracer TR("core.scoring.ScoreFunction.cxxtest");

// using declarations
using namespace core;
using namespace scoring;

///////////////////////////////////////////////////////////////////////////
/// @name ScoreFunctionUtilityTest
/// @brief: unified tests for associated score function utilities
/// Moved from ScoreFunction.cxxtest.hh.
///////////////////////////////////////////////////////////////////////////
class ScoreFunctionUtilityTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_find_weights_file() {
		// Local weights
		TS_ASSERT_EQUALS( find_weights_file("core/scoring/test",".wts"), "core/scoring/test.wts" );
		TS_ASSERT_EQUALS( find_weights_file("core/scoring/test.wts",".wts"), "core/scoring/test.wts" );
		// Local patch
		TS_ASSERT_EQUALS( find_weights_file("core/scoring/test",".wts_patch"), "core/scoring/test.wts_patch" );
		TS_ASSERT_EQUALS( find_weights_file("core/scoring/test.wts_patch",".wts_patch"), "core/scoring/test.wts_patch" );
		// Database weights
		TS_ASSERT_EQUALS( find_weights_file("pre_talaris_2013_standard",".wts"), basic::database::full_name( "scoring/weights/pre_talaris_2013_standard.wts" ) );
		TS_ASSERT_EQUALS( find_weights_file("pre_talaris_2013_standard.wts",".wts"), basic::database::full_name( "scoring/weights/pre_talaris_2013_standard.wts" ) );
		TS_ASSERT_EQUALS( find_weights_file("talaris2013.wts",".wts"), basic::database::full_name( "scoring/weights/talaris2013.wts" ) );
		TS_ASSERT_EQUALS( find_weights_file("talaris2013",".wts"), basic::database::full_name( "scoring/weights/talaris2013.wts" ) );
		TS_ASSERT_EQUALS( find_weights_file("talaris2014.wts",".wts"), basic::database::full_name( "scoring/weights/talaris2014.wts" ) );
		TS_ASSERT_EQUALS( find_weights_file("talaris2014",".wts"), basic::database::full_name( "scoring/weights/talaris2014.wts" ) );
		TS_ASSERT_EQUALS( find_weights_file("beta_nov15.wts",".wts"), basic::database::full_name( "scoring/weights/beta_nov15.wts" ) );
		TS_ASSERT_EQUALS( find_weights_file("beta_nov15",".wts"), basic::database::full_name( "scoring/weights/beta_nov15.wts" ) );
		TS_ASSERT_EQUALS( find_weights_file("beta_july15.wts",".wts"), basic::database::full_name( "scoring/weights/beta_july15.wts" ) );
		TS_ASSERT_EQUALS( find_weights_file("beta_july15",".wts"), basic::database::full_name( "scoring/weights/beta_july15.wts" ) );
		// Database patch
		TS_ASSERT_EQUALS( find_weights_file("score12",".wts_patch"), basic::database::full_name( "scoring/weights/score12.wts_patch" ) );
		TS_ASSERT_EQUALS( find_weights_file("score12.wts_patch",".wts_patch"), basic::database::full_name( "scoring/weights/score12.wts_patch" ) );
	}


}; //ScoreFunctionUtilityTest


///////////////////////////////////////////////////////////////////////////
/// @name ScoreFunctionFactoryTest
/// @brief: unified tests for ScoreFunctionFactory
///////////////////////////////////////////////////////////////////////////
class ScoreFunctionFactoryTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	/// @details the beta_15 family of scorefunctions requires both a weights file AND a command-line flag to activate loading of different SF parameters (like solvation parameters).  This tests that if the weights file is beta, the command line flags should have been set too.
	void test_beta_15_behavior() {

		using namespace core::scoring;

		core_init_with_additional_options(""); //init with NO beta params
		//These should utility_exit_with_message; I'm (SML) not sure how to capture or examine that it actually threw the right thing.
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING(ScoreFunctionFactory::validate_beta("beta_nov15.wts", basic::options::option));
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING(ScoreFunctionFactory::validate_beta("beta_july15.wts", basic::options::option));

		core_init_with_additional_options("-beta_nov15"); //init with beta_nov15 params
		//shouldn't throw
		TS_ASSERT_THROWS_NOTHING(ScoreFunctionFactory::validate_beta("beta_nov15.wts", basic::options::option));
		TS_ASSERT(ScoreFunctionFactory::validate_beta("beta_nov15.wts", basic::options::option));
		//This should utility_exit_with_message; I'm (SML) not sure how to capture or examine that it actually threw the right thing.
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING(ScoreFunctionFactory::validate_beta("beta_july15.wts", basic::options::option));

		core_init_with_additional_options("-beta_july15"); //init with NO beta params
		//This should utility_exit_with_message; I'm (SML) not sure how to capture or examine that it actually threw the right thing.
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING(ScoreFunctionFactory::validate_beta("beta_nov15.wts", basic::options::option));
		//shouldn't throw
		TS_ASSERT_THROWS_NOTHING(ScoreFunctionFactory::validate_beta("beta_july15.wts", basic::options::option));
		TS_ASSERT(ScoreFunctionFactory::validate_beta("beta_july15.wts", basic::options::option));
		return;
	}

	/// @details same test, but now without the .wts extension, since the ScoreFunction(Factory) machinery is agnostic
	void test_beta_15_behavior_noext() {

		using namespace core::scoring;

		core_init_with_additional_options(""); //init with NO beta params
		//These should utility_exit_with_message; I'm (SML) not sure how to capture or examine that it actually threw the right thing.
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING(ScoreFunctionFactory::validate_beta("beta_nov15", basic::options::option));
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING(ScoreFunctionFactory::validate_beta("beta_july15", basic::options::option));

		core_init_with_additional_options("-beta_nov15"); //init with beta_nov15 params
		//shouldn't throw
		TS_ASSERT_THROWS_NOTHING(ScoreFunctionFactory::validate_beta("beta_nov15", basic::options::option));
		TS_ASSERT(ScoreFunctionFactory::validate_beta("beta_nov15", basic::options::option));
		//This should utility_exit_with_message; I'm (SML) not sure how to capture or examine that it actually threw the right thing.
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING(ScoreFunctionFactory::validate_beta("beta_july15", basic::options::option));

		core_init_with_additional_options("-beta_july15"); //init with NO beta params
		//This should utility_exit_with_message; I'm (SML) not sure how to capture or examine that it actually threw the right thing.
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING(ScoreFunctionFactory::validate_beta("beta_nov15", basic::options::option));
		//shouldn't throw
		TS_ASSERT_THROWS_NOTHING(ScoreFunctionFactory::validate_beta("beta_july15", basic::options::option));
		TS_ASSERT(ScoreFunctionFactory::validate_beta("beta_july15", basic::options::option));
		return;
	}

	/// @brief this function checks the return values for validate_beta
	void test_beta_15_returns() {

		using namespace core::scoring;

		core_init_with_additional_options(""); //init with NO beta params

		TS_ASSERT(ScoreFunctionFactory::validate_beta("some_weights_file", basic::options::option));

		//test a short string; earlier code did str.length-4, which was bad
		//should just bail and return early with true
		//lit = little, because "short" abbreviates badly
		TS_ASSERT(ScoreFunctionFactory::validate_beta("lit", basic::options::option));

		//beta should return true but make noise doing it if mute isn't on
		//(not that the latter bit there is tested)
		TS_ASSERT(ScoreFunctionFactory::validate_beta("beta", basic::options::option));

		return;
	}

}; //ScoreFunctionFactoryTest
