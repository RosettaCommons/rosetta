// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/ScoreTest.cxxtest.hh
/// @brief  tests for container ScoreMover classe.
/// @author Monica Berrondo

// Test headers
#include <test/UMoverTest.hh>

// Unit headers
#include <protocols/moves/ScoreMover.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace protocols::moves;

///////////////////////////////////////////////////////////////////////////
/// @name ScoreTest
/// @brief: class for Score Mover with different scores unified testing
/// @author Monica Berrondo
///////////////////////////////////////////////////////////////////////////
class ScoreTest : public CxxTest::TestSuite, public test::UMoverTest {

public:
	void setUp() {
		test::UMoverTest::setUp();
	}

	/// @brief test score 12
	void test_Score12() {
		std::cout << "Start All Scoring tests" << "\n";
		core_init_with_additional_options( "-score:patch score12 -out:output" );
		one_mover_test(__FILE__, __LINE__, new ScoreMover,
						 "protocols/moves/test_in.pdb", "protocols/moves/score12.pdb",
						 0, "protocols/moves/score12.u", "protocols");
		std::cout << "End Scoring -score:patch score12 test" << "\n";
		std::cout << "End Scoring tests" << "\n";
		core_init_with_additional_options( "" );
	 }
};


