// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/MoversTest.cxxtest.hh
/// @brief  tests for container Movers classes.
/// @author Sergey Lyskov

// Test headers
#include <test/UMoverTest.hh>

// Unit headers
#include <protocols/moves/BackboneMover.hh>

//Auto Headers
#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace protocols::moves;

///////////////////////////////////////////////////////////////////////////
/// @name MoversTest
/// @brief: class for Movers unified testing
/// @author Sergey Lyskov
///////////////////////////////////////////////////////////////////////////
class MoversTest : public CxxTest::TestSuite, public test::UMoverTest {

public:
	void setUp() {
		test::UMoverTest::setUp();
	}

	void test_AllMovers() {
		TEST_MOVER(SmallMover, "protocols/moves/test_in.pdb", "protocols/moves/smallmoves_out.pdb");
		//one_mover_test(__FILE__, __LINE__, new SmallMover,
		//			   "protocols/moves/test_in.pdb", "protocols/moves/smallmoves_out.pdb",
		//			   "protocols/moves/smallmoves.u", "core protocols");
		std::cout << "End SmallMover test" << "\n";

		TEST_MOVER(ShearMover, "protocols/moves/test_in.pdb", "protocols/moves/shearmoves_out.pdb");
		std::cout << "End ShearMover test" << "\n";
	 }
};
