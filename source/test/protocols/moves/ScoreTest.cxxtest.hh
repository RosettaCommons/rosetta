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
#include <protocols/simple_moves/ScoreMover.hh>


//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/simple_moves/ScoreMover.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <vector>


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
		//std::cout << "Start All Scoring tests" << "\n";
		core_init_with_additional_options( "-score:patch score12 -out:output -restore_pre_talaris_2013_behavior" );
		one_mover_test(__FILE__, __LINE__, protocols::moves::MoverOP( new protocols::simple_moves::ScoreMover ),
						 "protocols/moves/test_in.pdb", "protocols/moves/score12.pdb",
						 0, "protocols/moves/score12.u", "protocols");
		//std::cout << "End Scoring -score:patch score12 test" << "\n";
		//std::cout << "End Scoring tests" << "\n";
		core_init_with_additional_options( "" );
	 }
};


