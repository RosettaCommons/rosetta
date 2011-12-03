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
#include <protocols/simple_moves/BackboneMover.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Range.fwd.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/TorsionID_Range.fwd.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/id/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
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
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <vector>


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
		TEST_MOVER( protocols::simple_moves::SmallMover, "protocols/moves/test_in.pdb", "protocols/moves/smallmoves_out.pdb");
		//one_mover_test(__FILE__, __LINE__, new SmallMover,
		//			   "protocols/moves/test_in.pdb", "protocols/moves/smallmoves_out.pdb",
		//			   "protocols/moves/smallmoves.u", "core protocols");
		std::cout << "End SmallMover test" << "\n";

		TEST_MOVER( protocols::simple_moves::ShearMover, "protocols/moves/test_in.pdb", "protocols/moves/shearmoves_out.pdb");
		std::cout << "End ShearMover test" << "\n";
	 }
};
