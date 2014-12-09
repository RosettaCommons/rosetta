// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/SequenceMapping.hh
/// @brief
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>
#include <test/UTracer.hh>

#include <core/sequence/util.hh>
// AUTO-REMOVED #include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/id/SequenceMapping.hh>

// AUTO-REMOVED #include <numeric/random/random.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/DerivedSequenceMapping.fwd.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
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
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <string>
#include <vector>



static basic::Tracer TR("test.core.sequence.SequenceMapping");


class SequenceMappingTests : public CxxTest::TestSuite
{
public:
	SequenceMappingTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

void test_transitive_mapping() {
	using namespace core::id;
	using namespace core::sequence;

	SequenceMapping map1, map2;
	map1.insert_aligned_residue_safe( 1, 1 );
	map1.insert_aligned_residue_safe( 2, 0 );
	map1.insert_aligned_residue_safe( 3, 2 );
	map1.insert_aligned_residue_safe( 4, 3 );
	map1.insert_aligned_residue_safe( 5, 5 );
	map1.insert_aligned_residue_safe( 6, 6 );

	map2.insert_aligned_residue_safe( 2, 3 );
	map2.insert_aligned_residue_safe( 3, 0 );
	map2.insert_aligned_residue_safe( 4, 4 );
	map2.insert_aligned_residue_safe( 5, 6 );
	map2.insert_aligned_residue_safe( 6, 7 );

	SequenceMapping trans_map = transitive_map( map1, map2 );
	TS_ASSERT( trans_map[1] == 0 );
	TS_ASSERT( trans_map[2] == 0 );
	TS_ASSERT( trans_map[3] == 3 );
	TS_ASSERT( trans_map[4] == 0 );
	TS_ASSERT( trans_map[5] == 6 );
	TS_ASSERT( trans_map[6] == 7 );
} // test_transitive_map

};
