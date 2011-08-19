// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
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
#include <basic/Tracer.hh>
#include <test/UTracer.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/id/SequenceMapping.hh>

#include <numeric/random/random.hh>

//Auto Headers
#include <utility/stream_util.hh>



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
