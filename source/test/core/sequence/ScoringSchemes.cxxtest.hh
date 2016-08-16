// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/sequence/ScoringSchemes.hh
/// @brief
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/UTracer.hh>
#include <core/types.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/PairScoringScheme.hh>
#include <core/sequence/PairScoringScheme.fwd.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/io/izstream.fwd.hh>
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


static basic::Tracer TR("test.core.sequence.ScoringSchemes");

class ScoringSchemeTests : public CxxTest::TestSuite
{
public:
	ScoringSchemeTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	void test_pair_scoring_scheme() {
		using core::Size;
		using core::Real;
		using namespace core::sequence;
		SequenceOP seq1( new Sequence( "ABCDEFGHIJ", "first",  1 ) );
		SequenceOP seq2( new Sequence( "----EFG-IJ", "second", 1 ) );

		core::Real const big_score( 537.5 );

		PairScoringSchemeOP ss( new PairScoringScheme );
		ss->add_scored_pair( 5, 5, big_score );
		TS_ASSERT( ss->score( seq1, seq2, 5, 5 ) == big_score );
	} // test_pair_scoring_scheme

	void test_simple_scoring_scheme() {
		using core::Size;
		using core::Real;
		using namespace core::sequence;

		SequenceOP seq1( new Sequence( "EFGILK", "first",  1 ) );
		SequenceOP seq2( new Sequence( "EFGIJX", "second", 1 ) );

		ScoringSchemeOP ss( new SimpleScoringScheme( 4, 1, -4, -1 ) );
		utility::vector1< core::Real > scores;
		scores.push_back(  4 );
		scores.push_back(  4 );
		scores.push_back(  4 );
		scores.push_back(  4 );
		scores.push_back(  1 );
		scores.push_back(  1 );

		for ( Size ii = 1; ii <= seq1->length() && ii <= seq2->length(); ++ii ) {
			TS_ASSERT( scores[ii] == ss->score( seq1, seq2, ii, ii ) );
		}
	} // test_pair_scoring_scheme

}; // ScoringSchemeTests
