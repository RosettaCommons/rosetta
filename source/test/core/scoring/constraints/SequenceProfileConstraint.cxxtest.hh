// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/SequenceProfileConstraint.hh.cxxtest.hh
/// @brief  test suite for SequenceProfileConstraint
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// basic headers
#include <basic/Tracer.hh>

// unit headers
#include <core/scoring/constraints/SequenceProfileConstraint.hh>

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/sequence/SequenceProfile.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class SequenceProfileConstraintTests : public CxxTest::TestSuite
{
public:


	void test_serialize_SequenceProfileConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::sequence;
		using namespace core::id;

		utility::vector1< utility::vector1< core::Real > > prof;
		utility::vector1< core::Real > line1, line2;
		for ( core::Size ii = 1; ii <= 20; ++ii ) { line1.push_back( 1 + 0.1*ii ); line2.push_back( 2 + 0.05 *ii ); }
		prof.push_back( line1 ); prof.push_back( line2 );
		SequenceProfileOP seqprof( new SequenceProfile( prof, "ID", "seqprof_1" ));
		ConstraintOP instance( new SequenceProfileConstraint( 2, seqprof ) ); // serialize this through a pointer to the base class

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		ConstraintOP instance2; // deserialize also through a pointer to the base class
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}

		// make sure the deserialized base class pointer points to a SequenceProfileConstraint
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< SequenceProfileConstraint > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};

