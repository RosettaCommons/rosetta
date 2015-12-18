// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/LocalCoordinateConstraint.hh.cxxtest.hh
/// @brief  test suite for LocalCoordinateConstraint
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// basic headers
#include <basic/Tracer.hh>

// unit headers
#include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class LocalCoordinateConstraintTests : public CxxTest::TestSuite
{

public:

	void test_serialize_LocalCoordinateConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		AtomID atid( 4, 5 );
		StubID stubid( AtomID( 2, 3 ), AtomID( 2, 3 ), AtomID( 3, 1 ), AtomID( 3, 2 ) );
		core::Vector target( 1.125, 2.5, -3.625 );
		LocalCoordinateConstraintOP instance( new LocalCoordinateConstraint( atid, stubid, target, some_func ) );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		LocalCoordinateConstraintOP instance2;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};

