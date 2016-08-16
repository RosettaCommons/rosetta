// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/DihedralPairConstraint.hh.cxxtest.hh
/// @brief  test suite for DihedralPairConstraint
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// basic headers
#include <basic/Tracer.hh>

// unit headers
#include <core/scoring/constraints/DihedralPairConstraint.hh>

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class DihedralPairConstraintTests : public CxxTest::TestSuite
{
public:


	void test_serialize_DihedralPairConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		AtomID at1a( 1, 1), at2a( 2, 1 ), at3a( 3, 1 ), at4a( 4, 1 );
		AtomID at1b( 1, 4), at2b( 2, 4 ), at3b( 3, 4 ), at4b( 4, 4 );
		// serialize this through a pointer to the base class
		ConstraintOP instance( new DihedralPairConstraint(at1a,at2a,at3a,at4a,at1b,at2b,at3b,at4b,some_func) );

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

		// make sure the deserialized base class pointer points to a DihedralPairConstraint
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< DihedralPairConstraint > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};

