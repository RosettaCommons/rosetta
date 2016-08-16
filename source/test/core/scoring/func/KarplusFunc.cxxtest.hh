// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/KarplusFunc.hh.cxxtest.hh
/// @brief  test suite for KarplusFunc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// basic headers
#include <basic/Tracer.hh>

// unit headers
#include <core/scoring/func/KarplusFunc.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class KarplusFuncTests : public CxxTest::TestSuite
{
public:


	void test_serialize_KarplusFunc() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::func;

		FuncOP instance( new KarplusFunc(1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7 ) ); // serialize this through a pointer to the base class

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		FuncOP instance2; // deserialize also through a pointer to the base class
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}

		// make sure the deserialized base class pointer points to a KarplusFunc
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< KarplusFunc > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};

