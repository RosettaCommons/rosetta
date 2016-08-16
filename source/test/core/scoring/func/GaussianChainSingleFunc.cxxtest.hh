// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/GaussianChainSingleFunc.hh.cxxtest.hh
/// @brief  test suite for GaussianChainSingleFunc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// basic headers
#include <basic/Tracer.hh>

// unit headers
#include <core/scoring/func/GaussianChainSingleFunc.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class GaussianChainSingleFuncTests : public CxxTest::TestSuite
{
public:


	void test_serialize_GaussianChainSingleFunc() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::func;

		FuncOP instance( new GaussianChainSingleFunc( 1.234, 2.345 ) ); // serialize this through a pointer to the base class

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

		// make sure the deserialized base class pointer points to a GaussianChainSingleFunc
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< GaussianChainSingleFunc > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};

