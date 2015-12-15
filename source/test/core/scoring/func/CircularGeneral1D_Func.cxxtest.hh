// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/CircularGeneral1D_Func.hh.cxxtest.hh
/// @brief  test suite for CircularGeneral1D_Func
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// basic headers
#include <basic/Tracer.hh>

// unit headers
#include <core/scoring/func/CircularGeneral1D_Func.hh>

#ifdef	SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class CircularGeneral1D_FuncTests : public CxxTest::TestSuite
{
public:


	void test_serialize_CircularGeneral1D_Func() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::func;

		ObjexxFCL::FArray1D< core::Real > func_values( 5 );
		for ( core::Size ii = 1; ii <= 5; ++ii ) { func_values( ii ) = 1.1 * ii; }
		FuncOP instance( new CircularGeneral1D_Func( func_values, 1.234, 2.345 ) ); // serialize this through a pointer to the base class

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

		// make sure the deserialized base class pointer points to a CircularGeneral1D_Func
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< CircularGeneral1D_Func > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};

