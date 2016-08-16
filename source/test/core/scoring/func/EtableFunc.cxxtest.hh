// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/EtableFunc.hh.cxxtest.hh
/// @brief  test suite for EtableFunc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// basic headers
#include <basic/Tracer.hh>

// unit headers
#include <core/scoring/func/EtableFunc.hh>

#ifdef SERIALIZATION

// C++ headers used in the serialization unit test
#include <string>
#include <sstream>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class EtableFuncTests : public CxxTest::TestSuite
{
public:


	void test_serialize_EtableFunc() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::func;

		//                       min max eleven values spanning the range w/ step size of 0.1
		std::string func_vals = "1.5 2.5 1.0 0.9 0.8 0.6 0.3 -0.1 -0.2 -0.2 -0.1 -0.05 0";
		std::istringstream func_iss( func_vals );

		FuncOP instance( new EtableFunc( 5.5, 6.6, 0.01 ) ); // these values are overwritten in the call to read_data
		instance->read_data( func_iss );

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

		// make sure the deserialized base class pointer points to a EtableFunc
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< EtableFunc > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};

