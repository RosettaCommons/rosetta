// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/func/SOGFunc.cxxtest.hh
/// @brief  test suite for SOGFunc function
/// @author James Thompson

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/func/SOGFunc_Impl.hh>
#include <core/types.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

using core::Size;
using core::Real;
using utility::vector1;
using namespace core::scoring::func;

class SOGFuncTests : public CxxTest::TestSuite {
public:
	void setUp() {
		core_init();
	}

	void test_func() {
		vector1< Real > means, sdevs, weights;
		means.push_back( 0.5 );
		sdevs.push_back( 1.0 );
		weights.push_back( 0.3333 );

		means.push_back( 1.5 );
		sdevs.push_back( 0.5 );
		weights.push_back( 0.3333 );

		means.push_back( 5.0 );
		sdevs.push_back( 2.0 );
		weights.push_back( 0.3334 );

		SOGFunc_ImplOP func( new SOGFunc_Impl(means, sdevs, weights) );
		Real const TOLERATED_ERROR = 0.001;
		Real const start = 2.0;
		Real const end = 20.0;
		Real const res = 0.5;

		Real func_values[] = {
			1.487, 2.4719, 3.01347, 2.96131, 2.83084, 2.74136, 2.71072, 2.74204,
			2.8358, 2.99205, 3.2108, 3.49205, 3.8358, 4.24205, 4.7108, 5.24205,
			5.8358, 6.49205, 7.2108, 7.99205, 8.8358, 9.74205, 10.7108, 11.742,
			12.8358, 13.992, 15.2108, 16.492, 17.8358, 19.242, 20.7108, 22.242,
			23.8358, 25.492, 27.2108, 28.992
			};

		Real dfunc_values[] = {
			1.642, 1.906, 0.248, -0.264, -0.231, -0.122, 0.000, 0.125, 0.250, 0.375,
			0.500, 0.625, 0.750, 0.875, 1.000, 1.125, 0.625, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0
			};

		const Real upper_bound_dist  = func->upper_bound();
		const Real upper_bound_score = func->upper_bound_score();
		const Size nsteps = Size((end - start) / res);

		for ( Size i = 0; i < nsteps; ++i ) {
			Real r = start + (i * res);

			Real actual_f = func->func(r);
			Real expect_f = (r >= upper_bound_dist) ? 0 : func_values[i] - upper_bound_score;
			TS_ASSERT_DELTA(actual_f, expect_f, TOLERATED_ERROR);

			Real actual_d = func->dfunc(r);
			Real expect_d = (r > upper_bound_dist) ? 0 : dfunc_values[i];
			TS_ASSERT_DELTA(actual_d, expect_d, TOLERATED_ERROR);
		}
	}

	void test_serialize_SOGFunc() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::func;

		vector1< Real > means, sdevs, weights;
		means.push_back( 0.5 );
		sdevs.push_back( 1.0 );
		weights.push_back( 0.3333 );

		means.push_back( 1.5 );
		sdevs.push_back( 0.5 );
		weights.push_back( 0.3333 );

		means.push_back( 5.0 );
		sdevs.push_back( 2.0 );
		weights.push_back( 0.3334 );

		FuncOP instance( new SOGFunc(means, sdevs, weights) ); // serialize this through a pointer to the base class

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

		// make sure the deserialized base class pointer points to a SOGFunc
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< SOGFunc > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};
