// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/GaussianChainGeneralFunc.hh.cxxtest.hh
/// @brief  test suite for GaussianChainGeneralFunc
/// @detailed adapted from src/apps/pilot/rhiju/gaussian_chain_func_test.cc
/// @author Rhiju Das
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <test/core/init_util.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// unit headers
#include <core/scoring/func/GaussianChainFunc.hh>
#include <core/scoring/func/GaussianChainGeneralFunc.hh>
#include <core/types.hh>

#include <utility/tools/make_vector1.hh>

#include <ObjexxFCL/format.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

static basic::Tracer TR("test.core.scoring.func.GaussianChainGeneralFunc");

using namespace utility;
using namespace utility::tools;
using namespace core;

class GaussianChainGeneralFuncTests : public CxxTest::TestSuite
{
public:

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	void test_serialize_GaussianChainGeneralFunc()
	{
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::func;

		utility::vector1< core::Real > other_dists( 5 );
		for ( core::Size ii = 1; ii  <= 5; ++ii ) other_dists.push_back( ii*1.1 );
		FuncOP instance( new GaussianChainGeneralFunc( 1.234, 2.345, other_dists ) ); // serialize this through a pointer to the base class

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

		// make sure the deserialized base class pointer points to a GaussianChainGeneralFunc
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< GaussianChainGeneralFunc > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}


	Real
	get_gaussian_chain_func( vector1< Distance > const & dists,
													 Real const & gaussian_variance,
													 Real const & loop_fixed_cost,
													 Size const idx,
													 bool const force_use_general )
	{
		using namespace core::scoring::func;
		Distance dist( dists[ idx ] );
		vector1< Distance > other_dists;
		for ( Size k = 1; k <= dists.size(); k++ ) {
			if ( k != idx ) other_dists.push_back( dists[ k ] );
		}
		FuncOP func;
		if ( force_use_general ) {
			func = FuncOP( new GaussianChainGeneralFunc( gaussian_variance, loop_fixed_cost, other_dists ) );
		} else {
			func = FuncOP( new GaussianChainFunc( gaussian_variance, loop_fixed_cost, other_dists ) );
		}

		// test derivs.
		Distance delta = 1.0e-6;
		//		TR << "Checking deriv " << ( func->func(dist + delta) - func->func(dist ) )/delta << " vs " << func->dfunc( dist ) << std::endl;

		TS_ASSERT_DELTA( ( func->func(dist + delta) - func->func(dist ) )/delta, func->dfunc( dist ), 1.0e-4 );

		return ( func->func( dist ) );
	}

	void test_GaussianChainGeneralFunc()
	{
		vector1< Real > all_dists = make_vector1(  13.7201, 9.0501, 10.7559 , 16.4339 , 9.93362 );
		Real const gaussian_variance = 20.0 * 20.0;
		Real loop_fixed_cost = 0.0;
		for ( Size n = 1; n <= all_dists.size(); n++ ) {

			vector1< Real > dists;
			for ( Size k = 1; k <= n; k++ ) dists.push_back( all_dists[ k ] );

			Real val1 = get_gaussian_chain_func( dists, gaussian_variance, loop_fixed_cost, 1, false );
			Real val2 = get_gaussian_chain_func( dists, gaussian_variance, loop_fixed_cost, std::min( Size(2), n ), false );
			Real val1_general = get_gaussian_chain_func( dists, gaussian_variance, loop_fixed_cost, 1, true );
			Real val2_general = get_gaussian_chain_func( dists, gaussian_variance, loop_fixed_cost, std::min( Size(2), n ), true );
			TS_ASSERT_DELTA( val1, val2, 1.0e-4 );
			TS_ASSERT_DELTA( val1, val1_general, 1.0e-4 );
			TS_ASSERT_DELTA( val1, val2_general, 1.0e-4 );
		}
	}

};

