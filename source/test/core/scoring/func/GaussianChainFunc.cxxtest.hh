// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/GaussianChainFunc.hh.cxxtest.hh
/// @brief  test suite for GaussianChainFunc
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
#include <core/types.hh>

#include <utility/tools/make_vector1.hh>

#include <ObjexxFCL/format.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

static basic::Tracer TR("test.core.scoring.func.GaussianChainFunc");

class GaussianChainFuncTests : public CxxTest::TestSuite
{
public:

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	void test_serialize_GaussianChainFunc()
	{
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::func;

		utility::vector1< core::Real > other_dists( 5 );
		for ( core::Size ii = 1; ii  <= 5; ++ii ) other_dists.push_back( ii*1.1 );
		FuncOP instance( new GaussianChainFunc( 1.234, 2.345, other_dists ) ); // serialize this through a pointer to the base class

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

		// make sure the deserialized base class pointer points to a GaussianChainFunc
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< GaussianChainFunc > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}


	void double_func_test( core::Distance const distance2, core::Real const &loop_fixed_cost,
		bool const & single_gaussian_approximation_should_be_good )
	{
		using namespace core;
		using namespace core::scoring::func;
		using ObjexxFCL::format::F;

		utility::vector1< Distance > other_distances;
		other_distances.push_back( distance2 );

		Real const gaussian_variance( 5.0 * 5.0 );
		FuncOP func2( new GaussianChainFunc( gaussian_variance, loop_fixed_cost, other_distances ) );
		FuncOP func1( new GaussianChainFunc( gaussian_variance, loop_fixed_cost ) );
		GaussianChainFuncOP func2_approx( new GaussianChainFunc( gaussian_variance, loop_fixed_cost, other_distances ) );
		func2_approx->set_force_combined_gaussian_approximation( true );

		// sweep through distance1.
		TR << "-- GAUSSIAN_CHAIN_DOUBLE_FUNC --" << std::endl;
		TR << "Fixing Distance2 at " << distance2 << std::endl;
		for ( Size i = 1; i <= 20; i++ ) {
			Distance const d = Real( i ) * 1.0;
			TR << "Distance1 " << F(6,2,d) << "   FUNC_EXACT " << F(10,5,func2->func( d )) << " FUNC_APPROX " << F(10,5,func2_approx->func( d )) << " SINGLE_FUNC " << F(10,5,func1->func( d )) << std::endl;

			TS_ASSERT_DIFFERS( func2->func( d ), func2_approx->func( d ) );
			TS_ASSERT_DIFFERS( func2->func( d ), func1->func( d ) );
			if ( single_gaussian_approximation_should_be_good ) {
				TS_ASSERT_DELTA( func2->func( d ), func2_approx->func( d ), 0.01 );
			}

		}
		TR << std::endl;

	}

	void test_GaussianDoubleChainFunc()
	{
		core::Real const loop_fixed_cost( basic::options::option[ basic::options::OptionKeys::score::loop_close::loop_fixed_cost ]() );
		double_func_test( 0.5, loop_fixed_cost, true /*approximation_should_be_good*/ );
		double_func_test( 15.0, loop_fixed_cost, false /*approximation_should_be_bad*/ );
	}


	void triple_func_test( core::Distance const distance2, core::Distance const distance3, core::Real const &loop_fixed_cost,
		bool const & single_gaussian_approximation_should_be_good  )
	{
		using namespace core;
		using namespace core::scoring::func;
		using ObjexxFCL::format::F;
		using utility::tools::make_vector1;

		utility::vector1< Distance > other_distances;
		other_distances.push_back( distance2 );
		other_distances.push_back( distance3 );

		Real const gaussian_variance( 5.0 * 5.0 );
		FuncOP func1( new GaussianChainFunc( gaussian_variance, loop_fixed_cost ) ); // single func
		FuncOP func2( new GaussianChainFunc( gaussian_variance, loop_fixed_cost, make_vector1( distance2 ) ) );
		FuncOP func3( new GaussianChainFunc( gaussian_variance, loop_fixed_cost, other_distances ) );
		GaussianChainFuncOP func3_approx( new GaussianChainFunc( gaussian_variance, loop_fixed_cost, other_distances ) );
		func3_approx->set_force_combined_gaussian_approximation( true );

		FuncOP func0( new GaussianChainFunc( gaussian_variance, loop_fixed_cost ) ); // single func

		// sweep through distance1.
		TR << "-- GAUSSIAN_CHAIN_TRIPLE_FUNC --" << std::endl;
		TR << "Fixing Distance2 at " << distance2 << ", Distance3 at " << distance3 << std::endl;
		for ( Size i = 1; i <= 20; i++ ) {
			Distance const d = Real( i ) * 1.0;
			TR << "Distance1 " << F(6,2,d) << "   FUNC_EXACT " << F(10,5,func3->func( d )) << " FUNC_APPROX " << F(10,5,func3_approx->func( d )) << " DOUBLE_FUNC " << F(10,5,func2->func( d )) << " SINGLE_FUNC " << F(10,5,func1->func( d ) ) << std::endl;

			TS_ASSERT_DIFFERS( func3->func( d ), func3_approx->func( d ) );
			TS_ASSERT_DIFFERS( func3->func( d ), func2->func( d ) );
			TS_ASSERT_DIFFERS( func3->func( d ), func1->func( d ) );
			if ( single_gaussian_approximation_should_be_good ) {
				TS_ASSERT_DELTA( func3->func( d ), func3_approx->func( d ), 0.01 );
			}

		}
		TR << std::endl;
	}


	void test_GaussianTripleChainFunc()
	{
		core::Real const loop_fixed_cost( basic::options::option[ basic::options::OptionKeys::score::loop_close::loop_fixed_cost ]() );
		triple_func_test( 0.5, 0.5, loop_fixed_cost, true /*approximation_should_be_good*/ );
		triple_func_test( 15.0, 0.5, loop_fixed_cost, false /*approximation_should_be_bad*/  );
		triple_func_test( 15.0, 15.0, loop_fixed_cost, false /*approximation_should_be_bad*/  );
	}

	void
	quadruple_func_test( core::Distance const distance2, core::Distance const distance3, core::Distance const distance4,
		core::Real const &loop_fixed_cost,
		bool const & single_gaussian_approximation_should_be_good  )
	{
		using namespace core;
		using namespace core::scoring::func;
		using ObjexxFCL::format::F;

		utility::vector1< Distance > other_distances;
		other_distances.push_back( distance2 );
		other_distances.push_back( distance3 );
		other_distances.push_back( distance4 );

		Real const gaussian_variance( 5.0 * 5.0 );
		// FuncOP func1( new GaussianChainFunc( gaussian_variance, loop_fixed_cost ) ); // single func
		// FuncOP func2( new GaussianChainFunc( gaussian_variance, loop_fixed_cost, make_vector1( distance2 ) ) );
		FuncOP func4( new GaussianChainFunc( gaussian_variance, loop_fixed_cost, other_distances ) );
		GaussianChainFuncOP func4_approx( new GaussianChainFunc( gaussian_variance, loop_fixed_cost, other_distances ) );
		func4_approx->set_force_combined_gaussian_approximation( true );

		FuncOP func0( new GaussianChainFunc( gaussian_variance, loop_fixed_cost ) ); // single func

		// sweep through distance1.
		TR << "-- GAUSSIAN_CHAIN_QUADRUPLE_FUNC --" << std::endl;
		TR << "Fixing Distance2 at " << distance2 << ", Distance3 at " << distance3  << ", Distance4 at " << distance4 << std::endl;
		for ( Size i = 1; i <= 20; i++ ) {
			Distance const d = Real( i ) * 1.0;
			TR << "Distance1 " << F(6,2,d) << "   FUNC_EXACT " << F(10,5,func4->func( d )) << " FUNC_APPROX " << F(10,5,func4_approx->func( d )) << std::endl; //<< " DOUBLE_FUNC " << F(10,5,func2->func( d )) << " SINGLE_FUNC " << F(10,5,func1->func( d ) ) << std::endl;

			TS_ASSERT_DIFFERS( func4->func( d ), func0->func( d ) );
			TS_ASSERT_DIFFERS( func4->func( d ), func4_approx->func( d ) );
			if ( single_gaussian_approximation_should_be_good ) {
				TS_ASSERT_DELTA( func4->func( d ), func4_approx->func( d ), 0.01 );
			}

		}
		TR << std::endl;
	}

	void test_GaussianQuadrupleChainFunc() {
		core::Real const loop_fixed_cost( basic::options::option[ basic::options::OptionKeys::score::loop_close::loop_fixed_cost ]() );
		quadruple_func_test( 0.01, 0.01, 0.01, loop_fixed_cost, true /*approximation_should_be_good*/ );
		quadruple_func_test( 15.0, 15.0, 15.0, loop_fixed_cost, false /*approximation_should_be_bad*/  );
	}

};

