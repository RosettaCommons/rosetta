// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/func/HarmonicFunc.cxxtest.hh
/// @brief  test suite for HarmonicFunc function
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/HarmonicFunc.fwd.hh>
#include <core/types.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

//Auto Headers


//static basic::Tracer TR("core.scoring.constraints.HarmonicFunc.cxxtest");

using namespace core;

class HarmonicFuncTests : public CxxTest::TestSuite {

public:
	HarmonicFuncTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {

	}


	///////////////////////////////////////////////////////////////////////////////
	// ------------------------------------------ //
	/// @brief simple test minimization
	void test_harmonic_func()
	{
		using namespace core::scoring::func;

		HarmonicFuncOP func( new HarmonicFunc( 5.16, 1.5 ) );

		float const TOLERATED_ERROR = 0.001;
		const core::Real start = 2;
		const core::Real end   = 20;
		const core::Real res   = 0.5;


		core::Real func_values[] = {
			4.43804, 3.14471, 2.0736, 1.22471, 0.598044, 0.1936, 0.0113778,
			0.0513778, 0.3136, 0.798044, 1.50471, 2.4336, 3.58471, 4.95804, 6.5536,
			8.37138, 10.4114, 12.6736, 15.158, 17.8647, 20.7936, 23.9447, 27.318,
			30.9136, 34.7314, 38.7714, 43.0336, 47.518, 52.2247, 57.1536, 62.3047,
			67.678, 73.2736, 79.0914, 85.1314, 91.3936
			};

		core::Real dfunc_values[] = {
			-2.80889, -2.36444, -1.92, -1.47556, -1.03111, -0.586667, -0.142222,
			0.302222, 0.746667, 1.19111, 1.63556, 2.08, 2.52444, 2.96889, 3.41333,
			3.85778, 4.30222, 4.74667, 5.19111, 5.63556, 6.08, 6.52444, 6.96889,
			7.41333, 7.85778, 8.30222, 8.74667, 9.19111, 9.63556, 10.08, 10.5244,
			10.9689, 11.4133, 11.8578, 12.3022, 12.7467
			};

		core::Size nsteps = core::Size( ( end - start ) / res );
		for ( core::Size i = 0; i < nsteps; ++i ) {
			core::Real r = start + (i * res);
			TS_ASSERT_DELTA( func->func(r),   func_values[i], TOLERATED_ERROR );
			TS_ASSERT_DELTA( func->dfunc(r), dfunc_values[i], TOLERATED_ERROR );
			TS_ASSERT_DELTA( func->dfunc(r), func->estimate_dfunc(r), TOLERATED_ERROR );
			// TR << r << ' ' << func->func(r) << ' ' << func->dfunc(r) << "\n";
		}
		//TR.flush();

	} // test_harmonic_func

	void test_serialize_HarmonicFunc() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::func;

		FuncOP instance( new HarmonicFunc( 5.16, 1.5 ) ); // serialize this through a pointer to the base class

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

		// make sure the deserialized base class pointer points to a HarmonicFunc
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< HarmonicFunc > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};
