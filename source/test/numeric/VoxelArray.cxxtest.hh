// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/VoxelArray.cxxtest.hh
/// @brief  numeric.VoxelArray.cxxtest: test suite for VoxelArray
/// @author Brian Coventry (bcov@uw.edu)

// Testing headers
#include <cxxtest/TestSuite.h>


// Project headers
#include <numeric/VoxelArray.hh>
#include <core/types.hh>

#include <numeric/random/random.hh>
#include <utility/vector0.hh>

// C/C++ headers
#include <iostream>
#include <string>
#include <sstream>

//#include <basic/Tracer.hh>

// static basic::Tracer TR( "numeric.VoxelArray.cxxtest.hh" );

class VoxelArrayTest : public CxxTest::TestSuite {
public:


	void test_1() {
		typedef utility::vector0< core::Real > F3;
		numeric::VoxelArray<float,float> a3( F3{-1,-2,-3}, F3{1,2,3}, F3{0.5, 0.5, 0.5} );
		TS_ASSERT_EQUALS(a3.shape()[0],5);
		TS_ASSERT_EQUALS(a3.shape()[1],9);
		TS_ASSERT_EQUALS(a3.shape()[2],13);
		a3[F3{0,0,0}] = 1;
		a3[F3{-1,1,2.3}] = 2;
		TS_ASSERT_EQUALS( a3[(F3{0,0,0})],1);
		TS_ASSERT_EQUALS( a3[(F3{0,0,0})] , 1 );
		TS_ASSERT_EQUALS( a3[(F3{-1,1,2.3})] , 2 );
	}

	void test_2() {
		typedef utility::vector0< core::Real > F3;
		numeric::VoxelArray<float,float> a3( F3{-1,-2,-3}, F3{1,2,3}, F3{0.49, 0.49, 0.49} );
		TS_ASSERT_EQUALS(a3.shape()[0],5);
		TS_ASSERT_EQUALS(a3.shape()[1],9);
		TS_ASSERT_EQUALS(a3.shape()[2],13);
		a3[F3{0,0,0}] = 1;
		a3[F3{-1,1,2.3}] = 2;
		TS_ASSERT_EQUALS(a3[(F3{0,0,0})],1);
		TS_ASSERT_EQUALS( a3[(F3{0,0,0})] , 1 );
		TS_ASSERT_EQUALS( a3[(F3{-1,1,2.3})] , 2 );
	}

	void test_3() {
		typedef utility::vector0< core::Real > F3;
		numeric::VoxelArray<double,double> a(F3{-6,-6,-6},F3{7,7,7},F3{1.6345,1.6345,1.6345});
		for ( size_t i = 0; i < a.num_elements(); ++i ) a.data()[i] = numeric::random::uniform();

		numeric::VoxelArray<double,double> b(F3{-6,-6,-6},F3{7,7,7},F3{1.6345,1.6345,1.6345});
		TS_ASSERT( !(a == b) );
		b = a;
		TS_ASSERT( a == b );
		for ( size_t i = 0; i < b.num_elements(); ++i ) b.data()[i] = 0;
		TS_ASSERT( !(a == b) );

		std::ostringstream oss;
		a.save( oss );
		std::istringstream iss(oss.str());
		numeric::VoxelArray<double,double> save_read_a;
		save_read_a.load( iss );

		TS_ASSERT( a == save_read_a );

	}


};
