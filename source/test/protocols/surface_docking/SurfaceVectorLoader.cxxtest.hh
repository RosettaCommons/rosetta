// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/surface_docking/SurfaceVectorLoader.cxxtest.hh
/// @brief test suite for protocols/surface_docking/SurfaceVectorLoader.cc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <protocols/surface_docking/SurfaceVectorLoader.hh>
#include <protocols/surface_docking/SurfaceParameters.hh>
#include <protocols/loops/LoopsFileOptions.hh> //this is crazy

// Utility headers


// Numeric headers

// C++ headers
#include <string>

using namespace protocols::surface_docking;

class SurfaceVectorLoaderTest : public CxxTest::TestSuite {

public:
	void setUp() {
		protocols_init();
	}

	// @brief test create_resource method
	void test_create_SurfaceParameters_from_SurfaceVectorLoader() {
		SurfaceVectorLoader loader;
		std::string surf_vec_file( "-28.188  33.809   3.201\n-46.450  30.380   6.325\n-31.129  24.953   4.699" );
		std::string surf_vec_file2( "-28.188  33.809   3.201\n-46.450  30.380   6.325\n-31.129  24.953   4.699\n" );
		std::istringstream lstream( surf_vec_file );
		std::istringstream lstream2( surf_vec_file2);

		protocols::loops::LoopsFileOptions opts;


		utility::pointer::ReferenceCountOP resource = loader.create_resource( opts, "unit_test", lstream );
		utility::pointer::ReferenceCountOP resource2 = loader.create_resource( opts, "unit_test", lstream2 );

		TS_ASSERT( resource ); // make sure a resource was returned

		SurfaceParametersOP spptr = utility::pointer::dynamic_pointer_cast< protocols::surface_docking::SurfaceParameters > ( resource );
		TS_ASSERT( spptr ); // make sure we're actually returned the correct type

		//SurfaceParameters const & sp( *spptr() ); // no copy ctor!
		//std::cout << sp << std::endl;
		/*TS_ASSERT( lfd.size() == 2 );
		TS_ASSERT( lfd[ 1 ].start_res().pose_index() == 1 );
		TS_ASSERT( lfd[ 1 ].end_res().pose_index() == 4 );
		TS_ASSERT( lfd[ 2 ].start_res().pose_index() == 5 );
		TS_ASSERT( lfd[ 2 ].end_res().pose_index() == 7 );
		*/
	}


};
