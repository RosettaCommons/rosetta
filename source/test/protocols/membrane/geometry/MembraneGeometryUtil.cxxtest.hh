// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 protocols/membrane/geometry/MembraneGeometryUtil.cxxtest.hh
/// @brief 	 Unit test for util functions
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/conformation/Residue.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>
#include <protocols/membrane/geometry/util.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Package Headers
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/conformation/membrane/Exceptions.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using namespace core;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace protocols::membrane::geometry;

class MembraneGeometryUtilTest : public CxxTest::TestSuite {
	
public: // test functions
    
    /// Test Setup Functions ////////
    
    /// @brief Setup Test
    void setUp(){
	
        // Initialize
        core_init();
	}
    
    /// @brief Standard Tear Down
    void tearDown() {}
	    
    ///// Test Methods /////////////

////////////////////////////////////////////////////////////////////////////////

	// TODO: compute_structure_based_membrane_position: replace it!

	// check vector for reasonable size
	void test_check_vector() {
		
		TS_TRACE("Test check vector");
				
		// define vectors and object
		Vector v1(1, 2, 1000);
		Vector v2(4, 5000, 3);

		// check for validity
		try {
			check_vector( v1 );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_error_message = "Unreasonable range for center or normal! Check your input vectors!";
			TS_ASSERT( expected_error_message == e.msg() );
		}
	}
	
	// average embeddings
	void test_average_embeddings() {
		
		TS_TRACE("Test average embeddings");
		
		// define vectors
		Vector v1(1, 2, 3);
		Vector v2(4, 5, 6);
		Vector v3(7, 8, 9);
		Vector v4(7, 5, 3);
		Vector avg_center(4, 5, 6);
//		Vector avg_normal(1.333, 0.666, 0);
		Vector avg_normal(6, 6, 6);
		
		// create embedding objects
		EmbeddingDefOP emb1( new EmbeddingDef( v1, v2 ) );
		EmbeddingDefOP emb2( new EmbeddingDef( v2, v3 ) );
		EmbeddingDefOP emb3( new EmbeddingDef( v3, v4 ) );
		
		// put them into a vector
		utility::vector1< EmbeddingDefOP > embeddings;
		embeddings.push_back( emb1 );
		embeddings.push_back( emb2 );
		embeddings.push_back( emb3 );
		
		// average them
		EmbeddingDefOP avg = average_embeddings( embeddings );
		
		// check for equality
		TS_ASSERT( position_equal_within_delta( avg->center(), avg_center, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( avg->normal(), avg_normal, 0.001 ) );
	}
	
		
////////////////////////////////////////////////////////////////////////////////

	/// @brief Position equal within delta (helper method)
	bool position_equal_within_delta( Vector a, Vector b, Real delta ) {
		
		TS_ASSERT_DELTA( a.x(), b.x(), delta );
		TS_ASSERT_DELTA( a.y(), b.y(), delta );
		TS_ASSERT_DELTA( a.z(), b.z(), delta );
		
		return true;
	}
};
