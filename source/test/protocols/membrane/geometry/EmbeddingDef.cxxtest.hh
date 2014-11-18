// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 protocols/membrane/geometry/EmbeddingDef.cxxtest.hh
/// @brief 	 Unit test for EmbeddingDef class
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/conformation/Residue.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
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

class EmbeddingDefTest : public CxxTest::TestSuite {
	
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

	// default constructor
	void test_default_constructor() {
		
		TS_TRACE("Test default constructor");
				
		// define vectors and object
		Vector center(0, 0, 0);
		Vector normal(0, 0, 1);
		EmbeddingDefOP embed( new EmbeddingDef() );

		// check positions
		TS_ASSERT( position_equal_within_delta( embed->center(), center, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed->normal(), normal, 0.001 ) );
	}
	
	// standard constructor
	void test_standard_constructor() {
	  
		TS_TRACE("Test constructor from center and normal");

		// define vectors and object
		Vector center(5, 4, 3);
		Vector normal(1, 5, 8);
		EmbeddingDefOP embed( new EmbeddingDef( center, normal ) );

		// check positions
		TS_ASSERT( position_equal_within_delta( embed->center(), center, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed->normal(), normal, 0.001 ) );
	}

	// from span
	void test_from_span() {
		
		TS_TRACE("Test from_span function");
		
		// read in pose
		core::pose::PoseOP pose = core::import_pose::pose_from_pdb("protocols/membrane/geometry/1AFO_AB.pdb");
		
		// create object
		Size res1(15);
		Size res2(32);
		EmbeddingDefOP embed( new EmbeddingDef( pose, res1, res2 ) );
		
		// define start, end, center and normal
		Vector start(-0.97, -1.864, -12.281);
		Vector end(-1.246, -7.692, 13.217);
		Vector center(-1.108, -4.778, 0.468);
		Vector normal(-0.0105517, -0.222808, 0.974805);
		
		// check positions
		TS_ASSERT( position_equal_within_delta( pose->residue( res1 ).atom( 2 ).xyz(), start, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( pose->residue( res2 ).atom( 2 ).xyz(), end, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed->center(), center, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed->normal(), normal, 0.001 ) );
		
	}

	// translate by
	void test_translate_by() {

		TS_TRACE("Test translate function");

		// define vectors and object
		Vector center(5, 4, 3);
		Vector normal(1, 5, 8);
		Vector add2center(2, 3, 4);
		Vector add2normal(8, 4, 1);
		EmbeddingDefOP embed( new EmbeddingDef( center, normal ) );
		EmbeddingDefOP translation( new EmbeddingDef( add2center, add2normal ) );
		
		// define new center and normal
		Vector new_center(7, 7, 7);
		Vector new_normal(9, 9, 9);

		embed->translate_by( translation );

		// check positions
		TS_ASSERT( position_equal_within_delta( embed->center(), new_center, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed->normal(), new_normal, 0.001 ) );
	}
	
	// equals
	void test_equals() {

		TS_TRACE("Test equals function");

		// define vectors and object
		Vector center(5, 4, 3);
		Vector normal(1, 5, 8);
		Vector center1(2, 3, 4);
		Vector normal1(8, 4, 1);
		EmbeddingDefOP embed( new EmbeddingDef( center, normal ) );
		EmbeddingDefOP embed1( new EmbeddingDef( center, normal ) );
		EmbeddingDefOP embed2( new EmbeddingDef( center1, normal1 ) );

		// check positions
		TS_ASSERT( embed->equals( *embed1 ) );
		TS_ASSERT( ! embed->equals( *embed2 ) );
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
