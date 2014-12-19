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
#include <protocols/membrane/AddMembraneMover.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Package Headers
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/conformation/membrane/Exceptions.hh>
#include <core/conformation/membrane/types.hh>
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

	// test constructor from jumpnum
	void test_constructor_from_jumpnum () {
		
		TS_TRACE("Test compute_structure_based_membrane_position");
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;
		
		// 1AFO
		// read in pose and create topology object
		TS_TRACE("1AFO");
		Pose pose1;
		core::import_pose::pose_from_pdb( pose1, "protocols/membrane/1AFO_AB.pdb" );

		// create membrane pose
		Vector center( mem_center );
		Vector normal( mem_normal );
		AddMembraneMoverOP addmem1( new AddMembraneMover( center, normal, "protocols/membrane/geometry/1AFO_AB.span" ) );
		addmem1->apply( pose1 );

		// define vectors and object
		Vector center1(0.53225, 0.361, 0.095);
		Vector normal1(12.7367, 7.68036, -1.94611);
		
		// compute embedding
		EmbeddingDefOP embed1( compute_structure_based_membrane_position( pose1 ) );
		
		// compare
		TS_ASSERT( position_equal_within_delta( embed1->center(), center1, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed1->normal(), normal1, 0.001 ) );
		
		// 1BL8
		TS_TRACE("1BL8");
		Pose pose2;
		core::import_pose::pose_from_pdb( pose2, "protocols/membrane/geometry/1BL8_.pdb" );
		AddMembraneMoverOP addmem2( new AddMembraneMover( center, normal, "protocols/membrane/geometry/1BL8__tr.span" ) );
		addmem2->apply(pose2);
		Vector center2(73.9421, 26.7549, 24.4493);
		Vector normal2(-5.7604, 0.605734, -13.8366);
		EmbeddingDefOP embed2( compute_structure_based_membrane_position( pose2 ) );
		TS_ASSERT( position_equal_within_delta( embed2->center(), center2, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed2->normal(), normal2, 0.001 ) );

	}

	////////////////////////////////////////////////////////////////////////////////
	
	// test constructor from jumpnum and axis
	void test_constructor_from_jumpnum () {
		
		TS_TRACE("Test constructor from jumpnum and axis");
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;
		
	}
	
	////////////////////////////////////////////////////////////////////////////////

	// test constructor from jumpnum and angle
	void test_constructor_from_jumpnum_and_angle () {
		
		TS_TRACE("Test constructor from jumpnum and angle");
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;
		
	}
	
	////////////////////////////////////////////////////////////////////////////////

	// test constructor from jumpnum, axis and angle
	void test_constructor_from_jumpnum_axis_angle () {
		
		TS_TRACE("Test constructor from jumpnum, axis and angle");
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;
		
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
