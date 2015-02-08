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
#include <core/conformation/membrane/SpanningTopology.hh>
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

	// compute_structure_based_embedding
	void test_compute_structure_based_membrane_position() {
		
		TS_TRACE("Test compute_structure_based_embedding");
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;
		
		// 1AFO
		// read in pose and create topology object
		TS_TRACE("1AFO");
		Pose pose1;
		core::import_pose::pose_from_pdb( pose1, "protocols/membrane/geometry/1AFO_.pdb" );

		// create membrane pose
		AddMembraneMoverOP addmem1( new AddMembraneMover( "protocols/membrane/geometry/1AFO__tr.span" ) );
		addmem1->apply( pose1 );

		// define vectors and object
		Vector center1(0.53225, 0.361, 0.095);
		Vector normal1(-12.7367, -7.68036, 1.94611);
		
		// compute embedding
		EmbeddingDefOP embed1( compute_structure_based_embedding( pose1 ) );
		
		// compare
		TS_ASSERT( position_equal_within_delta( embed1->center(), center1, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed1->normal(), normal1, 0.001 ) );
		
		// 1BL8
		TS_TRACE("1BL8");
		Pose pose2;
		core::import_pose::pose_from_pdb( pose2, "protocols/membrane/geometry/1BL8_.pdb" );
		AddMembraneMoverOP addmem2( new AddMembraneMover( "protocols/membrane/geometry/1BL8__tr.span" ) );
		addmem2->apply(pose2);
		Vector center2(73.9421, 26.7549, 24.4493);
		Vector normal2(5.7604, -0.605734, 13.8366);
		EmbeddingDefOP embed2( compute_structure_based_embedding( pose2 ) );
		TS_ASSERT( position_equal_within_delta( embed2->center(), center2, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed2->normal(), normal2, 0.001 ) );

		// 1QJP - beta-barrel
		TS_TRACE("1QJP");
		Pose pose3;
		core::import_pose::pose_from_pdb( pose3, "protocols/membrane/geometry/1QJP_.pdb" );
		AddMembraneMoverOP addmem3( new AddMembraneMover( "protocols/membrane/geometry/1QJP__tr.span" ) );
		addmem3->apply(pose3);
		Vector center3(31.2161, 16.9685, 37.6119);
		Vector normal3(13.1689, -7.07507, 1.23442);
		EmbeddingDefOP embed3( compute_structure_based_embedding( pose3 ) );
		TS_ASSERT( position_equal_within_delta( embed3->center(), center3, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed3->normal(), normal3, 0.001 ) );

		// 2BS2
		TS_TRACE("2BS2");
		Pose pose4;
		core::import_pose::pose_from_pdb( pose4, "protocols/membrane/geometry/2BS2_CF.pdb" );
		AddMembraneMoverOP addmem4( new AddMembraneMover( "protocols/membrane/geometry/2BS2_CF_tr.span" ) );
		addmem4->apply(pose4);
		Vector center4(21.4326, 6.0464, -41.0573);
		Vector normal4(0.0900585, 0.176022, 14.9987);
		EmbeddingDefOP embed4( compute_structure_based_embedding( pose4 ) );
		TS_ASSERT( position_equal_within_delta( embed4->center(), center4, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed4->normal(), normal4, 0.001 ) );

		// 2MPN
		TS_TRACE("2MPN");
		Pose pose5;
		core::import_pose::pose_from_pdb( pose5, "protocols/membrane/geometry/2MPN_.pdb" );
		AddMembraneMoverOP addmem5( new AddMembraneMover( "protocols/membrane/geometry/2MPN__tr.span" ) );
		addmem5->apply(pose5);
		Vector center5(0.3645, 3.66025, 41.345);
		Vector normal5(0.104341, 14.9499, 1.221);
		EmbeddingDefOP embed5( compute_structure_based_embedding( pose5 ) );
		TS_ASSERT( position_equal_within_delta( embed5->center(), center5, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed5->normal(), normal5, 0.001 ) );

		// 2OAR
		TS_TRACE("2OAR");
		Pose pose6;
		core::import_pose::pose_from_pdb( pose6, "protocols/membrane/geometry/2OAR_.pdb" );
		AddMembraneMoverOP addmem6( new AddMembraneMover( "protocols/membrane/geometry/2OAR__tr.span" ) );
		addmem6->apply(pose6);
		Vector center6(18.8453, 122.117, 1.079);
		Vector normal6(-11.3264, -9.83365, 0.106647);
		EmbeddingDefOP embed6( compute_structure_based_embedding( pose6 ) );
		TS_ASSERT( position_equal_within_delta( embed6->center(), center6, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed6->normal(), normal6, 0.001 ) );

		// 2UUH
		TS_TRACE("2UUH");
		Pose pose7;
		core::import_pose::pose_from_pdb( pose7, "protocols/membrane/geometry/2UUH__tr.pdb" );
		AddMembraneMoverOP addmem7( new AddMembraneMover( "protocols/membrane/geometry/2UUH__tr.span" ) );
		addmem7->apply(pose7);
		Vector center7(-0.000166667, -0.000125, 0.295625);
		Vector normal7(1.19923e-05, 2.70232e-05, 15);
		EmbeddingDefOP embed7( compute_structure_based_embedding( pose7 ) );
		TS_ASSERT( position_equal_within_delta( embed7->center(), center7, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed7->normal(), normal7, 0.001 ) );

		// 3PXO
		TS_TRACE("3PXO");
		Pose pose8;
		core::import_pose::pose_from_pdb( pose8, "protocols/membrane/geometry/3PXO_.pdb" );
		AddMembraneMoverOP addmem8( new AddMembraneMover( "protocols/membrane/geometry/3PXO__tr.span" ) );
		addmem8->apply(pose8);
		Vector center8(-36.1201, -7.59636, 37.6713);
		Vector normal8(-14.793, -2.47196, 0.237567);
		EmbeddingDefOP embed8( compute_structure_based_embedding( pose8 ) );
		TS_ASSERT( position_equal_within_delta( embed8->center(), center8, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed8->normal(), normal8, 0.001 ) );
	}

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
		Vector avg_normal(8.66, 8.66, 8.66);
		
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
	
	// TODO: average antiparallel embeddings!!!
	
	// split topology by jump
	void test_split_topology_by_jump() {

		TS_TRACE("Test average embeddings");
		using namespace core::conformation::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose, pose_up, pose_down;
		core::import_pose::pose_from_pdb( pose, "protocols/membrane/geometry/1AFO_.pdb" );
		
		SpanningTopology topo = SpanningTopology( "protocols/membrane/geometry/1AFO__tr.span" );
		SpanningTopology topo_up, topo_down;
		
		// call function
		split_topology_by_jump( pose, 1, topo, pose_up, pose_down, topo_up, topo_down );
		
		// test
		TS_ASSERT_EQUALS( topo_up.span(1)->start(), 15 );
		TS_ASSERT_EQUALS( topo_up.span(1)->end(), 31 );
		TS_ASSERT_EQUALS( topo_down.span(1)->start(), 15 );
		TS_ASSERT_EQUALS( topo_down.span(1)->end(), 33 );
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
