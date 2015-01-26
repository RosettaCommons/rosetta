// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/membrane/SpanningTopology.cxxtest.hh
///
/// @brief 	 Unit test for SpanningTopology class
/// @details Test object for transmembrane spans
///
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <protocols/membrane/geometry/util.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Package Headers
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/conformation/membrane/Exceptions.hh>
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

class SpanningTopologyTest : public CxxTest::TestSuite {
	
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

	// create object from PDB
	void test_from_pdb(){

		TS_TRACE("Test constructor from structure");
		using namespace core::pose;
		using namespace core::import_pose;

		// Load in pose from pdb
        PoseOP pose( new Pose() );
        pose_from_pdb( *pose, "core/conformation/membrane/1AFO_AB.pdb" );
		
		// get info from pose
		std::pair< utility::vector1< Real >, utility::vector1< Size > > pose_info( get_chain_and_z( pose ));
		zcoord_ = pose_info.first;
		chainID_ = pose_info.second;
		Real thickness = 12.5;
		
		// create object
		topo_from_pdb_ = SpanningTopologyOP( new SpanningTopology( zcoord_, chainID_, thickness ) );

		// check topology
		TS_ASSERT_EQUALS( topo_from_pdb_->nspans(), 2 );
		TS_ASSERT_EQUALS( topo_from_pdb_->span(1)->start(), 15 );
		TS_ASSERT_EQUALS( topo_from_pdb_->span(1)->end(),   31 );
		TS_ASSERT_EQUALS( topo_from_pdb_->span(2)->start(), 55 );
		TS_ASSERT_EQUALS( topo_from_pdb_->span(2)->end(),   73 );
		TS_TRACE("Finished testing constructor from structure.");
	}

////////////////////////////////////////////////////////////////////////////////
	
	// create object from spanfile
	void test_from_spanfile(){
		TS_TRACE("Test constructor from spanfile");

		// create object
		topo_from_spanfile_ = SpanningTopologyOP( new SpanningTopology( "core/conformation/membrane/1AFO_AB.span", 80 ) );
		
		// check topology
		TS_ASSERT_EQUALS( topo_from_spanfile_->nspans(), 2 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->span(1)->start(), 15 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->span(1)->end(),   31 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->span(2)->start(), 55 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->span(2)->end(),   73 );
		TS_TRACE("Finished testing constructor from spanfile.");
	}

////////////////////////////////////////////////////////////////////////////////
	
	// is residue in membrane?
	void test_res_in_membrane(){
		TS_TRACE("Test whether residue is in the membrane");
		
		// check specific residues for membrane location
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span( 8), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span(14), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span(15), 1 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span(31), 1 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span(32), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span(40), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span(41), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span(54), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span(55), 1 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span(73), 1 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span(74), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile_->in_span(80), 0 );
		TS_TRACE("Finished whether residue is in the membrane.");
	}

////////////////////////////////////////////////////////////////////////////////
	
	// does span cross the membrane?
	void test_spanning_membrane(){
		TS_TRACE("Test whether Span crosses the membrane");
		
		// check specific spans for first helix
		SpanOP span1( new Span(  1, 10 ) );
		TS_ASSERT_EQUALS( topo_from_spanfile_->spanning( zcoord_, span1 ), 0 );
		
		SpanOP span2( new Span( 10, 20 ) );
		TS_ASSERT_EQUALS( topo_from_spanfile_->spanning( zcoord_, span2 ), 0 );

		SpanOP span3( new Span( 15, 27 ) );
		TS_ASSERT_EQUALS( topo_from_spanfile_->spanning( zcoord_, span3 ), 1 );

		SpanOP span4( new Span( 27, 40 ) );
		TS_ASSERT_EQUALS( topo_from_spanfile_->spanning( zcoord_, span4 ), 0 );
		
		// check specific spans for second helix
		SpanOP span5( new Span( 40, 41 ) );
		TS_ASSERT_EQUALS( topo_from_spanfile_->spanning( zcoord_, span5 ), 1 );
		
		SpanOP span6( new Span( 41, 56 ) );
		TS_ASSERT_EQUALS( topo_from_spanfile_->spanning( zcoord_, span6 ), 0 );
		
		SpanOP span7( new Span( 57, 68 ) );
		TS_ASSERT_EQUALS( topo_from_spanfile_->spanning( zcoord_, span7 ), 1 );
		
		SpanOP span8( new Span( 66, 80 ) );
		TS_ASSERT_EQUALS( topo_from_spanfile_->spanning( zcoord_, span8 ), 0 );
		
		TS_TRACE("Finished testing whether span crosses the membrane.");
	}

////////////////////////////////////////////////////////////////////////////////
	
	// is the spanning topology valid?
	void test_is_valid(){
		TS_TRACE("Testing whether spans are valid");
		using namespace core::conformation::membrane;
		
		// test if span too short
		SpanningTopology topo1( SpanningTopology( "core/conformation/membrane/1AFO_AB_invalid_too-short.span", 80 ) );
		TS_ASSERT_EQUALS( topo1.is_valid(), true );
		TS_TRACE("...topo1 done.");
		
		// test if span too long
		SpanningTopology topo2( SpanningTopology( "core/conformation/membrane/1AFO_AB_invalid_too-long.span", 80 ) );
		TS_ASSERT_EQUALS( topo2.is_valid(), true );
		TS_TRACE("...topo2 done.");

		// test whether span start > end
		try {
			SpanningTopology topo3( SpanningTopology( "core/conformation/membrane/1AFO_AB_invalid_order.span", 80 ) );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_error_message = "SpanningTopology invalid: check your span file!";
			TS_ASSERT( expected_error_message == e.msg() );
		}
		TS_TRACE("...topo3 done.");

		// test increasing order of spans in topology
		try {
			SpanningTopology topo4( SpanningTopology( "core/conformation/membrane/1AFO_AB_invalid_spanorder.span", 80 ) );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_error_message = "SpanningTopology invalid: check your span file!";
			TS_ASSERT( expected_error_message == e.msg() );
		}
		TS_TRACE("...topo4 done.");
		
}

////////////////////////////////////////////////////////////////////////////////

private: // data

	// zcoord of pose
	utility::vector1< Real > zcoord_;
	
	// chainID of the pose
	utility::vector1< Size > chainID_;

	// topology from pdb
	SpanningTopologyOP topo_from_pdb_;
	
	// topology from spanfile
	SpanningTopologyOP topo_from_spanfile_;

};
