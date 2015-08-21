// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/SpanningTopology.cxxtest.hh
///
/// @brief   Unit test for SpanningTopology class
/// @details Test object for transmembrane spans
///
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <protocols/membrane/util.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

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
using namespace protocols::membrane;

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
	void test_from_pdb1(){

		TS_TRACE("Test constructor from structure for 1AFO");

		SpanningTopologyOP topo_from_pdb( topology_from_pdb( "core/conformation/membrane/1AFO_AB.pdb") );

		// check topology
		TS_ASSERT_EQUALS( topo_from_pdb->nspans(), 2 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(1)->start(), 13 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(1)->end(),   33 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(2)->start(), 53 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(2)->end(),   76 );
		TS_TRACE("Finished testing constructor from structure.");
	}

	////////////////////////////////////////////////////////////////////////////////

	// create object from PDB
	void test_from_pdb2(){

		TS_TRACE("Test constructor from structure for 2LEG");

		SpanningTopologyOP topo_from_pdb( topology_from_pdb( "core/conformation/membrane/2LEG_B_tr.pdb") );

		// check topology
		TS_ASSERT_EQUALS( topo_from_pdb->nspans(), 4 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(1)->start(), 2 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(1)->end(),   18 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(2)->start(), 34 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(2)->end(),   49 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(3)->start(), 59 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(3)->end(),   75 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(4)->start(), 118 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(4)->end(),   133 );
		TS_TRACE("Finished testing constructor from structure.");
	}

	////////////////////////////////////////////////////////////////////////////////

	// create object from PDB
	void test_from_pdb3(){

		TS_TRACE("Test constructor from structure for 3MP7");

		SpanningTopologyOP topo_from_pdb( topology_from_pdb( "core/conformation/membrane/3MP7__tr.pdb") );

		// check topology
		TS_ASSERT_EQUALS( topo_from_pdb->nspans(), 11 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(1)->start(), 5 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(1)->end(),   20 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(2)->start(), 54 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(2)->end(),   62 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(3)->start(), 67 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(3)->end(),   83 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(4)->start(), 93 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(4)->end(),   114 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(5)->start(), 122 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(5)->end(),   139 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(6)->start(), 184 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(6)->end(),   201 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(7)->start(), 218 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(7)->end(),   237 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(8)->start(), 276 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(8)->end(),   297 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(9)->start(), 328 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(9)->end(),   344 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(10)->start(), 351 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(10)->end(),   361 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(11)->start(), 413 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(11)->end(),   438 );
		TS_TRACE("Finished testing constructor from structure.");
	}

	////////////////////////////////////////////////////////////////////////////////

	// create object from PDB
	void test_from_pdb4(){

		TS_TRACE("Test constructor from structure for 3UKM");

		SpanningTopologyOP topo_from_pdb( topology_from_pdb( "core/conformation/membrane/3UKM__tr.pdb") );

		// check topology
		TS_ASSERT_EQUALS( topo_from_pdb->nspans(), 12 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(1)->start(), 4 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(1)->end(),   25 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(2)->start(), 81 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(2)->end(),   92 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(3)->start(), 106 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(3)->end(),   134 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(4)->start(), 155 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(4)->end(),   175 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(5)->start(), 183 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(5)->end(),   194 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(6)->start(), 215 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(6)->end(),   239 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(7)->start(), 256 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(7)->end(),   277 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(8)->start(), 333 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(8)->end(),   344 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(9)->start(), 358 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(9)->end(),   386 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(10)->start(), 406 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(10)->end(),   427 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(11)->start(), 435 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(11)->end(),   446 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(12)->start(), 467 );
		TS_ASSERT_EQUALS( topo_from_pdb->span(12)->end(),   491 );
		TS_TRACE("Finished testing constructor from structure.");
	}

	////////////////////////////////////////////////////////////////////////////////

	// create object from spanfile
	void test_from_spanfile(){
		TS_TRACE("Test constructor from spanfile");

		// create object
		SpanningTopologyOP topo_from_spanfile( new SpanningTopology( "core/conformation/membrane/1AFO_AB.span", 80 ) );

		// check topology
		TS_ASSERT_EQUALS( topo_from_spanfile->nspans(), 2 );
		TS_ASSERT_EQUALS( topo_from_spanfile->span(1)->start(), 15 );
		TS_ASSERT_EQUALS( topo_from_spanfile->span(1)->end(),   31 );
		TS_ASSERT_EQUALS( topo_from_spanfile->span(2)->start(), 55 );
		TS_ASSERT_EQUALS( topo_from_spanfile->span(2)->end(),   73 );
		TS_TRACE("Finished testing constructor from spanfile.");
	}

	////////////////////////////////////////////////////////////////////////////////

	// is residue in membrane?
	void test_res_in_membrane(){
		TS_TRACE("Test whether residue is in the membrane");

		// create object
		SpanningTopologyOP topo_from_spanfile( new SpanningTopology( "core/conformation/membrane/1AFO_AB.span", 80 ) );

		// check specific residues for membrane location
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span( 8), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span(14), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span(15), 1 );
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span(31), 1 );
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span(32), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span(40), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span(41), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span(54), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span(55), 1 );
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span(73), 1 );
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span(74), 0 );
		TS_ASSERT_EQUALS( topo_from_spanfile->in_span(80), 0 );
		TS_TRACE("Finished whether residue is in the membrane.");
	}

	////////////////////////////////////////////////////////////////////////////////

	// does span cross the membrane?
	void test_spanning_membrane(){
		TS_TRACE("Test whether Span crosses the membrane");

		using namespace core::pose;
		using namespace core::import_pose;

		SpanningTopologyOP topo_from_pdb( topology_from_pdb( "core/conformation/membrane/1AFO_AB.pdb" ) );

		// Load in pose from pdb
		Pose pose;
		pose_from_pdb( pose, "core/conformation/membrane/1AFO_AB.pdb" );

		// get info from pose
		std::pair< utility::vector1< Real >, utility::vector1< Size > > pose_info( get_chain_and_z( pose ));
		utility::vector1< Real > zcoord = pose_info.first;

		// check specific spans for first helix
		SpanOP span1( new Span(  1, 10 ) );
		TS_ASSERT_EQUALS( topo_from_pdb->spanning( zcoord, *span1 ), 0 );

		SpanOP span2( new Span( 10, 20 ) );
		TS_ASSERT_EQUALS( topo_from_pdb->spanning( zcoord, *span2 ), 0 );

		SpanOP span3( new Span( 15, 27 ) );
		TS_ASSERT_EQUALS( topo_from_pdb->spanning( zcoord, *span3 ), 1 );

		SpanOP span4( new Span( 27, 40 ) );
		TS_ASSERT_EQUALS( topo_from_pdb->spanning( zcoord, *span4 ), 0 );

		// check specific spans for second helix
		SpanOP span5( new Span( 40, 41 ) );
		TS_ASSERT_EQUALS( topo_from_pdb->spanning( zcoord, *span5 ), 1 );

		SpanOP span6( new Span( 41, 56 ) );
		TS_ASSERT_EQUALS( topo_from_pdb->spanning( zcoord, *span6 ), 0 );

		SpanOP span7( new Span( 57, 68 ) );
		TS_ASSERT_EQUALS( topo_from_pdb->spanning( zcoord, *span7 ), 1 );

		SpanOP span8( new Span( 66, 80 ) );
		TS_ASSERT_EQUALS( topo_from_pdb->spanning( zcoord, *span8 ), 0 );

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

	// concatenate topology
	void test_concatenate_topology(){
		TS_TRACE("Testing concatenate topology");
		using namespace core::conformation::membrane;

		// create topology objects
		SpanningTopology topo1( SpanningTopology( "core/conformation/membrane/pmp22_withSS_PPM_opt.span" ) );
		SpanningTopology topo2( SpanningTopology( "core/conformation/membrane/caln_tm_finalB.span" ) );

		// concatenate them
		topo1.concatenate_topology( topo2 );

		// test topology
		TS_ASSERT_EQUALS( topo1.span(1)->start(), 4 );
		TS_ASSERT_EQUALS( topo1.span(1)->end(),   24 );
		TS_ASSERT_EQUALS( topo1.span(2)->start(), 64 );
		TS_ASSERT_EQUALS( topo1.span(2)->end(),   84 );
		TS_ASSERT_EQUALS( topo1.span(3)->start(), 98 );
		TS_ASSERT_EQUALS( topo1.span(3)->end(),   116 );
		TS_ASSERT_EQUALS( topo1.span(4)->start(), 135 );
		TS_ASSERT_EQUALS( topo1.span(4)->end(),   154 );
		TS_ASSERT_EQUALS( topo1.span(5)->start(), 167 );
		TS_ASSERT_EQUALS( topo1.span(5)->end(),   186 );

	}

	////////////////////////////////////////////////////////////////////////////////

	SpanningTopologyOP topology_from_pdb( std::string pdbfile ) {

		using namespace core::pose;
		using namespace core::import_pose;

		// Load in pose from pdb
		Pose pose;
		pose_from_pdb( pose, pdbfile );

		// get info from pose
		std::pair< utility::vector1< Real >, utility::vector1< Size > > pose_info( get_chain_and_z( pose ));
		utility::vector1< Real > zcoord = pose_info.first;
		utility::vector1< Size > chainID = pose_info.second;
		utility::vector1< char > secstruct = get_secstruct( pose );
		Real thickness = 15;

		// create object
		SpanningTopologyOP topo( new SpanningTopology( zcoord, chainID, secstruct, thickness ) );

		return topo;
	}

};
