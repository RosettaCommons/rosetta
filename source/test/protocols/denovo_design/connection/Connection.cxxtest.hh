// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/denovo_design/Connection.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::connection::Connection
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
// this should be the only protocols file that includes this file
#include <test/protocols/denovo_design/test_utils.cc>

// Unit headers
#include <protocols/denovo_design/connection/Connection.hh>

// Project headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>

// Protocol headers
#include <protocols/matdes/SymDofMover.hh>
#include <protocols/moves/DsspMover.hh>

// Core headers
#include <core/kinematics/FoldTree.hh>
#include <core/io/pdb/file_data.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// C++ headers
static thread_local basic::Tracer TR( "protocols.denovo_design.Connection.cxxtest" );

// --------------- Test Class --------------- //
class ConnectionTests : public CxxTest::TestSuite {
public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init();

		// set preserve header always
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_api() {
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		using namespace protocols::denovo_design::connection;
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/test_pdbcomp.pdb" );
		StructureDataOP sd = StructureData::create_from_pose( input_pose, "UnitTest" );

		BridgeTomponents conn;
		conn.set_id( "bridge" );

		// choose some lengths
		conn.set_cut_resis( "3" );
		conn.set_lengths( "3:4" );
		for ( Connection::MotifList::const_iterator m = conn.motifs().begin(); m != conn.motifs().end(); ++m ) {
			TS_ASSERT( (m->len == core::Size( 3 )) || (m->len == core::Size( 4 )) );
		}

		// choose some motifs
		conn.set_motifs( "1LX-3EB-2LG-3EB-1LX,4LX,2LX-10HA-1LB-1LA" );
		TS_ASSERT_EQUALS( conn.motifs().size(), 3 );
		TS_ASSERT_EQUALS( conn.motifs()[ 1 ].len, 10 );
		TS_ASSERT_EQUALS( conn.motifs()[ 2 ].len, 4 );
		TS_ASSERT_EQUALS( conn.motifs()[ 3 ].len, 14 );
		TS_ASSERT_EQUALS( conn.motifs()[ 1 ].ss, "LEEELLEEEL" );
		TS_ASSERT_EQUALS( conn.motifs()[ 2 ].ss, "LLLL" );
		TS_ASSERT_EQUALS( conn.motifs()[ 3 ].ss, "LLHHHHHHHHHHLL" );
		TS_ASSERT_EQUALS( conn.motifs()[ 1 ].abego, "XBBBGGBBBX" );
		TS_ASSERT_EQUALS( conn.motifs()[ 2 ].abego, "XXXX" );
		TS_ASSERT_EQUALS( conn.motifs()[ 3 ].abego, "XXAAAAAAAAAABA" );

		// find available termini for connections
		utility::vector1< std::string > upper =
			conn.find_available_upper_termini( *sd );
		utility::vector1< std::string > lower =
			conn.find_available_lower_termini( *sd );

		TR << *sd << std::endl;
		TS_ASSERT_EQUALS( upper.size(), 5 );
		TS_ASSERT_EQUALS( lower.size(), 5 );

		TS_ASSERT_EQUALS( upper[ 1 ], "rot_comp.rot_comp2.catalytic.1" );
		TS_ASSERT_EQUALS( upper[ 2 ], "sheet1.s2" );
		TS_ASSERT_EQUALS( upper[ 3 ], "sheet1.s3" );
		TS_ASSERT_EQUALS( upper[ 4 ], "sheet1.s4" );
		TS_ASSERT_EQUALS( upper[ 5 ], "rot_comp.rot_comp2.catalytic.2" );

		TS_ASSERT_EQUALS( std::count( lower.begin(), lower.end(), "sheet1.s1" ), 1 );
		TS_ASSERT_EQUALS( std::count( lower.begin(), lower.end(), "rot_comp.rot_comp2.catalytic.2" ), 1 );
		TS_ASSERT_EQUALS( std::count( lower.begin(), lower.end(), "sheet1.s2" ), 1 );
		TS_ASSERT_EQUALS( std::count( lower.begin(), lower.end(), "sheet1.s3" ), 1 );
		TS_ASSERT_EQUALS( std::count( lower.begin(), lower.end(), "sheet1.s4" ), 1 );

		// when set set ids, these should change
		conn.set_comp1_ids( "rot_comp.rot_comp2.catalytic.1,sheet1.s2" );
		conn.set_comp2_ids( "rot_comp.rot_comp2.catalytic.1,sheet1.s2" );
		upper = conn.find_available_upper_termini( *sd );
		lower = conn.find_available_lower_termini( *sd );
		TS_ASSERT_EQUALS( upper.size(), 2 );
		TS_ASSERT_EQUALS( lower.size(), 1 );
		TS_ASSERT_EQUALS( std::count( upper.begin(), upper.end(), "rot_comp.rot_comp2.catalytic.1" ), 1 );
		TS_ASSERT_EQUALS( std::count( upper.begin(), upper.end(), "sheet1.s2" ), 1 );
		TS_ASSERT_EQUALS( std::count( lower.begin(), lower.end(), "sheet1.s2" ), 1 );

		// test user-specified chains
		conn.set_comp1_ids( "" );
		conn.set_comp2_ids( "" );
		conn.set_user_chain1( 1 );
		conn.set_user_chain2( 1 );

		// this should fail because we're connecting chain 1 to itself
		StructureDataOP sd2 = sd->clone();
		protocols::moves::MoverStatus status = conn.setup_from_random( *sd2, core::Real( 0.0000 ) );
		TS_ASSERT_EQUALS( status, protocols::moves::FAIL_DO_NOT_RETRY );

		// should succeed after enabling cyclic connections
		sd2 = sd->clone();
		conn.set_allow_cyclic( true );
		status = conn.setup_from_random( *sd2, core::Real( 0.000 ) );
		conn.set_allow_cyclic( false );
		TS_ASSERT_EQUALS( status, protocols::moves::MS_SUCCESS );
		TS_ASSERT_EQUALS( conn.build_len( *sd2 ), 10 );
		TS_ASSERT_EQUALS( conn.build_ss( *sd2 ), "LEEELLEEEL" );
		TS_ASSERT_EQUALS( conn.build_abego( *sd2 ), "XBBBGGBBBX" );
		TS_ASSERT_EQUALS( conn.lower_segment_id( *sd2 ), "rot_comp.rot_comp2.catalytic.1" );
		TS_ASSERT_EQUALS( conn.comp1_lower( *sd2 ), "sheet1.s1" );
		TS_ASSERT_EQUALS( conn.upper_segment_id( *sd2 ), "sheet1.s1" );
		TS_ASSERT_EQUALS( conn.comp2_upper( *sd2 ), "rot_comp.rot_comp2.catalytic.1" );

		// this should succeed
		sd2 = sd->clone();
		conn.set_user_chain1( 2 );
		status = conn.setup_from_random( *sd2, core::Real( 0.0000 ) );
		TS_ASSERT_EQUALS( status, protocols::moves::MS_SUCCESS );
		TS_ASSERT_EQUALS( conn.build_len( *sd2 ), 10 );
		TS_ASSERT_EQUALS( conn.build_ss( *sd2 ), "LEEELLEEEL" );
		TS_ASSERT_EQUALS( conn.build_abego( *sd2 ), "XBBBGGBBBX" );
		TS_ASSERT_EQUALS( conn.lower_segment_id( *sd2 ), "sheet1.s2" );
		TS_ASSERT_EQUALS( conn.comp1_lower( *sd2 ), "sheet1.s2" );
		TS_ASSERT_EQUALS( conn.upper_segment_id( *sd2 ), "sheet1.s1" );
		TS_ASSERT_EQUALS( conn.comp2_upper( *sd2 ), "rot_comp.rot_comp2.catalytic.1" );

		// test building loop residues
		conn.setup( *sd2 );
		TS_ASSERT( sd2 );
		TS_ASSERT( sd2->pose() );
		TS_ASSERT( sd2->has_segment( "bridge" ) );

		// look at connectivity
		TS_ASSERT_EQUALS( sd2->segment( "bridge" ).lower_segment(), "sheet1.s2" );
		TS_ASSERT_EQUALS( sd2->segment( "bridge" ).upper_segment(), "bridge_1" );
		TS_ASSERT_EQUALS( sd2->segment( "bridge_1" ).upper_segment(), "sheet1.s1" );

		// make sure residues were moved properly
		TS_ASSERT_EQUALS( sd2->segment( "sheet1.s2" ).cterm_resi(), sd2->segment( "sheet1.s2" ).stop() );
		TS_ASSERT_EQUALS( sd2->segment( "bridge" ).nterm_resi(), sd2->segment( "bridge" ).start() );
		TS_ASSERT_EQUALS( sd2->segment( "bridge" ).nterm_resi(), 1 + sd2->segment( "sheet1.s2" ).cterm_resi() );
		TS_ASSERT_EQUALS( sd2->segment( "bridge" ).cterm_resi() + 1, sd2->segment( "bridge_1" ).nterm_resi() );
		TS_ASSERT_EQUALS( sd2->segment( "bridge_1" ).cterm_resi() + 1, sd2->segment( "sheet1.s1" ).nterm_resi() );
		TS_ASSERT_EQUALS( sd2->segment( "bridge_1" ).cterm_resi(), sd2->segment( "bridge_1" ).stop() );
		TS_ASSERT_EQUALS( sd2->segment( "sheet1.s1" ).nterm_resi(), sd2->segment( "sheet1.s1" ).start() );
		TS_ASSERT_EQUALS( conn.build_left( *sd2 ), sd2->segment( "bridge" ).nterm_resi() );
		TS_ASSERT_EQUALS( conn.build_right( *sd2 ), sd2->segment( "bridge_1" ).cterm_resi() );

		// make sure residue were added: 10 residue loop
		TS_ASSERT_EQUALS( sd2->pose()->total_residue(), sd->pose()->total_residue() + conn.build_len( *sd2 ) );

		// make sure the right amount of residues are added on either side of the loop
		TS_ASSERT_EQUALS( sd2->segment( conn.loop_lower( *sd2 ) ).stop() - sd2->segment( conn.loop_lower( *sd2 ) ).start() + 1, conn.cut_resi( *sd2 ) );
		TS_ASSERT_EQUALS( sd2->segment( conn.loop_upper( *sd2 ) ).stop() - sd2->segment( conn.loop_upper( *sd2 ) ).start() + 1, conn.build_len( *sd2 ) - conn.cut_resi( *sd2 ) );

		// make sure nothing moved
		check_unwanted_movement( *sd, *sd2 );

		// test overlap
		sd2 = sd->clone();
		conn.set_overlap( 4 );
		status = conn.setup_from_random( *sd2, core::Real( 0.000 ) );
		conn.setup( *sd2 );
		TS_ASSERT_EQUALS( status, protocols::moves::MS_SUCCESS );
		TS_ASSERT_EQUALS( conn.build_left( *sd2 ), sd2->segment( conn.loop_lower( *sd2 ) ).nterm_resi() - 4 );
		TS_ASSERT_EQUALS( conn.build_right( *sd2 ), sd2->segment( conn.loop_upper( *sd2 ) ).cterm_resi() + 4 );
	}
};

