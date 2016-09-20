// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/denovo_design/ConsensusLoopDesignOperation.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::components::ConsensusLoopDesignOperation
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/components/Segment.hh>

// Protocol headers
#include <protocols/denovo_design/util.hh>

// Core headers

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// C++ headers

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.SegmentTests.cxxtest" );

// --------------- Test Class --------------- //
class SegmentTests : public CxxTest::TestSuite {
public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init();

		// set preserve header always for "tomponent"
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_segment_lookups() {
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		Segment s(
			"s", // id
			"LEEELLEEEL", // ss
			"XBBBGGBBBX", // abego
			false, //lower terminus included
			false ); //upper terminus included
		s.set_safe( 5 );

		// pretend this segment starts at residue 20
		s.set_pose_start( core::Size( 20 ) );

		TS_ASSERT_EQUALS( s.lower(), 20 );
		TS_ASSERT_EQUALS( s.upper(), 29 );
		TS_ASSERT_EQUALS( s.start(), 21 );
		TS_ASSERT_EQUALS( s.stop(), 28 );
		TS_ASSERT_EQUALS( s.safe(), 25 );
		TS_ASSERT_EQUALS( s.cutpoint(), 0 );
		TS_ASSERT_EQUALS( s.elem_length(), 8 );
		TS_ASSERT_EQUALS( s.length(), 10 );
		TS_ASSERT_EQUALS( s.lower_padding(), 1 );
		TS_ASSERT_EQUALS( s.upper_padding(), 1 );
		TS_ASSERT( !s.nterm_included() );
		TS_ASSERT( !s.cterm_included() );

		// test contains()
		TS_ASSERT( ! s.contains( 19 ) );
		TS_ASSERT( s.contains( 20 ) );
		TS_ASSERT( s.contains( 29 ) );
		TS_ASSERT( ! s.contains( 31 ) );

		// test delete nterm padding
		Segment t2 = s;
		t2.delete_lower_padding();
		TR << "W leading  " << s << std::endl;
		TR << "No leading " << t2 << std::endl;
		TS_ASSERT_EQUALS( t2.lower_padding(), 0 );
		TS_ASSERT_EQUALS( t2.upper_padding(), 1 );
		TS_ASSERT_EQUALS( s.lower(), t2.lower() );
		TS_ASSERT_EQUALS( t2.start() + 1, s.start() );
		TS_ASSERT_EQUALS( t2.stop() + 1, s.stop() );
		TS_ASSERT_EQUALS( t2.upper() + 1, s.upper() );
		TS_ASSERT_EQUALS( t2.length() + 1, s.length() );
		TS_ASSERT_EQUALS( t2.elem_length(), s.elem_length() );
		TS_ASSERT_EQUALS( "L" + t2.ss(), s.ss() );
		TS_ASSERT_EQUALS( "X" + t2.abego(), s.abego() );
		TS_ASSERT_EQUALS( t2.id(), s.id() );
		TS_ASSERT( !t2.contains( 29 ) );

		// test delete cterm padding
		t2.delete_upper_padding();
		TR << "No trailing " << t2 << std::endl;
		TS_ASSERT_EQUALS( t2.lower_padding(), 0 );
		TS_ASSERT_EQUALS( t2.upper_padding(), 0 );
		TS_ASSERT_EQUALS( s.lower(), t2.lower() );
		TS_ASSERT_EQUALS( t2.start() + 1, s.start() );
		TS_ASSERT_EQUALS( t2.stop() + 1, s.stop() );
		TS_ASSERT_EQUALS( t2.upper() + 2, s.upper() );
		TS_ASSERT_EQUALS( t2.length() + 2, s.length() );
		TS_ASSERT_EQUALS( t2.elem_length(), s.elem_length() );
		TS_ASSERT_EQUALS( "L" + t2.ss() + "L", s.ss() );
		TS_ASSERT_EQUALS( "X" + t2.abego() + "X", s.abego() );
		TS_ASSERT( !t2.contains( 28 ) );

		// test engulf nterm padding
		//t2 = s;
		//t2.engulf_lower_padding();
		//TS_ASSERT_EQUALS( t2.lower_padding(), 0 );
		//TS_ASSERT_EQUALS( t2.upper_padding(), 1 );
		//TS_ASSERT_EQUALS( s.lower(), t2.lower() );
		//TS_ASSERT_EQUALS( t2.start() + 1, s.start() );
		//TS_ASSERT_EQUALS( t2.stop(), s.stop() );
		//TS_ASSERT_EQUALS( t2.upper(), s.upper() );
		//TS_ASSERT_EQUALS( t2.length(), s.length() );
		//TS_ASSERT_EQUALS( t2.elem_length(), s.elem_length() + 1 );
		//TS_ASSERT_EQUALS( t2.ss(), s.ss() );
		//TS_ASSERT_EQUALS( t2.abego(), s.abego() );
		//TS_ASSERT( t2.contains( 29 ) );

		//// test engulf cterm padding
		//t2.engulf_upper_padding();
		//TS_ASSERT_EQUALS( t2.lower_padding(), 0 );
		//TS_ASSERT_EQUALS( t2.upper_padding(), 0 );
		//TS_ASSERT_EQUALS( s.lower(), t2.lower() );
		//TS_ASSERT_EQUALS( t2.start() + 1, s.start() );
		//TS_ASSERT_EQUALS( t2.stop(), s.stop() + 1 );
		//TS_ASSERT_EQUALS( t2.upper(), s.upper() );
		//TS_ASSERT_EQUALS( t2.length(), s.length() );
		//TS_ASSERT_EQUALS( t2.elem_length(), s.elem_length() + 2 );
		//TS_ASSERT_EQUALS( t2.ss(), s.ss() );
		//TS_ASSERT_EQUALS( t2.abego(), s.abego() );
		//TS_ASSERT( t2.contains( 29 ) );

		// connections
		t2 = s;
		TS_ASSERT( t2.has_free_lower_terminus() );
		TS_ASSERT( t2.has_free_upper_terminus() );
		t2.set_lower_segment( "othersegment" );
		TS_ASSERT( !t2.has_free_lower_terminus() );
		TS_ASSERT( t2.has_free_upper_terminus() );
		t2.set_upper_segment( "othersegment2" );
		TS_ASSERT( !t2.has_free_lower_terminus() );
		TS_ASSERT( !t2.has_free_upper_terminus() );
		TS_ASSERT_EQUALS( t2.upper_segment(), "othersegment2" );
		TS_ASSERT_EQUALS( t2.lower_segment(), "othersegment" );

		t2 = s;
		t2.delete_residue( 2 );
		TS_ASSERT_EQUALS( t2.length() + 1, s.length() );
		TS_ASSERT_EQUALS( t2.lower_padding(), 1 );
		TS_ASSERT_EQUALS( t2.upper_padding(), 1 );
		TS_ASSERT_EQUALS( s.lower(), t2.lower() );
		TS_ASSERT_EQUALS( t2.start(), s.start() );
		TS_ASSERT_EQUALS( t2.stop() + 1, s.stop() );
		TS_ASSERT_EQUALS( t2.upper() + 1, s.upper() );
		TS_ASSERT_EQUALS( t2.length() + 1, s.length() );
		TS_ASSERT_EQUALS( t2.elem_length() + 1, s.elem_length() );
		TS_ASSERT_EQUALS( t2.ss(), s.ss()[0] + s.ss().substr( 2,std::string::npos ) );
		TS_ASSERT_EQUALS( t2.abego().size() + 1, s.abego().size() );
	}

	void test_cutpoint_residues() {
		using namespace protocols::denovo_design::components;

		Segment s1( "s1", "LL", "XX", true, true );
		TS_ASSERT_EQUALS( s1.cutpoint(), 0 );
		TS_ASSERT_EQUALS( s1.n_residues_before_cutpoint(), 0 );
		TS_ASSERT_EQUALS( s1.n_residues_after_cutpoint(), 2 );

		s1.set_cutpoint( 1 );
		TS_ASSERT_EQUALS( s1.cutpoint(), 1 );
		TS_ASSERT_EQUALS( s1.n_residues_before_cutpoint(), 1 );
		TS_ASSERT_EQUALS( s1.n_residues_after_cutpoint(), 1 );

		s1.set_cutpoint( 2 );
		TS_ASSERT_EQUALS( s1.cutpoint(), 2 );
		TS_ASSERT_EQUALS( s1.n_residues_before_cutpoint(), 2 );
		TS_ASSERT_EQUALS( s1.n_residues_after_cutpoint(), 0 );

		Segment s2( "s2", "L", "X", true, true );
		TS_ASSERT_EQUALS( s2.cutpoint(), 0 );
		TS_ASSERT_EQUALS( s2.n_residues_before_cutpoint(), 0 );
		TS_ASSERT_EQUALS( s2.n_residues_after_cutpoint(), 1 );

		s2.set_cutpoint( 1 );
		TS_ASSERT_EQUALS( s2.cutpoint(), 1 );
		TS_ASSERT_EQUALS( s2.n_residues_before_cutpoint(), 1 );
		TS_ASSERT_EQUALS( s2.n_residues_after_cutpoint(), 0 );

		// should be invalid
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( s2.set_cutpoint( 2 ) );

		Segment s3( "s3", "", "", true, true );
		TS_ASSERT_EQUALS( s3.cutpoint(), 0 );
		TS_ASSERT_EQUALS( s3.n_residues_before_cutpoint(), 0 );
		TS_ASSERT_EQUALS( s3.n_residues_after_cutpoint(), 0 );
	}

	void test_numbering() {
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;

		Segment s1( "s1", "LHHHHHHHHL", "XAAAAAAAAO", false, false );
		TR << s1.start_local() << " " << s1.lower_padding() <<  " " << s1.lower() << " " << s1.start() << std::endl;
		TS_ASSERT_EQUALS( s1.id(), "s1" );
		TS_ASSERT_EQUALS( s1.movable_group(), 1 );
		TS_ASSERT_EQUALS( s1.start_local(), 2 );
		TS_ASSERT_EQUALS( s1.stop_local(), 9 );
		TS_ASSERT_EQUALS( s1.start(), 2 );
		TS_ASSERT_EQUALS( s1.stop(), 9 );
		TS_ASSERT_EQUALS( s1.lower(), 1 );
		TS_ASSERT_EQUALS( s1.upper(), 10 );
		TS_ASSERT_EQUALS( s1.length(), 10 );
		TS_ASSERT_EQUALS( s1.elem_length(), 8 );
		TS_ASSERT_EQUALS( s1.segment_to_pose( 1 ), s1.start() );
		TS_ASSERT_EQUALS( s1.segment_to_pose( 8 ), s1.stop() );

		// now let's say it starts at residue 20
		Segment s2( s1 );
		s2.set_pose_start( 20 );
		TS_ASSERT_EQUALS( s2.id(), s1.id() );
		TS_ASSERT_EQUALS( s2.movable_group(), 1 );
		TS_ASSERT_EQUALS( s2.start_local(), s1.start_local() );
		TS_ASSERT_EQUALS( s2.stop_local(), s1.stop_local() );
		TS_ASSERT_EQUALS( s2.lower(), 20 );
		TS_ASSERT_EQUALS( s2.upper(), 29 );
		TS_ASSERT_EQUALS( s2.start(), 21 );
		TS_ASSERT_EQUALS( s2.stop(), 28 );
		TS_ASSERT_EQUALS( s2.length(), 10 );
		TS_ASSERT_EQUALS( s2.elem_length(), 8 );
		TS_ASSERT_EQUALS( s2.segment_to_pose( 1 ), s2.start() );
		TS_ASSERT_EQUALS( s2.segment_to_pose( 8 ), s2.stop() );

		// test n-terminus included
		Segment s3( "s3", "LHHHHHHHHL", "XAAAAAAAAO", true, false );
		TS_ASSERT_EQUALS( s3.id(), "s3" );
		TS_ASSERT_EQUALS( s3.start_local(), 1 );
		TS_ASSERT_EQUALS( s3.stop_local(), 9 );
		TS_ASSERT_EQUALS( s3.start(), 1 );
		TS_ASSERT_EQUALS( s3.lower(), s3.start() );
		TS_ASSERT_EQUALS( s3.stop(), 9 );
		TS_ASSERT_EQUALS( s3.lower(), 1 );
		TS_ASSERT_EQUALS( s3.upper(), 10 );
		TS_ASSERT_EQUALS( s3.length(), 10 );
		TS_ASSERT_EQUALS( s3.elem_length(), 9 );
		TS_ASSERT_EQUALS( s3.segment_to_pose( 1 ), s3.start() );
		TS_ASSERT_EQUALS( s3.segment_to_pose( 9 ), s3.stop() );

		Segment s4( s3 );
		s4.set_pose_start( 20 );
		TS_ASSERT_EQUALS( s4.id(), s3.id() );
		TS_ASSERT_EQUALS( s4.start_local(), s4.start_local() );
		TS_ASSERT_EQUALS( s4.stop_local(), s4.stop_local() );
		TS_ASSERT_EQUALS( s4.lower(), s4.start() );
		TS_ASSERT_EQUALS( s4.lower(), 20 );
		TS_ASSERT_EQUALS( s4.upper(), 29 );
		TS_ASSERT_EQUALS( s4.start(), 20 );
		TS_ASSERT_EQUALS( s4.stop(), 28 );
		TS_ASSERT_EQUALS( s4.length(), 10 );
		TS_ASSERT_EQUALS( s4.elem_length(), 9 );
		TS_ASSERT_EQUALS( s4.segment_to_pose( 1 ), s4.start() );
		TS_ASSERT_EQUALS( s4.segment_to_pose( 9 ), s4.stop() );

	}

	void test_motif_str_parse() {
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;

		Segment s1( "s1" );
		s1.parse_motif( "10HA-1LB-1LG-1LA" );
		TS_ASSERT_EQUALS( s1.ss(), "LHHHHHHHHHHLLLL" );
		TS_ASSERT_EQUALS( s1.abego(), "XAAAAAAAAAABGAX" );
		TS_ASSERT_EQUALS( s1.length(), 15 );
		TS_ASSERT_EQUALS( s1.elem_length(), 13 );
		TR << "SS=" << s1.ss() << " Abego=" << s1.abego() << std::endl;

		Segment s2( "s2" );
		s2.parse_motif( "10:HA-2LB-3:LG-2HA" );
		TS_ASSERT_EQUALS( s2.ss(), "LHHHHHHHHHHLLLLLHHL" );
		TS_ASSERT_EQUALS( s2.abego(), "XAAAAAAAAAABBGGGAAX" );

		Segment s3( "s3" );
		s3.parse_motif( " \n" );
		TS_ASSERT_EQUALS( s3.ss(), "LL" );
		TS_ASSERT_EQUALS( s3.abego(), "XX" );
	}
};
