// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

static thread_local basic::Tracer TR( "protocols.denovo_design.SegmentTests.cxxtest" );

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
				10, //pose len
				5,  //safe res
				0,  //cutpoint
				1,  //movable group
				false, //loop
				false, //lower terminus included
				false, //upper terminus included
				"", //lower segment
				"", //upper segment
				"LEEELLEEEL", // ss
				protocols::denovo_design::abego_vector( "XBBBGGBBBX" ) );

		// pretend this segment starts at residue 20
		s.set_pose_start( core::Size( 20 ) );
		TS_ASSERT_EQUALS( s.nterm_resi(), 20 );
		TS_ASSERT_EQUALS( s.cterm_resi(), 29 );
		TS_ASSERT_EQUALS( s.safe(), 24 );
		TS_ASSERT_EQUALS( s.cutpoint(), 0 );
		TS_ASSERT_EQUALS( s.elem_length(), 8 );
		TS_ASSERT_EQUALS( s.length(), 10 );
		TS_ASSERT_EQUALS( s.nterm_pad(), 1 );
		TS_ASSERT_EQUALS( s.cterm_pad(), 1 );
		TS_ASSERT( !s.nterm_included() );
		TS_ASSERT( !s.cterm_included() );

		// test contains()
		TS_ASSERT( ! s.contains( 19 ) );
		TS_ASSERT( s.contains( 20 ) );
		TS_ASSERT( s.contains( 29 ) );
		TS_ASSERT( ! s.contains( 31 ) );

		// test delete nterm padding
		Segment t2 = s;
		t2.delete_leading_residues();
		TS_ASSERT_EQUALS( t2.nterm_pad(), 0 );
		TS_ASSERT_EQUALS( t2.cterm_pad(), 1 );
		TS_ASSERT_EQUALS( s.nterm_resi(), t2.nterm_resi() );
		TS_ASSERT_EQUALS( t2.start() + 1, s.start() );
		TS_ASSERT_EQUALS( t2.stop() + 1, s.stop() );
		TS_ASSERT_EQUALS( t2.cterm_resi() + 1, s.cterm_resi() );
		TS_ASSERT_EQUALS( t2.length() + 1, s.length() );
		TS_ASSERT_EQUALS( t2.elem_length(), s.elem_length() );
		TS_ASSERT_EQUALS( "L" + t2.ss(), s.ss() );
		TS_ASSERT_EQUALS( "X" + abego_str( t2.abego() ), abego_str( s.abego() ) );
		TS_ASSERT( !t2.contains( 29 ) );

		// test delete cterm padding
		t2.delete_trailing_residues();
		TS_ASSERT_EQUALS( t2.nterm_pad(), 0 );
		TS_ASSERT_EQUALS( t2.cterm_pad(), 0 );
		TS_ASSERT_EQUALS( s.nterm_resi(), t2.nterm_resi() );
		TS_ASSERT_EQUALS( t2.start() + 1, s.start() );
		TS_ASSERT_EQUALS( t2.stop() + 1, s.stop() );
		TS_ASSERT_EQUALS( t2.cterm_resi() + 2, s.cterm_resi() );
		TS_ASSERT_EQUALS( t2.length() + 2, s.length() );
		TS_ASSERT_EQUALS( t2.elem_length(), s.elem_length() );
		TS_ASSERT_EQUALS( "L" + t2.ss() + "L", s.ss() );
		TS_ASSERT_EQUALS( "X" + abego_str( t2.abego() ) + "X", abego_str( s.abego() ) );
		TS_ASSERT( !t2.contains( 28 ) );

		// test engulf nterm padding
		t2 = s;
		t2.engulf_leading_residues();
		TS_ASSERT_EQUALS( t2.nterm_pad(), 0 );
		TS_ASSERT_EQUALS( t2.cterm_pad(), 1 );
		TS_ASSERT_EQUALS( s.nterm_resi(), t2.nterm_resi() );
		TS_ASSERT_EQUALS( t2.start() + 1, s.start() );
		TS_ASSERT_EQUALS( t2.stop(), s.stop() );
		TS_ASSERT_EQUALS( t2.cterm_resi(), s.cterm_resi() );
		TS_ASSERT_EQUALS( t2.length(), s.length() );
		TS_ASSERT_EQUALS( t2.elem_length(), s.elem_length() + 1 );
		TS_ASSERT_EQUALS( t2.ss(), s.ss() );
		TS_ASSERT_EQUALS( abego_str( t2.abego() ), abego_str( s.abego() ) );
		TS_ASSERT( t2.contains( 29 ) );

		// test engulf cterm padding
		t2.engulf_trailing_residues();
		TS_ASSERT_EQUALS( t2.nterm_pad(), 0 );
		TS_ASSERT_EQUALS( t2.cterm_pad(), 0 );
		TS_ASSERT_EQUALS( s.nterm_resi(), t2.nterm_resi() );
		TS_ASSERT_EQUALS( t2.start() + 1, s.start() );
		TS_ASSERT_EQUALS( t2.stop(), s.stop() + 1 );
		TS_ASSERT_EQUALS( t2.cterm_resi(), s.cterm_resi() );
		TS_ASSERT_EQUALS( t2.length(), s.length() );
		TS_ASSERT_EQUALS( t2.elem_length(), s.elem_length() + 2 );
		TS_ASSERT_EQUALS( t2.ss(), s.ss() );
		TS_ASSERT_EQUALS( abego_str( t2.abego() ), abego_str( s.abego() ) );
		TS_ASSERT( t2.contains( 29 ) );

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
		t2.delete_residues( 2, 2 );
		TS_ASSERT_EQUALS( t2.length() + 1, s.length() );
		TS_ASSERT_EQUALS( t2.nterm_pad(), 1 );
		TS_ASSERT_EQUALS( t2.cterm_pad(), 1 );
		TS_ASSERT_EQUALS( s.nterm_resi(), t2.nterm_resi() );
		TS_ASSERT_EQUALS( t2.start(), s.start() );
		TS_ASSERT_EQUALS( t2.stop() + 1, s.stop() );
		TS_ASSERT_EQUALS( t2.cterm_resi() + 1, s.cterm_resi() );
		TS_ASSERT_EQUALS( t2.length() + 1, s.length() );
		TS_ASSERT_EQUALS( t2.elem_length() + 1, s.elem_length() );
		TS_ASSERT_EQUALS( t2.ss(), s.ss()[0] + s.ss().substr( 2,std::string::npos ) );
		TS_ASSERT_EQUALS( t2.abego().size() + 1, s.abego().size() );
	}
};
