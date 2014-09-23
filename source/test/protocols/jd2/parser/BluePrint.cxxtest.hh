// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/parser/BluePrint.cxxtest.hh
/// @brief  test suite for BluePrint
/// @author Nobuyasu Koga

// Test headers
#include <cxxtest/TestSuite.h>

#include <core/types.hh>
#include <test/core/init_util.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("test.protocols.jd2.parser:BluePrint");

// --------------- Test Class --------------- //

class BluePrintTests : public CxxTest::TestSuite {


	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef core::kinematics::MoveMap MoveMap;
	typedef core::kinematics::MoveMapOP MoveMapOP;

	
public:


	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	// ------------- Helper Functions ------------- //

	void test_blueprint() 
	{
		using protocols::jd2::parser::BluePrint;
		BluePrint bbskel( "protocols/jd2/parser/test.bbskel" );
		
		TS_ASSERT( bbskel.total_residue() == 14 );
		TS_ASSERT( bbskel.total_residue_wolig() == 13 );
		TS_ASSERT( bbskel.secstruct() == "LEEELLLHHHHDDL" );
		TS_ASSERT( bbskel.sequence() == "VVVGGILFVVVVVX" );

		Size lig( 0 );
		for( Size ii=1; ii<=bbskel.total_residue(); ii++ ) {
			if( bbskel.extra( ii ) == "LIGAND" ) lig ++;								
		}
		TS_ASSERT( bbskel.total_residue() == ( bbskel.total_residue_wolig() + lig ) );
		
		TS_ASSERT( bbskel.strand_pairings() == "1-3.A.2;1-4.A.0;2-3.P.0" );
		TS_ASSERT( bbskel.helix_pairings() == "1-2.A;3-4.P" );
		TS_ASSERT( bbskel.hss_triplets() == "1,1-2;2,2-3");
		
		TS_ASSERT( bbskel.insertion( 1 ) == "I_like_play" );
		TS_ASSERT( bbskel.insertion( 2 ) == "the_piano" );

		TS_ASSERT( bbskel.abego( 1 ) == "X" );
		TS_ASSERT( bbskel.abego( 2 ) == "B" );
		TS_ASSERT( bbskel.abego( 4 ) == "SPZYD" );
		TS_ASSERT( bbskel.abego( 9 ) == "A" );
		
		TS_ASSERT( bbskel.buildtype( 1 ) == '.' );
		TS_ASSERT( bbskel.buildtype( 2 ) == 'W' );
		TS_ASSERT( bbskel.buildtype( 3 ) == 'F' );
		TS_ASSERT( bbskel.buildtype( 4 ) == 'P' );
		TS_ASSERT( bbskel.buildtype( 5 ) == 'R' );		
		TS_ASSERT( bbskel.buildtype( 6 ) == 'C' );		
		TS_ASSERT( bbskel.buildtype( 7 ) == 'I' );		
		
		TS_ASSERT( bbskel.resnum_map( 6 ) == 9 );
		TS_ASSERT( bbskel.resnum_map( 7 ) == 10 );
		TS_ASSERT( bbskel.resnum_map( 8 ) == 11 );
		TS_ASSERT( bbskel.resnum_map( 9 ) == 12 );
		TS_ASSERT( bbskel.resnum_map( 13 ) == 0 );
		TS_ASSERT( bbskel.resnum_map( 14 ) == 0 );
		
		MoveMapOP mm( new MoveMap );
		bbskel.set_movemap( mm );

		// .
		TS_ASSERT( mm->get_bb( 1 ) );
		TS_ASSERT( mm->get_chi( 1 ) );
		// W
		TS_ASSERT( mm->get_bb( 2 ) );
		TS_ASSERT( !mm->get_chi( 2 ) );
		// F
		TS_ASSERT( !mm->get_bb( 3 ) );
		TS_ASSERT( !mm->get_chi( 3 ) );
		// P		
		TS_ASSERT( !mm->get_bb( 4 ) );
		TS_ASSERT( mm->get_chi( 4 ) );
		// R		
		TS_ASSERT( mm->get_bb( 5 ) );
		TS_ASSERT( mm->get_chi( 5 ) );
			
		
	}// end test_blueprint

};//end class
