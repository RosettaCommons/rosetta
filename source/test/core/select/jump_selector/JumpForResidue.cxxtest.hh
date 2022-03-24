// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/jump_selector/JumpForResidue.cxxtest.hh
/// @brief  test suite for core::select::jump_selector::JumpForResidue
/// @author Jack Maguire, jackmaguire1444@gmail.com


// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <core/select/jump_selector/JumpForResidue.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/select/residue_selector/ChainSelector.hh>

// Utility headers

// Basic headers

// C++ headers
#include <string>

#include <core/init_util.hh> // AUTO IWYU For core_init

using namespace core::select::jump_selector;
using namespace core::select::residue_selector;

class JumpForResidueTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	core::pose::Pose
	create_n_chain_pose( int n ){
		core::pose::Pose pose;

		std::string seq = "AAA";
		for ( int chain = 2; chain <= n; ++chain ) {
			seq += "/AAA";
		}
		core::pose::make_pose_from_sequence( pose, seq, core::chemical::FA_STANDARD );
		TS_ASSERT_EQUALS( pose.num_jump(), n-1 );
		return pose;
	}


	void test_on_ala_chains() {
		JumpForResidue js;
		core::pose::Pose const three_chains =
			create_n_chain_pose( 3 );

		// ONE CHAIN AT A TIME

		{ //select chain 2
			ChainSelector c2( 2 );
			js.set_residue_selector( c2.clone() );
			auto const sele = js.apply( three_chains );

			TS_ASSERT( sele.size() == 2 );
			TS_ASSERT( sele[ 1 ] );
			TS_ASSERT( not sele[ 2 ] );

			utility::vector1< core::Size > jumps = js.selection_jumps( three_chains );
			TS_ASSERT_EQUALS( jumps.size(), 1 );
			TS_ASSERT_EQUALS( jumps.front(), 1 );
		}

		{ //select chain 3
			ChainSelector c3( 3 );
			js.set_residue_selector( c3.clone() );
			auto const sele = js.apply( three_chains );

			TS_ASSERT( sele.size() == 2 );
			TS_ASSERT( not sele[ 1 ] );
			TS_ASSERT( sele[ 2 ] );

			utility::vector1< core::Size > jumps = js.selection_jumps( three_chains );
			TS_ASSERT_EQUALS( jumps.size(), 1 );
			TS_ASSERT_EQUALS( jumps.front(), 2 );
		}

		/*{ //select chain 4
		ChainSelector c4( 4 );
		js.set_residue_selector( c4.clone() );
		auto const sele = js.apply( three_chains );

		TS_ASSERT( sele.size() == 2 );
		TS_ASSERT( not sele[ 1 ] );
		TS_ASSERT( not sele[ 2 ] );
		}*/

		/*{ //select chain 1
		ChainSelector c1( 1 );
		js.set_residue_selector( c1.clone() );
		TS_ASSERT_THROWS( js.apply( three_chains ) );
		}*/

		// MULTIPLE CHAINS

		/*{ //select chains 2 and 3 : DISALLOWED
		js.set_allow_multiple_results( false );
		ChainSelector c23( "2,3" );
		js.set_residue_selector( c23.clone() );
		TS_ASSERT_THROWS( js.apply( three_chains ) );
		}*/

		{ //select chains 2 and 3 : ALLOWED
			js.set_allow_multiple_results( true );
			ChainSelector c23( "2,3" );
			js.set_residue_selector( c23.clone() );
			auto const sele = js.apply( three_chains );

			TS_ASSERT( sele.size() == 2 );
			TS_ASSERT( sele[ 1 ] );
			TS_ASSERT( sele[ 2 ] );

			utility::vector1< core::Size > jumps = js.selection_jumps( three_chains );
			TS_ASSERT_EQUALS( jumps.size(), 2 );
			TS_ASSERT_EQUALS( jumps.front(), 1 );
			TS_ASSERT_EQUALS( jumps.back(), 2 );
		}

		// OTHER

		/*{ //nullptr
		js.set_residue_selector( nullptr );
		TS_ASSERT_THROWS( js.apply( three_chains ) );
		}*/

	}


};
