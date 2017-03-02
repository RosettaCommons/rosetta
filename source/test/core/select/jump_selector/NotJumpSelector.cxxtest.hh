// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/jump_selector/NotJumpSelector.cxxtest.hh
/// @brief test suite for core::select::jump_selector::NotJumpSelector
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/jump_selector/DummySelectors.hh>

// Package headers
#include <core/select/jump_selector/NotJumpSelector.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::select::jump_selector;

class NotJumpSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	// Create a silly fold tree with as many jumps as there are residues.
	void fragment_pose_ft( core::pose::Pose & pose ) {
		using namespace core::kinematics;
		FoldTree ft;
		for ( core::Size ii = 2; ii <= pose.total_residue(); ++ii ) {
			ft.add_edge( 1, ii, ii-1 );
		}
		pose.fold_tree( ft );
	}

	// @brief make sure that the not selector negates the odd selector
	void test_not_jump_selector_odd_lots_of_jumps() {
		JumpSelectorOP odd_js( new OddJumpSelector );
		JumpSelectorOP not_js( new NotJumpSelector( odd_js ) );
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		fragment_pose_ft( trpcage );
		JumpSubset subset = not_js->apply( trpcage );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 2 == 0) );
		}
	}

	void test_not_jump_selector_odd_no_jumps() {
		JumpSelectorOP odd_js( new OddJumpSelector );
		JumpSelectorOP not_js( new NotJumpSelector( odd_js ) );
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		JumpSubset subset = not_js->apply( trpcage );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 2 == 0) );
		}
	}

	/// @brief Test NotJumpSelector::parse_my_tag
	void test_NotJumpSelector_parse_my_tag() {
		std::string tag_string = "<Not name=not_js selector=odd/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		JumpSelectorOP odd_js( new OddJumpSelector );
		dm.add( "JumpSelector", "odd", odd_js );

		JumpSelectorOP not_js( new NotJumpSelector );
		try {
			not_js->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		JumpSubset subset = not_js->apply( trpcage );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 2 == 0) );
		}
	}


	/// @brief Test that an excpetion is thrown if the NotJumpSelector is ever initialized
	/// from parse_my_tag where no JumpSelector has been provided.
	void test_NotJumpSelector_parse_my_tag_no_provided_selectors() {
		std::string tag_string = "<Not name=not_js />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		JumpSelectorOP not_js( new NotJumpSelector );
		try {
			not_js->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected = "No JumpSelector given to the NotJumpSelector; NotJumpSelector requires a JumpSelector as input\n";
			TS_ASSERT_EQUALS( e.msg(), expected );
		}
	}

	/// @brief Test than an exception is thrown if the AndJumpSelector is initialized
	/// from parse_my_tag where the JumpSelectors it requests are not in the datamap
	void test_NotJumpSelector_parse_my_tag_selectors_not_in_datamap() {
		std::string tag_string = "<Not name=not_js selector=odd/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		JumpSelectorOP not_js( new NotJumpSelector );
		try {
			not_js->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected = "Failed to find JumpSelector named 'odd' from the Datamap from NotJumpSelector::parse_my_tag.\nERROR: Could not find JumpSelector and name odd in Datamap\n";
			TS_ASSERT_EQUALS( e.msg(), expected );
		}
	}


};
