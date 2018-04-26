// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/jump_selector/OrJumpSelector.cxxtest.hh
/// @brief  test suite for core::select::jump_selector::OrJumpSelector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/jump_selector/DummySelectors.hh>

// Package headers
#include <core/select/jump_selector/OrJumpSelector.hh>

// Project headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::select::jump_selector;

class OrJumpSelectorTests : public CxxTest::TestSuite {

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

	// @brief make sure that when we register a residue selector, we can later get it back
	void test_and_jump_selector_odd_1mod5_lots_of_jumps() {
		JumpSelectorOP odd_js( new OddJumpSelector );
		JumpSelectorOP one_mod_5_js( new XModYJumpSelector( 1, 5 ) );
		JumpSelectorOP and_js( new OrJumpSelector( odd_js, one_mod_5_js ) );
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		fragment_pose_ft( trpcage );

		JumpSubset subset = and_js->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.num_jump() );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 5 == 1 || ii % 2 == 1) );
		}
	}

	// Same as above, except, no jumps at all
	void test_and_jump_selector_odd_1mod5_no_jumps() {
		JumpSelectorOP odd_js( new OddJumpSelector );
		JumpSelectorOP one_mod_5_js( new XModYJumpSelector( 1, 5 ) );
		JumpSelectorOP and_js( new OrJumpSelector( odd_js, one_mod_5_js ) );
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		JumpSubset subset = and_js->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.num_jump() );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 5 == 1 || ii % 2 == 1) );
		}
	}

	/// @brief Test OrJumpSelector::parse_my_tag
	void test_OrJumpSelector_parse_my_tag() {
		std::string tag_string = "<And name=and_js selectors=odd,one_mod_five/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		JumpSelectorOP odd_js( new OddJumpSelector );
		JumpSelectorOP one_mod_5_js( new XModYJumpSelector( 1, 5 ) );
		dm.add( "JumpSelector", "odd", odd_js );
		dm.add( "JumpSelector", "one_mod_five", one_mod_5_js );

		JumpSelectorOP and_js( new OrJumpSelector );
		try {
			and_js->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		fragment_pose_ft( trpcage );

		JumpSubset subset = and_js->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.num_jump() );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 5 == 1 || ii % 2 == 1) );
		}
	}

	/// @brief Test that an excpetion is thrown if the OrJumpSelector is ever initialized
	/// from parse_my_tag where no JumpSelectors have been provided.
	void test_OrJumpSelector_parse_my_tag_no_provided_selectors() {
		std::string tag_string = "<And name=and_js/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		JumpSelectorOP and_js( new OrJumpSelector );
		try {
			and_js->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch (utility::excn::Exception & e ) {
			std::string expected = "No JumpSelectors given to the OrJumpSelector; OrJumpSelector requires at least one JumpSelector as input\n";
			TS_ASSERT( e.msg().find(expected) != std::string::npos );
		}

	}

	/// @brief Test than an exception is thrown if the OrJumpSelector is initialized
	/// from parse_my_tag where the JumpSelectors it requests are not in the datamap
	void test_OrJumpSelector_parse_my_tag_selectors_not_in_datamap() {
		std::string tag_string = "<And name=and_js selectors=odd,one_mod_five/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		JumpSelectorOP and_js( new OrJumpSelector );
		try {
			and_js->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch (utility::excn::Exception & e ) {
			std::string expected = "Failed to find JumpSelector named 'odd' from the Datamap from OrJumpSelector::parse_my_tag."; //\nERROR: Could not find JumpSelector and name odd in Datamap\n
			TS_ASSERT( e.msg().find(expected) != std::string::npos );
		}
	}

	void test_OrJumpSelector_parse_subtag_v1() {
		std::string tag_string = "<And name=and_js>\n\t<JumpIndex jump=\"2\" />\n\t<JumpIndex jump=\"2\" />\n</And>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		JumpSelectorOP and_js( new OrJumpSelector );
		try {
			and_js->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e ) {
			TS_ASSERT( false );
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		fragment_pose_ft( trpcage );
		JumpSubset subset = and_js->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.num_jump() );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == ( ii == 2 ) );
		}
	}

	void test_OrJumpSelector_parse_subtag_v2() {
		std::string tag_string = "<And name=and_js>\n\t<JumpIndex jump=\"2\" />\n\t<JumpIndex jump=\"3\" />\n</And>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		JumpSelectorOP and_js( new OrJumpSelector );
		try {
			and_js->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e ) {
			TS_ASSERT( false );
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		fragment_pose_ft( trpcage );
		JumpSubset subset = and_js->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.num_jump() );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == ( ii==2 || ii==3 ) );
		}
	}

	void test_OrJumpSelector_fail_parse_subtag() {
		std::string tag_string = "<And name=and_js>\n\t<JumpIndex jump=\"2\"/>\n\t<Bogus />\n</And>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		JumpSelectorOP and_js( new OrJumpSelector );
		try {
			and_js->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // parsing should fail here
		} catch (utility::excn::Exception & e ) {
			std::string err_msg =  "No JumpSelectorCreator with the name 'Bogus' has been registered with the JumpSelectorFactory";
			TS_ASSERT( e.msg().find(err_msg) != std::string::npos );
		}
	}


};
