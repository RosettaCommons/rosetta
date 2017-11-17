// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/jump_selector/InterchainJumpSelector.cxxtest.hh
/// @brief  test suite for core::select::jump_selector::InterchainJumpSelector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/jump_selector/DummySelectors.hh>

// Package headers
#include <core/select/jump_selector/InterchainJumpSelector.hh>

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

class InterchainJumpSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void create_two_chain_pose( core::pose::Pose & pose ) {
		core::pose::Pose second_pose( pose );
		pose.append_pose_by_jump( second_pose, 1 );
		core::kinematics::FoldTree ft;
		ft.add_edge( 1, second_pose.size() + 1, 1 );
		for ( core::Size ii = 2; ii <= second_pose.size(); ++ii ) {
			ft.add_edge( 1, ii, ii );
		}
		for ( core::Size ii = 2; ii <= second_pose.size(); ++ii ) {
			ft.add_edge( second_pose.size()+1, second_pose.size() + ii, second_pose.size()+ii-1 );
		}
		pose.fold_tree( ft );
	}

	// @brief make sure that when we register a residue selector, we can later get it back
	void test_interchain_jump_selector() {
		JumpSelectorOP interchain_js( new InterchainJumpSelector );
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		create_two_chain_pose( trpcage );

		JumpSubset subset = interchain_js->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.num_jump() );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii == 1 ) );
		}
	}

	/// @brief Test InterchainJumpSelector::parse_my_tag
	void test_InterchainJumpSelector_parse_my_tag() {
		std::string tag_string = "<Interchain/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		JumpSelectorOP interchain_js( new InterchainJumpSelector );
		try {
			interchain_js->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		create_two_chain_pose( trpcage );
		JumpSubset subset = interchain_js->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.num_jump() );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii == 1) );
		}
	}

};
