// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/jump_selector/JumpIndexSelector.cxxtest.hh
/// @brief  test suite for core::select::jump_selector::JumpIndexSelector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/jump_selector/DummySelectors.hh>

// Package headers
#include <core/select/jump_selector/JumpIndexSelector.hh>

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

class JumpIndexSelectorTests : public CxxTest::TestSuite {

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
	void test_jump_index_jump_selector() {
		JumpIndexSelectorOP index_js( new JumpIndexSelector( 5 ) );
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		fragment_pose_ft( trpcage );

		JumpSubset subset = index_js->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.num_jump() );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii == 5) );
		}
	}

	/// @brief Test JumpIndexSelector::parse_my_tag
	void test_JumpIndexSelector_parse_my_tag() {
		std::string tag_string = "<JumpIndex name=\"ind\" jump=\"5\"/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		JumpSelectorOP index_js( new JumpIndexSelector );
		try {
			index_js->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		fragment_pose_ft( trpcage );

		JumpSubset subset = index_js->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.num_jump() );
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii == 5) );
		}
	}

};
