// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pack/task/residue_selector/JumpUpstreamSelector.cxxtest.hh
/// @brief test suite for core::pack::task::residue_selector::JumpUpstreamSelector
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>

// Package headers
#include <core/pack/task/residue_selector/JumpUpstreamSelector.hh>

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

using namespace core::pack::task::residue_selector;


class JumpUpstreamSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_JumpUpstreamSelector_parse_my_tag() {
		std::string tag_string = "<JumpUp name=jump_u_rs jump=1/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP jump_u_rs( new JumpUpstreamSelector );

		try {
			jump_u_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::kinematics::FoldTree ft(trpcage.fold_tree());
		ft.new_jump(3, 7, 5);
		trpcage.fold_tree(ft);

		ResidueSubset subset = jump_u_rs->apply( trpcage );

		TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( !subset[ ii ] || ii < 6 );
		}
	}

	// make sure we fail if no selection string is provided
	void test_NeighbohoodResidueSelector_fail_no_resnums() {
		std::string tag_string = "<JumpUp name=jump_u_rs />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP jump_u_rs( new JumpUpstreamSelector );
		try {
			jump_u_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); //parsing should fail!
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			// std::cerr << "Exception (fail_no_resnums): " << e.msg();
			std::string expected_err = "Failed to access required option 'jump' from JumpUpstreamSelector::parse_my_tag.\nOption jump not found.\n";
			TS_ASSERT( e.msg() == expected_err);
		}
	}

	// everything else is done by assertion
};
