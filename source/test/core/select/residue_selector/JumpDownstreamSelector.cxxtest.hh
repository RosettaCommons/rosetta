// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/JumpDownstreamSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::JumpDownstreamSelector
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>

// Package headers
#include <core/select/residue_selector/JumpDownstreamSelector.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>

#include <core/init_util.hh> // AUTO IWYU For core_init
#include <iostream> // AUTO IWYU For operator<<, basic_ostream, endl, cerr

using namespace core::select::residue_selector;
static basic::Tracer TR( "tests.core.select.residue_selector.JumpDownstreamSelector" );

class JumpDownstreamSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_JumpDownstreamSelector_parse_my_tag() {
		std::string tag_string = "<JumpDown name=jump_d_rs jump=1/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		JumpDownstreamSelector jump_d_rs;

		try {
			jump_d_rs.parse_my_tag( tag, dm );

		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::kinematics::FoldTree ft(trpcage.fold_tree());
		ft.new_jump(3, 7, 5);
		trpcage.fold_tree(ft);

		ResidueSubset subset = jump_d_rs.apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( !subset[ ii ] || ii >= 6 );
		}
	}

	// make sure we fail if no selection string is provided
	void test_NeighbohoodResidueSelector_fail_no_resnums() {
		std::string tag_string = "<JumpDown name=jump_d_rs />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		JumpDownstreamSelector jump_d_rs;
		TS_ASSERT_THROWS_ANYTHING( jump_d_rs.parse_my_tag( tag, dm ) );
	}

	///@author Jack Maguire
	void test_JumpDownstreamSelector_simple_two_chains() {
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "AAA/GGG", "fa_standard" );
		JumpDownstreamSelector const selector( 1 );
		{
			utility::vector1< bool > const selection = selector.apply( pose );
			TS_ASSERT( !selection[1] );
			TS_ASSERT( !selection[2] );
			TS_ASSERT( !selection[3] );
			TS_ASSERT( selection[4] );
			TS_ASSERT( selection[5] );
			TS_ASSERT( selection[6] );
		}

		//Re-root
		core::kinematics::FoldTree ft = pose.fold_tree();
		bool const reorder_success = ft.reorder( 4 );
		TS_ASSERT( reorder_success );
		pose.fold_tree( ft );
		for ( core::kinematics::Edge const & e : pose.fold_tree() ) {
			TR << e << std::endl;
		}
		{
			utility::vector1< bool > const selection = selector.apply( pose );
			TS_ASSERT( selection[1] );
			TS_ASSERT( selection[2] );
			TS_ASSERT( selection[3] );
			TS_ASSERT( !selection[4] );
			TS_ASSERT( !selection[5] );
			TS_ASSERT( !selection[6] );
		}
	}

	// everything else is done by assertion
};
