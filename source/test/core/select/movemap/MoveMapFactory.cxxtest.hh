// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/movemap/MoveMapFactory.cxxtest.hh
/// @brief  test suite for core::select::movemap::MoveMapFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>
#include <test/core/select/jump_selector/DummySelectors.hh>

// Package headers
#include <core/select/movemap/MoveMapFactory.hh>

// Project headers
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::select::movemap;

class MoveMapFactoryTests : public CxxTest::TestSuite {

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

	void test_MoveMapFactory_parse_my_tag_v1() {
		std::string tag_string =
			"<MoveMapFactory name=\"dummy\">\n"
			"  <Backbone residue_selector=\"odd\"/>\n"
			"  <Jumps jump_selector=\"one_mod_five\"/>\n"
			"</MoveMapFactory>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		ResidueSelectorOP odd_rs( new OddResidueSelector );
		JumpSelectorOP one_mod_5_js( new XModYJumpSelector( 1, 5 ) );
		basic::datacache::DataMap dm;
		dm.add( "ResidueSelector", "odd", odd_rs );
		dm.add( "JumpSelector", "one_mod_five", one_mod_5_js );

		MoveMapFactoryOP mmf( new MoveMapFactory );
		try {
			mmf->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		fragment_pose_ft( trpcage );

		core::kinematics::MoveMapOP mm = mmf->create_movemap_from_pose( trpcage );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			TS_ASSERT( mm->get_bb(ii) == (ii%2 == 1) );
			TS_ASSERT( ! mm->get_chi(ii) );
		}
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( mm->get_jump(ii) == ( ii%5 == 1 ) );
		}

	}

	void test_MoveMapFactory_parse_my_tag_v2() {
		std::string tag_string =
			"<MoveMapFactory name=\"dummy\" jumps=\"true\">\n"
			"  <Backbone residue_selector=\"odd\"/>\n"
			"  <Jumps jump_selector=\"one_mod_five\" enable=\"false\"/>\n"
			"</MoveMapFactory>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		ResidueSelectorOP odd_rs( new OddResidueSelector );
		JumpSelectorOP one_mod_5_js( new XModYJumpSelector( 1, 5 ) );
		basic::datacache::DataMap dm;
		dm.add( "ResidueSelector", "odd", odd_rs );
		dm.add( "JumpSelector", "one_mod_five", one_mod_5_js );

		MoveMapFactoryOP mmf( new MoveMapFactory );
		try {
			mmf->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		fragment_pose_ft( trpcage );

		core::kinematics::MoveMapOP mm = mmf->create_movemap_from_pose( trpcage );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			TS_ASSERT( mm->get_bb(ii) == (ii%2 == 1) );
			TS_ASSERT( ! mm->get_chi(ii) );
		}
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( mm->get_jump(ii) == ( ii%5 != 1 ) );
		}

	}

	void test_MoveMapFactory_parse_my_tag_v3() {
		std::string tag_string =
			"<MoveMapFactory name=\"dummy\" bb=\"true\" jumps=\"true\">\n"
			"  <Backbone residue_selector=\"odd\" enable=\"false\"/>\n"
			"  <Chi residue_selector=\"odd\"/>\n"
			"  <Jumps jump_selector=\"one_mod_five\" enable=\"false\"/>\n"
			"</MoveMapFactory>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		ResidueSelectorOP odd_rs( new OddResidueSelector );
		JumpSelectorOP one_mod_5_js( new XModYJumpSelector( 1, 5 ) );
		basic::datacache::DataMap dm;
		dm.add( "ResidueSelector", "odd", odd_rs );
		dm.add( "JumpSelector", "one_mod_five", one_mod_5_js );

		MoveMapFactoryOP mmf( new MoveMapFactory );
		try {
			mmf->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		fragment_pose_ft( trpcage );

		core::kinematics::MoveMapOP mm = mmf->create_movemap_from_pose( trpcage );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			TS_ASSERT( mm->get_bb(ii) != (ii%2 == 1) );
			TS_ASSERT( mm->get_chi(ii) == (ii%2 == 1 ) );
		}
		for ( core::Size ii = 1; ii <= trpcage.num_jump(); ++ii ) {
			TS_ASSERT( mm->get_jump(ii) == ( ii%5 != 1 ) );
		}

	}

};
