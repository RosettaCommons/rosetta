// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/SetupMetalsMoverTests.cxxtest.hh
/// @brief  Test suite for SetupMetalsMover, which applies -auto_setup_metals in mover form
/// @author Sharon Guffy

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/simple_moves/SetupMetalsMover.hh>
// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

class SetupMetalsMoverTests: public CxxTest::TestSuite {

private:
	core::pose::PoseOP ref_pose_;

public:

	// Shared initialization goes here.
	void setUp() {
		//Extra options are provided so we can test grabbing default values from command line
		protocols_init_with_additional_options( "-in:metals_distance_constraint_multiplier 2.0 -in:metals_angle_constraint_multiplier 0.5" );
		ref_pose_ = core::import_pose::pose_from_file( "core/util/2c9v_stripped.pdb" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_default_options(){
		core::pose::PoseOP pose_ = ref_pose_->clone();
		//This should behave exactly like the -in:auto_setup_metals flag
		protocols::simple_moves::SetupMetalsMover new_mover;
		//Check all default values
		TS_ASSERT_EQUALS( new_mover.get_remove_hydrogens(), true );
		TS_ASSERT_EQUALS( new_mover.get_metal_resnums_string(), "" );
		TS_ASSERT_EQUALS( new_mover.get_metal_selector(), nullptr );

		//The following should have been set from command line values in the constructor
		TS_ASSERT_DELTA( new_mover.get_metals_detection_LJ_multiplier(), 1.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_distance_constraint_multiplier(), 2.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_angle_constraint_multiplier(), 0.5, 1e-6 );

		//Apply the mover
		new_mover.apply( *pose_ );
		//Check that all bonds were created
		TS_ASSERT( pose_->residue(46).is_bonded(154) );
		TS_ASSERT( pose_->residue(48).is_bonded(154) );
		TS_ASSERT( pose_->residue(63).is_bonded(154) );
		TS_ASSERT( pose_->residue(120).is_bonded(154) );

		TS_ASSERT( pose_->residue(63).is_bonded(155) );
		TS_ASSERT( pose_->residue(71).is_bonded(155) );
		TS_ASSERT( pose_->residue(80).is_bonded(155) );
		TS_ASSERT( pose_->residue(83).is_bonded(155) );

		TS_ASSERT( pose_->residue(201).is_bonded(309) );
		TS_ASSERT( pose_->residue(203).is_bonded(309) );
		TS_ASSERT( pose_->residue(218).is_bonded(309) );
		TS_ASSERT( pose_->residue(275).is_bonded(309) );

		TS_ASSERT( pose_->residue(218).is_bonded(310) );
		TS_ASSERT( pose_->residue(226).is_bonded(310) );
		TS_ASSERT( pose_->residue(235).is_bonded(310) );
		TS_ASSERT( pose_->residue(238).is_bonded(310) );
		//Check that the pose has the proper number of constraints (should have 18 constraints per metal)
		TS_ASSERT_EQUALS( pose_->constraint_set()->get_all_constraints().size(), 72 );
	}

	void test_parse_tag_all_metals(){
		core::pose::PoseOP pose_ = ref_pose_->clone();
		protocols::simple_moves::SetupMetalsMover new_mover;
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss;
		ss << "<SetupMetalsMover name=\"metal\" remove_hydrogens=\"false\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" />" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		TS_TRACE( "Tag with no selector or resnums" );
		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		protocols::moves::Movers_map mm;
		TS_ASSERT_THROWS_NOTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		//Check all default values
		TS_ASSERT_EQUALS( new_mover.get_remove_hydrogens(), false );
		TS_ASSERT_EQUALS( new_mover.get_metal_resnums_string(), "" );
		TS_ASSERT_EQUALS( new_mover.get_metal_selector(), nullptr );

		//The following should have been set from command line values in the constructor
		TS_ASSERT_DELTA( new_mover.get_metals_detection_LJ_multiplier(), 1.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_distance_constraint_multiplier(), 2.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_angle_constraint_multiplier(), 5.0, 1e-6 );
		//Apply the mover
		new_mover.apply( *pose_ );
		//Check that all bonds were created
		TS_ASSERT( pose_->residue(46).is_bonded(154) );
		TS_ASSERT( pose_->residue(48).is_bonded(154) );
		TS_ASSERT( pose_->residue(63).is_bonded(154) );
		TS_ASSERT( pose_->residue(120).is_bonded(154) );

		TS_ASSERT( pose_->residue(63).is_bonded(155) );
		TS_ASSERT( pose_->residue(71).is_bonded(155) );
		TS_ASSERT( pose_->residue(80).is_bonded(155) );
		TS_ASSERT( pose_->residue(83).is_bonded(155) );

		TS_ASSERT( pose_->residue(201).is_bonded(309) );
		TS_ASSERT( pose_->residue(203).is_bonded(309) );
		TS_ASSERT( pose_->residue(218).is_bonded(309) );
		TS_ASSERT( pose_->residue(275).is_bonded(309) );

		TS_ASSERT( pose_->residue(218).is_bonded(310) );
		TS_ASSERT( pose_->residue(226).is_bonded(310) );
		TS_ASSERT( pose_->residue(235).is_bonded(310) );
		TS_ASSERT( pose_->residue(238).is_bonded(310) );
		//Check that the pose has the proper number of constraints (should have 16 atom pair constraints and 16 angle constraints)
		//It also has coordinate constraints on virtual atoms--for each metal, there will be 1 + (4 - binder#) of these (4 on the 1st, 3 on the second, 2 on the 3rd, 1 on the 4th for a total of 10 per metal)
		//This comes out to 18 constraints per metal (4 distance, 4 angle, 10 coordinate) = 72 total
		TS_ASSERT_EQUALS( pose_->constraint_set()->get_all_constraints().size(), 72 );

	}
	void test_parse_invalid_tag(){
		protocols::simple_moves::SetupMetalsMover new_mover;
		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		protocols::moves::Movers_map mm;
		core::pose::PoseOP pose_ = ref_pose_->clone();

		/*
		//Try to provide both a selector and a resnum string
		std::stringstream ss_string_and_selector;
		ss_string_and_selector << "<SetupMetalsMover name=\"metal\" metal_resnums=\"154,155\" remove_hydrogens=\"false\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\"><Index resnums=\"154,155\" /></SetupMetalsMover>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss_string_and_selector );
		TS_TRACE( "Tag with embedded selector and resnums" );
		TS_ASSERT_THROWS_ANYTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		*/
		/*
		//Try to provide two embedded selectors
		std::stringstream ss_two_embedded;
		ss_two_embedded << "<SetupMetalsMover name=\"metal\" remove_hydrogens=\"false\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\"><Index resnums=\"154,155\" /><Chain chains=\"A\" /></SetupMetalsMover>" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_two_embedded );
		TS_TRACE( "Tag with two embedded selectors" );
		TS_ASSERT_THROWS_ANYTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		*/
		//Provide a selector that isn't in the data map
		std::stringstream ss_bad_selector;
		ss_bad_selector << "<SetupMetalsMover name=\"metal\" remove_hydrogens=\"false\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_residue_selector=\"dummy\" />" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss_bad_selector );
		TS_TRACE( "Tag with undefined residue selector" );
		TS_ASSERT_THROWS_ANYTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );

		//Add selector to the data map
		core::select::residue_selector::ResidueIndexSelectorOP dummy( new core::select::residue_selector::ResidueIndexSelector );
		dummy->set_index( "154,155" );
		dm.add( "ResidueSelector", "dummy", dummy );
		//Confirm that the above tag would work once the selector was added
		std::stringstream ss_good_selector;
		ss_good_selector << "<SetupMetalsMover name=\"metal\" remove_hydrogens=\"false\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_residue_selector=\"dummy\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_good_selector );
		TS_TRACE( "Tag with defined residue selector" );
		TS_ASSERT_THROWS_NOTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		/*
		//Try to provide a selector string and an embedded selector
		std::stringstream ss_selector_two_ways;
		ss_selector_two_ways << "<SetupMetalsMover name=\"metal\" remove_hydrogens=\"false\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_residue_selector=\"dummy\" ><Index resnums=\"154,155\" /></SetupMetalsMover>" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_selector_two_ways );
		TS_TRACE( "Tag with selector two ways" );
		TS_ASSERT_THROWS_ANYTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		*/
		//Try to provide a selector string and a resnum string
		std::stringstream ss_select_string_and_resnum_string;
		ss_select_string_and_resnum_string << "<SetupMetalsMover name=\"metal\" metal_resnums=\"154,155\" remove_hydrogens=\"false\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_residue_selector=\"dummy\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_select_string_and_resnum_string );
		TS_TRACE( "Tag with selector string and resnum string" );
		TS_ASSERT_THROWS_ANYTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
	}

	void test_parse_tag_named_selector(){
		//This should only set up metal ions specified by the selector
		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		protocols::moves::Movers_map mm;
		core::pose::PoseOP pose_ = ref_pose_->clone();
		//Add the selectors to the data map
		//Names will be all, chain, zinc, and none
		core::select::residue_selector::ResidueIndexSelectorOP all( new core::select::residue_selector::ResidueIndexSelector );
		all->set_index( "154,155,309,310" );
		dm.add( "ResidueSelector", "all", all );
		core::select::residue_selector::ResidueIndexSelectorOP none( new core::select::residue_selector::ResidueIndexSelector );
		none->set_index( "1,2,3" );
		dm.add( "ResidueSelector", "none", none );
		core::select::residue_selector::ResidueNameSelectorOP zinc( new core::select::residue_selector::ResidueNameSelector );
		zinc->set_residue_name3( " ZN" );
		dm.add( "ResidueSelector", "zinc", zinc );
		core::select::residue_selector::ChainSelectorOP chain( new core::select::residue_selector::ChainSelector );
		utility::vector1< std::string > chain_ids;
		chain_ids.push_back( "A" );
		chain->set_chain_strings( chain_ids );
		dm.add( "ResidueSelector", "chain", chain );

		//Select all metals
		protocols::simple_moves::SetupMetalsMover new_mover;
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_all;
		ss_all << "<SetupMetalsMover name=\"metal\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_residue_selector=\"all\" />" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss_all );
		TS_TRACE( "Tag with selector to select all metals" );
		TS_ASSERT_THROWS_NOTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		//Check all default values
		TS_ASSERT_EQUALS( new_mover.get_remove_hydrogens(), true );
		TS_ASSERT_EQUALS( new_mover.get_metal_resnums_string(), "" );

		//The following should have been set from command line values in the constructor
		TS_ASSERT_DELTA( new_mover.get_metals_detection_LJ_multiplier(), 1.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_distance_constraint_multiplier(), 2.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_angle_constraint_multiplier(), 5.0, 1e-6 );
		//Apply the mover
		new_mover.apply( *pose_ );
		//Check that all bonds were created
		TS_ASSERT( pose_->residue(46).is_bonded(154) );
		TS_ASSERT( pose_->residue(48).is_bonded(154) );
		TS_ASSERT( pose_->residue(63).is_bonded(154) );
		TS_ASSERT( pose_->residue(120).is_bonded(154) );

		TS_ASSERT( pose_->residue(63).is_bonded(155) );
		TS_ASSERT( pose_->residue(71).is_bonded(155) );
		TS_ASSERT( pose_->residue(80).is_bonded(155) );
		TS_ASSERT( pose_->residue(83).is_bonded(155) );

		TS_ASSERT( pose_->residue(201).is_bonded(309) );
		TS_ASSERT( pose_->residue(203).is_bonded(309) );
		TS_ASSERT( pose_->residue(218).is_bonded(309) );
		TS_ASSERT( pose_->residue(275).is_bonded(309) );

		TS_ASSERT( pose_->residue(218).is_bonded(310) );
		TS_ASSERT( pose_->residue(226).is_bonded(310) );
		TS_ASSERT( pose_->residue(235).is_bonded(310) );
		TS_ASSERT( pose_->residue(238).is_bonded(310) );
		//Check that the pose has the proper number of constraints (should have 18 per setup metal)
		TS_ASSERT_EQUALS( pose_->constraint_set()->get_all_constraints().size(), 72 );

		//Select metals only on chain A
		pose_ = ref_pose_->clone();
		new_mover = protocols::simple_moves::SetupMetalsMover();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_a;
		ss_a << "<SetupMetalsMover name=\"metal\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_residue_selector=\"chain\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_a );
		TS_TRACE( "Tag with selector to select chain A" );
		TS_ASSERT_THROWS_NOTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		//Check all default values
		TS_ASSERT_EQUALS( new_mover.get_remove_hydrogens(), true );
		TS_ASSERT_EQUALS( new_mover.get_metal_resnums_string(), "" );

		//The following should have been set from command line values in the constructor
		TS_ASSERT_DELTA( new_mover.get_metals_detection_LJ_multiplier(), 1.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_distance_constraint_multiplier(), 2.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_angle_constraint_multiplier(), 5.0, 1e-6 );
		//Apply the mover
		new_mover.apply( *pose_ );
		//Check that all bonds were created
		TS_ASSERT( pose_->residue(46).is_bonded(154) );
		TS_ASSERT( pose_->residue(48).is_bonded(154) );
		TS_ASSERT( pose_->residue(63).is_bonded(154) );
		TS_ASSERT( pose_->residue(120).is_bonded(154) );

		TS_ASSERT( pose_->residue(63).is_bonded(155) );
		TS_ASSERT( pose_->residue(71).is_bonded(155) );
		TS_ASSERT( pose_->residue(80).is_bonded(155) );
		TS_ASSERT( pose_->residue(83).is_bonded(155) );

		TS_ASSERT( !pose_->residue(201).is_bonded(309) );
		TS_ASSERT( !pose_->residue(203).is_bonded(309) );
		TS_ASSERT( !pose_->residue(218).is_bonded(309) );
		TS_ASSERT( !pose_->residue(275).is_bonded(309) );
		TS_ASSERT( !pose_->residue(218).is_bonded(310) );
		TS_ASSERT( !pose_->residue(226).is_bonded(310) );
		TS_ASSERT( !pose_->residue(235).is_bonded(310) );
		TS_ASSERT( !pose_->residue(238).is_bonded(310) );
		//Check that the pose has the proper number of constraints (should have 18 per setup metal)
		TS_ASSERT_EQUALS( pose_->constraint_set()->get_all_constraints().size(), 36 );

		//Select only zinc ions
		pose_ = ref_pose_->clone();
		new_mover = protocols::simple_moves::SetupMetalsMover();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_zn;
		ss_zn << "<SetupMetalsMover name=\"metal\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_residue_selector=\"zinc\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_zn );
		TS_TRACE( "Tag with selector to select zinc ions" );
		TS_ASSERT_THROWS_NOTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		//Check all default values
		TS_ASSERT_EQUALS( new_mover.get_remove_hydrogens(), true );
		TS_ASSERT_EQUALS( new_mover.get_metal_resnums_string(), "" );

		//The following should have been set from command line values in the constructor
		TS_ASSERT_DELTA( new_mover.get_metals_detection_LJ_multiplier(), 1.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_distance_constraint_multiplier(), 2.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_angle_constraint_multiplier(), 5.0, 1e-6 );
		//Apply the mover
		new_mover.apply( *pose_ );
		//Check that all bonds were created
		TS_ASSERT( !pose_->residue(46).is_bonded(154) );
		TS_ASSERT( !pose_->residue(48).is_bonded(154) );
		TS_ASSERT( !pose_->residue(63).is_bonded(154) );
		TS_ASSERT( !pose_->residue(120).is_bonded(154) );

		TS_ASSERT( pose_->residue(63).is_bonded(155) );
		TS_ASSERT( pose_->residue(71).is_bonded(155) );
		TS_ASSERT( pose_->residue(80).is_bonded(155) );
		TS_ASSERT( pose_->residue(83).is_bonded(155) );

		TS_ASSERT( !pose_->residue(201).is_bonded(309) );
		TS_ASSERT( !pose_->residue(203).is_bonded(309) );
		TS_ASSERT( !pose_->residue(218).is_bonded(309) );
		TS_ASSERT( !pose_->residue(275).is_bonded(309) );

		TS_ASSERT( pose_->residue(218).is_bonded(310) );
		TS_ASSERT( pose_->residue(226).is_bonded(310) );
		TS_ASSERT( pose_->residue(235).is_bonded(310) );
		TS_ASSERT( pose_->residue(238).is_bonded(310) );
		//Check that the pose has the proper number of constraints (should have 18 per setup metal)
		TS_ASSERT_EQUALS( pose_->constraint_set()->get_all_constraints().size(), 36 );

		//Select none of the metals (none of the metals are in the selection)
		pose_ = ref_pose_->clone();
		new_mover = protocols::simple_moves::SetupMetalsMover();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_none;
		ss_none << "<SetupMetalsMover name=\"metal\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_residue_selector=\"none\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_none );
		TS_TRACE( "Tag with selector that should not include any of the metal ions" );
		TS_ASSERT_THROWS_NOTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		//Check all default values
		TS_ASSERT_EQUALS( new_mover.get_remove_hydrogens(), true );
		TS_ASSERT_EQUALS( new_mover.get_metal_resnums_string(), "" );

		//The following should have been set from command line values in the constructor
		TS_ASSERT_DELTA( new_mover.get_metals_detection_LJ_multiplier(), 1.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_distance_constraint_multiplier(), 2.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_angle_constraint_multiplier(), 5.0, 1e-6 );
		//Apply the mover
		new_mover.apply( *pose_ );
		//Check that all bonds were created
		TS_ASSERT( !pose_->residue(46).is_bonded(154) );
		TS_ASSERT( !pose_->residue(48).is_bonded(154) );
		TS_ASSERT( !pose_->residue(63).is_bonded(154) );
		TS_ASSERT( !pose_->residue(120).is_bonded(154) );

		TS_ASSERT( !pose_->residue(63).is_bonded(155) );
		TS_ASSERT( !pose_->residue(71).is_bonded(155) );
		TS_ASSERT( !pose_->residue(80).is_bonded(155) );
		TS_ASSERT( !pose_->residue(83).is_bonded(155) );

		TS_ASSERT( !pose_->residue(201).is_bonded(309) );
		TS_ASSERT( !pose_->residue(203).is_bonded(309) );
		TS_ASSERT( !pose_->residue(218).is_bonded(309) );
		TS_ASSERT( !pose_->residue(275).is_bonded(309) );

		TS_ASSERT( !pose_->residue(218).is_bonded(310) );
		TS_ASSERT( !pose_->residue(226).is_bonded(310) );
		TS_ASSERT( !pose_->residue(235).is_bonded(310) );
		TS_ASSERT( !pose_->residue(238).is_bonded(310) );
		//Check that the pose has the proper number of constraints (should have 18 per setup metal)
		TS_ASSERT_EQUALS( pose_->constraint_set()->get_all_constraints().size(), 0 );
	}
	void test_parse_tag_resnum_list(){
		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		protocols::moves::Movers_map mm;
		core::pose::PoseOP pose_ = ref_pose_->clone();
		//Select all metals
		protocols::simple_moves::SetupMetalsMover new_mover;
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_all;
		ss_all << "<SetupMetalsMover name=\"metal\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_resnums=\"154,155,309,310\" />" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss_all );
		TS_TRACE( "Tag with resnums to select all metals" );
		TS_ASSERT_THROWS_NOTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		//Check all default values
		TS_ASSERT_EQUALS( new_mover.get_remove_hydrogens(), true );
		TS_ASSERT_EQUALS( new_mover.get_metal_selector(), nullptr );

		//The following should have been set from command line values in the constructor
		TS_ASSERT_DELTA( new_mover.get_metals_detection_LJ_multiplier(), 1.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_distance_constraint_multiplier(), 2.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_angle_constraint_multiplier(), 5.0, 1e-6 );
		//Apply the mover
		new_mover.apply( *pose_ );
		//Check that all bonds were created
		TS_ASSERT( pose_->residue(46).is_bonded(154) );
		TS_ASSERT( pose_->residue(48).is_bonded(154) );
		TS_ASSERT( pose_->residue(63).is_bonded(154) );
		TS_ASSERT( pose_->residue(120).is_bonded(154) );

		TS_ASSERT( pose_->residue(63).is_bonded(155) );
		TS_ASSERT( pose_->residue(71).is_bonded(155) );
		TS_ASSERT( pose_->residue(80).is_bonded(155) );
		TS_ASSERT( pose_->residue(83).is_bonded(155) );

		TS_ASSERT( pose_->residue(201).is_bonded(309) );
		TS_ASSERT( pose_->residue(203).is_bonded(309) );
		TS_ASSERT( pose_->residue(218).is_bonded(309) );
		TS_ASSERT( pose_->residue(275).is_bonded(309) );

		TS_ASSERT( pose_->residue(218).is_bonded(310) );
		TS_ASSERT( pose_->residue(226).is_bonded(310) );
		TS_ASSERT( pose_->residue(235).is_bonded(310) );
		TS_ASSERT( pose_->residue(238).is_bonded(310) );
		//Check that the pose has the proper number of constraints (should have 18 per setup metal)
		TS_ASSERT_EQUALS( pose_->constraint_set()->get_all_constraints().size(), 72 );

		//Select metals only on chain A
		pose_ = ref_pose_->clone();
		new_mover = protocols::simple_moves::SetupMetalsMover();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_a;
		ss_a << "<SetupMetalsMover name=\"metal\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_resnums=\"1-155\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_a );
		TS_TRACE( "Tag with resnums to select chain A" );
		TS_ASSERT_THROWS_NOTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		//Check all default values
		TS_ASSERT_EQUALS( new_mover.get_remove_hydrogens(), true );
		TS_ASSERT_EQUALS( new_mover.get_metal_resnums_string(), "1-155" );

		//The following should have been set from command line values in the constructor
		TS_ASSERT_DELTA( new_mover.get_metals_detection_LJ_multiplier(), 1.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_distance_constraint_multiplier(), 2.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_angle_constraint_multiplier(), 5.0, 1e-6 );
		//Apply the mover
		new_mover.apply( *pose_ );
		//Check that all bonds were created
		TS_ASSERT( pose_->residue(46).is_bonded(154) );
		TS_ASSERT( pose_->residue(48).is_bonded(154) );
		TS_ASSERT( pose_->residue(63).is_bonded(154) );
		TS_ASSERT( pose_->residue(120).is_bonded(154) );

		TS_ASSERT( pose_->residue(63).is_bonded(155) );
		TS_ASSERT( pose_->residue(71).is_bonded(155) );
		TS_ASSERT( pose_->residue(80).is_bonded(155) );
		TS_ASSERT( pose_->residue(83).is_bonded(155) );

		TS_ASSERT( !pose_->residue(201).is_bonded(309) );
		TS_ASSERT( !pose_->residue(203).is_bonded(309) );
		TS_ASSERT( !pose_->residue(218).is_bonded(309) );
		TS_ASSERT( !pose_->residue(275).is_bonded(309) );
		TS_ASSERT( !pose_->residue(218).is_bonded(310) );
		TS_ASSERT( !pose_->residue(226).is_bonded(310) );
		TS_ASSERT( !pose_->residue(235).is_bonded(310) );
		TS_ASSERT( !pose_->residue(238).is_bonded(310) );
		//Check that the pose has the proper number of constraints (should have 18 per setup metal)
		TS_ASSERT_EQUALS( pose_->constraint_set()->get_all_constraints().size(), 36 );

		//Select only zinc ions
		pose_ = ref_pose_->clone();
		new_mover = protocols::simple_moves::SetupMetalsMover();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_zn;
		ss_zn << "<SetupMetalsMover name=\"metal\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_resnums=\"155,310\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_zn );
		TS_TRACE( "Tag with resnums to select zinc ions" );
		TS_ASSERT_THROWS_NOTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		//Check all default values
		TS_ASSERT_EQUALS( new_mover.get_remove_hydrogens(), true );
		TS_ASSERT_EQUALS( new_mover.get_metal_selector(), nullptr );

		//The following should have been set from command line values in the constructor
		TS_ASSERT_DELTA( new_mover.get_metals_detection_LJ_multiplier(), 1.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_distance_constraint_multiplier(), 2.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_angle_constraint_multiplier(), 5.0, 1e-6 );
		//Apply the mover
		new_mover.apply( *pose_ );
		//Check that all bonds were created
		TS_ASSERT( !pose_->residue(46).is_bonded(154) );
		TS_ASSERT( !pose_->residue(48).is_bonded(154) );
		TS_ASSERT( !pose_->residue(63).is_bonded(154) );
		TS_ASSERT( !pose_->residue(120).is_bonded(154) );

		TS_ASSERT( pose_->residue(63).is_bonded(155) );
		TS_ASSERT( pose_->residue(71).is_bonded(155) );
		TS_ASSERT( pose_->residue(80).is_bonded(155) );
		TS_ASSERT( pose_->residue(83).is_bonded(155) );

		TS_ASSERT( !pose_->residue(201).is_bonded(309) );
		TS_ASSERT( !pose_->residue(203).is_bonded(309) );
		TS_ASSERT( !pose_->residue(218).is_bonded(309) );
		TS_ASSERT( !pose_->residue(275).is_bonded(309) );

		TS_ASSERT( pose_->residue(218).is_bonded(310) );
		TS_ASSERT( pose_->residue(226).is_bonded(310) );
		TS_ASSERT( pose_->residue(235).is_bonded(310) );
		TS_ASSERT( pose_->residue(238).is_bonded(310) );
		//Check that the pose has the proper number of constraints (should have 18 per setup metal)
		TS_ASSERT_EQUALS( pose_->constraint_set()->get_all_constraints().size(), 36 );

		//Select none of the metals (none of the metals are in the selection)
		pose_ = ref_pose_->clone();
		new_mover = protocols::simple_moves::SetupMetalsMover();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_none;
		ss_none << "<SetupMetalsMover name=\"metal\" metals_detection_LJ_multiplier=\"1.0\" metals_angle_constraint_multiplier=\"5.0\" metal_resnums=\"1,2,3\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_none );
		TS_TRACE( "Tag with selection that should not include any of the metal ions" );
		TS_ASSERT_THROWS_NOTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );
		//Check all default values
		TS_ASSERT_EQUALS( new_mover.get_remove_hydrogens(), true );
		TS_ASSERT_EQUALS( new_mover.get_metal_selector(), nullptr );

		//The following should have been set from command line values in the constructor
		TS_ASSERT_DELTA( new_mover.get_metals_detection_LJ_multiplier(), 1.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_distance_constraint_multiplier(), 2.0, 1e-6 );
		TS_ASSERT_DELTA( new_mover.get_metals_angle_constraint_multiplier(), 5.0, 1e-6 );
		//Apply the mover
		new_mover.apply( *pose_ );
		//Check that all bonds were created
		TS_ASSERT( !pose_->residue(46).is_bonded(154) );
		TS_ASSERT( !pose_->residue(48).is_bonded(154) );
		TS_ASSERT( !pose_->residue(63).is_bonded(154) );
		TS_ASSERT( !pose_->residue(120).is_bonded(154) );

		TS_ASSERT( !pose_->residue(63).is_bonded(155) );
		TS_ASSERT( !pose_->residue(71).is_bonded(155) );
		TS_ASSERT( !pose_->residue(80).is_bonded(155) );
		TS_ASSERT( !pose_->residue(83).is_bonded(155) );

		TS_ASSERT( !pose_->residue(201).is_bonded(309) );
		TS_ASSERT( !pose_->residue(203).is_bonded(309) );
		TS_ASSERT( !pose_->residue(218).is_bonded(309) );
		TS_ASSERT( !pose_->residue(275).is_bonded(309) );

		TS_ASSERT( !pose_->residue(218).is_bonded(310) );
		TS_ASSERT( !pose_->residue(226).is_bonded(310) );
		TS_ASSERT( !pose_->residue(235).is_bonded(310) );
		TS_ASSERT( !pose_->residue(238).is_bonded(310) );
		//Check that the pose has the proper number of constraints (should have 18 per setup metal)
		TS_ASSERT_EQUALS( pose_->constraint_set()->get_all_constraints().size(), 0 );

	}

	void test_contact_selection(){
		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		protocols::moves::Movers_map mm;
		core::pose::PoseOP pose_ = ref_pose_->clone();
		//Select all metals
		protocols::simple_moves::SetupMetalsMover new_mover;
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_all;
		ss_all << "<SetupMetalsMover name=\"metal\" metal_resnums=\"154,155,309,310\" metals_detection_LJ_multiplier=\"1.0\" add_constraints=\"false\" contact_resnums=\"46,71,218\" />" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss_all );
		TS_ASSERT_THROWS_NOTHING( new_mover.parse_my_tag( tag, dm, fm, mm, *pose_ ) );

		//Apply the mover
		new_mover.apply( *pose_ );
		//Check that all bonds were created
		TS_ASSERT( pose_->residue(46).is_bonded(154) );
		TS_ASSERT( !pose_->residue(48).is_bonded(154) );
		TS_ASSERT( !pose_->residue(63).is_bonded(154) );
		TS_ASSERT( !pose_->residue(120).is_bonded(154) );

		TS_ASSERT( !pose_->residue(63).is_bonded(155) );
		TS_ASSERT( pose_->residue(71).is_bonded(155) );
		TS_ASSERT( !pose_->residue(80).is_bonded(155) );
		TS_ASSERT( !pose_->residue(83).is_bonded(155) );

		TS_ASSERT( !pose_->residue(201).is_bonded(309) );
		TS_ASSERT( !pose_->residue(203).is_bonded(309) );
		TS_ASSERT( pose_->residue(218).is_bonded(309) );
		TS_ASSERT( !pose_->residue(275).is_bonded(309) );

		TS_ASSERT( pose_->residue(218).is_bonded(310) );
		TS_ASSERT( !pose_->residue(226).is_bonded(310) );
		TS_ASSERT( !pose_->residue(235).is_bonded(310) );
		TS_ASSERT( !pose_->residue(238).is_bonded(310) );
		//Check that the pose has the proper number of constraints (should have 18 per setup metal)
		TS_ASSERT_EQUALS( pose_->constraint_set()->get_all_constraints().size(), 0 );


	}

};
