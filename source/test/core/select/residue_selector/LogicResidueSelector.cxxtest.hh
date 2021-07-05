// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/select/residue_selector/LogicResidueSelector.cxxtest.hh
/// @brief  Testing the core::select::residue_selector::LogicResidueSelector
/// @author Frances Chu (francesc345@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <core/select/residue_selector/LogicResidueSelector.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

#include <iostream> // AUTO IWYU For operator<<, basic_ostream, endl, cerr, ostream

static basic::Tracer TR("LogicResidueSelector");

using namespace core::select::residue_selector;


class LogicResidueSelectorTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();
	}


	/// @brief Test LogicResidueSelector::parse_my_tag when "selector" is given
	void test_logic_selector_parse_my_tag_selector() {
		std::string tag_string = "<LogicResidueSelector name=logic_rs selector=\"odd OR one_mod_five\" />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		ResidueSelectorOP one_mod_5_rs( new XModYResidueSelector( 1, 5 ));
		dm.add( "ResidueSelector", "odd", odd_rs );
		dm.add( "ResidueSelector", "one_mod_five", one_mod_5_rs );

		ResidueSelectorOP logic_rs( new LogicResidueSelector );
		try {
			logic_rs->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = logic_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 5 == 1 || ii % 2 == 1) );
		}
	}

	/// @brief Test LogicResidueSelector::parse_my_tag when "selectors" is given
	void test_logic_selector_parse_my_tag_selectors() {
		std::string tag_string = "<LogicResidueSelector name=logic_rs selectors=\"odd OR one_mod_five\" />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		ResidueSelectorOP one_mod_5_rs( new XModYResidueSelector( 1, 5 ));
		dm.add( "ResidueSelector", "odd", odd_rs );
		dm.add( "ResidueSelector", "one_mod_five", one_mod_5_rs );

		ResidueSelectorOP logic_rs( new LogicResidueSelector );
		try {
			logic_rs->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = logic_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 5 == 1 || ii % 2 == 1) );
		}
	}

	/// @brief Test LogicResidueSelector::parse_my_tag when "residue_selector" is given
	void test_logic_selector_parse_my_tag_residue_selector() {
		std::string tag_string = "<LogicResidueSelector name=logic_rs residue_selector=\"odd OR one_mod_five\" />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		ResidueSelectorOP one_mod_5_rs( new XModYResidueSelector( 1, 5 ));
		dm.add( "ResidueSelector", "odd", odd_rs );
		dm.add( "ResidueSelector", "one_mod_five", one_mod_5_rs );

		ResidueSelectorOP logic_rs( new LogicResidueSelector );
		try {
			logic_rs->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = logic_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 5 == 1 || ii % 2 == 1) );
		}
	}

	/// @brief Test LogicResidueSelector::parse_my_tag when "residue_selectors" is given
	void test_logic_selector_parse_my_tag_residue_selectors() {
		std::string tag_string = "<LogicResidueSelector name=logic_rs residue_selectors=\"odd OR one_mod_five\" />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		ResidueSelectorOP one_mod_5_rs( new XModYResidueSelector( 1, 5 ));
		dm.add( "ResidueSelector", "odd", odd_rs );
		dm.add( "ResidueSelector", "one_mod_five", one_mod_5_rs );

		ResidueSelectorOP logic_rs( new LogicResidueSelector );
		try {
			logic_rs->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = logic_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 5 == 1 || ii % 2 == 1) );
		}
	}

	/// @brief Test that an exception is thrown if the LogicResidueSelector is ever initialized
	/// from parse_my_tag where no ResidueSelectors have been provided.
	void test_logic_selector_parse_my_tag_more_than_one_provided_selectors() {
		std::string tag_string = "<LogicResidueSelector name=logic_rs selector=\"odd OR one_mod_five\" residue_selector=\"odd OR one_mod_five\"/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP logic_rs( new LogicResidueSelector );
		try {
			logic_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch (utility::excn::Exception & e ) {
			std::string expected = "LogicResidueSelector must be given exactly one of the following residue selector attributes: \"selector\", \"selectors\", \"residue_selector\", or \"residue_selectors\". The LogicResidueSelector performs the same function no matter the attribute used.\n";
			TS_ASSERT( e.msg().find(expected) != std::string::npos );
		}

	}

	void test_logic_selector_parse_my_tag_no_name_option() {
		std::string tag_string = "<LogicResidueSelector selector=\"odd OR one_mod_five\" residue_selector=\"odd OR one_mod_five\"/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP logic_rs( new LogicResidueSelector );
		try {
			logic_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch (utility::excn::Exception & e ) {
			std::string expected = "LogicResidueSelector must be given exactly one of the following residue selector attributes: \"selector\", \"selectors\", \"residue_selector\", or \"residue_selectors\". The LogicResidueSelector performs the same function no matter the attribute used.\n";
			TS_ASSERT( e.msg().find(expected) != std::string::npos );
		}
	}



	/// @brief Test that an exception is thrown if the LogicResidueSelector is ever initialized
	/// from parse_my_tag where no ResidueSelectors have been provided.
	void test_logic_selector_parse_my_tag_no_provided_selectors() {
		std::string tag_string = "<LogicResidueSelector name=logic_rs />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP logic_rs( new LogicResidueSelector );
		try {
			logic_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch (utility::excn::Exception & e ) {
			std::string expected = "LogicResidueSelector must be given exactly one of the following residue selector attributes: \"selector\", \"selectors\", \"residue_selector\", or \"residue_selectors\". The LogicResidueSelector performs the same function no matter the attribute used.\n";
			TS_ASSERT( e.msg().find(expected) != std::string::npos );
		}

	}

	void test_logic_selector_parse_my_tag_comma_in_selectors() {
		std::string tag_string = "<LogicResidueSelector name=logic_rs residue_selectors=\"odd,one_mod_five\"/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP logic_rs( new LogicResidueSelector );
		try {
			logic_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch (utility::excn::Exception & e ) {
			std::string expected = "LogicResidueSelector does not take a comma separated list of residue selectors. Residue selectors should be separated by boolean(s) \"and/or/not\" with spaces on both sides.";
			TS_ASSERT( e.msg().find(expected) != std::string::npos );
		}

	}

	void tearDown(){

	}




};
