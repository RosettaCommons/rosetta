// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pack/task/residue_selector/AndResidueSelector.cxxtest.hh
/// @brief  test suite for core::pack::task::residue_selector::AndResidueSelector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/pack/task/residue_selector/DummySelectors.hh>

// Package headers
#include <core/pack/task/residue_selector/AndResidueSelector.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::pack::task::residue_selector;

class AndResidueSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	// @brief make sure that when we register a residue selector, we can later get it back
	void test_and_residue_selector_odd_1mod5() {
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		ResidueSelectorOP one_mod_5_rs( new XModYResidueSelector( 1, 5 ) );
		ResidueSelectorOP and_rs( new AndResidueSelector( odd_rs, one_mod_5_rs ) );
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		ResidueSubset subset = and_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			TS_ASSERT_EQUALS( subset[ ii ], (ii % 5 == 1 && ii % 2 == 1) );
		}
	}

	/// @brief Test AndResidueSelector::parse_my_tag
	void test_AndResidueSelector_parse_my_tag() {
		std::string tag_string = "<And name=and_rs selectors=odd,one_mod_five/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		ResidueSelectorOP one_mod_5_rs( new XModYResidueSelector( 1, 5 ) );
		dm.add( "ResidueSelector", "odd", odd_rs );
		dm.add( "ResidueSelector", "one_mod_five", one_mod_5_rs );

		ResidueSelectorOP and_rs( new AndResidueSelector );
		try {
			and_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = and_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			TS_ASSERT_EQUALS( subset[ ii ], (ii % 5 == 1 && ii % 2 == 1) );
		}
	}

	/// @brief Test that an excpetion is thrown if the AndResidueSelector is ever initialized
	/// from parse_my_tag where no ResidueSelectors have been provided.
	void test_AndResidueSelector_parse_my_tag_no_provided_selectors() {
		std::string tag_string = "<And name=and_rs/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP and_rs( new AndResidueSelector );
		try {
			and_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected = "No ResidueSelectors given to the AndResidueSelector; AndResidueSelector requires at least one ResidueSelector as input\n";
			TS_ASSERT_EQUALS( e.msg(), expected );
		}

	}

	/// @brief Test than an exception is thrown if the AndResidueSelector is initialized
	/// from parse_my_tag where the ResidueSelectors it requests are not in the datamap
	void test_AndResidueSelector_parse_my_tag_selectors_not_in_datamap() {
		std::string tag_string = "<And name=and_rs selectors=odd,one_mod_five/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP and_rs( new AndResidueSelector );
		try {
			and_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected = "Failed to find ResidueSelector named 'odd' from the Datamap from AndResidueSelector::parse_my_tag.\nERROR: Could not find ResidueSelector and name odd in Datamap\n";
			TS_ASSERT_EQUALS( e.msg(), expected );
		}
	}

	void test_AndResidueSelector_parse_subtag() {
		std::string tag_string = "<And name=and_rs>\n\t<Index resnums=2-4 />\n\t<Index resnums=3-5 />\n</And>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP and_rs( new AndResidueSelector );
		try {
			and_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			TS_ASSERT( false );
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = and_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			TS_ASSERT_EQUALS( subset[ ii ], (ii == 3 || ii == 4) );
		}
	}

	void test_AndResidueSelector_fail_parse_subtag() {
		std::string tag_string = "<And name=and_rs>\n\t<Index resnums=2-4 />\n\t<Bogus />\n</And>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP and_rs( new AndResidueSelector );
		try {
			and_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // parsing should fail here
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string err_msg =  "No ResidueSelectorCreator with the name 'Bogus' has been registered with the ResidueSelectorFactory";
			TS_ASSERT( e.msg() == err_msg );
		}
	}


};
