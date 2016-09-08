// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/NotResidueSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::NotResidueSelector
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <core/select/residue_selector/NotResidueSelector.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::select::residue_selector;

class NotResidueSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	// @brief make sure that the not selector negates the odd selector
	void test_not_residue_selector_odd() {
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		ResidueSelectorOP not_rs( new NotResidueSelector( odd_rs ) );
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = not_rs->apply( trpcage );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 2 == 0) );
		}
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_NotResidueSelector_parse_my_tag() {
		std::string tag_string = "<Not name=not_rs selector=odd/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		dm.add( "ResidueSelector", "odd", odd_rs );

		ResidueSelectorOP not_rs( new NotResidueSelector );
		try {
			not_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = not_rs->apply( trpcage );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == (ii % 2 == 0) );
		}
	}


	/// @brief Test that an excpetion is thrown if the NotResidueSelector is ever initialized
	/// from parse_my_tag where no ResidueSelector has been provided.
	void test_NotResidueSelector_parse_my_tag_no_provided_selectors() {
		std::string tag_string = "<Not name=not_rs />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP not_rs( new NotResidueSelector );
		try {
			not_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected = "No ResidueSelector given to the NotResidueSelector; NotResidueSelector requires a ResidueSelector as input\n";
			TS_ASSERT_EQUALS( e.msg(), expected );
		}
	}

	/// @brief Test than an exception is thrown if the AndResidueSelector is initialized
	/// from parse_my_tag where the ResidueSelectors it requests are not in the datamap
	void test_NotResidueSelector_parse_my_tag_selectors_not_in_datamap() {
		std::string tag_string = "<Not name=not_rs selector=odd/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP not_rs( new NotResidueSelector );
		try {
			not_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected = "Failed to find ResidueSelector named 'odd' from the Datamap from NotResidueSelector::parse_my_tag.\nERROR: Could not find ResidueSelector and name odd in Datamap\n";
			TS_ASSERT_EQUALS( e.msg(), expected );
		}
	}


};
