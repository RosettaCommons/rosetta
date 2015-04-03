// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pack/task/residue_selector/NeighborhoodResidueSelector.cxxtest.hh
/// @brief test suite for core::pack::task::residue_selector::NeighborhoodResidueSelector
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/pack/task/residue_selector/DummySelectors.hh>

// Package headers
#include <core/pack/task/residue_selector/NeighborhoodResidueSelector.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::pack::task::residue_selector;


class NeighborhoodResidueSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_NeighborhoodResidueSelector_parse_my_tag_selector() {
		std::string tag_string = "<Neighborhood name=neighbor_rs selector=odd distance=5.2/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		dm.add( "ResidueSelector", "odd", odd_rs );

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		try {
			neighbor_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = neighbor_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );

		// check the result
		// 1. generate fake focus
		std::set< core::Size > testFocus;
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ii += 2 ) {
			testFocus.insert( ii );
		}
		// test
		TS_ASSERT( check_calculation( trpcage, subset, testFocus, 5.2 ) );
	}

	// make sure to fail if a tag specifies both a selector and resnums to parse
	void test_NeighborhoodResidueSelector_fail_resnum_and_selector_set() {
		std::string tag_string = "<Neighborhood name=neighbor_rs selector=odd resnums=2,4,6 distance=5.2/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		dm.add( "ResidueSelector", "odd", odd_rs );

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		try {
			neighbor_rs->parse_my_tag( tag, dm );
			TS_ASSERT ( false );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_err =  "NeighborhoodResidueSelector takes EITHER 'selector' OR 'resnum' options, not both!\n";
			TS_ASSERT ( e.msg() == expected_err );
		}
	}


	void test_NeighborhoodResidueSelector_parse_my_tag_str() {
		std::string tag_string = "<Neighborhood name=neighbor_rs resnums=2,3,5 distance=5.2/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		try {
			neighbor_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = neighbor_rs->apply( trpcage );

		std::set< core::Size > testFocus;
		testFocus.insert(2);
		testFocus.insert(3);
		testFocus.insert(5);
		TS_ASSERT( check_calculation( trpcage, subset, testFocus, 5.2 ) );

	}

	// make sure we fail if neither selector nor focus string are provided
	void test_NeighbohoodResidueSelector_fail_no_focus() {
		std::string tag_string = "<Neighborhood name=neighbor_rs distance=5.2/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		try {
			neighbor_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); //parsing should fail!
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_err = "Failed to access option 'resnums' from NeighborhoodResidueSelector::parse_my_tag.\nOption resnums not found.\n";
			TS_ASSERT( e.msg() == expected_err);
		}
	}

	// make sure we fail if focus residues are out of range
	void test_NeighborhoodResidueSelector_fail_focus_out_of_range() {
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		std::set< core::Size > bad_focus;
		bad_focus.insert(1);
		bad_focus.insert(trpcage.total_residue() + 1);

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector(bad_focus, 10.0) );

		ResidueSubset subset( trpcage.total_residue(), false );
		try {
			neighbor_rs->apply( trpcage );
			TS_ASSERT( false );
		} catch( utility::excn::EXCN_Msg_Exception e) {
			std::string expected_err = "Residue 21 not found in pose!\n";
			TS_ASSERT( e.msg() == expected_err);
		}
	}

	/// @brief Test than an exception is thrown if the AndResidueSelector is initialized
	/// from parse_my_tag where the ResidueSelectors it requests are not in the datamap
	void test_NeighborhoodResidueSelector_parse_my_tag_selector_not_in_datamap() {
		std::string tag_string = "<Neighborhood name=neighbor_rs selector=odd/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		try {
			neighbor_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected = "Failed to find ResidueSelector named 'odd' from the Datamap from NeighborhoodResidueSelector::parse_my_tag.\nERROR: Could not find ResidueSelector and name odd in Datamap\n";
			TS_ASSERT_EQUALS( e.msg(), expected );
		}
	}

	// desired behavior is that the most recent call to set_focus or set_focus_selector
	// determines which source of focus residues is used
	void test_NeighborhoodResidueSelector_use_last_provided_source_of_focus() {
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		std::set< core::Size > focus_set;
		focus_set.insert(2);
		focus_set.insert(3);
		NeighborhoodResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector(focus_set, 5.0) );
		ResidueSelectorOP odd_rs( new OddResidueSelector );

		ResidueSubset subset( trpcage.total_residue(), false );
		TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );

		std::set< core::Size > testFocus_odd;
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ii += 2 ) {
			testFocus_odd.insert( ii );
		}

		try {
			subset = neighbor_rs->apply( trpcage );
			TS_ASSERT( check_calculation( trpcage, subset, focus_set, 5.0 ) );

			neighbor_rs->set_focus_selector( odd_rs );
			subset = neighbor_rs->apply( trpcage );
			TS_ASSERT( check_calculation( trpcage, subset, testFocus_odd, 5.0 ) );

			neighbor_rs->set_focus( focus_set );
			subset = neighbor_rs->apply( trpcage );
			TS_ASSERT( check_calculation( trpcage, subset, focus_set, 5.0 ) );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception! " << e.msg();
			TS_ASSERT( false );
		}
	}

	//make sure to fail if ResidueSelector subtags are provided as well as the resnums option
	void test_NeighborhoodSelector_fail_subtag_and_resnums() {
		std::string tag_string = "<Neighborhood name=neighbor_rs resnums=4-8>\n\t<Index resnums=2-3 />\n</Neighborhood>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		try {
			neighbor_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected = "NeighborhoodResidueSelector takes EITHER a 'resnums' tag or a selector subtag, not both!\n";
			TS_ASSERT_EQUALS( e.msg(), expected );
		}
	}

	//make sure to fail if ResidueSelector subtags are provided as well as a selector option
	void test_NeighborhoodSelector_fail_subtag_and_selector() {
		std::string tag_string = "<Neighborhood name=neighbor_rs selector=odd>\n\t<Index resnums=2-3 />\n</Neighborhood>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		try {
			neighbor_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected = "NeighborhoodResidueSelector can only have one ResidueSelector loaded!\n";
			TS_ASSERT_EQUALS( e.msg(), expected );
		}
	}


	bool
	check_calculation( core::pose::Pose const & pose,
		ResidueSubset const & subset,
		std::set< core::Size > const & focus,
		core::Real distance) {

		ResidueSubset ctrl_subset(subset.size(), false);

		core::Real const dst_squared = distance * distance;
		for( core::Size ii = 1; ii < subset.size() ; ++ii )
		{
			if( focus.find( ii ) != focus.end() ) {
				ctrl_subset[ ii ] = true;
				continue;
			}
			core::conformation::Residue const & r1( pose.residue( ii ) );
			for( std::set< core::Size >::const_iterator it = focus.begin();
				it != focus.end(); ++it )
			{
				if( *it == 0 || *it > pose.total_residue())
						  return false;

				core::conformation::Residue const & r2( pose.residue( *it ) );
				core::Real const d_sq( r1.xyz( r1.nbr_atom() ).distance_squared( r2.xyz( r2.nbr_atom() ) ) );
				if(d_sq <= dst_squared)
					ctrl_subset[ ii ] = true;
			} // focus set
			if(ctrl_subset[ ii ] != subset[ ii ])
					  return false;
		} // subset

		// no mismatches found
		return true;
}

};
