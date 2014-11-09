// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/InterGroupInterfaceByVectorSelector.cxxtest.hh
/// @brief  test suite for core::pack::task::residue_selector::InterGroupInterfaceByVectorSelector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/pack/task/residue_selector/DummySelectors.hh>

// Package headers
#include <core/pack/task/residue_selector/InterGroupInterfaceByVectorSelector.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pack/task/operation/util/interface_vector_calculate.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::pack::task::residue_selector;


class ResRangeSelector : public ResidueSelector {
public:
	ResRangeSelector( Size lower, Size upper ) : lower_( lower ), upper_( upper ) {}

	virtual
	ResidueSubset
	apply( core::pose::Pose const & pose ) const
	{
		ResidueSubset subset( pose.total_residue(), false );
		for ( core::Size ii = lower_; ii <= upper_; ++ii ) subset[ ii ] = true;
		return subset;
	}

	virtual std::string get_name() const { return "ResRange"; }
private:
	core::Size lower_, upper_;
};


class InterGroupInterfaceByVectorSelectorTests : public CxxTest::TestSuite {
public:
	typedef std::pair< std::set< core::Size >, std::set< core::Size > > two_sets;
public:

	void setUp() {
		core_init();
	}

	two_sets
	two_test_in_groups() {
		// group1
		std::set< core::Size > grp1;
		for ( core::Size ii = 1; ii <= 20; ++ii ) grp1.insert( ii );

		// group2
		std::set< core::Size > grp2;
		for ( core::Size ii = 45; ii <= 66; ++ii ) grp2.insert( ii );
		return std::make_pair( grp1, grp2 );
	}

	ResidueSubset
	gold_result( core::pose::Pose const & pose ) {
		two_sets both_sets = two_test_in_groups();
		return core::pack::task::operation::util::calc_interacting_vector(
			pose, both_sets.first, both_sets.second, 11.0, 5.5, 75.0, 9.0 );
	}


	/// @brief make sure the igibv selector properly identifies the interface residues in test_in.pdb
	/// where the interface is between the two groups grp1= 49-69 and grp2= 95-116 both on chain A.
	void test_interface_vector_selector_from_std_set_on_test_in_pdb() {
		InterGroupInterfaceByVectorSelectorOP igibv_rs( new InterGroupInterfaceByVectorSelector );

		two_sets both_sets = two_test_in_groups();

		igibv_rs->group1_set( both_sets.first );
		igibv_rs->group2_set( both_sets.second );

		core::pose::Pose pose = create_test_in_pdb_pose();
		ResidueSubset subset = igibv_rs->apply( pose );

		ResidueSubset ground_truth = gold_result( pose );
		TS_ASSERT_EQUALS( subset.size(), ground_truth.size() );
		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			//std::cout << ii << " subset[" << ii << "] = " << subset[ ii ] << " ground_truth[ " << ii << " ] = " << ground_truth[ ii ] << std::endl;
			TS_ASSERT( subset[ ii ] == ground_truth[ ii ] );
		}
	}

	/// @brief make sure the igibv selector properly identifies the interface residues in test_in.pdb
	/// where the interface is between the two groups grp1= 49-69 and grp2= 95-116 both on chain A.
	void test_interface_vector_selector_from_resselector_on_test_in_pdb() {
		InterGroupInterfaceByVectorSelectorOP igibv_rs( new InterGroupInterfaceByVectorSelector );

		igibv_rs->group1_selector( ResidueSelectorCOP( new ResRangeSelector( 1, 20 ) ) );
		igibv_rs->group2_selector( ResidueSelectorCOP( new ResRangeSelector( 45, 66 ) ) );

		core::pose::Pose pose = create_test_in_pdb_pose();
		ResidueSubset subset = igibv_rs->apply( pose );

		ResidueSubset ground_truth = gold_result( pose );

		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			TS_ASSERT( subset[ ii ] == ground_truth[ ii ] );
		}
	}


	/// @brief Test InterGroupInterfaceByVectorSelector::parse_my_tag
	void test_InterGroupInterfaceByVectorSelector_parse_my_tag() {
		std::string tag_string = "<InterfaceByVector name=int_rs grp1_selector=res1to20 grp2_selector=res45to66 />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP range1_rs( new ResRangeSelector( 1, 20 ) );
		ResidueSelectorOP range2_rs( new ResRangeSelector( 45, 66 ) );
		dm.add( "ResidueSelector", "res1to20", range1_rs );
		dm.add( "ResidueSelector", "res45to66", range2_rs );

		ResidueSelectorOP igibv_rs( new InterGroupInterfaceByVectorSelector );
		try {
			igibv_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
			return;
		}

		core::pose::Pose pose = create_test_in_pdb_pose();
		ResidueSubset subset =  igibv_rs->apply( pose );

		ResidueSubset ground_truth = gold_result( pose );

		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			TS_ASSERT( subset[ ii ] == ground_truth[ ii ] );
		}
	}


	/// @brief Test that an excpetion is thrown if the InterGroupInterfaceByVectorSelector is ever initialized
	/// from parse_my_tag where no ResidueSelectors have been provided.
	void test_InterGroupInterfaceByVectorSelector_parse_my_tag_no_provided_grp1_data() {
		std::string tag_string = "<InterfaceByVector name=int_rs grp2_selector=res45to66 />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP range1_rs( new ResRangeSelector( 1, 20 ) );
		ResidueSelectorOP range2_rs( new ResRangeSelector( 45, 66 ) );
		dm.add( "ResidueSelector", "res1to20", range1_rs );
		dm.add( "ResidueSelector", "res45to66", range2_rs );

		ResidueSelectorOP igibv_rs( new InterGroupInterfaceByVectorSelector );
		try {
			igibv_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should not succeed
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			//std::cerr << "Exception!" << e.msg() << std::endl;
			std::string expected_error = "InterGroupInterfaceByVectorSelector::parse_my_tag requires either grp1_selector or grp1_residues to be specified\n";
			TS_ASSERT_EQUALS( e.msg(), expected_error );
			return;
		}

	}

	void test_InterGroupInterfaceByVectorSelector_parse_my_tag_no_provided_grp2_data() {
		std::string tag_string = "<InterfaceByVector name=int_rs grp1_selector=res45to66 />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP range1_rs( new ResRangeSelector( 1, 20 ) );
		ResidueSelectorOP range2_rs( new ResRangeSelector( 45, 66 ) );
		dm.add( "ResidueSelector", "res1to20", range1_rs );
		dm.add( "ResidueSelector", "res45to66", range2_rs );

		ResidueSelectorOP igibv_rs( new InterGroupInterfaceByVectorSelector );
		try {
			igibv_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should not succeed
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			//std::cerr << "Exception!" << e.msg() << std::endl;
			std::string expected_error = "InterGroupInterfaceByVectorSelector::parse_my_tag requires either grp2_selector or grp2_residues to be specified\n";
			TS_ASSERT_EQUALS( e.msg(), expected_error );
		}

	}

	/// @brief Test than an exception is thrown if the InterGroupInterfaceByVectorSelector is initialized
	/// from parse_my_tag where the ResidueSelectors it requests are not in the datamap
	void test_InterGroupInterfaceByVectorSelector_parse_my_tag_selectors_not_in_datamap() {
		std::string tag_string = "<InterfaceByVector name=int_rs grp1_selector=bogus grp2_selector=res45to66 />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP igibv_rs( new InterGroupInterfaceByVectorSelector );
		try {
			igibv_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should not succeed
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
	//		std::cerr << "Exception!" << e.msg() << std::endl;
			std::string expected_error = "Failed to find ResidueSelector named 'bogus' from the Datamap from InterGroupInterfaceByVectorSelector::parse_my_tag\nERROR: Could not find ResidueSelector and name bogus in Datamap\n";
			TS_ASSERT_EQUALS( e.msg(), expected_error );
			return;
		}

	}

	/// @brief Test than an exception is thrown if the InterGroupInterfaceByVectorSelector is initialized
	/// from parse_my_tag where the ResidueSelectors it requests are not in the datamap
	void test_InterGroupInterfaceByVectorSelector_parse_my_tag_fail_only_one_subselector() {
		std::string tag_string = "<InterfaceByVector name=int_rs grp1_selector=bogus>\n\t<Index resnames=12-17 />\n</InterfaceByVector>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
	
		ResidueSelectorOP igibv_rs( new InterGroupInterfaceByVectorSelector );
		try {
			igibv_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should not succeed
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
		//	std::cerr << "Exception!" << e.msg() << std::endl;
			std::string expected_error = "InterGroupInterfaceByVectorSelector takes either two or zero subtags to specify residue groups!\n";
			TS_ASSERT_EQUALS( e.msg(), expected_error );
			return;
		}
}



};
