// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ChainSelector.cxxtest.hh
/// @brief  test suite for core::select::residue_selector::ChainSelector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>

// Package headers
#include <core/select/residue_selector/ChainSelector.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

#include <core/conformation/Residue.hh> // AUTO IWYU For Pose::Residue Residue
#include <core/init_util.hh> // AUTO IWYU For core_init

using namespace core::select::residue_selector;


class ChainSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief test that we select all the residues in trpcage
	void test_chain_selector_chain_A_select_all() {
		ChainSelectorOP chain_rs( new ChainSelector );
		utility::vector1< std::string > chain_strings;
		chain_strings.push_back( "A" );
		chain_rs->set_chain_strings( chain_strings );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = chain_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == true );
		}
	}

	/// @brief hack trpcage to set half its residues to chain B
	/// test that we select only the chain A residues
	void test_chain_selector_chain_A_select_half() {
		ChainSelectorOP chain_rs( new ChainSelector );
		utility::vector1< std::string > chain_strings; chain_strings.push_back( "A" );
		chain_rs->set_chain_strings( chain_strings );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		for ( core::Size ii = 11; ii <= trpcage.size(); ++ii ) {
			trpcage.pdb_info()->chain( ii, "B" );
		}
		ResidueSubset subset = chain_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == ( ii <= 10 ) );
		}

		utility::vector1< core::Size > const positions = chain_rs->selection_positions( trpcage );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			auto const iter = std::find( positions.begin(), positions.end(), ii );
			bool const present = ( iter != positions.end() );
			if ( ii <= 10 ) {
				TS_ASSERT( present );
			} else {
				TS_ASSERT( ! present );
			}
		}

	}

	/// @brief hack trpcage to set half its residues to chain B
	/// test that when we give two chains to the ChainSelector, we get
	/// all the residues in trpcage
	void test_chain_selector_chain_AandB_select_all() {
		ChainSelectorOP chain_rs( new ChainSelector );
		utility::vector1< std::string > chain_strings;
		chain_strings.push_back( "A" );
		chain_strings.push_back( "B" );
		chain_rs->set_chain_strings( chain_strings );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		for ( core::Size ii = 11; ii <= trpcage.size(); ++ii ) {
			trpcage.pdb_info()->chain( ii, "B" );
		}
		ResidueSubset subset = chain_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == true );
		}
	}

	/// @brief test that we select all the residues in trpcage
	void test_chain_selector_chain_1_select_all() {
		ChainSelectorOP chain_rs( new ChainSelector );
		utility::vector1< std::string > chain_strings; chain_strings.push_back( "1" );
		chain_rs->set_chain_strings( chain_strings );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = chain_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == true );
		}
	}

	/// @brief hack together a two-chained trpcage to give it a chain B.
	/// Test that we select only the chain A residues
	void test_chain_selector_chain_1_select_half() {
		ChainSelectorOP chain_rs( new ChainSelector );
		utility::vector1< std::string > chain_strings; chain_strings.push_back( "1" );
		chain_rs->set_chain_strings( chain_strings );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::pose::Pose trpcage2 = trpcage;

		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			core::conformation::ResidueOP iiclone = trpcage.residue( ii ).clone();
			iiclone->chain( 'B' );
			if ( ii == 1 ) {
				trpcage2.append_residue_by_jump( *iiclone, 1, "", "", true );
			} else {
				trpcage2.append_residue_by_bond( *iiclone );
			}
		}
		ResidueSubset subset = chain_rs->apply( trpcage2 );
		for ( core::Size ii = 1; ii <= trpcage2.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == ( ii <= trpcage.size() ) );
		}
	}

	/// @brief hack trpcage to set half its residues to chain B
	/// test that when we give two chains to the ChainSelector, we get
	/// all the residues in trpcage
	void test_chain_selector_chain_1and2_select_all() {
		ChainSelectorOP chain_rs( new ChainSelector );
		utility::vector1< std::string > chain_strings;
		chain_strings.push_back( "1" );
		chain_strings.push_back( "2" );
		chain_rs->set_chain_strings( chain_strings );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::pose::Pose trpcage2 = trpcage;

		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			core::conformation::ResidueOP iiclone = trpcage.residue( ii ).clone();
			iiclone->chain( 'B' );
			if ( ii == 1 ) {
				trpcage2.append_residue_by_jump( *iiclone, 1, "", "", true );
			} else {
				trpcage2.append_residue_by_bond( *iiclone );
			}
		}

		ResidueSubset subset = chain_rs->apply( trpcage2 );
		for ( core::Size ii = 1; ii <= trpcage2.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == true );
		}
	}

	/// @brief Test ChainSelector::parse_my_tag
	void test_ChainSelector_parse_my_tag() {
		std::string tag_string = "<Chain name=chain_rs chains=A/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP chain_rs( new ChainSelector );
		try {
			chain_rs->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
			return;
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = chain_rs->apply( trpcage );
		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == true );
		}
	}

	/// @brief Test ChainSelector::parse_my_tag two chains
	void test_ChainSelector_parse_my_tag_w_two_chains() {
		std::string tag_string = "<Chain name=chain_rs chains=\"1,2\"/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP chain_rs( new ChainSelector );
		try {
			chain_rs->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
			return;
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::pose::Pose trpcage2 = trpcage;

		for ( core::Size ii = 1; ii <= trpcage.size(); ++ii ) {
			core::conformation::ResidueOP iiclone = trpcage.residue( ii ).clone();
			iiclone->chain( 'B' );
			if ( ii == 1 ) {
				trpcage2.append_residue_by_jump( *iiclone, 1, "", "", true );
			} else {
				trpcage2.append_residue_by_bond( *iiclone );
			}
		}

		ResidueSubset subset = chain_rs->apply( trpcage2 );
		for ( core::Size ii = 1; ii <= trpcage2.size(); ++ii ) {
			TS_ASSERT( subset[ ii ] == true );
		}
	}


	/// @brief Test that an excpetion is thrown if the ChainSelector is ever initialized
	/// from parse_my_tag where no "chains" option was provided.
	void test_ChainSelector_parse_my_tag_no_provided_chains() {
		std::string tag_string = "<Chain name=chain_A_rs />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP chain_rs( new ChainSelector );
		try {
			chain_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); // this parsing should fail
		} catch (utility::excn::Exception & e ) {
			std::string expected = "ChainSelector::parse_my_tag was not able to find the required option 'chains' in the input Tag"; // \nOption chains not found.\n
			TS_ASSERT( e.msg().find(expected) != std::string::npos );
		}

	}


};
