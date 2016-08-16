// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/ResidueIndexSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::ResidueIndexSelector
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::select::residue_selector;


class ResidueIndexSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_ResidueIndexSelector_parse_my_tag() {
		std::string tag_string = "<Index name=index_rs resnums=5,7-8/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP index_rs( new ResidueIndexSelector );
		try {
			index_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = index_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );

		// test
		std::set < core::Size > acceptTrue;
		acceptTrue.insert( 5 );
		acceptTrue.insert( 7 );
		acceptTrue.insert( 8 );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( !subset[ ii ] || ( acceptTrue.find( ii ) != acceptTrue.end() ) )
				}
				}

				void test_ResidueIndexSelector_parse_my_tag_mixed_str() {
				std::string tag_string = "<Index name=index_rs resnums=2,3-4A />";
			std::stringstream ss( tag_string );
			utility::tag::TagOP tag( new utility::tag::Tag() );
			tag->read( ss );
			basic::datacache::DataMap dm;

			ResidueSelectorOP index_rs( new ResidueIndexSelector );
			try {
				index_rs->parse_my_tag( tag, dm );
			} catch ( utility::excn::EXCN_Msg_Exception e ) {
				std::cerr << "Exception!" << e.msg() << std::endl;
				TS_ASSERT( false ); // this parsing should succeed
			}

			core::pose::Pose trpcage = create_trpcage_ideal_pose();
			ResidueSubset subset = index_rs->apply( trpcage );

			TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );
			std::set< core::Size > acceptTrue;
			acceptTrue.insert(2);
			acceptTrue.insert( trpcage.pdb_info()->pdb2pose('A', 3) );
			acceptTrue.insert( trpcage.pdb_info()->pdb2pose('A', 4) );
			for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
				//if ( subset[ii] ) std::cerr << ii << "\n";
				TS_ASSERT( !subset[ ii ] || ( acceptTrue.find( ii ) != acceptTrue.end() ) )
					}
					}

					// make sure we fail if no selection string is provided
					void test_NeighbohoodResidueSelector_fail_no_resnums() {
					std::string tag_string = "<Index name=index_rs />";
				std::stringstream ss( tag_string );
				utility::tag::TagOP tag( new utility::tag::Tag() );
				tag->read( ss );
				basic::datacache::DataMap dm;

				ResidueSelectorOP index_rs( new ResidueIndexSelector );
				try {
					index_rs->parse_my_tag( tag, dm );
					TS_ASSERT( false ); //parsing should fail!
				} catch ( utility::excn::EXCN_Msg_Exception e ) {
					//std::cerr << "Exception (fail_no_resnums): " << e.msg();
					std::string expected_err = "Failed to access required option 'resnums' from ResidueIndexSelector::parse_my_tag.\nOption resnums not found.\n";
					TS_ASSERT_EQUALS( e.msg(), expected_err);
				}
			}

			// make sure we fail if indexed residues are out of range
			void test_ResidueIndexSelector_fail_index_out_of_range() {
				core::pose::Pose trpcage = create_trpcage_ideal_pose();
			std::stringstream bad_index;

			bad_index << "2,4," << trpcage.total_residue() + 1;

			ResidueSelectorOP index_rs( new ResidueIndexSelector( bad_index.str() ) );

			try {
				index_rs->apply( trpcage );
				TS_ASSERT( false );
			} catch( utility::excn::EXCN_Msg_Exception e) {
				//std::cerr << "Exception (fail_index_out_of_range): '" << e.msg() << "'";
				std::string expected_err = "Residue 21 not found in pose!\n";
				TS_ASSERT( e.msg() == expected_err);
			}
		}


		// cannot test since bad chain does not throw but exits internally
		void dont_test_ResidueIndexSelector_fail_chain_out_of_range() {
			core::pose::Pose trpcage = create_trpcage_ideal_pose();
		std::string bad_index = "2,4,7D";

		ResidueSelectorOP index_rs( new ResidueIndexSelector( bad_index ) );

		try {
			index_rs->apply( trpcage );
			TS_ASSERT( false );
		} catch( utility::excn::EXCN_Msg_Exception e) {
			//std::cerr << "Exception (fail_chain_out_of_range): " << e.msg();
			std::string expected_err = "Residue 0 not found in pose!\n";
			TS_ASSERT( e.msg() == expected_err);
		}
	}

};
