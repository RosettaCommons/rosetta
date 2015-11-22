// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/select/residue_selector/ResidueNameSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::ResidueNameSelector
/// @author Tom Linsky (tlinsky@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <core/select/residue_selector/ResidueNameSelector.hh>

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


class ResidueNameSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_ResidueNameSelector_parse_my_tag() {
		std::stringstream ss;
		ss << "<ResidueName name=allala residue_names=ASP,SER:CtermProteinFull,LEU,CYS />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP name_rs( new ResidueNameSelector );
		try {
			name_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = name_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );

		// test
		std::set < core::Size > acceptTrue;
		acceptTrue.insert( 9 );
		acceptTrue.insert( 20 );
		acceptTrue.insert( 2 );
		acceptTrue.insert( 7 );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( !subset[ ii ] || ( acceptTrue.find( ii ) != acceptTrue.end() ) );
		}
	}
	/// @brief Test NotResidueSelector::parse_my_tag
	void test_ResidueNameSelector_name3() {
		std::stringstream ss;
		ss << "<ResidueName name=allala residue_name3=ASP,SER,CYS />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP name_rs( new ResidueNameSelector );
		try {
			name_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = name_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );

		// test
		std::set < core::Size > acceptTrue;
		acceptTrue.insert( 9 );
		acceptTrue.insert( 13 );
		acceptTrue.insert( 14 );
		acceptTrue.insert( 20 );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( !subset[ ii ] || ( acceptTrue.find( ii ) != acceptTrue.end() ) );
		}
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_ResidueNameSelector_name_and_name3() {
		std::stringstream ss;
		ss << "<ResidueName name=allala residue_names=SER:CtermProteinFull residue_name3=ASP,LEU,CYS />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP name_rs( new ResidueNameSelector );
		try {
			name_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = name_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );

		// test
		std::set < core::Size > acceptTrue;
		acceptTrue.insert( 2 );
		acceptTrue.insert( 7 );
		acceptTrue.insert( 9 );
		acceptTrue.insert( 20 );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( !subset[ ii ] || ( acceptTrue.find( ii ) != acceptTrue.end() ) );
		}
	}

	// make sure we fail if no selection string is provided
	void test_NeighbohoodResidueSelector_fail_no_resnums() {
		std::stringstream ss;
		ss << "<ResidueName name=I_should_fail />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP name_rs( new ResidueNameSelector );
		try {
			name_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); //parsing should fail!
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			// should fail with a message
			TS_ASSERT( e.msg() != "" );
		}
	}

	// make sure we fail if indexed residues are out of range
	void test_ResidueNameSelector_fail_bad_name() {
		std::string names = "ALA,HDA";
		/// borrow trp cage code from ResidueIndexSelector Unit test
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		ResidueSelectorOP name_rs( new ResidueNameSelector( names ) );

		try {
			name_rs->apply( trpcage );
			TS_ASSERT( false );
		} catch( utility::excn::EXCN_Msg_Exception e) {
			TS_ASSERT( e.msg() != "" );
		}
	}


	// unknown variant should throw user error
	void test_ResidueNameSelector_fail_bad_variant() {
		// The T in NTermProteinFull should not be capitalized
		std::string names = "ALA:NTermProteinFull";
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		ResidueSelectorOP name_rs( new ResidueNameSelector( names ) );

		try {
			name_rs->apply( trpcage );
			TS_ASSERT( false );
		} catch( utility::excn::EXCN_Msg_Exception e) {
			TS_ASSERT( e.msg() != "" );
		}
	}

};
