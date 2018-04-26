// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/SecondaryStructureSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::SecondaryStructureSelector
/// @author Tom Linsky (tlinsky@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <core/select/residue_selector/SecondaryStructureSelector.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/dssp/Dssp.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::select::residue_selector;


class SecondaryStructureSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_SecondaryStructureSelector_parse_my_tag() {
		std::stringstream ss;
		ss << "<SecondaryStructure name=\"allloops\" ss=\"L\" />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		SecondaryStructureSelectorOP rs( new SecondaryStructureSelector );
		try {
			rs->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::scoring::dssp::Dssp dssp( trpcage );
		std::string const secstruct = dssp.get_dssp_secstruct();
		for ( core::Size i=1, endi=trpcage.size(); i<=endi; ++i ) {
			trpcage.set_secstruct( i, secstruct[ i - 1 ] );
		}

		ResidueSubset subset = rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		// test
		std::set < core::Size > acceptTrue;
		acceptTrue.insert( 8 );
		acceptTrue.insert( 9 );
		acceptTrue.insert( 10 );
		acceptTrue.insert( 15 );
		acceptTrue.insert( 16 );
		acceptTrue.insert( 17 );
		acceptTrue.insert( 18 );
		acceptTrue.insert( 19 );
		acceptTrue.insert( 20 );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( ( !subset[ ii ] && ( acceptTrue.find( ii ) == acceptTrue.end() ) ) ||
				( subset[ ii ] && ( acceptTrue.find( ii ) != acceptTrue.end() ) ) );
		}

		// now include terminal "loop"
		rs->set_include_terminal_loops( true );
		subset = rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		acceptTrue.insert( 1 );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( ( !subset[ ii ] && ( acceptTrue.find( ii ) == acceptTrue.end() ) ) ||
				( subset[ ii ] && ( acceptTrue.find( ii ) != acceptTrue.end() ) ) );
		}
	}

	void test_SecondaryStructureSelector_overlap() {
		std::stringstream ss;
		ss << "<SecondaryStructure name=allala ss=\"HE\" overlap=\"1\" />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP rs( new SecondaryStructureSelector );
		try {
			rs->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::scoring::dssp::Dssp dssp( trpcage );
		std::string const secstruct = dssp.get_dssp_secstruct();
		for ( core::Size i=1, endi=trpcage.size(); i<=endi; ++i ) {
			trpcage.set_secstruct( i, secstruct[ i - 1 ] );
		}

		ResidueSubset subset = rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		// test
		std::set < core::Size > acceptTrue;
		acceptTrue.insert( 1 );
		acceptTrue.insert( 2 );
		acceptTrue.insert( 3 );
		acceptTrue.insert( 4 );
		acceptTrue.insert( 5 );
		acceptTrue.insert( 6 );
		acceptTrue.insert( 7 );
		acceptTrue.insert( 8 );
		acceptTrue.insert( 10 );
		acceptTrue.insert( 11 );
		acceptTrue.insert( 12 );
		acceptTrue.insert( 13 );
		acceptTrue.insert( 14 );
		acceptTrue.insert( 15 );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( ( !subset[ ii ] && ( acceptTrue.find( ii ) == acceptTrue.end() ) ) ||
				( subset[ ii ] && ( acceptTrue.find( ii ) != acceptTrue.end() ) ) );
		}
	}

	// make sure we fail if no selection string is provided
	void test_SecondaryStructureSelector_fail_no_ss() {
		std::stringstream ss;
		ss << "<SecondaryStructure name=I_should_fail />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP name_rs( new SecondaryStructureSelector );
		try {
			name_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); //parsing should fail!
		} catch (utility::excn::Exception & e ) {
			// should fail with a message
			TS_ASSERT( e.msg() != "" );
		}
	}

	// make sure we fail if no selection string is provided
	void test_SecondaryStructureSelector_fail_no_ss_C() {
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSelectorOP rs( new SecondaryStructureSelector );
		try {
			ResidueSubset subset = rs->apply( trpcage );
			TS_ASSERT( false ); //parsing should fail!
		} catch (utility::excn::Exception & e ) {
			// should fail with a message
			TS_ASSERT( e.msg() != "" );
		}
	}

	// make sure we fail if indexed residues are out of range
	void test_SecondaryStructureSelector_fail_bad_ss() {
		std::string names = "DL";
		/// borrow trp cage code from ResidueIndexSelector Unit test
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		try {
			ResidueSelectorOP name_rs( new SecondaryStructureSelector( names ) );
			name_rs->apply( trpcage );
			TS_ASSERT( false );
		} catch (utility::excn::Exception & e) {
			TS_ASSERT( e.msg() != "" );
		}
	}

};
