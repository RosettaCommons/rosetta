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
#include <core/select/residue_selector/SelectorLogicParser.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/dssp/Dssp.hh>

// Utility headers
#include <utility/backtrace.hh>
//#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/memory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::select::residue_selector;


class SelectorLogicParserTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_SelectorParser_or_combination() {

		auto res5( utility::pointer::make_shared< ResidueIndexSelector >( 5 ) );
		auto res10( utility::pointer::make_shared< ResidueIndexSelector >( 10 ) );

		basic::datacache::DataMap dm;
		dm.add( "ResidueSelector", "res5", res5 );
		dm.add( "ResidueSelector", "res10", res10 );

		ResidueSelectorOP rs;
		SelectorLogicParser parser;
		try {
			rs = parser.parse_string_to_residue_selector( dm, "res5 OR res10" );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		ResidueSubset subset = rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT_EQUALS( bool(subset[ ii ]), ii == 5 || ii == 10  );
		}
	}

	void test_SelectorParser_and_combination() {

		auto res5to15( utility::pointer::make_shared< ResidueIndexSelector >( "5-15" ) );
		auto res10to20( utility::pointer::make_shared< ResidueIndexSelector >( "10-20" ) );

		basic::datacache::DataMap dm;
		dm.add( "ResidueSelector", "res5to15", res5to15 );
		dm.add( "ResidueSelector", "res10to20", res10to20 );

		ResidueSelectorOP rs;
		SelectorLogicParser parser;
		try {
			rs = parser.parse_string_to_residue_selector( dm, "res5to15 AND res10to20" );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		ResidueSubset subset = rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT_EQUALS( bool(subset[ ii ]), ii >= 10 && ii <= 15  );
		}
	}

	void test_SelectorParser_not_and_combination() {

		auto res5to15( utility::pointer::make_shared< ResidueIndexSelector >( "5-15" ) );
		auto res10to20( utility::pointer::make_shared< ResidueIndexSelector >( "10-20" ) );

		basic::datacache::DataMap dm;
		dm.add( "ResidueSelector", "res5to15", res5to15 );
		dm.add( "ResidueSelector", "res10to20", res10to20 );

		ResidueSelectorOP rs;
		SelectorLogicParser parser;
		try {
			rs = parser.parse_string_to_residue_selector( dm, "res5to15 AND ! res10to20" );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		ResidueSubset subset = rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT_EQUALS( bool(subset[ ii ]), ii >= 5 && ii < 10  );
		}
	}

	void test_SelectorParser_not_or_combination() {

		auto res5to15( utility::pointer::make_shared< ResidueIndexSelector >( "5-15" ) );
		auto res10to20( utility::pointer::make_shared< ResidueIndexSelector >( "10-20" ) );

		basic::datacache::DataMap dm;
		dm.add( "ResidueSelector", "res5to15", res5to15 );
		dm.add( "ResidueSelector", "res10to20", res10to20 );

		ResidueSelectorOP rs;
		SelectorLogicParser parser;
		try {
			rs = parser.parse_string_to_residue_selector( dm, "res5to15 OR ! res10to20" );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		ResidueSubset subset = rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT_EQUALS( bool(subset[ ii ]), ii >= 1 && ii <= 15  );
		}
	}


	void test_SelectorParser_not_or_combination2() {

		auto res5to15( utility::pointer::make_shared< ResidueIndexSelector >( "5-15" ) );
		auto res10to20( utility::pointer::make_shared< ResidueIndexSelector >( "10-20" ) );

		basic::datacache::DataMap dm;
		dm.add( "ResidueSelector", "res5to15", res5to15 );
		dm.add( "ResidueSelector", "res10to20", res10to20 );

		ResidueSelectorOP rs;
		SelectorLogicParser parser;
		try {
			// make sure that ! "sticks to" res5to15 and not to the OR of the two RSs.
			rs = parser.parse_string_to_residue_selector( dm, "! res5to15 OR res10to20" );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		ResidueSubset subset = rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT_EQUALS( bool(subset[ ii ]), ( ii >= 1 && ii < 5 ) || ( ii >= 10 )  );
		}
	}

	void test_SelectorParser_missing_residue_selector() {

		auto res5to15( utility::pointer::make_shared< ResidueIndexSelector >( "5-15" ) );
		//auto res10to20( utility::pointer::make_shared< ResidueIndexSelector >( "10-20" ) );

		basic::datacache::DataMap dm;
		dm.add( "ResidueSelector", "res5to15", res5to15 );
		//dm.add( "ResidueSelector", "res10to20", res10to20 );

		ResidueSelectorOP rs;
		SelectorLogicParser parser;
		try {
			// make sure that ! "sticks to" res5to15 and not to the OR of the two RSs.
			set_throw_on_next_assertion_failure();
			rs = parser.parse_string_to_residue_selector( dm, "! res5to15 OR res10to20" );
			TS_ASSERT( false );
		} catch (utility::excn::Exception & e ) {
			//std::cerr << "Exception!" << e.msg() << std::endl;
			std::string gold_standard1 =
				"Failed to tokenize string logically combining ResidueSelectors.\n"
				"\"! res5to15 OR res10to20\"\n"
				"Allowed names for ResidueSelectors are:\n"
				"   res5to15\n"
				"Error message from Scanner:\n";

			std::string gold_standard2 =
				"ERROR: Error in scan_identifier: Failed to find input_string in either variables_ or functions_ maps\n"
				"res10to20\n";
			TS_ASSERT_STRING_CONTAINS( e.msg(), gold_standard1 );
			TS_ASSERT_STRING_CONTAINS( e.msg(), gold_standard2 );
		}

	}

};
