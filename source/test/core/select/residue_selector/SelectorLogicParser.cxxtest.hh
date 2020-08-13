// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/SelectorLogicParser.cxxtest.hh
/// @brief test suite for core::select::residue_selector::SelectorLogicParser
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) //Selection logic within all RS selector.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>
#include <core/select/residue_selector/util.hh>
#include <test/core/select/residue_selector/utilities_for_testing.hh>

// Package headers
#include <core/select/residue_selector/SelectorLogicParser.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>


// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/parser/DataLoader.hh>

// Utility headers
#include <utility/backtrace.hh>
//#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/memory.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

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

		// repeat with a lower-case "or"
		try {
			rs = parser.parse_string_to_residue_selector( dm, "res5 or res10" );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		subset = rs->apply( trpcage );
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

		// repeat with lower-case "and"
		try {
			rs = parser.parse_string_to_residue_selector( dm, "res5to15 and res10to20" );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		subset = rs->apply( trpcage );
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

		try {
			rs = parser.parse_string_to_residue_selector( dm, "res5to15 and not res10to20" );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		subset = rs->apply( trpcage );
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
			rs = parser.parse_string_to_residue_selector( dm, "res5to15 OR NOT res10to20" );
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

		try {
			rs = parser.parse_string_to_residue_selector( dm, "res5to15 or not res10to20" );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		subset = rs->apply( trpcage );
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

		// repeat with lower-case "not" and "or"
		try {
			// make sure that ! "sticks to" res5to15 and not to the OR of the two RSs.
			rs = parser.parse_string_to_residue_selector( dm, "not res5to15 or res10to20" );
		} catch (utility::excn::Exception & e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		subset = rs->apply( trpcage );
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
			rs = parser.parse_string_to_residue_selector( dm, "not res5to15 OR res10to20" );
			TS_ASSERT( false );
		} catch (utility::excn::Exception & e ) {
			//std::cerr << "Exception!" << e.msg() << std::endl;
			std::string gold_standard1 =
				"Failed to tokenize string logically combining ResidueSelectors.\n"
				"\"not res5to15 OR res10to20\"\n"
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

	void test_combine_logic_within_all_selections(){
		//JAB - this tests the addition of logic parsing within get_residue_selector/parse_residue_selector
		// It is not exhaustive, but tests the base get/parsing functions used by most selectors

		using namespace core::pose;

		using utility::tag::TagOP;

		basic::datacache::DataMap dm;
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		auto first( utility::pointer::make_shared< ResidueIndexSelector >( "1" ) );
		auto last( utility::pointer::make_shared< ResidueIndexSelector >( "10" ) );

		auto first3( utility::pointer::make_shared< ResidueIndexSelector >( "1-3" ) );

		dm.add("ResidueSelector", "first3", first3);
		dm.add("ResidueSelector", "last",  last);
		dm.add("ResidueSelector", "first",  first);

		//Or

		ResidueSelectorCOP rs = get_residue_selector("first or last", dm);


		ResidueSubset subset = rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		utility::vector1< bool > or_correct(trpcage.size(), false);
		or_correct[1] = true;
		or_correct[10] = true;
		compare_bool_vector(subset, or_correct);

		//And + Not

		ResidueSelectorCOP rs2 = get_residue_selector("first and not last", dm);
		subset = rs2->apply(trpcage);
		utility::vector1< bool > and_correct(trpcage.size(), false);
		and_correct[1] = true;
		compare_bool_vector(subset, and_correct);

		//Regular Selector

		ResidueSelectorCOP rs3 = get_residue_selector("first", dm);

		subset = rs3->apply(trpcage);
		compare_bool_vector(subset, and_correct);

		//With tag Parsing
		utility::tag::TagOP tag = utility::tag::Tag::create( "<True name=\"true_sel\" selector=\"first or last\"/>\n" );
		ResidueSelectorCOP rs4 = parse_residue_selector(tag, dm, "selector");
		and_correct[10] = true;

		subset = rs4->apply(trpcage);
		compare_bool_vector(subset, and_correct);
	}

	void test_SelectorParser_error_handling_with_illegal_varnames() {
		auto threonines( utility::pointer::make_shared< ResidueNameSelector >( "THR", true ));
		auto not_thr( utility::pointer::make_shared< NotResidueSelector >( threonines ) );

		basic::datacache::DataMap dm1, dm2, dm3;
		dm1.add( "ResidueSelector", "T", threonines );
		dm1.add( "ResidueSelector", "noT", not_thr );
		ResidueSelectorOP rs;
		SelectorLogicParser parser;
		try {
			set_throw_on_next_assertion_failure();
			rs = parser.parse_string_to_residue_selector( dm1, "T or noT" );
			TS_ASSERT( false ); // the line above should throw
		} catch (utility::excn::Exception & e ) {
			std::string gold_standard1 =
				"Cannot proceed because of the residue selector named 'noT'."
				" This selector should be renamed.";
			std::string gold_standard2 =
				"Illegal variable named 'noT' ('NOT') conflicts with the boolean"
				" operator of the same name.";
			TS_ASSERT_STRING_CONTAINS( e.msg(), gold_standard1 );
			TS_ASSERT_STRING_CONTAINS( e.msg(), gold_standard2 );
		}


		dm2.add( "ResidueSelector", "T", threonines );
		dm2.add( "ResidueSelector", "oR", not_thr );
		try {
			set_throw_on_next_assertion_failure();
			rs = parser.parse_string_to_residue_selector( dm2, "T or oR" );
			TS_ASSERT( false ); // the line above should throw
		} catch (utility::excn::Exception & e ) {
			std::string gold_standard1 =
				"Cannot proceed because of the residue selector named 'oR'."
				" This selector should be renamed.";
			std::string gold_standard2 =
				"Illegal variable named 'oR' ('OR') conflicts with the boolean"
				" operator of the same name.";
			TS_ASSERT_STRING_CONTAINS( e.msg(), gold_standard1 );
			TS_ASSERT_STRING_CONTAINS( e.msg(), gold_standard2 );
		}

		dm3.add( "ResidueSelector", "T", threonines );
		dm3.add( "ResidueSelector", "aNd", not_thr );
		try {
			set_throw_on_next_assertion_failure();
			rs = parser.parse_string_to_residue_selector( dm3, "T or aNd" );
			TS_ASSERT( false ); // the line above should throw
		} catch (utility::excn::Exception & e ) {
			std::string gold_standard1 =
				"Cannot proceed because of the residue selector named 'aNd'."
				" This selector should be renamed.";
			std::string gold_standard2 =
				"Illegal variable named 'aNd' ('AND') conflicts with the boolean"
				" operator of the same name.";
			TS_ASSERT_STRING_CONTAINS( e.msg(), gold_standard1 );
			TS_ASSERT_STRING_CONTAINS( e.msg(), gold_standard2 );
		}
	}


};
