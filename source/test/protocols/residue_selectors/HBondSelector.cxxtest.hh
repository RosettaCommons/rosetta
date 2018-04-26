// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/BondedResidueSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::BondedResidueSelector
/// @author Sharon Guffy (guffy@email.unc.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <protocols/residue_selectors/HBondSelector.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace protocols::residue_selectors;

static basic::Tracer TR( "test.core.select.HBondSelectorTests" );

class HBondSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_HBondSelector_parse_my_tag() {
		std::stringstream ss;
		ss << "<HBond name=\"hbond\" />"; //No attributes are required, so this should work
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		HBondSelectorOP rs( new HBondSelector );

		TS_TRACE( "Tag with no options" );
		TS_ASSERT_THROWS_ANYTHING( rs->parse_my_tag( tag, dm ) );

		core::scoring::ScoreFunctionOP dummy_scorefxn( new core::scoring::ScoreFunction );
		dm.add( "scorefxns", "dummy", dummy_scorefxn );

		std::stringstream ss_scorefxn;
		ss_scorefxn << "<HBond name=\"hbond\" scorefxn=\"dummy\" />"; //No attributes are required, so this should work
		tag->read( ss_scorefxn );
		TS_TRACE( "Tag with dummy scorefxn" );
		TS_ASSERT_THROWS_NOTHING( rs->parse_my_tag( tag, dm ) );
		//Make sure defaults are set correctly
		TS_TRACE( "Testing default values" );
		TS_ASSERT_EQUALS( rs->get_hbond_energy_cutoff(), -0.5 );
		TS_ASSERT_EQUALS( rs->get_include_bb_bb(), false );
		//The score function should be the default
		//It should have input_set_defined_ false, use_input_set_selector_ true (and is a TrueResidueSelector)
		TS_ASSERT_EQUALS( rs->get_input_set_defined(), false );
		TS_ASSERT( rs->get_use_input_set_selector() );
		TS_ASSERT_EQUALS( rs->get_input_set_selector()->get_name(), "True" );


		std::stringstream ss_set_options;
		ss_set_options << "<HBond name=\"hbond\" scorefxn=\"dummy\" hbond_energy_cutoff=\"-1\" include_bb_bb=\"true\" />";
		TS_TRACE( "Tag with options but no selector(s)" );
		tag->read( ss_set_options );
		TS_ASSERT_THROWS_NOTHING( rs->parse_my_tag( tag, dm ) );
		TS_ASSERT_EQUALS( rs->get_hbond_energy_cutoff(), -1 );
		TS_ASSERT( rs->get_include_bb_bb() );


		std::stringstream ss_two_selectors;
		ss_two_selectors << "<HBond name=\"hbond\" scorefxn=\"dummy\" residue_selector=\"dummy\"><Index name=\"index\" resnums=\"2,3\" /></HBond>";
		tag->read( ss_two_selectors );
		TS_ASSERT_THROWS( rs->parse_my_tag( tag, dm ), utility::excn::Exception & );

		std::stringstream ssgood_resnums;
		ssgood_resnums << "<HBond name=\"hbond\" scorefxn=\"dummy\" resnums=\"2,3\" />";
		tag->read( ssgood_resnums );
		TS_ASSERT_THROWS_NOTHING( rs->parse_my_tag( tag, dm ) );
		TS_ASSERT_EQUALS( rs->get_input_set_str(), "2,3" );
		TS_ASSERT_EQUALS( rs->get_use_input_set_selector(), false );
		TS_ASSERT_EQUALS( rs->get_input_set_defined(), true );


		std::stringstream ss_good_selector;
		ss_good_selector << "<HBond name=\"hbond\" scorefxn=\"dummy\" ><Index name=\"index\" resnums=\"2,3\" /></HBond>";
		tag->read( ss_good_selector );
		TS_ASSERT_THROWS_NOTHING( rs->parse_my_tag( tag, dm ) );
		TS_ASSERT_EQUALS( rs->get_use_input_set_selector(), true );
		TS_ASSERT_EQUALS( rs->get_input_set_defined(), true );
		TS_ASSERT_EQUALS( rs->get_input_set_selector()->get_name(), "Index" );

		std::stringstream ss_undefined_selector;
		ss_undefined_selector << "<HBond name=\"hbond\" scorefxn=\"dummy\" residue_selector=\"dummy\" />";
		tag->read( ss_undefined_selector );
		TS_ASSERT_THROWS_ANYTHING( rs->parse_my_tag( tag, dm ) );

		core::select::residue_selector::ResidueIndexSelectorOP dummy( new core::select::residue_selector::ResidueIndexSelector );
		dummy->set_index( "2,3" );
		dm.add( "ResidueSelector", "dummy", dummy );
		TS_ASSERT_THROWS_NOTHING( rs->parse_my_tag( tag, dm ) );

	}

};
