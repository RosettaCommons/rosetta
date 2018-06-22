// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/InterfaceHydrophobicResidueContactsFilter.cxxtest.hh
/// @brief test suite for protocols::simple_moves::InterfaceHydrophobicResidueContactsFilter
/// @author Brian Coventry (bcov@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <protocols/simple_filters/InterfaceHydrophobicResidueContactsFilter.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/ChainSelector.hh>

#include <protocols/rosetta_scripts/XmlObjects.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>

using namespace protocols::simple_filters;

static basic::Tracer TR("protocols.simple_filters.InterfaceHydrophobicResidueContactsFilter");


class InterfaceHydrophobicResidueContactsFilterTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		core::import_pose::pose_from_file( pose_1af0_ab, "protocols/membrane/1AFO_AB.pdb" );
	}

	void test_filter( InterfaceHydrophobicResidueContactsFilterOP filter, core::Size exp_num_contacts, bool exp_pass ) {
		core::Size num_contacts = filter->compute( pose_1af0_ab );
		bool pass = filter->apply( pose_1af0_ab );

		TR << "Expected " << exp_num_contacts << " contacts. Found " << num_contacts << std::endl;
		TR << "Expected " << ( exp_pass ? "pass" : "fail" ) << " . Found " << ( pass ? "pass" : "fail" ) << std::endl;

		TS_ASSERT( num_contacts == exp_num_contacts );
		TS_ASSERT( exp_pass == pass );
	}


	void test_InterfaceHydrophobicResidueContactsFilter_from_xml_defaults() {
		TR << "test_InterfaceHydrophobicResidueContactsFilter_from_xml_defaults" << std::endl;

		std::string xml =
			"    <SCOREFXNS> "
			"        <ScoreFunction name=\"sfxn\" weights=\"ref2015_soft\" /> "
			"    </SCOREFXNS> "
			"\t<RESIDUE_SELECTORS> "
			"        <Chain name=\"chainA\" chains=\"A\"/> "
			"        <Chain name=\"chainB\" chains=\"B\"/> "
			"\t</RESIDUE_SELECTORS> "
			"    <FILTERS> "
			"        <InterfaceHydrophobicResidueContacts name=\"hydrophobic_residue_contacts\"  "
			"        \ttarget_selector=\"chainA\" binder_selector=\"chainB\" scorefxn=\"sfxn\" /> "
			"    </FILTERS> "
			;

		protocols::rosetta_scripts::XmlObjectsCOP objs = protocols::rosetta_scripts::XmlObjects::create_from_string( xml );
		InterfaceHydrophobicResidueContactsFilterOP filter = std::dynamic_pointer_cast<InterfaceHydrophobicResidueContactsFilter>(
			objs->get_filter("hydrophobic_residue_contacts") );

		test_filter( filter, 6, true );
	}

	void test_InterfaceHydrophobicResidueContactsFilter_from_xml_not_defaults() {
		TR << "test_InterfaceHydrophobicResidueContactsFilter_from_xml_not_defaults" << std::endl;

		std::string xml =
			"   <SCOREFXNS> "
			"       <ScoreFunction name=\"sfxn\" weights=\"ref2015_soft\" /> "
			"   </SCOREFXNS> "
			"\t<RESIDUE_SELECTORS> "
			"       <Chain name=\"chainA\" chains=\"A\"/> "
			"       <Chain name=\"chainB\" chains=\"B\"/> "
			"\t</RESIDUE_SELECTORS> "
			"   <FILTERS> "
			"       <InterfaceHydrophobicResidueContacts name=\"hydrophobic_residue_contacts\"  "
			"        \ttarget_selector=\"chainA\" binder_selector=\"chainB\" scorefxn=\"sfxn\"  "
			"        \tapolar_res=\"VAL\" score_cut=\"-0.001\" threshold=\"1\" "
			"        \t/> "
			"   </FILTERS> "
			;

		protocols::rosetta_scripts::XmlObjectsCOP objs = protocols::rosetta_scripts::XmlObjects::create_from_string( xml );
		InterfaceHydrophobicResidueContactsFilterOP filter = std::dynamic_pointer_cast<InterfaceHydrophobicResidueContactsFilter>(
			objs->get_filter("hydrophobic_residue_contacts") );

		test_filter( filter, 2, true );
	}

	void test_InterfaceHydrophobicResidueContactsFilter_from_constructor_not_defaults() {
		TR << "test_InterfaceHydrophobicResidueContactsFilter_from_constructor_not_defaults" << std::endl;

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("ref2015_soft");
		core::select::residue_selector::ChainSelectorOP target( new core::select::residue_selector::ChainSelector("A"));
		core::select::residue_selector::ChainSelectorOP binder( new core::select::residue_selector::ChainSelector("B"));

		InterfaceHydrophobicResidueContactsFilterOP filter( new InterfaceHydrophobicResidueContactsFilter(
			1, target, binder, scorefxn, -0.001, "VAL" ));

		test_filter( filter, 2, true );
	}

private:
	core::pose::Pose pose_1af0_ab;

};
