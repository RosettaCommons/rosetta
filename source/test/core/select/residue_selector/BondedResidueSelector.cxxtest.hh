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
#include <core/select/residue_selector/BondedResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Project headers

//This isn't ideal, but it's for testing purposes

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::select::residue_selector;

static THREAD_LOCAL basic::Tracer TR( "test.core.select.BondedResidueSelectorTests" );

class BondedSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_BondedResidueSelector_parse_my_tag() {
		std::stringstream ss;
		ss << "<Bonded name=\"bonded\" resnums=\"2,3\" residue_selector=\"dummy\" />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		BondedResidueSelectorOP rs( new BondedResidueSelector );
		TS_TRACE( "Resnums and selector" );
		TS_ASSERT_THROWS( rs->parse_my_tag( tag, dm ), utility::excn::EXCN_Msg_Exception );

		std::stringstream ss_two_selectors;
		ss_two_selectors << "<Bonded name=\"bonded\" residue_selector=\"dummy\"><Index name=\"index\" resnums=\"2,3\" /></Bonded>";
		tag->read( ss_two_selectors );
		TS_TRACE( "Two selectors" );
		TS_ASSERT_THROWS( rs->parse_my_tag( tag, dm ), utility::excn::EXCN_Msg_Exception );

		std::stringstream ssgood_resnums;
		ssgood_resnums << "<Bonded name=\"bonded\" resnums=\"2,3\" />";
		tag->read( ssgood_resnums );
		TS_TRACE( "Resnums only" );
		TS_ASSERT_THROWS_NOTHING( rs->parse_my_tag( tag, dm ) );

		std::stringstream ss_good_selector;
		ss_good_selector << "<Bonded name=\"bonded\" ><Index name=\"index\" resnums=\"2,3\" /></Bonded>";
		tag->read( ss_good_selector );
		TS_TRACE( "Selector subtag" );
		TS_ASSERT_THROWS_NOTHING( rs->parse_my_tag( tag, dm ) );

		std::stringstream ss_undefined_selector;
		ss_undefined_selector << "<Bonded name=\"bonded\" residue_selector=\"dummy\" />";
		tag->read( ss_undefined_selector );
		TS_TRACE( "Undefined selector" );
		TS_ASSERT_THROWS_ANYTHING( rs->parse_my_tag( tag, dm ) )

			core::select::residue_selector::ResidueIndexSelectorOP dummy( new core::select::residue_selector::ResidueIndexSelector );
		dummy->set_index( "2,3" );
		dm.add( "ResidueSelector", "dummy", dummy );
		TS_TRACE( "Undefined selector after defining the selector" );
		std::stringstream ss_defined_selector;
		ss_defined_selector << "<Bonded name=\"bonded\" residue_selector=\"dummy\" />";
		tag->read( ss_defined_selector );
		TS_ASSERT_THROWS_NOTHING( rs->parse_my_tag( tag, dm ) );
	}

	void test_detects_adjacent_residues(){
		BondedResidueSelectorOP rs( new BondedResidueSelector );
		rs->set_input_set( "2,3" );
		//Make a very simple pose--structure is irrelevant in this case
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "AAAAAAAAA", core::chemical::FA_STANDARD );

		core::select::residue_selector::ResidueSubset result = rs->apply( pose );
		TS_ASSERT( result[ 1 ] );
		TS_ASSERT( result[ 4 ] );
		TS_ASSERT( !result[ 5 ] );
		//The input residues are also included
		TS_ASSERT( result[ 2 ] );
		TS_ASSERT( result[ 3 ] );
	}

	void test_detects_declared_chemical_bonds(){
		BondedResidueSelectorOP rs( new BondedResidueSelector );
		rs->set_input_set( "1,3" );
		//Make a very simple pose--structure is irrelevant in this case
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "AAAAAA", core::chemical::FA_STANDARD );
		//Remove variant types from ends
		core::chemical::ResidueTypeCOP res1_type = pose.residue(1).type().get_base_type_cop();
		core::pose::replace_pose_residue_copying_existing_coordinates( pose, 1, *res1_type );
		core::chemical::ResidueTypeCOP res6_type = pose.residue(6).type().get_base_type_cop();
		core::pose::replace_pose_residue_copying_existing_coordinates( pose, 6, *res6_type );

		pose.conformation().declare_chemical_bond( 1, "N", 6, "C" );
		core::select::residue_selector::ResidueSubset result = rs->apply( pose );
		//These should still work
		TS_ASSERT( result[ 1 ] );
		TS_ASSERT( result[ 2 ] );
		TS_ASSERT( result[ 3 ] );
		TS_ASSERT( result[ 4 ] );
		TS_ASSERT( result[ 6 ] );
		TS_ASSERT( !result[ 5 ] );
	}
};
