// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	 ResidueProperties.cxxtest.hh
/// @brief   Test suite for ResidueProperties
/// @author  Labonte <JWLabonte@jhu.edu>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>

// Utility header
#include <utility/excn/EXCN_Base.hh>

// C++ header
#include <map>


using namespace core::chemical;
using namespace utility::excn;


class ResiduePropertiesTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init();
		ResidueTypeCAP residue_type;
		test_properties_ = new ResidueProperties( residue_type );
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	// Confirm that property-related methods function properly.
	void test_properties()
	{
		TS_TRACE( "Testing property-related methods of ResidueProperties." );
		TS_ASSERT( ! test_properties_->has_property( LIPID ) );
		TS_ASSERT( ! test_properties_->has_property( "LIPID" ) );
		TS_ASSERT( ! test_properties_->has_property( "FAT" ) );

		test_properties_->set_property( LIPID, true );
		TS_ASSERT_THROWS_EQUALS(
				test_properties_->set_property( "FAT", true ),
				EXCN_Base const & e,
				e.msg().substr( e.msg().find( "ERROR: " ) ),
				"ERROR: Rosetta does not recognize the property: FAT; "
				"has it been added to general_properties.list?\n\n" );

		TS_ASSERT( test_properties_->has_property( LIPID ) );
		TS_ASSERT( test_properties_->has_property( "LIPID") );
		TS_ASSERT( ! test_properties_->has_property( "FAT") );

		TS_ASSERT_EQUALS( test_properties_->get_list_of_properties()[ 1 ], "LIPID" );

		TS_ASSERT_EQUALS( test_properties_->get_list_of_properties().size(), 1 );
		test_properties_->set_property( "LIPID", true );
		TS_ASSERT_EQUALS( test_properties_->get_list_of_properties().size(), 1 );
		test_properties_->set_property( LIPID, false );
		TS_ASSERT_EQUALS( test_properties_->get_list_of_properties().size(), 0 );
	}

	// Confirm that variant-related methods function properly.
	void test_variant_types()
	{
		TS_TRACE( "Testing variant-related methods of ResidueProperties." );
		TS_ASSERT( ! test_properties_->is_variant_type( SC_FRAGMENT ) );
		TS_ASSERT( ! test_properties_->is_variant_type( "SC_FRAGMENT" ) );
		TS_ASSERT( ! test_properties_->is_variant_type( "BIZARRO" ) );

		test_properties_->set_variant_type( SC_FRAGMENT, true );
		TS_ASSERT_THROWS_EQUALS(
				test_properties_->set_variant_type( "BIZARRO", true ),
				EXCN_Base const & e,
				e.msg().substr( e.msg().find( "ERROR: " ) ),
				"ERROR: Rosetta does not recognize the variant: BIZARRO; "
				"has it been added to variant_types.list?\n\n" );
		test_properties_->enable_custom_variant_types();
		TS_ASSERT( test_properties_->has_custom_variant_types() );
		test_properties_->set_variant_type( "WHACKY", true );
		test_properties_->set_variant_type( "BIZARRO", true );
		test_properties_->set_variant_type( "WHACKY", false );
		TS_ASSERT_THROWS_EQUALS(
				test_properties_->set_variant_type( "GNARLY", false ),
				EXCN_Base const & e,
				e.msg().substr( e.msg().find( "ERROR: " ) ),
				"ERROR: Rosetta does not recognize the custom variant: GNARLY\n\n" );

		TS_ASSERT( test_properties_->is_variant_type( SC_FRAGMENT ) );
		TS_ASSERT( test_properties_->is_variant_type( "SC_FRAGMENT" ) );
		TS_ASSERT( test_properties_->is_variant_type( "BIZARRO" ) );

		TS_ASSERT_EQUALS( test_properties_->get_list_of_variants()[ 1 ], "SC_FRAGMENT");
		TS_ASSERT_EQUALS( test_properties_->get_list_of_variants()[ 2 ], "BIZARRO");
		TS_ASSERT_EQUALS( test_properties_->get_list_of_variants().size(), 2 );
		test_properties_->set_variant_type( SPECIAL_ROT, true );
		TS_ASSERT_EQUALS( test_properties_->get_list_of_variants().size(), 3 );
		test_properties_->set_variant_type( SC_FRAGMENT, false );
		TS_ASSERT_EQUALS( test_properties_->get_list_of_variants().size(), 2 );
		test_properties_->set_variant_type( "BIZARRO", true );
		TS_ASSERT_EQUALS( test_properties_->get_list_of_variants().size(), 2 );
	}


private:  // Private data /////////////////////////////////////////////////////
	core::chemical::ResiduePropertiesOP test_properties_;

};  // class ResiduePropertiesTests
