// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  ResidueProperties.cxxtest.hh
/// @brief   Test suite for ResidueProperties
/// @author  Labonte <JWLabonte@jhu.edu>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>

// Package headers
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/MMAtomType.hh>

// Utility header
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// C++ header
#include <map>

static basic::Tracer TR("core.chemical.ResidueProperties.cxxtest");

using namespace core::chemical;
using namespace utility::excn;


class ResiduePropertiesTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init();

		ChemicalManager * manager( ChemicalManager::get_instance() );
		AtomTypeSetCOP atom_types = manager->atom_type_set( FA_STANDARD );
		ElementSetCOP element_types = manager->element_set( "default" );
		MMAtomTypeSetCOP mm_atom_types = manager->mm_atom_type_set( FA_STANDARD );
		orbitals::OrbitalTypeSetCOP orbital_types = manager->orbital_type_set( FA_STANDARD );
		residue_type_ = ResidueTypeOP( new ResidueType( atom_types, element_types, mm_atom_types, orbital_types ) );
		residue_type_->name( "test_residue" );

		test_properties_ = core::chemical::ResiduePropertiesOP( new ResidueProperties( residue_type_.get() ) );
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////

	void exception_message_matches( utility::excn::Exception const & e, std::string const & expected_output )
	{
		//std::string msg = e.msg();
		//std::vector< std::string > msg_lines = utility::split_by_newlines( msg );
		//TS_ASSERT_EQUALS( msg_lines.size(), 3 );
		//TS_ASSERT_EQUALS( msg_lines[1], expected_output );
		//std::cout << "_____________" << e.msg() << std::endl;
		TS_ASSERT( e.msg().find(expected_output) != std::string::npos );
	}

	// Confirm that property-related methods function properly.
	void test_properties()
	{
		TR << "Testing property-related methods of ResidueProperties."  << std::endl;
		TS_ASSERT( ! test_properties_->has_property( LIPID ) );
		TS_ASSERT( ! test_properties_->has_property( "LIPID" ) );
		TS_ASSERT( ! test_properties_->has_property( "FAT" ) );

		test_properties_->set_property( LIPID, true );
		set_throw_on_next_assertion_failure();
		try {
			test_properties_->set_property( "FAT", true );
		} catch (Exception const & e )  {
			exception_message_matches( e,
				"Rosetta does not recognize the property: FAT; has it been added to general_properties.list?" );
		}

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
		TR << "Testing variant-related methods of ResidueProperties."  << std::endl;
		TS_ASSERT( ! test_properties_->is_variant_type( SC_FRAGMENT ) );
		TS_ASSERT( ! test_properties_->is_variant_type( "SC_FRAGMENT" ) );
		TS_ASSERT( ! test_properties_->is_variant_type( "BIZARRO" ) );

		test_properties_->set_variant_type( SC_FRAGMENT, true );

		set_throw_on_next_assertion_failure();
		try {
			test_properties_->set_variant_type( "BIZARRO", true );
		} catch (Exception const & e )  {
			exception_message_matches( e,
				"Rosetta does not recognize the variant: BIZARRO; "
				"has it been added to variant_types.list?" );
		}

		test_properties_->enable_custom_variant_types();
		TS_ASSERT( test_properties_->has_custom_variant_types() );
		test_properties_->set_variant_type( "WHACKY", true );
		test_properties_->set_variant_type( "BIZARRO", true );
		test_properties_->set_variant_type( "WHACKY", false );

		set_throw_on_next_assertion_failure();
		try {
			test_properties_->set_variant_type( "GNARLY", false );
		} catch (Exception const & e )  {
			exception_message_matches( e,
				"Rosetta does not recognize the custom variant GNARLY in test_residue" );
		}

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
	core::chemical::ResidueTypeOP residue_type_;
	core::chemical::ResiduePropertiesOP test_properties_;

};  // class ResiduePropertiesTests
