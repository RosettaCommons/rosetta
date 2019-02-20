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

		TR << "Setting up ResiduePropertiesTests." << std::endl;
		TR.flush();

		ChemicalManager * manager( ChemicalManager::get_instance() );
		AtomTypeSetCOP atom_types = manager->atom_type_set( FA_STANDARD );
		ElementSetCOP element_types = manager->element_set( "default" );
		MMAtomTypeSetCOP mm_atom_types = manager->mm_atom_type_set( FA_STANDARD );
		orbitals::OrbitalTypeSetCOP orbital_types = manager->orbital_type_set( FA_STANDARD );
		residue_type_ = utility::pointer::make_shared< ResidueType >( atom_types, element_types, mm_atom_types, orbital_types );
		residue_type_->name( "test_residue" );

		test_properties_ = utility::pointer::make_shared< ResidueProperties >( residue_type_.get() );

		TR << "Finished setting up ResiduePropertiesTests." << std::endl;
		TR.flush();
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
		TS_ASSERT( e.msg().find(expected_output) != std::string::npos );
	}

	// Confirm that property-related methods function properly.
	void test_properties()
	{
		TR << "Testing property-related methods of ResidueProperties."  << std::endl;
		TR << "Ensuring that LIPID property is absent (by enum)." << std::endl; TR.flush();
		TS_ASSERT( ! test_properties_->has_property( LIPID ) );
		TR << "Ensuring that \"LIPID\" property is absent (by string)." << std::endl; TR.flush();
		TS_ASSERT( ! test_properties_->has_property( "LIPID" ) );
		TR << "Ensuring that \"FAT\" property is absent (by string)." << std::endl; TR.flush();
		TS_ASSERT( ! test_properties_->has_property( "FAT" ) );


		TR << "Setting property LIPID to true (by enum)." << std::endl; TR.flush();
		test_properties_->set_property( LIPID, true );
		set_throw_on_next_assertion_failure();
		try {
			TR << "Trying to set property \"FAT\" to true (by string).  This should throw an exception." << std::endl; TR.flush();
			test_properties_->set_property( "FAT", true );
		} catch (Exception const & e )  {
			exception_message_matches( e,
				"ERROR: Error in core::chemical::ResidueProperties::set_property(): Rosetta does not recognize the property \"FAT\".  Has it been added to the \"general_properties.list\" file?");
		}

		TR << "Confirming that property LIPID is added (by enum)." << std::endl; TR.flush();
		TS_ASSERT( test_properties_->has_property( LIPID ) );
		TR << "Confirming that property \"LIPID\" is added (by string)." << std::endl; TR.flush();
		TS_ASSERT( test_properties_->has_property( "LIPID") );
		TR << "Confirming that property \"FAT\" is still absent (by string)." << std::endl; TR.flush();
		TS_ASSERT( ! test_properties_->has_property( "FAT") );

		TR << "Getting list of properties, and confirming that \"LIPID\" is the only entry." << std::endl; TR.flush();
		TS_ASSERT_EQUALS( test_properties_->get_list_of_properties()[ 1 ], "LIPID" );
		TS_ASSERT_EQUALS( test_properties_->get_list_of_properties().size(), 1 );

		TR << "Setting \"LIPID\" to true again (by string), and confirming that \"LIPID\" is still the only entry." << std::endl; TR.flush();
		test_properties_->set_property( "LIPID", true );
		TS_ASSERT_EQUALS( test_properties_->get_list_of_properties().size(), 1 );

		TR << "Setting LIPID to false again (by enum), and confirming that the properties list is now empty." << std::endl; TR.flush();
		test_properties_->set_property( LIPID, false );
		TS_ASSERT_EQUALS( test_properties_->get_list_of_properties().size(), 0 );

		TR << "Test completed!" << std::endl; TR.flush();
	}

	/// @brief Ensure that the get_property_from_string() and get_string_from_property() functions work.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void test_property_enums() {
		using namespace core::chemical;
		TS_ASSERT_EQUALS( ResidueProperties::get_property_from_string( "BOGUS_PROPERTY_THAT_DOESN'T_EXIST" ), NO_PROPERTY );
		TS_ASSERT_EQUALS( ResidueProperties::get_property_from_string( "AROMATIC" ), AROMATIC );
		TS_ASSERT_EQUALS( ResidueProperties::get_property_from_string( "D_AA" ), D_AA );
		TS_ASSERT_EQUALS( ResidueProperties::get_property_from_string( "POLAR" ), POLAR );
		TS_ASSERT_EQUALS( ResidueProperties::get_property_from_string( "CHARGED" ), CHARGED );

		TS_ASSERT_EQUALS( ResidueProperties::get_string_from_property( HYDROPHOBIC ), "HYDROPHOBIC" );
		TS_ASSERT_EQUALS( ResidueProperties::get_string_from_property( ALIPHATIC ), "ALIPHATIC" );
		TS_ASSERT_EQUALS( ResidueProperties::get_string_from_property( L_AA ), "L_AA" );
		TS_ASSERT_EQUALS( ResidueProperties::get_string_from_property( CANONICAL_AA ), "CANONICAL_AA" );
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
