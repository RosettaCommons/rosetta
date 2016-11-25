// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/ResidueSelectorFactory.cxxtest.hh
/// @brief test suite for basic::resource_manager::ResourceOptionsFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreator.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/XMLSchemaValidation.hh>
#include <utility/excn/Exceptions.hh>

// Test utility headers
#include <test/util/schema_utilities.hh>

// C++ headers
#include <string>

static THREAD_LOCAL basic::Tracer TR("core.select.residue_selector.ResidueSelectorFactory.cxxtest.hh");

using namespace utility::tag;
using namespace core::select::residue_selector;

class DummyResidueSelector : public ResidueSelector {
public:
	DummyResidueSelector() {}

	DummyResidueSelector( DummyResidueSelector const & ) {}

	ResidueSelectorOP clone() const { return ResidueSelectorOP( new DummyResidueSelector(*this) ); }

	virtual
	ResidueSubset apply( core::pose::Pose const & ) const
	{
		return ResidueSubset( 10, true );
	}

	virtual std::string get_name() const { return "DummyResidueSelector"; }

};

class DummyResidueSelectorCreator : public ResidueSelectorCreator {
public:
	virtual
	ResidueSelectorOP
	create_residue_selector() const {
		return ResidueSelectorOP( new DummyResidueSelector );
	}

	virtual
	std::string keyname() const { return "DummyResidueSelector"; }

	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
		core::select::residue_selector::xsd_type_definition_w_attributes( xsd, keyname(), "testing 1 2 3", utility::tag::AttributeList() );
	}
};


class ResidueSelectorFactoryTests : public CxxTest::TestSuite {

public:

	void setUp() { core_init(); }

	// @brief make sure that when we register a residue selector, we can later get it back
	void test_register_one_creator_with_ResidueSelectorFactory() {
		ResidueSelectorFactory * factory = ResidueSelectorFactory::get_instance();
		factory->set_throw_on_double_registration();
		factory->factory_register( ResidueSelectorCreatorOP( new DummyResidueSelectorCreator ) );
		basic::datacache::DataMap dm;
		ResidueSelectorOP selector = factory->new_residue_selector( "DummyResidueSelector", 0, dm );
		TS_ASSERT( selector.get() ); // make sure we got back a non-null pointer
		DummyResidueSelector * dselector = dynamic_cast< DummyResidueSelector * > ( selector.get() );
		TS_ASSERT( dselector ); // make sure we got back the right residue selector kind
	}

	void test_all_residue_selectors_define_valid_xsds() {
		XMLSchemaDefinition xsd;
		ResidueSelectorFactory::get_instance()->define_residue_selector_xml_schema( xsd );

		XMLSchemaModelGroupOP rosetta_scripts_seq( new XMLSchemaModelGroup( xsmgt_sequence ));

		XMLSchemaComplexTypeOP rosetta_scripts_type( new XMLSchemaComplexType );
		rosetta_scripts_type->set_model_group( rosetta_scripts_seq );

		XMLSchemaElementOP rosetta_scripts_element( new XMLSchemaElement );
		rosetta_scripts_element->name( "ROSETTASCRIPTS" );
		rosetta_scripts_element->element_type_def( rosetta_scripts_type );

		xsd.add_top_level_element( *rosetta_scripts_element );

		//std::cout << "xsd for residue selector:\n" << xsd.full_definition() << std::endl;
		XMLValidationOutput output = test_if_schema_is_valid( xsd.full_definition() );
		TS_ASSERT( output.valid() );
		if ( ! output.valid() ) {
			std::cout << output.error_messages() << std::endl;
		}
	}


	/// @brief make sure that if a ResidueSelector with the name of an already-registered
	/// ResidueSelector gets registered with the factory, that the factory throws an exception
	void test_register_one_creator_twice_with_ResidueSelectorFactory() {
		try {
			ResidueSelectorFactory * factory = ResidueSelectorFactory::get_instance();
			factory->set_throw_on_double_registration();
			factory->factory_register( ResidueSelectorCreatorOP( new DummyResidueSelectorCreator ) );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_err_msg = "Factory Name Conflict: Two or more ResidueSelectorCreators registered with the name DummyResidueSelector";
			TS_ASSERT( expected_err_msg == e.msg() );

		}

	}

	// @brief Make sure that the factory will throw a meaningful error when asked for a
	// residue selector that has not been defined.
	void test_throw_on_unregistered_selector_name_ResidueSelectorFactory() {
		try {
			ResidueSelectorFactory * factory = ResidueSelectorFactory::get_instance();
			basic::datacache::DataMap dm;
			ResidueSelectorOP selector = factory->new_residue_selector( "DummyResidueSelector2", 0, dm );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_err_msg = "No ResidueSelectorCreator with the name 'DummyResidueSelector2' has been registered with the ResidueSelectorFactory";
			TS_ASSERT( expected_err_msg == e.msg() );
		}
	}


	void test_output_residue_selector_xsds() {
		utility::tag::XMLSchemaDefinition xsd;
		ResidueSelectorFactory * factory = ResidueSelectorFactory::get_instance();
		try {
			// make sure the RSF doesn't find errors with the xml schema it produces
			factory->define_residue_selector_xml_schema( xsd );
			TS_ASSERT( true );
		} catch ( ... ) {
			TS_ASSERT( false );
		}
		TR << "XSD: " << std::endl;
		TR << xsd.full_definition() << std::endl;
	}

	void test_ResidueSelectorFactory_all_ResidueSelector_complexTypes_have_descriptions() {
		ensure_all_cts_for_creators_have_documentation_strings(
			ResidueSelectorFactory::get_instance()->creator_map(),
			"ResidueSelector",
			& complex_type_name_for_residue_selector );
	}

	void test_ResidueSelectorFactory_all_ResidueSelector_all_attributes_have_descriptions() {

		ResidueSelectorFactory * rsf = ResidueSelectorFactory::get_instance();
		std::map< std::string, ResidueSelectorCreatorOP > const & creator_map = rsf->creator_map();

		for ( auto iter : creator_map ) {
			XMLSchemaDefinition xsd;
			iter.second->provide_xml_schema( xsd );
			std::string full_def = xsd.full_definition();
			//std::cout << "full def: " << iter.first << "\n" << full_def << std::endl;
			TagCOP tag( Tag::create( full_def ) );

			recurse_through_subtags_for_attribute_descriptions( tag, iter.first );
		}
	}


};
