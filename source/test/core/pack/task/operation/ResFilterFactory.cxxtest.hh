// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/select/residue_selector/OperateOnResidueSubset.cxxtest.hh
/// @brief  test suite for core::select::residue_selector::OperateOnResidueSubset
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <core/pack/task/operation/ResFilterFactory.hh>
#include <core/pack/task/operation/ResFilter.hh>
#include <core/pack/task/operation/ResFilterCreator.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/XMLSchemaValidation.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>

using namespace core::pack::task::operation;
using namespace utility::tag;

class DummyFilterCreator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const { return ResFilterOP( 0 ); }
	virtual std::string keyname() const { return "DummyResFilter"; }
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
	{
		res_filter_schema_w_attributes( xsd, keyname(), AttributeList() );
	}

};

class DummyFilter2Creator : public ResFilterCreator {
public:
	virtual ResFilterOP create_res_filter() const { return ResFilterOP( 0 ); }
	virtual std::string keyname() const { return "DummyResFilter2"; }
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
	{
		XMLSchemaComplexType ct;

		// this is the error: the ComplexType name needs to come from
		// the complex_type_name_for_res_filter function
		ct.name( keyname() );

		xsd.add_top_level_element( ct );
	}

};

class ResFilterFactoryTests : public CxxTest::TestSuite {
private:
	bool dummy_initialized_;
public:

	// ctor is called exactly once per unit test suite execution (and once per execution of core.test, for that matter)
	// register the dummy creator here to avoid it happening multiple times.
	ResFilterFactoryTests() : dummy_initialized_( false ) {}

	void setUp() {
		if ( !dummy_initialized_ ) {
			dummy_initialized_ = true;
			ResFilterFactory::get_instance()->factory_register( ResFilterCreatorOP( new DummyFilterCreator ));
		}
		core_init();
	}

	void test_res_filter_factory_xsd_creation() {
		XMLSchemaDefinition xsd;
		ResFilterFactory::get_instance()->define_res_filter_xml_schema( xsd );
		TS_ASSERT( xsd.has_top_level_element( complex_type_name_for_res_filter( "DummyResFilter" )));
		//std::cout << "Res filter factory XSD:\n" << xsd.full_definition() << std::endl;
	}

	void test_all_res_filters_define_valid_xsds() {
		XMLSchemaDefinition xsd;
		ResFilterFactory::get_instance()->define_res_filter_xml_schema( xsd );

		XMLSchemaModelGroupOP rosetta_scripts_seq( new XMLSchemaModelGroup( xsmgt_sequence ));

		XMLSchemaComplexTypeOP rosetta_scripts_type( new XMLSchemaComplexType );
		rosetta_scripts_type->set_model_group( rosetta_scripts_seq );

		XMLSchemaElementOP rosetta_scripts_element( new XMLSchemaElement );
		rosetta_scripts_element->name( "ROSETTASCRIPTS" );
		rosetta_scripts_element->element_type_def( rosetta_scripts_type );

		xsd.add_top_level_element( *rosetta_scripts_element );

		//std::cout << "xsd for res filters:\n" << xsd.full_definition() << std::endl;
		XMLValidationOutput output = test_if_schema_is_valid( xsd.full_definition() );
		TS_ASSERT( output.valid() );
		if ( ! output.valid() ) {
			std::cout << output.error_messages() << std::endl;
		}
	}


	void test_res_filter_factory_w_bad_xsd() {
		ResFilterFactory::get_instance()->factory_register( ResFilterCreatorOP( new DummyFilter2Creator ));

		try {
			XMLSchemaDefinition xsd;
			ResFilterFactory::get_instance()->define_res_filter_xml_schema( xsd );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::string expected_message =
				"Could not generate an XML Schema for ResFilters from ResFilterFactory; offending class"
	      " must call core::pack::task::operation::complex_type_name_for_res_filter when defining"
	      " its XML Schema\ndefine_xml_schema_group: failed to detect a complex type of name \"" +
				complex_type_name_for_res_filter( "DummyResFilter2" ) + "\" for \"DummyResFilter2\"\n";

			TS_ASSERT_EQUALS( e.msg(), expected_message );

			//std::cout << e.msg() << "\n" << expected_message << std::endl;
		}
	}

};
