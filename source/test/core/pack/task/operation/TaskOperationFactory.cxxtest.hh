// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/OperateOnResidueSubset.cxxtest.hh
/// @brief  test suite for core::select::residue_selector::OperateOnResidueSubset
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperationCreator.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/XMLSchemaValidation.hh>
#include <utility/excn/Exceptions.hh>

// Test utility headers
#include <test/util/schema_utilities.hh>

// C++ headers
#include <string>

using namespace core::pack::task::operation;
using namespace utility::tag;

class DummyTaskOpCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const { return TaskOperationOP( 0 ); }
	virtual std::string keyname() const { return "DummyTaskOperation"; }
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
	{
		task_op_schema_w_attributes( xsd, keyname(), AttributeList() );
	}

};

class DummyTaskOp2Creator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const { return TaskOperationOP( 0 ); }
	virtual std::string keyname() const { return "DummyTaskOperation2"; }
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
	{
		XMLSchemaComplexType ct;

		// this is the error: the ComplexType name needs to come from
		// the complex_type_name_for_task_op function
		ct.name( keyname() );

		xsd.add_top_level_element( ct );
	}

};

class TaskOperationFactoryTests : public CxxTest::TestSuite {
private:
	bool dummy_initialized_;
public:

	// ctor is called exactly once per unit test suite execution (and once per execution of core.test, for that matter)
	// register the dummy creator here to avoid it happening multiple times.
	TaskOperationFactoryTests() : dummy_initialized_( false ) {}

	void setUp() {
		if ( !dummy_initialized_ ) {
			dummy_initialized_ = true;
			TaskOperationFactory::get_instance()->factory_register( TaskOperationCreatorOP( new DummyTaskOpCreator ));
		}
		core_init();
	}

	void test_all_task_operations_define_valid_xsds() {
		XMLSchemaDefinition xsd;
		TaskOperationFactory::get_instance()->define_task_op_xml_schema( xsd );

		XMLSchemaComplexTypeOP rosetta_scripts_type( new XMLSchemaComplexType );
		XMLSchemaModelGroupOP rosetta_scripts_seq( new XMLSchemaModelGroup( xsmgt_sequence ));
		XMLSchemaElementOP rosetta_scripts_element( new XMLSchemaElement );

		XMLSchemaElementOP task_operations_element( new XMLSchemaElement );
		//XMLSchemaModelGroupOP task_operations_choice( new XMLSchemaModelGroup( xsmgt_choice ));
		XMLSchemaComplexTypeOP task_operations_type( new XMLSchemaComplexType );

		XMLSchemaModelGroupOP group_ref( new XMLSchemaModelGroup( xsmgt_group ) );
		group_ref->group_name( "task_operation" );
		group_ref->min_occurs( 0 );
		group_ref->max_occurs( xsminmax_unbounded );
		task_operations_type->name( "TASK_OPERATIONS_Type" );
		task_operations_type->set_model_group( group_ref );

		xsd.add_top_level_element( *task_operations_type );

		task_operations_element->name( "TASK_OPERATIONS" );
		task_operations_element->type_name( "TASK_OPERATIONS_Type" );
		task_operations_element->min_occurs( 0 );

		rosetta_scripts_seq->append_particle( task_operations_element );
		rosetta_scripts_type->set_model_group( rosetta_scripts_seq );

		rosetta_scripts_element->name( "ROSETTASCRIPTS" );
		rosetta_scripts_element->element_type_def( rosetta_scripts_type );

		xsd.add_top_level_element( *rosetta_scripts_element );

		//std::cout << "xsd for task operations:\n" << xsd.full_definition() << std::endl;
		XMLValidationOutput output = test_if_schema_is_valid( xsd.full_definition() );
		TS_ASSERT( output.valid() );
		if ( ! output.valid() ) {
			std::cout << output.error_messages() << std::endl;
		}
	}

	void test_task_op_factory_xsd_creation() {
		XMLSchemaDefinition xsd;
		TaskOperationFactory::get_instance()->define_task_op_xml_schema( xsd );
		TS_ASSERT( xsd.has_top_level_element( complex_type_name_for_task_op( "DummyTaskOperation" )));
		//std::cout << "Res filter factory XSD:\n" << xsd.full_definition() << std::endl;
	}


	void test_TaskOperationFactory_all_TaskOperation_complexTypes_have_descriptions() {
		ensure_all_cts_for_creators_have_documentation_strings(
			TaskOperationFactory::get_instance()->creator_map(),
			"TaskOperation",
			& complex_type_name_for_task_op );
	}

	void test_TaskOperationFactory_all_TaskOperation_all_attributes_have_descriptions() {
		TaskOperationFactory * mf = TaskOperationFactory::get_instance();
		TaskOperationFactory::TaskOperationCreatorMap const & creator_map = mf->creator_map();

		for ( auto iter : creator_map ) {
			XMLSchemaDefinition xsd;
			iter.second->provide_xml_schema( xsd );
			std::string full_def = xsd.full_definition();
			//std::cout << "full def: " << iter.first << "\n" << full_def << std::endl;
			TagCOP tag( Tag::create( full_def ) );
			recurse_through_subtags_for_attribute_descriptions( tag, iter.first );
		}
	}

	/// @brief From here forward, the TaskOperationFactory is "corrupted" in the sense that the DummyTaskOp2Creator
	/// is always going to cause an exception to be thrown.
	/// This should be the last test in this test suite
	void test_task_op_factory_w_bad_xsd() {
		TaskOperationFactory::get_instance()->factory_register( TaskOperationCreatorOP( new DummyTaskOp2Creator ));
		try {
			XMLSchemaDefinition xsd;
			TaskOperationFactory::get_instance()->define_task_op_xml_schema( xsd );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::string expected_message =
				"Could not generate an XML Schema for TaskOperations from TaskOperationFactory\ndefine_xml_schema_group: failed to detect a complex type of name \"" +
				complex_type_name_for_task_op( "DummyTaskOperation2" ) + "\" for \"DummyTaskOperation2\"\n";

			TS_ASSERT_EQUALS( e.msg(), expected_message );

			//std::cout << e.msg() << "\n" << expected_message << std::endl;
		}
	}


};
