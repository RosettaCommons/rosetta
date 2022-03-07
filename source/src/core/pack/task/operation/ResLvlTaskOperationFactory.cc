// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResLvlTaskOperationFactory.cc
/// @brief
/// @author ashworth

#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>

#ifdef WIN32
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/ResFilter.hh>

#endif

#include <utility/exit.hh> // runtime_assert, utility_exit_with_message

#include <core/pack/task/operation/ResLvlTaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreator.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/xml_schema_group_initialization.hh>

// Basic headers
#include <basic/citation_manager/CitationManager.hh>

// C++ headers
#include <iostream>
#include <sstream>

// Boost headers

namespace core {
namespace pack {
namespace task {
namespace operation {

void
ResLvlTaskOperationFactory::factory_register( ResLvlTaskOperationCreatorOP creator )
{
	add_creator( creator );
}

ResLvlTaskOperationFactory::ResLvlTaskOperationFactory() = default;

ResLvlTaskOperationFactory::~ResLvlTaskOperationFactory()= default;

/// @brief add a ResLvlTaskOperation prototype, using its default type name as the map key
void
ResLvlTaskOperationFactory::add_creator( ResLvlTaskOperationCreatorOP rltoc )
{
	runtime_assert( rltoc != nullptr );
	rltoc_map_[ rltoc->keyname() ] = rltoc;
}

bool ResLvlTaskOperationFactory::has_type( std::string const & type ) const
{
	return ( rltoc_map_.find( type ) != rltoc_map_.end() );
}

/// @brief Get the XML schema for a given ResLvlTaskOperation.
/// @details Throws an error if the ResLvlTaskOperation is unknown to Rosetta.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
ResLvlTaskOperationFactory::provide_xml_schema(
	std::string const &task_operation_name,
	utility::tag::XMLSchemaDefinition & xsd
) const {
	auto iter( rltoc_map_.find( task_operation_name ) );
	if ( iter != rltoc_map_.end() ) {
		iter->second->provide_xml_schema( xsd );
	} else {

		utility_exit_with_message( task_operation_name + " is not known to the ResLvlTaskOperationFactory. Was its ResLvlTaskOperationCreator class registered at initialization?" );
	}
}

/// @brief Get a human-readable listing of the citations for a given residue level task operation, by taskop name.
/// @details Returns an empty string if there are no citations.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
std::string
ResLvlTaskOperationFactory::get_citation_humanreadable(
	std::string const & taskop_name
) const {
	using namespace basic::citation_manager;
	CitationCollectionList citations;
	std::string const errmsg( "Error in ResLvlTaskOperationFactory::get_citation_humanreadable(): " );

	auto const iter( rltoc_map_.find( taskop_name ) );
	runtime_assert_string_msg( iter != rltoc_map_.end(), errmsg + taskop_name + " is not known to the ResLvlTaskOperationFactory. Was its ResLvlTaskOperationCreator class registered at initialization?" );
	if ( iter->second->skip_citation_info_in_schema() ) return "";
	ResLvlTaskOperationOP rlto( iter->second->create_res_level_task_operation() );
	runtime_assert_string_msg( rlto != nullptr, errmsg + "Could not instantiate " + taskop_name + "!" );
	rlto->provide_citation_info(citations);
	if ( citations.empty() ) return "";
	std::ostringstream ss;
	ss << "References and author information for the " << taskop_name << " residue level task operation:" << std::endl;
	ss << std::endl;
	basic::citation_manager::CitationManager::get_instance()->write_all_citations_and_unpublished_author_info_from_list_to_stream( citations, ss );
	return ss.str();
}

/// @brief return new ResLvlTaskOperation by key lookup in rltoc_map_ (new ResLvlTaskOperation parses Tag if provided)
ResLvlTaskOperationOP
ResLvlTaskOperationFactory::newRLTO( std::string const & type ) const
{
	auto iter( rltoc_map_.find( type ) );
	if ( iter != rltoc_map_.end() ) {
		ResLvlTaskOperationOP rlto( iter->second->create_res_level_task_operation() );
		return rlto;
	} else {
		utility_exit_with_message( type + " is not known to the ResLvlTaskOperationFactory. Was its ResLvlTaskOperationCreator class registered at initialization?" );
		return nullptr;
	}
}

ResLvlTaskOperationOP
ResLvlTaskOperationFactory::newRLTO(
	std::string const & type,
	TagCOP tag
) const
{
	auto iter( rltoc_map_.find( type ) );
	if ( iter != rltoc_map_.end() ) {
		ResLvlTaskOperationOP rlto( iter->second->create_res_level_task_operation() );
		rlto->parse_tag( tag );
		return rlto;
	} else {
		utility_exit_with_message( type + " is not known to the ResLvlTaskOperationFactory. Was its ResLvlTaskOperationCreator class registered at initialization?" );
		return nullptr;
	}
}

/// @details By convention, the named assigned to each of the complexTypes for ResLvlTaskOperations should be
/// what is returned by the function "complex_type_name_for_res_lvl_task_op" (declared in
/// core/pack/task/operation/task_op_schemas.hh) when given the argument returned by that ResLvlTaskOperation's
/// ResLvlTaskOperationCreator's keyname() function. So long as the writing of XML schema for your ResLvlTaskOperation
/// is accomplished by calling the functions in core/select/res_lvl_task_operations/task_op_schemas.hh, then
/// this should happen automatically.
void ResLvlTaskOperationFactory::define_res_lvl_task_op_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	try {
		utility::tag::define_xml_schema_group(
			rltoc_map_,
			res_lvl_task_op_xml_schema_group_name(),
			& complex_type_name_for_res_lvl_task_op,
			xsd );
	} catch ( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Could not generate an XML Schema for ResLvlTaskOperations from ResLvlTaskOperationFactory; offending class"
			" must call core::pack::task::operation::complex_type_name_for_res_lvl_task_op when defining"
			" its XML Schema\n" + e.msg() );
	}


}

std::string ResLvlTaskOperationFactory::res_lvl_task_op_xml_schema_group_name()
{
	return "res_lvl_task_op";
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
