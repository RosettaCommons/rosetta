// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResFilterFactory.cc
/// @brief
/// @author ashworth

// Unit headers
#include <core/pack/task/operation/ResFilterFactory.hh>
#include <core/pack/task/operation/ResFilterCreator.hh>

// Package headers
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/operation/ResFilter.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace core {
namespace pack {
namespace task {
namespace operation {

ResFilterFactory::ResFilterFactory() = default;

ResFilterFactory::~ResFilterFactory()= default;

void
ResFilterFactory::factory_register( ResFilterCreatorOP creator )
{
	add_creator( creator );
}

/// @brief add a ResFilter prototype, using its default type name as the map key
void
ResFilterFactory::add_creator( ResFilterCreatorOP creator )
{
	runtime_assert( creator != nullptr );
	filter_creator_map_[ creator->keyname() ] = creator;
}

bool ResFilterFactory::has_type( std::string const & type ) const
{
	return ( filter_creator_map_.find( type ) != filter_creator_map_.end() );
}

ResFilterOP
ResFilterFactory::newResFilter(
	std::string const & type
) const
{
	return newResFilter( type, TagCOP( TagOP( new Tag )));
}

ResFilterOP
ResFilterFactory::newResFilter(
	std::string const & type,
	TagCOP tag
) const
{
	auto iter( filter_creator_map_.find( type ) );
	if ( iter != filter_creator_map_.end() ) {
		ResFilterOP filter( iter->second->create_res_filter() );
		// parse tag if tag pointer is pointing to one
		if ( tag.get() != nullptr ) filter->parse_tag( tag );
		return filter;
	} else {
		utility_exit_with_message( type + " is not known to the ResFilterFactory. Was its ResFilterCreator class registered at initialization?" );
		return nullptr;
	}
}

/// @details By convention, the named assigned to each of the complexTypes for ResFilter s should be
/// what is returned by the function "complex_type_name_for_res_lvl_task_op" (declared in
/// core/pack/task/operation/task_op_schemas.hh) when given the argument returned by that ResFilter's
/// ResFilterCreator's keyname() function. So long as the writing of XML schema for your ResFilter
/// is accomplished by calling the functions in core/select/res_lvl_task_operations/task_op_schemas.hh, then
/// this should happen automatically.
void ResFilterFactory::define_res_filter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	try {
		utility::tag::define_xml_schema_group(
			filter_creator_map_,
			res_filter_xml_schema_group_name(),
			& complex_type_name_for_res_filter,
			xsd );
	} catch ( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Could not generate an XML Schema for ResFilters from ResFilterFactory; offending class"
			" must call core::pack::task::operation::complex_type_name_for_res_filter when defining"
			" its XML Schema\n" + e.msg() );
	}
}

std::string ResFilterFactory::res_filter_xml_schema_group_name()
{
	return "res_filter";
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
