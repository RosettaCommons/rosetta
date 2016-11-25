// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/task_op_schemas.cc
/// @brief  Utility functions for creating XML Schema definitions for TaskOperations
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/pack/task/operation/task_op_schemas.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace task {
namespace operation {

std::string
complex_type_name_for_task_op( std::string const & task_op_name )
{
	return "to_" + task_op_name + "Type";
}

std::string
complex_type_name_for_res_lvl_task_op( std::string const & task_op_name )
{
	return "rlto_" + task_op_name + "Type";
}

std::string
complex_type_name_for_res_filter( std::string const & res_filter_name )
{
	return "resfilter_" + res_filter_name + "Type";
}


void
task_op_schema_empty(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & task_op_name,
	std::string const & description
)
{
	using namespace utility::tag;
	XMLSchemaComplexType task_op_def;
	task_op_def.name( complex_type_name_for_task_op( task_op_name ));
	task_op_def.description( description ) ;
	append_name_and_attributes_to_complex_type( AttributeList(), task_op_def );
	xsd.add_top_level_element( task_op_def );
}

void
task_op_schema_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & task_op_name,
	utility::tag::AttributeList const & attributes,
	std::string const & description
)
{
	using namespace utility::tag;
	XMLSchemaComplexType task_op_def;
	task_op_def.name( complex_type_name_for_task_op( task_op_name ));
	task_op_def.description ( description ) ;
	append_name_and_attributes_to_complex_type( attributes, task_op_def );
	xsd.add_top_level_element( task_op_def );
}

void
res_lvl_task_op_schema_empty(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & task_op_name,
	std::string const & description
)
{
	using namespace utility::tag;
	XMLSchemaComplexType task_op_def;
	task_op_def.name( complex_type_name_for_res_lvl_task_op( task_op_name ));
	task_op_def.description( description ) ;
	xsd.add_top_level_element( task_op_def );
}

void
res_lvl_task_op_schema_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & task_op_name,
	utility::tag::AttributeList const & attributes,
	std::string const & description
)
{
	using namespace utility::tag;
	XMLSchemaComplexType task_op_def;
	task_op_def.name( complex_type_name_for_res_lvl_task_op( task_op_name ));
	task_op_def.description ( description ) ;
	for ( auto const & attribute : attributes ) {
		task_op_def.add_attribute( attribute );
	}
	xsd.add_top_level_element( task_op_def );
}

void
res_filter_schema_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & res_filter_name,
	utility::tag::AttributeList const & attributes,
	std::string const & description
)
{
	using namespace utility::tag;
	XMLSchemaComplexType res_filter_def;
	res_filter_def.name( complex_type_name_for_res_filter( res_filter_name ));
	res_filter_def.description ( description ) ;
	for ( auto const & attribute : attributes ) {
		res_filter_def.add_attribute( attribute );
	}
	xsd.add_top_level_element( res_filter_def );
}

}
}
}
}
