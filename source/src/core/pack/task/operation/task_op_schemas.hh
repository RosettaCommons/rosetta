// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/task_op_schemas.hh
/// @brief  Utility functions for creating XML Schema definitions for TaskOperations
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_task_operation_task_op_schemas_hh
#define INCLUDED_core_pack_task_operation_task_op_schemas_hh

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace task {
namespace operation {

std::string
complex_type_name_for_task_op( std::string const & task_op_name );

std::string
complex_type_name_for_res_lvl_task_op( std::string const & res_lvl_task_op_name );

std::string
complex_type_name_for_res_filter( std::string const & res_filter_name );

void
task_op_schema_empty(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & task_op_name,
	std::string const & description = "XRW TO DO"
);

void
task_op_schema_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & task_op_name,
	utility::tag::AttributeList const & attributes,
	std::string const & description = "XRW TO DO"
);

void
res_lvl_task_op_schema_empty(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & task_op_name,
	std::string const & description = "XRW TO DO"
);

void
res_lvl_task_op_schema_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & task_op_name,
	utility::tag::AttributeList const & attributes,
	std::string const & description = "XRW TO DO"
);

void
res_filter_schema_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & res_filter_name,
	utility::tag::AttributeList const & attributes,
	std::string const & description = "XRW TO DO"
);

}
}
}
}

#endif
