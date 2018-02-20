// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pack/task/operation/parsing_utilities.hh
/// @brief  Utility functions for extracting packer-related objects out of the DataMap
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_task_operation_parsing_utilities_HH
#define INCLUDED_core_pack_task_operation_parsing_utilities_HH

// Unit headers

// Project Headers
#include <core/types.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utillity Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>
#include <list>


namespace core {
namespace pack {
namespace task {
namespace operation {

///////////////////////////////////////////////////////////
////////Residue Level Task Operations /////////////////////


std::string const &
rlto_datamap_category();

void
parse_residue_level_task_operations(
	utility::tag::TagCOP const & ,
	basic::datacache::DataMap &,
	std::list< core::pack::task::operation::ResLvlTaskOperationOP > &
);

void
parse_residue_level_task_operations(
	utility::tag::TagCOP const & ,
	basic::datacache::DataMap &,
	std::list< core::pack::task::operation::ResLvlTaskOperationCOP > &
);

///////////////////// Attributes /////////////////////////

void
attributes_for_parse_residue_level_operations(
	utility::tag::AttributeList & attributes
);

}
}
}
}

#endif
