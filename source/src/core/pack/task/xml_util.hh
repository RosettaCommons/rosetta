// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/xml_util.hh
/// @brief Task operation/factory utilities useful in RosettaScripts.
/// @author Sarel Fleishman (sarelf@u.washington.edu)
///         Jacob Corn (jecorn@u.washington.edu)
///         Rocco Moretti (rmoretti@u.washington.edu)
///         Eva-Maria Strauch (evas01@uw.edu)

#ifndef INCLUDED_core_pack_task_rosetta_scripts_HH
#define INCLUDED_core_pack_task_rosetta_scripts_HH

// Unit headers

// Project Headers
#include <core/types.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace task {

static const std::string TASK_OPERATIONS_TAG( "task_operations" );
static const std::string PACKER_PALETTE_TAG( "packer_palette" );
static const std::string TASK_FACTORY_TAG( "task_factory" );

/// @brief Return a list of the task operations referenced in the given tag.
utility::vector1< core::pack::task::operation::TaskOperationOP >
get_task_operations(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data
);

/// @brief Construct a TaskFactory from the task operations referenced in the
/// given tag.
core::pack::task::TaskFactoryOP
parse_task_operations(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data
);

/// @brief Construct a task factory from the task operations referenced in the
/// given comma-separated list of names.
core::pack::task::TaskFactoryOP
parse_task_operations(
	std::string const & task_list,
	basic::datacache::DataMap const & data
);

/// @brief Construct a TaskFactory by adding the task operations referenced in
/// the given tag to a task factory already present in the data map.
/// @details This allows the transfer of whole task factories on the data map.
/// This way a "base" task factory can be created, transferred on the data map,
/// and individual mover's specific task operations can be added on top.
core::pack::task::TaskFactoryOP
parse_task_operations(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap /*const*/ & data,
	core::pack::task::TaskFactoryOP & task_factory
);

/// @brief Append the 'task_operation' attribute.
/// @details "description" can be used to specify for what the TaskOperations are being used for.
void
attributes_for_parse_task_operations(
	utility::tag::AttributeList & attributes,
	std::string const & description = "",
	std::string const & packerpalette_description = ""
);

/// @brief Append the 'task_operation' and 'task_factory' attributes.
void
attributes_for_parse_task_operations_w_factory(
	utility::tag::AttributeList & attributes,
	std::string const & used_for_descr = ""
);

/// @brief Parse a "packer_palette" tag from XML, and add a PackerPalette to a TaskFactory.
/// @details Added as part of the 2016 Chemical XRW (eXtreme Rosetta Workshop), New Brunswick, NJ,
/// 25 January 2016 to 1 February 2016.
/// @param[in] tag Const owning pointer to an XML tag that may or may not contain a single "packer_palette=xxx" string.
/// @param[in] data Const reference to the data map, from which defined PackerPalettes will be sought.
/// @param[out] task_factory Non-const owning pointer to the task_factory that will have a PackerPalette added to it.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @note Moved to core on 7 Nov. 2018 by VKM, vmulligan@flatironinstitute.org.
void
set_up_packer_palette(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	core::pack::task::TaskFactoryOP & task_factory
);


}
}
}

#endif
