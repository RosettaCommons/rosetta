// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/movemap/util.hh
/// @brief  Utility functions for the MoveMapFactory class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_movemap_util_HH
#define INCLUDED_core_select_movemap_util_HH

// Core headers
#include <core/select/movemap/MoveMapFactory.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>
#include <list>
#include <algorithm>

namespace core {
namespace select {
namespace movemap {

/// @brief Category in the DataMap in which the utility functions below
/// will go hunting for a particular MoveMapFactory
std::string
movemap_factory_category();

std::string
default_movemap_factory_attribute_name();

/// @brief returns a movemap factory given a tag and datamap
MoveMapFactoryOP
parse_movemap_factory(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data
);

/// @brief returns a movemap factory given a tag and datamap
MoveMapFactoryOP
parse_movemap_factory(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	std::string const & attribute_name
);

/// @brief Companion function for parse_movemap_factory.
/// This uses the default movemap factory attribute name
/// If no documentation string is provided, a default documentation
/// string is used.
void
attributes_for_parse_movemap_factory_default_attr_name(
	utility::tag::AttributeList & attlist,
	std::string const & documentation_string = ""
);

/// @brief Companion function for parse_movemap where a non-default attribute name
/// is used (e.g. you might want two move-map factories given by attributes "mmf1" and "mmf2")
/// If no documentation string is provided, a default documentation
/// string is used.
void
attributes_for_parse_movemap_factory(
	utility::tag::AttributeList & attlist,
	std::string const & attribute_name,
	std::string const & documentation_string = ""
);

/// @brief Companion function for parse_movemap to be used when it is unacceptible
/// for the parse_movemap function to return a null pointer
/// If no documentation string is provided, a default documentation
/// string is used.
void
attributes_for_parse_movemap_factory_when_required_default_name(
	utility::tag::AttributeList & attlist,
	std::string const & documentation_string = ""
);

/// @brief Companion function for parse_movemap to be used when it is unacceptible
/// for the parse_movemap function to return a null pointer
/// If no documentation string is provided, a default documentation
/// string is used.
void
attributes_for_parse_movemap_factory_when_required(
	utility::tag::AttributeList & attlist,
	std::string const & attribute_name,
	std::string const & documentation_string = ""
);

/// @brief returns a MoveMapFactory given the factory's name and the datamap holding it
/// @details Looks for the factory in the datamap
///          Returns a ptr to it
/// @throws utility::excn::EXCN_Msg_Exception if selector is not found in datamap
MoveMapFactoryOP
get_movemap_factory( std::string const & factory_name, basic::datacache::DataMap const & data );


}
}
}

#endif
