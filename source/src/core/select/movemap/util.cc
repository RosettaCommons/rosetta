// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/movemap/util.cc
/// @brief  Utility functions for the move map factory class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/movemap/util.hh>

// Package headers
#include <core/select/movemap/MoveMapFactory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <sstream>

static basic::Tracer TR( "core.select.movemap.util" );

namespace core {
namespace select {
namespace movemap {

std::string
movemap_factory_category()
{
	return "MoveMapFactory";
}

std::string
default_movemap_factory_attribute_name() {
	return "movemap_factory";
}


/// @details Looks for the attribute with the name given by
/// default_movemap_factory_attribute_name() in the input tag;
///    If that option isn't found, returns a NULL ptr
///    If that option is found, it calls get_movemap_factory()
MoveMapFactoryOP
parse_movemap_factory(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data
)
{
	return parse_movemap_factory( tag, data, default_movemap_factory_attribute_name() );
}


/// @details Looks for the attribute with name attrname
/// in tag.
///    If that option isn't found, it returns a NULL ptr
///    If that option is found, it calls get_movemap_factory()
MoveMapFactoryOP
parse_movemap_factory(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	std::string const & attribute_name
)
{
	std::string const factoryname = tag->getOption< std::string >( attribute_name, "" );
	if ( factoryname.empty() ) {
		return MoveMapFactoryOP();
	}
	return get_movemap_factory( factoryname, data );
}

/// @brief Companion function for parse_movemap_factory.
/// This uses the default movemap factory attribute name
void
attributes_for_parse_movemap_factory_default_attr_name(
	utility::tag::AttributeList & attlist,
	std::string const & documentation_string
)
{
	attributes_for_parse_movemap_factory(attlist, default_movemap_factory_attribute_name(), documentation_string);
}

/// @brief Companion function for parse_movemap where a non-default attribute name
/// is used (e.g. you might want two move-map factories given by attributes "mmf1" and "mmf2")
void
attributes_for_parse_movemap_factory(
	utility::tag::AttributeList & attlist,
	std::string const & attribute_name,
	std::string const & documentation_string
)
{
	using namespace utility::tag;
	attlist + XMLSchemaAttribute( attribute_name, xs_string, documentation_string == "" ? "The name of the already defined MoveMapFactory that will be used by this object" : documentation_string );

}

/// @brief Companion function for parse_movemap to be used when it is unacceptible
/// for the parse_movemap function to return a null pointer
void
attributes_for_parse_movemap_factory_when_required_default_name(
	utility::tag::AttributeList & attlist,
	std::string const & documentation_string
)
{
	attributes_for_parse_movemap_factory_when_required( attlist,
		default_movemap_factory_attribute_name(), documentation_string );

}

/// @brief Companion function for parse_movemap to be used when it is unacceptible
/// for the parse_movemap function to return a null pointer
void
attributes_for_parse_movemap_factory_when_required(
	utility::tag::AttributeList & attlist,
	std::string const & attribute_name,
	std::string const & documentation_string
)
{
	using namespace utility::tag;
	attlist + XMLSchemaAttribute::required_attribute( attribute_name, xs_string, documentation_string == "" ? "The name of the already-defined MoveMapFactory that will be used by this object" : documentation_string );

}

MoveMapFactoryOP
get_movemap_factory( std::string const & factory_name, basic::datacache::DataMap const & data )
{
	core::select::movemap::MoveMapFactoryOP factory;
	try {
		factory = data.get_ptr< core::select::movemap::MoveMapFactory >( movemap_factory_category(), factory_name );
	} catch ( utility::excn::Exception & e ) {
		std::stringstream error_msg;
		error_msg << "Failed to find MoveMapFactory named '" << factory_name << "' in the DataMap.\n";
		error_msg << e.msg();
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
	}
	debug_assert( factory );
	TR << "Found MoveMapFactory " << factory_name << std::endl;
	return factory;
}


}
}
}
