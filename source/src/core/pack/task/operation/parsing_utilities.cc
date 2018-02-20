// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pack/task/operation/parsing_utilities.cc
/// @brief  Utility functions for extracting packer-related objects out of the DataMap
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/task/operation/parsing_utilities.hh>

// Project Headers
#include <core/pack/task/operation/ResLvlTaskOperation.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utillity Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

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
rlto_datamap_category() {
	// Note that construction of static data is threadsafe in C++11
	static std::string const cat( "residue_level_task_operations" );
	return cat;
}

void
parse_residue_level_task_operations(
	utility::tag::TagCOP const & tag,
	basic::datacache::DataMap & dm,
	std::list< core::pack::task::operation::ResLvlTaskOperationOP > & rltos
)
{
	if ( ! tag->hasOption( "residue_level_operations" ) ) return;
	std::string const rlto_cs_list = tag->getOption< std::string >( "residue_level_operations" );
	utility::vector1< std::string > rlto_names_vector = utility::string_split( rlto_cs_list, ',' );
	for ( auto const & rlto_name : rlto_names_vector ) {
		if ( dm.has( rlto_datamap_category(), rlto_name ) ) {
			rltos.push_back( dm.get_ptr< ResLvlTaskOperation >( rlto_datamap_category(), rlto_name ) );
		} else {
			std::ostringstream oss;
			oss << "Could not locate residue-level-task operation named \"" << rlto_name << "\";\n";
			oss << "The residue-level-task operations that have been defined are:\n";
			for ( auto const & decl_rlto_pair : dm[ rlto_datamap_category() ] ) {
				oss << "   " << decl_rlto_pair.first << "\n";
			}
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}
	}
}

void
parse_residue_level_task_operations(
	utility::tag::TagCOP const & tag,
	basic::datacache::DataMap & dm,
	std::list< core::pack::task::operation::ResLvlTaskOperationCOP > & rltos
)
{
	if ( ! tag->hasOption( "residue_level_operations" ) ) return;
	std::string const rlto_cs_list = tag->getOption< std::string >( "residue_level_operations" );
	utility::vector1< std::string > rlto_names_vector = utility::string_split( rlto_cs_list, ',' );
	for ( auto const & rlto_name : rlto_names_vector ) {
		if ( dm.has( rlto_datamap_category(), rlto_name ) ) {
			rltos.push_back( dm.get_ptr< ResLvlTaskOperation >( rlto_datamap_category(), rlto_name ) );
		} else {
			std::ostringstream oss;
			oss << "Could not locate residue-level-task operation named \"" << rlto_name << "\";\n";
			oss << "The residue-level-task operations that have been defined are:\n";
			for ( auto const & decl_rlto_pair : dm[ rlto_datamap_category() ] ) {
				oss << "   " << decl_rlto_pair.first << "\n";
			}
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}
	}
}

///////////////////// Attributes /////////////////////////

void
attributes_for_parse_residue_level_operations(
	utility::tag::AttributeList & attributes
)
{
	using namespace utility::tag;
	attributes +
		XMLSchemaAttribute( "residue_level_operations", xs_string, "A comma-separated list of"
		" residue-level-task operations that will be retrieved from the DataMap." );
}

}
}
}
}

