// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    utility/xsd_util/util.hh
/// @brief   General utility functions for parsing an XSD without needing access to Rosetta's higher
/// libraries.  Useful for UI development.
/// @author  Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_utility_xsd_util_util_hh
#define INCLUDED_utility_xsd_util_util_hh

#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>

// LibXML includes
#include <libxml/parser.h>
#include <libxml/tree.h>

// Platform includes
#include <platform/types.hh>

namespace utility {
namespace xsd_util {

/// @brief Given an xsd for RosettaScripts, get all of the names of mover elements.
/// @details The mover_names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void get_mover_names_from_xsd( std::string const &xsd, utility::vector1< std::string > &mover_names_out );

/// @brief Given an xsd for RosettaScripts, get all of the names of filter elements.
/// @details The filter_names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void get_filter_names_from_xsd( std::string const &xsd, utility::vector1< std::string > &filter_names_out );

/// @brief Given an xsd for RosettaScripts, get all of the names of task operation elements.
/// @details The taskop_names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void get_task_operation_names_from_xsd( std::string const &xsd, utility::vector1< std::string > &taskop_names_out );

/// @brief Given an xsd for RosettaScripts, get all of the names of residue-level task operation elements.
/// @details The rltaskop_names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void get_residue_level_task_operation_names_from_xsd( std::string const &xsd, utility::vector1< std::string > &rltaskop_names_out );

/// @brief Given an xsd for RosettaScripts, get all of the names of residue selector elements.
/// @details The residue_selector_names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void get_residue_selector_names_from_xsd( std::string const &xsd, utility::vector1< std::string > &residue_selector_names_out );

/// @brief Given an xsd for RosettaScripts, get all of the names of a particular type of Rosetta module.
/// @details The names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void get_module_names_from_xsd( std::string const &xsd, std::string const &module_name, utility::vector1< std::string > &names_out );

/// @brief Given a node in an XML tree, get the value of an option.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string get_node_option( xmlNode* node, std::string const &optionname );

/// @brief Go through an XSD and generate a human-readable description of a particular Rosetta module.
/// @details Calls generate_human_readable_recursive.
/// @note If component_name and component_type are provided, then summary information is only returned for that particular
/// object type.  For example, component_name="DeclareBond" and component_type="mover" would return information on the DeclareBond
/// mover.  Also note, this function uses raw pointers, unfortunately -- the libxml2 functions that I'm calling require it.
/// As such, no premature return statements should be added prior to the xmlFreeDoc and xmlCleanupParser statements at the end of the function.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string generate_human_readable_summary( std::string const &xsd, std::string const &component_name, std::string const &component_type );

/// @brief Go through an XSD XML and generate a human-readable version, recursively.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void generate_human_readable_recursive( xmlNode* rootnode, std::stringstream &description, std::stringstream &tags, std::stringstream &options, platform::Size const level, std::string const &complextype, std::string const &tag_name_to_print  );

/// @brief Given a node in an XML tree representing a tag, get the string that starts the tag.  (e.g. "DeclareBond" in "<DeclareBond .... />").
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string get_tag_name( xmlNode* node );

/// @brief Given a node in an XML tree representing a tag, get all options in that tag, write their descriptions to "options",
/// and write their names and types to "tags".
/// @details Heads the section in "options" with the tag name.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void output_all_tag_options( xmlNode* node, std::string const &tagname, platform::Size const level, std::stringstream &tags, std::stringstream &options);

/// @brief Given a documentation node in an XML tree, parse out the documentation and print it.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void generate_human_readable_documentation( xmlNode* doc_node, std::stringstream &outstream );

/// @brief Given an XSD simple type (e.g. "rosetta_bool"), return a more human-readable type (e.g. "bool").
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string get_type_name( std::string const& xsd_type );

} //xsd_util
} //xsd_util

#endif
