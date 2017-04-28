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

// C++ headers
#include <tuple>

namespace utility {
namespace xsd_util {

/// @brief Enum for the types of Rosetta modules covered by the XSD that a UI would want to be able to
/// provide an interface for.
/// @details If you add to this list, update the get_module_type_as_it_appears_in_xsd() function and the
/// get_specific_module_name_as_it_appears_in_xsd() function.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
enum RosettaModuleType {
	RMT_mover=1, //Keep this first
	RMT_filter,
	RMT_task_operation,
	RMT_res_lvl_task_operation,
	RMT_residue_selector,
	RMT_sub_element,
	RMT_invalid, //Keep this second-to-last
	RMT_end_of_list = RMT_invalid //Keep this last.
};

/// @brief Enum for the types of options (boolean, task operation, string, integer, etc.) that a Rosetta module might have.
/// @details This is provided for the convenience of user interface code that will want to know whether a particular option is
/// a boolean, an integer, etc.  It should NOT be confused with the enums in XMLSchemaGeneration.hh and XMLSchemaGeneration.fwd.hh,
/// which are a more comprehensive and detailed list of the possible types of data.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
enum RosettaModuleOptionType {
	RMOT_rosetta_bool=1, //Keep this first
	RMOT_integer,
	RMOT_nonnegative_integer,
	RMOT_real,
	RMOT_task_operation, //A single task operation
	RMOT_task_operation_list, //A list of task operations
	RMOT_residue_selector, //A single residue selector
	RMOT_residue_selector_list, //A list of residue selectors
	RMOT_general_string,
	RMOT_invalid, //Keep this second-to-last
	RMOT_end_of_list = RMOT_invalid
};

/// @brief Given the name of a tag (tag_name) and a list of options (as pairs of option name, option value), generate
/// the tag.
/// @details If terminal_slash is true, a slash is added at the end of the tag.  If triangular_brackets is true, everything is enclosed in
/// triangular brackets.
void generate_tag_given_options( std::stringstream &stringstream_out, std::string const &tag_name, utility::vector1< std::pair< std::string, std::string > > const &option_list, bool const terminal_slash, bool const triangular_brackets );

/// @brief Given a type of module (mover, filter, etc.), get its name as it appears in the RosettaScripts xsd (e.g. "mover", "filter", etc.).
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string get_module_type_as_it_appears_in_xsd( RosettaModuleType const module_type );

/// @brief Given a name of a specific module (e.g. "FastDesign") of a particular type (e.g. mover), return the name
/// as it appears in the XSD (e.g. "mover_FastDesign_type").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string get_specific_module_name_as_it_appears_in_xsd( std::string const &module_name, RosettaModuleType const module_type );

/// @brief Given an option type that appears in the xsd (e.g. "xs:string", "rosetta_bool", etc.), convert it to an enum.
/// @details Note that for certain string types, this function looks at the specific name to determine whether this is a task operation/residue selector list.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
RosettaModuleOptionType get_option_type_from_name( std::string const &name, std::string const &type );

/// @brief Given an xsd for RosettaScripts and the name of a module of a given type, get its
/// first-level child options (as a vector of pairs of <RosettaModuleOptionType, option name>).
/// @details The option_list_out vector is cleared and repopulated by this operation.  The tuple is (type, name, description).
/// If module_type is utility::xsd_util::RMT_sub_element, then we search for the sub-element complextype by its complextype name.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
get_rosetta_module_options_from_xsd(
	std::string const &xsd,
	RosettaModuleType const module_type,
	std::string const &module_name,
	std::string &module_description_out,
	utility::vector1< std::tuple< RosettaModuleOptionType, std::string, std::string > > &option_list_out,
	utility::vector1< std::tuple< std::string, std::string, platform::Size, platform::Size > > &sub_element_information_out
);

/// @brief Given an xsd for Rosettascripts and the name of a particular ComplexType element or sub-element, get a
/// pointer to its XMLnode.
/// @details The module_node pointer is set by this operation.  Returns "false" if not found.  No recusion limit.
/// @note Yes, this accepts a pointer to a pointer.  The address stored in the pointer is modified by this function.
bool recursively_find_module_in_xsd( xmlNode* const parent_node, std::string const & module_name, xmlNode** module_node );

/// @brief Given a "xs:choice" element in an XSD, get the minOccurs and maxOccurs options and store them in unsigned integers.
/// @details Stores 0 for "unbounded" or for anything that can't be interpreted as a nonnegative integer.
/// @param[in] choice_node A pointer ot an xmlNode representing an <xs:choice ... /> tag.
/// @param[out] choice_min The value parsed from minOccurs.
/// @param[out] choice_max The value parsed from maxOccurs.
void get_choice_min_and_max( xmlNode* choice_node, platform::Size & choice_min, platform::Size & choice_max );

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
