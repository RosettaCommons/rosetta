// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    utility/xsd_util/util.cc
/// @brief   General utility functions for parsing an XSD without needing access to Rosetta's higher
/// libraries.  Useful for UI development.
/// @author  Vikram K. Mulligan (vmullig@uw.edu).


// Unit header
#include <utility/xsd_util/util.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>

namespace utility {
namespace xsd_util {

#define MAX_HUMAN_READABLE_RECURSION_LEVELS 10 //Maximum number of levels of nested tags when generating human-readable summary of an XSD for a mover, filter, etc.

/// @brief Given the name of a tag (tag_name) and a list of options (as pairs of option name, option value), generate
/// the tag.
/// @details If terminal_slash is true, a slash is added at the end of the tag.  If triangular_brackets is true, everything is enclosed in
/// triangular brackets.
void
generate_tag_given_options(
	std::stringstream &stringstream_out,
	std::string const &tag_name,
	utility::vector1< std::pair< std::string, std::string > > const &option_list,
	bool const terminal_slash,
	bool const triangular_brackets
) {
	if ( triangular_brackets ) {
		stringstream_out << "<";
	}
	stringstream_out << tag_name;
	for ( platform::Size i=1, imax=option_list.size(); i<=imax; ++i ) {
		stringstream_out << " " << option_list[i].first << "=\"" << option_list[i].second << "\"";
	}
	if ( terminal_slash ) stringstream_out << " /";
	if ( triangular_brackets ) stringstream_out << ">";
}


/// @brief Given a type of module, get its name as it appears in the RosettaScripts xsd.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
get_module_type_as_it_appears_in_xsd(
	RosettaModuleType const module_type
) {
	switch(module_type) {
	case RMT_mover :
		return "mover";
	case RMT_filter :
		return "filter";
	case RMT_task_operation :
		return "task_operation";
	case RMT_res_lvl_task_operation :
		return "res_lvl_task_op";
	case RMT_residue_selector :
		return "residue_selector";
	case RMT_sub_element : //Not actually something that appears in the xsd, but used to refer to sub-elements of movers, filters, etc.
		break;
	default :
		break;
	}
	return "INVALID";
}

/// @brief Given a name of a specific module (e.g. "FastDesign") of a particular type (e.g. mover), return the name
/// as it appears in the XSD (e.g. "mover_FastDesign_type").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
get_specific_module_name_as_it_appears_in_xsd(
	std::string const &module_name,
	RosettaModuleType const module_type
) {
	switch(module_type) {
	case RMT_mover :
		return "mover_" + module_name + "_type";
	case RMT_filter :
		return "filter_" + module_name + "_type";
	case RMT_task_operation :
		return "to_" + module_name + "_type";
	case RMT_res_lvl_task_operation :
		return "rlto_" + module_name + "_type";
	case RMT_residue_selector :
		return "rs_" + module_name + "_type";
	case RMT_sub_element : //Not actually something that appears in the xsd, but used to refer to sub-elements of movers, filters, etc.
		break;
	default :
		break;
	}
	return "INVALID";
}

/// @brief Given an option type that appears in the xsd (e.g. "xs:string", "rosetta_bool", etc.), convert it to an enum.
/// @details Note that for certain string types, this function looks at the specific name to determine whether this is a task operation/residue selector list.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
RosettaModuleOptionType
get_option_type_from_name(
	std::string const &/*name*/,
	std::string const &type
) {
	if ( type == "xs:string" ) {
		return RMOT_general_string; //TODO: distinguish residue selectors and task operations from general strings.
	}
	if ( type == "rosetta_bool" ) {
		return RMOT_rosetta_bool;
	}
	if ( type == "xs:integer" ) {
		return RMOT_integer;
	}
	if ( type == "non_negative_integer" ) {
		return RMOT_nonnegative_integer;
	}
	if ( type == "real" ) {
		return RMOT_real;
	}

	return RMOT_general_string;
}

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
) {
	option_list_out.clear();

	std::string const xsd_module_name( module_type == RMT_sub_element ?  module_name : get_specific_module_name_as_it_appears_in_xsd(module_name, module_type) );

	xmlDoc* doc( xmlReadMemory( xsd.c_str(), xsd.length()+1, nullptr, nullptr, 0 ) ); //Must be deleted later!  Unfortunately, libxml2 likes raw pointers.

	xmlNode* rootnode( xmlDocGetRootElement(doc) );
	xmlNode* module_node( nullptr );

	//std::cout << "module_name=" << xsd_module_name << std::endl; //DELETE ME

	runtime_assert_string_msg(
		recursively_find_module_in_xsd( rootnode, xsd_module_name, &module_node ),
		"Error in utility::xsd_util::get_rosetta_module_options_from_xsd(): Could not find element with the name \"" + module_name + "\" in the XSD."
	);

	//First, get the description (annotation node):
	std::stringstream description;
	for ( xmlNode* annotation_node = module_node->children; annotation_node!=nullptr; annotation_node = annotation_node->next ) {
		if ( !strcmp(reinterpret_cast<const char*>(annotation_node->name),"annotation") ) {
			generate_human_readable_documentation( annotation_node, description );
		}
	} //looping through annotation node
	module_description_out = description.str();
	//Next, find the attribute nodes for each option that can be set:
	for ( xmlNode* attribute_node = module_node->children; attribute_node!=nullptr; attribute_node = attribute_node->next ) { //Loop through first-level sub-nodes
		if ( //Find the first-level options:
				attribute_node->type != XML_ELEMENT_NODE ||
				strcmp( reinterpret_cast<const char*>(attribute_node->name), "attribute" )
				) {
			continue;
		}
		std::string const optionname( get_node_option(attribute_node, "name") );
		RosettaModuleOptionType const optiontype( get_option_type_from_name( optionname, get_node_option(attribute_node, "type") ) );
		std::stringstream optiondescription;
		generate_human_readable_documentation( attribute_node, optiondescription );
		option_list_out.push_back( std::tuple<RosettaModuleOptionType, std::string, std::string>( optiontype, optionname, optiondescription.str() ) );
	} //Looping through attribute nodes

	//Get allowed sub-elements.
	for ( xmlNode* element_node = module_node->children; element_node != nullptr; element_node = element_node->next ) {
		if ( element_node->type == XML_ELEMENT_NODE ) {
			if ( !strcmp(reinterpret_cast<const char*>(element_node->name),"element") ) {
				std::string element_type_name( get_node_option( element_node, "type" ) );
				std::string const element_name( get_node_option( element_node, "name" ) );
				if ( element_type_name.empty() ) element_type_name = element_name;
				sub_element_information_out.push_back( std::tuple< std::string, std::string, platform::Size, platform::Size >( element_name, element_type_name, 0, 0 ) );
			} else if ( !strcmp(reinterpret_cast<const char*>(element_node->name),"choice") ) {
				platform::Size choice_min, choice_max;
				get_choice_min_and_max( element_node, choice_min, choice_max );
				for ( xmlNode* element_node2 = element_node->children; element_node2 != nullptr; element_node2 = element_node2->next ) {
					if ( element_node2->type == XML_ELEMENT_NODE && !strcmp(reinterpret_cast<const char*>(element_node2->name),"element") ) {
						std::string element_type_name( get_node_option( element_node2, "type" ) );
						std::string const element_name( get_node_option( element_node2, "name" ) );
						if ( element_type_name.empty() ) element_type_name = element_name;
						//std::cout << "Adding element_name=" << element_name << ", element_type=" << element_type_name << std::endl; //DELETE ME
						sub_element_information_out.push_back( std::tuple< std::string, std::string, platform::Size, platform::Size >( element_name, element_type_name, choice_min, choice_max ) );
					}
				}
			}
		}
	}

	//Deleting:
	xmlFreeDoc(doc); doc=nullptr;
	xmlCleanupParser();
}

/// @brief Given an xsd for Rosettascripts and the name of a particular ComplexType element or sub-element, get a
/// pointer to its XMLnode.
/// @details The module_node pointer is set by this operation.  Returns "false" if not found.  No recusion limit.
/// @note Yes, this accepts a pointer to a pointer.  The address stored in the pointer is modified by this function.
bool
recursively_find_module_in_xsd(
	xmlNode* const parent_node,
	std::string const & module_name,
	xmlNode** module_node
) {
	if ( parent_node->type == XML_ELEMENT_NODE ) {
		if ( !strcmp(reinterpret_cast<const char*>(parent_node->name),"complexType") && !get_node_option(parent_node, "name").compare( module_name ) ) {
			*module_node = parent_node;
			//std::cout << "Found " << module_name << " complexType." << std::endl; //DELETE ME.
			return true;
		}
		if ( !strcmp(reinterpret_cast<const char*>(parent_node->name),"element") && !get_node_option(parent_node, "name").compare( module_name ) ) {
			//std::cout << "Found " << module_name << " element." << std::endl; //DELETE ME.
			for ( xmlNode* sub_node(parent_node->children); sub_node != nullptr; sub_node = sub_node->next ) {
				if ( !strcmp(reinterpret_cast<const char*>(sub_node->name),"complexType") ) {
					*module_node = sub_node;
					//std::cout << "Found sub-node for complexType!" << std::endl; //DELETE ME.
					return true;
				}
			}
			return false;
		}
		for ( xmlNode* sub_node( parent_node->children ); sub_node != nullptr; sub_node = sub_node->next ) {
			if ( recursively_find_module_in_xsd( sub_node, module_name, module_node ) ) return true;
		}
	}
	return false;
}

/// @brief Given a "xs:choice" element in an XSD, get the minOccurs and maxOccurs options and store them in unsigned integers.
/// @details Stores 0 for "unbounded" or for anything that can't be interpreted as a nonnegative integer.
/// @param[in] choice_node A pointer ot an xmlNode representing an <xs:choice ... /> tag.
/// @param[out] choice_min The value parsed from minOccurs.
/// @param[out] choice_max The value parsed from maxOccurs.
void
get_choice_min_and_max(
	xmlNode* choice_node,
	platform::Size & choice_min,
	platform::Size & choice_max
) {
	std::string const choice_min_string( get_node_option( choice_node, "minOccurs" ) );
	std::string const choice_max_string( get_node_option( choice_node, "maxOccurs" ) );

	//Process choice_min:
	std::stringstream ssmin( choice_min_string );
	signed long int minval; //Deliberately signed.
	ssmin >> minval;
	if ( ssmin.fail() || minval < 0 ) {
		choice_min = 0;
	} else {
		choice_min = static_cast< platform::Size >(minval);
	}

	//Process choice_max:
	if ( !choice_max_string.compare("unbounded") ) {
		choice_max = 0;
	} else {
		std::stringstream ssmax( choice_max_string );
		signed long int maxval; //Deliberately signed.
		ssmax >> maxval;
		if ( ssmax.fail() || maxval < 0 ) {
			choice_max = 0;
		} else {
			choice_max = static_cast< platform::Size >(maxval);
		}
	}
}


/// @brief Given an xsd for RosettaScripts, get all of the names of mover elements.
/// @details The mover_names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
get_mover_names_from_xsd(
	std::string const &xsd,
	utility::vector1< std::string > &mover_names_out
) {
	get_module_names_from_xsd(xsd, "mover", mover_names_out);
}

/// @brief Given an xsd for RosettaScripts, get all of the names of filter elements.
/// @details The filter_names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
get_filter_names_from_xsd(
	std::string const &xsd,
	utility::vector1< std::string > &filter_names_out
) {
	get_module_names_from_xsd(xsd, "filter", filter_names_out);
}

/// @brief Given an xsd for RosettaScripts, get all of the names of task operation elements.
/// @details The taskop_names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
get_task_operation_names_from_xsd(
	std::string const &xsd,
	utility::vector1< std::string > &taskop_names_out
) {
	get_module_names_from_xsd( xsd, "task_operation", taskop_names_out );
}

/// @brief Given an xsd for RosettaScripts, get all of the names of residue-level task operation elements.
/// @details The rltaskop_names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
get_residue_level_task_operation_names_from_xsd(
	std::string const &xsd,
	utility::vector1< std::string > &rltaskop_names_out
) {
	get_module_names_from_xsd( xsd, "res_lvl_task_op", rltaskop_names_out );
}

/// @brief Given an xsd for RosettaScripts, get all of the names of residue selector elements.
/// @details The residue_selector_names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
get_residue_selector_names_from_xsd(
	std::string const &xsd,
	utility::vector1< std::string > &residue_selector_names_out
) {
	get_module_names_from_xsd( xsd, "residue_selector", residue_selector_names_out );
}

/// @brief Given an xsd for RosettaScripts, get all of the names of a particular type of Rosetta module.
/// @details The names_out vector is cleared and repopulated by this operation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
get_module_names_from_xsd(
	std::string const &xsd,
	std::string const &module_name,
	utility::vector1< std::string > &names_out
) {
	names_out.clear();

	xmlDoc* doc( xmlReadMemory( xsd.c_str(), xsd.length()+1, nullptr, nullptr, 0 ) ); //Must be deleted later!  Unfortunately, libxml2 likes raw pointers.

	xmlNode* rootnode( xmlDocGetRootElement(doc) );
	for ( xmlNode* group_node = rootnode->children; group_node!=nullptr; group_node = group_node->next ) {
		if ( //Find the list of the module type
				group_node->type != XML_ELEMENT_NODE ||
				strcmp( reinterpret_cast<const char*>(group_node->name), "group" ) ||
				get_node_option(group_node, "name").compare( module_name )
				) {
			continue;
		}
		for ( xmlNode* choice_node = group_node->children; choice_node!=nullptr; choice_node = choice_node->next ) { //Loop through first-level sub-nodes
			if ( //Find the xs:choice node
					choice_node->type != XML_ELEMENT_NODE ||
					strcmp( reinterpret_cast<const char*>(choice_node->name), "choice" )
					) {
				continue;
			}
			for ( xmlNode* module_node = choice_node->children; module_node!=nullptr; module_node = module_node->next ) {
				if ( //Find the module list
						module_node->type != XML_ELEMENT_NODE ||
						strcmp( reinterpret_cast<const char*>(module_node->name), "element")
						) {
					continue;
				}
				names_out.push_back( get_node_option(module_node, "name") );
			} //Looping through module nodes
		} //Looping through choice nodes
	} //Looping through group nodes

	//Deleting:
	xmlFreeDoc(doc); doc=nullptr;
	xmlCleanupParser();
}

/// @brief Given a node in an XML tree, get the value of an option.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string
get_node_option(
	xmlNode* node,
	std::string const &optionname
) {
	for ( xmlAttr* curoption = node->properties; curoption != nullptr; curoption = curoption->next ) {
		if ( curoption->type ==  XML_ATTRIBUTE_NODE && !strcmp(reinterpret_cast<const char*>(curoption->name),optionname.c_str()) ) {
			debug_assert( curoption->children != nullptr );
			debug_assert( curoption->children->type == XML_TEXT_NODE );
			return std::string(reinterpret_cast<const char*>(curoption->children->content) );
		}
	}
	return "";
}

/// @brief Go through an XSD and generate a human-readable description of a particular Rosetta module.
/// @details Calls generate_human_readable_recursive.
/// @note If component_name and component_type are provided, then summary information is only returned for that particular
/// object type.  For example, component_name="DeclareBond" and component_type="mover" would return information on the DeclareBond
/// mover.  Also note, this function uses raw pointers, unfortunately -- the libxml2 functions that I'm calling require it.
/// As such, no premature return statements should be added prior to the xmlFreeDoc and xmlCleanupParser statements at the end of the function.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
generate_human_readable_summary(
	std::string const &xsd,
	std::string const &component_name,
	std::string const &component_type
) {
	bool const has_name_and_type( !component_name.empty() && !component_type.empty() );
	std::string const name_and_type( has_name_and_type ? component_type + "_" + component_name + "_type" : "" );

	std::stringstream description(""), tags(""), options("");
	description << "DESCRIPTION:\n\n";
	tags << "USAGE:\n\n";
	options << "OPTIONS:\n\n";

	xmlDoc* doc( xmlReadMemory( xsd.c_str(), xsd.length()+1, nullptr, nullptr, 0 ) );

	generate_human_readable_recursive( xmlDocGetRootElement(doc), description, tags, options, 1, name_and_type, "" );

	//Deleting:
	xmlFreeDoc(doc); doc=nullptr;
	xmlCleanupParser();

	return description.str() + "\n" + tags.str() + "\n" + options.str();
}


/// @brief Go through an XSD XML and generate a human-readable version, recursively.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
generate_human_readable_recursive(
	xmlNode* rootnode,
	std::stringstream &description,
	std::stringstream &tags,
	std::stringstream &options,
	platform::Size const level,
	std::string const &complextype,
	std::string const &tag_name_to_print
) {

	if ( level > MAX_HUMAN_READABLE_RECURSION_LEVELS ) return;

	for ( xmlNode* cplxtype_node = rootnode->children; cplxtype_node!=nullptr; cplxtype_node = cplxtype_node->next ) {
		if ( cplxtype_node->type == XML_ELEMENT_NODE && !strcmp(reinterpret_cast<const char*>(cplxtype_node->name),"complexType") ) {

			if ( complextype.empty() ||  get_node_option(cplxtype_node, "name") == complextype ) {
				for ( platform::Size i=1; i<level; ++i ) { tags << "\t"; }

				tags << "<";
				std::string const tagname( tag_name_to_print.empty() ? get_tag_name(cplxtype_node) : tag_name_to_print );
				tags << tagname;
				output_all_tag_options(cplxtype_node, tagname, level, tags, options);
				tags << ">\n";

				//If this is level 1, generate the documentation information.
				if ( level==1 ) {
					for ( xmlNode* annotation_node = cplxtype_node->children; annotation_node!=nullptr; annotation_node = annotation_node->next ) {
						if ( !strcmp(reinterpret_cast<const char*>(annotation_node->name),"annotation") ) {
							generate_human_readable_documentation( annotation_node, description );
						}
					}
				}

				//Get sub-tags.
				for ( xmlNode* element_node = cplxtype_node->children; element_node != nullptr; element_node = element_node->next ) {
					if ( element_node->type == XML_ELEMENT_NODE ) {
						if ( !strcmp(reinterpret_cast<const char*>(element_node->name),"element") ) {
							std::string type_name( get_node_option( element_node, "type" ) );
							generate_human_readable_recursive( type_name.empty() ? element_node : rootnode, description, tags, options, level+1, type_name, type_name.empty() ? get_node_option( element_node, "name" ) : "" );
						} else if ( !strcmp(reinterpret_cast<const char*>(element_node->name),"choice") ) {
							for ( xmlNode* element_node2 = element_node->children; element_node2 != nullptr; element_node2 = element_node2->next ) {
								if ( element_node2->type == XML_ELEMENT_NODE && !strcmp(reinterpret_cast<const char*>(element_node2->name),"element") ) {
									std::string type_name( get_node_option( element_node2, "type" ) );
									generate_human_readable_recursive( type_name.empty() ? element_node2 : rootnode, description, tags, options, level+1, type_name, type_name.empty() ? get_node_option( element_node2, "name" ) : "" );
								}
							}
						}
					}
				}

				for ( platform::Size i=1; i<level; ++i ) { tags << "\t"; }
				tags << "</" << tagname << ">\n";
			}
		}
	}
}

/// @brief Given a node in an XML tree representing a tag, get the string that starts the tag.  (e.g. "DeclareBond" in "<DeclareBond .... />").
/// @details From a string of format "*_name_type" or "*_name", extracts out "name".
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string
get_tag_name(
	xmlNode* node
) {
	std::string const extendedname( get_node_option(node, "name") );
	utility::vector1< std::string > const splitstr( string_split(extendedname, '_' ) );
	debug_assert( splitstr.size() > 1 );
	if ( splitstr[ splitstr.size() ].compare( "type" ) ) return splitstr[ splitstr.size() ];
	return splitstr[ splitstr.size() - 1 ];
}

/// @brief Given a node in an XML tree representing a tag, get all options in that tag, write their descriptions to "options",
/// and write their names and types to "tags".
/// @details Heads the section in "options" with the tag name.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
output_all_tag_options(
	xmlNode* node,
	std::string const &tagname,
	platform::Size const level,
	std::stringstream &tags,
	std::stringstream &options
) {
	options << (level > 1 ? "\n" : "" ) << "\"" << tagname << "\" ";
	for ( platform::Size i=1; i<level; ++i ) options << "sub-";
	options << "tag:";

	if ( level > 1 ) {
		options << " ";
		generate_human_readable_documentation(node, options);
		options << "\n";
	} else {
		options << "\n\n";
	}

	for ( xmlNode* subnode = node->children; subnode != nullptr; subnode = subnode->next ) { //Loop through sub-nodes.
		if ( subnode->type == XML_ELEMENT_NODE && !strcmp(reinterpret_cast<const char*>(subnode->name),"attribute") ) {
			std::string const name( get_node_option( subnode, "name" ) );
			std::string const type( get_type_name( get_node_option( subnode, "type" ) ) );
			std::string const def( get_node_option(subnode, "default") );
			tags << " " << name << "=(" << type << (def.empty() ? "" : ",\"" + def + "\"" ) << ")";

			options << "\t" << name << " (" << type << (def.empty() ? "" : ",\"" + def + "\"" ) << "):  ";
			generate_human_readable_documentation( subnode, options );
		}
	}
}

/// @brief Given a documentation node in an XML tree, parse out the documentation and print it.
void
generate_human_readable_documentation(
	xmlNode* parent_node,
	std::stringstream &outstream
) {
	std::stringstream outstream2;

	for ( xmlNode* doc_node = parent_node->children; doc_node != nullptr; doc_node = doc_node->next ) {
		if ( doc_node->type == XML_ELEMENT_NODE && !strcmp(reinterpret_cast<const char*>(doc_node->name),"documentation") ) {
			for ( xmlNode* cur_node = doc_node->children; cur_node != nullptr; cur_node = cur_node->next ) {
				if ( cur_node->type == XML_TEXT_NODE ) {
					outstream2 << cur_node->content;
					break;
				}
			}
		} else if ( doc_node->type == XML_ELEMENT_NODE && !strcmp(reinterpret_cast<const char*>(doc_node->name),"annotation") ) {
			generate_human_readable_documentation( doc_node, outstream );
		}
	}

	outstream <<  utility::strip( outstream2.str(), " \t\n" ) << "\n" ;
}

/// @brief Given an XSD simple type (e.g. "rosetta_bool"), return a more human-readable type (e.g. "bool").
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string
get_type_name(
	std::string const& xsd_type
) {
	if ( !xsd_type.compare( "rosetta_bool" ) ) { return "bool"; }
	else if (
			!xsd_type.compare( "xs:string" ) ||
			!xsd_type.compare( "task_operation" ) ||
			!xsd_type.compare( "task_operation_comma_separated_list" ) ||
			!xsd_type.compare( "pose_cached_task_operation" )
			) { return "string"; }
	else if ( !xsd_type.compare( "xs:integer" ) ) { return "int"; }
	else if ( !xsd_type.compare( "xs:decimal" ) ) { return "real"; }
	else if ( !xsd_type.compare( "non_negative_integer" ) ) { return "int"; }
	return xsd_type;
}

}  // namespace io
}  // namespace utility
