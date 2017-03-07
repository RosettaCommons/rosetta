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
	else if ( !xsd_type.compare( "xs:string" ) ) { return "string"; }
	else if ( !xsd_type.compare( "xs:integer" ) ) { return "int"; }
	else if ( !xsd_type.compare( "xs:decimal" ) ) { return "real"; }
	else if ( !xsd_type.compare( "non_negative_integer" ) ) { return "int"; }
	return xsd_type;
}

}  // namespace io
}  // namespace utility
