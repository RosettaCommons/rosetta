// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/utility/tag/XMLSchemaValidation.cc
/// @brief  functions and classes needed to validate an XML file against a schema
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <utility/tag/XMLSchemaValidation.hh>

// Utility headers
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>

// LibXML2 headers
#include <libxml/xmlschemas.h>
#include <libxml/xmlstring.h>
#include <libxml/tree.h>

#include <platform/types.hh>

// Utility headers
//#include <utility/pointer/ReferenceCount.hh>

// Boost headers
//#include <boost/function.hpp>

// C++ headers
//#include <list>
//#include <vector>
//#include <map>
//#include <iosfwd>
#include <sstream>


namespace utility {
namespace tag {

std::string element_type_name(
	xmlElementType element_type
)
{
	switch ( element_type ) {
	case XML_ELEMENT_NODE : return "XML_ELEMENT_NODE";
	case XML_ATTRIBUTE_NODE : return "XML_ATTRIBUTE_NODE";
	case XML_TEXT_NODE : return "XML_TEXT_NODE";
	case XML_CDATA_SECTION_NODE : return "XML_CDATA_SECTION_NODE";
	case XML_ENTITY_REF_NODE : return "XML_ENTITY_REF_NODE";
	case XML_ENTITY_NODE : return "XML_ENTITY_NODE";
	case XML_PI_NODE : return "XML_PI_NODE";
	case XML_COMMENT_NODE : return "XML_COMMENT_NODE";
	case XML_DOCUMENT_NODE : return "XML_DOCUMENT_NODE";
	case XML_DOCUMENT_TYPE_NODE : return "XML_DOCUMENT_TYPE_NODE";
	case XML_DOCUMENT_FRAG_NODE : return "XML_DOCUMENT_FRAG_NODE";
	case XML_NOTATION_NODE : return "XML_NOTATION_NODE";
	case XML_HTML_DOCUMENT_NODE : return "XML_HTML_DOCUMENT_NODE";
	case XML_DTD_NODE : return "XML_DTD_NODE";
	case XML_ELEMENT_DECL : return "XML_ELEMENT_DECL";
	case XML_ATTRIBUTE_DECL : return "XML_ATTRIBUTE_DECL";
	case XML_ENTITY_DECL : return "XML_ENTITY_DECL";
	case XML_NAMESPACE_DECL : return "XML_NAMESPACE_DECL";
	case XML_XINCLUDE_START : return "XML_XINCLUDE_START";
	case XML_XINCLUDE_END : return "XML_XINCLUDE_END";
	case XML_DOCB_DOCUMENT_NODE : return "XML_DOCB_DOCUMENT_NODE";
	}
	return "(Unrecognized XML element type!)";
}


int
length_of_message( const char * message, va_list args )
{
	return vsnprintf(NULL, 0, message, args );
}

char *
make_message( int length, const char * message, va_list args ) {
	char * formatted_message = new char[ length+1 ];
	vsnprintf( formatted_message, length+1, message, args );
	//std::cout << "message: " << formatted_message << std::endl;
	return formatted_message;
}

void handle_structured_xml_error( void * ctxt, xmlErrorPtr error ) {
	XMLErrorHandler * handler = static_cast< XMLErrorHandler * > ( ctxt );
	std::ostringstream oss;
	if ( error->node ) {
		xmlNode * node( ( xmlNode * ) error->node );
		//oss << "Node " << node->content << " on line " << node->line << " of type " << element_type_name( node->type ) << "\n";
		oss << "From line " << xmlGetLineNo(node) << ":\n";
	}
	oss << "Error: " << error->message << "\n";
	//std::cout << "Error message from handle_structured_xml_error: " << oss.str() << " and node? " << error->node << " on line? " << error->line << "\n" << std::endl;
	handler->handle_xml_error( oss.str(), error->line );
}

void handle_xml_error( void * ctxt, char const * message, ... )
{
	va_list args;
	va_start( args, message );
	int length = length_of_message( message, args );
	va_end( args );

	va_start( args, message );
	char * formatted_message = make_message( length, message, args );
	va_end( args );

	std::string formatted_string( formatted_message );
	delete [] formatted_message;

	XMLErrorHandler * handler = static_cast< XMLErrorHandler * > ( ctxt );
	handler->handle_xml_error( formatted_string, 0 );
}

void handle_xml_warning( void * ctxt, char const * message, ... )
{
	va_list args;
	va_start( args, message );
	int length = length_of_message( message, args );
	va_end( args );

	va_start( args, message );
	char * formatted_message = make_message( length, message, args );
	va_end( args );

	std::string formatted_string( formatted_message );
	delete [] formatted_message;

	XMLErrorHandler * handler = static_cast< XMLErrorHandler * > ( ctxt );
	handler->handle_xml_warning( formatted_string, 0 );
}


void print_node_stack( utility::vector1< xmlNode * > const & node_stack, std::ostringstream & oss ) {
	for ( platform::Size ii = 1; ii <= node_stack.size(); ++ii ) {
		xmlNode * iinode = node_stack[ ii ];
		oss << "Node " << iinode->content << " on line " << iinode->line << " of type " << element_type_name( iinode->type ) << "\n";
	}
}

void XMLErrorHandler::set_file_contents( std::string const & file_contents )
{
	file_lines_ = string_split( file_contents, '\n' );
}

void XMLErrorHandler::handle_xml_error( std::string const & message, int line )
{

	//std::ostringstream oss;
	//oss << "Error detected\n";
	//print_node_stack( node_stack_, oss );
	//oss << message << "\n";
	//error_list_.push_back( oss.str() );
	error_list_.push_back( message + lines_near_error( line ) );
}

void XMLErrorHandler::handle_xml_warning( std::string const & message, int line )
{
	//std::ostringstream oss;
	//oss << "Warning detected\n";
	//print_node_stack( node_stack_, oss );
	//oss << message << "\n";
	//warning_list_.push_back( oss.str() );
	warning_list_.push_back( message + lines_near_error( line ) );
}

std::list< std::string > const & XMLErrorHandler::errors() const
{
	return error_list_;
}
std::list< std::string > const & XMLErrorHandler::warnings() const
{
	return warning_list_;
}

void XMLErrorHandler::push_node( xmlNode * node )
{
	node_stack_.push_back( node );
}

void XMLErrorHandler::pop_node() {
	node_stack_.pop_back();
}

std::string
XMLErrorHandler::lines_near_error( int line ) const
{
	using ObjexxFCL::format::RJ;

	std::string extra_message = "";
	if ( line > 0 && line <= (int) file_lines_.size() ) {
		std::ostringstream extra_ss;
		int line_begin = std::max( line - 5, 1 );
		int line_end   = std::min( line + 5, (int) file_lines_.size() );
		int last_line_width = utility::to_string( line_end ).size();
		for ( int ii = line_begin; ii <= line_end; ++ii ) {
			extra_ss << RJ( last_line_width, ii ) << ": " << file_lines_[ ii ] << "\n";
		}
		extra_message = extra_ss.str();
	}
	return extra_message;
}

XMLValidationOutput::XMLValidationOutput() : valid_( true ) {}

void XMLValidationOutput::valid( bool setting ) { valid_ = setting; }
bool XMLValidationOutput::valid() const { return valid_; }

std::list< std::string > const & XMLValidationOutput::errors() const { return error_list_; }
std::list< std::string > const & XMLValidationOutput::warnings() const { return warning_list_; }

void XMLValidationOutput::errors( std::list< std::string > const & error_list ) { error_list_ = error_list; }
void XMLValidationOutput::warnings( std::list< std::string > const & warning_list ) { warning_list_ = warning_list; }


std::string concat_stringlist( std::list< std::string > const & strings ) {
	std::ostringstream oss;
	for ( std::list< std::string >::const_iterator iter = strings.begin(); iter != strings.end(); ++iter ) {
		oss << *iter;
	}
	return oss.str();
}

std::string XMLValidationOutput::error_messages() const
{
	return concat_stringlist( error_list_ );
}

std::string XMLValidationOutput::warning_messages() const
{
	return concat_stringlist( warning_list_ );
}

int
recurse_through_tree(
	xmlSchemaValidCtxtPtr ctxt,
	xmlNodePtr node,
	XMLErrorHandler handler
)
{
	if ( node == 0 ) return 0;
	handler.push_node( node );
	int error = xmlSchemaValidateOneElement( ctxt, node );
	xmlNodePtr child = xmlFirstElementChild( node );
	while ( child ) {
		int error_in_child = recurse_through_tree( ctxt, child, handler );
		child = xmlNextElementSibling( child );
		if ( error == 0 && error_in_child != 0 ) error = error_in_child;
	}
	handler.pop_node();
	return error;
}


XMLValidationOutput
validate_xml_against_xsd(
	std::string const & xml_string,
	std::string const & xsd_string
)
{
	XMLValidationOutput output;

	xmlLineNumbersDefault( 1 );
	XMLErrorHandler handler1;
	handler1.set_file_contents( xsd_string );
	xmlSetStructuredErrorFunc( & handler1, handle_structured_xml_error );

	//xmlChar * xsd_xmlchar_string = xmlCharStrdup( xsd_string.c_str() );
	//xmlChar * xml_version = xmlCharStrdup( "1.0" );
	//xmlDoc * xsd_doc = xmlParseDoc( xsd_xmlchar_string );

	xmlDoc * xsd_doc = xmlParseMemory( xsd_string.c_str(), xsd_string.size() );

	xmlSchemaParserCtxtPtr schema_parser_context = xmlSchemaNewDocParserCtxt( xsd_doc );
	//xmlSchemaSetParserErrors( schema_parser_context, handle_xml_error, handle_xml_warning, &handler );
	xmlSchemaPtr schema = xmlSchemaParse( schema_parser_context );

	//std::cout << "Parsed the schema" << std::endl;
	xmlSchemaValidCtxtPtr schema_validator = xmlSchemaNewValidCtxt( schema );
	//std::cout << "Created schema validator" << std::endl;

	//xmlSchemaSetValidErrors( schema_validator, handle_xml_error, handle_xml_warning, &handler );

	xmlChar * xml_input_xmlchar = xmlCharStrdup( xml_string.c_str() );
	xmlDoc * xml_doc = xmlParseDoc( xml_input_xmlchar );

	//std::cout << "Validating XML document" << std::endl;
	XMLErrorHandler handler2;
	handler2.set_file_contents( xml_string );
	xmlSetStructuredErrorFunc( & handler2, handle_structured_xml_error );

	int validation_output = xmlSchemaValidateDoc( schema_validator, xml_doc );
	//int first_error_code = recurse_through_tree( schema_validator, xmlDocGetRootElement( xml_doc ), handler );

	// clean up
	//free( xsd_xmlchar_string );
	free( xml_input_xmlchar );
	xmlSchemaFree( schema );
	xmlSchemaFreeParserCtxt( schema_parser_context );
	xmlSchemaFreeValidCtxt( schema_validator );
	xmlFreeDoc( xsd_doc );

	//output.valid( first_error_code == 0 );
	output.valid( validation_output == 0 );
	output.errors( handler2.errors() );
	output.warnings( handler2.warnings() );

	return output;
}

XMLValidationOutput
test_if_schema_is_valid(
	std::string const & xsd_string
)
{
	XMLValidationOutput output;

	xmlLineNumbersDefault( 1 );
	XMLErrorHandler handler;
	handler.set_file_contents( xsd_string );
	xmlSetStructuredErrorFunc( & handler, handle_structured_xml_error );

	xmlDoc * xsd_doc = xmlParseMemory( xsd_string.c_str(), xsd_string.size() );

	xmlSchemaParserCtxtPtr schema_parser_context = xmlSchemaNewDocParserCtxt( xsd_doc );
	//xmlSchemaSetParserErrors( schema_parser_context, handle_xml_error, handle_xml_warning, &handler );
	xmlSchemaPtr schema = xmlSchemaParse( schema_parser_context );

	xmlSchemaFree( schema );
	xmlSchemaFreeParserCtxt( schema_parser_context );
	xmlFreeDoc( xsd_doc );


	output.valid( handler.errors().empty() );
	output.errors( handler.errors() );
	output.warnings( handler.warnings() );

	return output;
}



}
}
