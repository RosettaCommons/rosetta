// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/utility/tag/XMLSchemaValidation.hh
/// @brief  functions and classes needed to validate an XML file against a schema
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_tag_XMLSchemaValidation_HH
#define INCLUDED_utility_tag_XMLSchemaValidation_HH

// Unit headers
#include <utility/tag/XMLSchemaValidation.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// LibXML headers
#include <libxml/tree.h>
#include <libxml/xmlerror.h>

// Boost headers
//#include <boost/function.hpp>

// C++ headers
#include <list>
#include <string>
//#include <vector>
//#include <map>
//#include <iosfwd>
//#include <sstream>

namespace utility {
namespace tag {

void handle_xml_error( void * ctxt, char const * message, ... );
void handle_xml_warning( void * ctxt, char const * message, ... );
void handle_structured_xml_error( void * ctxt, xmlErrorPtr error );

class XMLErrorHandler {
public:
	typedef std::list< std::string > Errors;
	typedef std::list< std::string > Warnings;

	void set_file_contents( std::string const & file_contents );

	void handle_xml_error( std::string const & message, int line );
	void handle_xml_warning( std::string const & message, int line );

	std::list< std::string > const & errors() const;
	std::list< std::string > const & warnings() const;

	void push_node( xmlNode * node );
	void pop_node();

private:
	std::string lines_near_error( int line ) const;

private:
	std::list< std::string > error_list_;
	std::list< std::string > warning_list_;
	utility::vector1< xmlNode * > node_stack_;
	utility::vector1< std::string > file_lines_;
};


class XMLValidationOutput
{
public:
	XMLValidationOutput();

	void valid( bool );
	bool valid() const;

	std::list< std::string > const & errors() const;
	std::list< std::string > const & warnings() const;

	void errors( std::list< std::string > const & error_list );
	void warnings( std::list< std::string > const & warning_list );

	std::string error_messages() const;
	std::string warning_messages() const;

private:
	bool valid_;
	std::list< std::string > error_list_;
	std::list< std::string > warning_list_;
};

XMLValidationOutput
validate_xml_against_xsd(
	std::string const & xml,
	std::string const & xsd
);

XMLValidationOutput
test_if_schema_is_valid(
	std::string const & xsd
);


}
}

#endif

