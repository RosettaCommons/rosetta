// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/utility/tag/XMLSchemaGeneration.fwd.hh
/// @brief  forward declaration of the classes used to define an XML Schema
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_tag_XMLSchemaGeneration_HH
#define INCLUDED_utility_tag_XMLSchemaGeneration_HH

// Unit headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>
#include <vector>
#include <map>
#include <iosfwd>
#include <sstream>

namespace utility {
namespace tag {


class XMLSchemaType : public utility::pointer::ReferenceCount {
public:
	XMLSchemaType();
	XMLSchemaType( XMLSchemaDataType setting );
	XMLSchemaType( std::string const & custom_type );
	XMLSchemaType( char const * custom_type );
	void type( XMLSchemaDataType setting );
	void custom_type_name( std::string const & setting );

	std::string type_name() const;

private:
	XMLSchemaDataType type_;
	std::string custom_type_name_;
};

class XMLSchemaAttribute : public utility::pointer::ReferenceCount {
public:
	XMLSchemaAttribute();
	XMLSchemaAttribute( std::string const & name, XMLSchemaType type, bool is_required=false );
	XMLSchemaAttribute( std::string const & name, XMLSchemaType type, std::string const & default_value );
	XMLSchemaAttribute( std::string const & name, XMLSchemaType type, char const * default_value );
	void name( std::string const & setting );
	void type( XMLSchemaType setting );
	void default_value( std::string const & setting );
	void is_required( bool setting );

	XMLSchemaType const & type() const;
	void write_definition( int indentation, std::ostream & os ) const;

private:
	std::string name_;
	XMLSchemaType type_;
	std::string default_value_;
	bool is_required_;
};

class XMLSchemaRestriction : public utility::pointer::ReferenceCount {
public:
	XMLSchemaRestriction();
	void name( std::string const & setting );
	void base_type( XMLSchemaType setting );
	void add_restriction( XMLSchemaRestrictionType type, std::string const & value );

	void write_definition( int indentation, std::ostream & os ) const;

private:
	std::string name_;
	XMLSchemaType base_type_;
	std::list< std::pair< XMLSchemaRestrictionType, std::string > > restrictions_;
};

class XMLSchemaComplexType : public utility::pointer::ReferenceCount {
public:
	XMLSchemaComplexType();
	void name( std::string const & setting );
	void type( XMLSchemaComplexTypeType setting );
	void reference_type( std::string const & setting ); // used if this is simply a reference to an existing type
	void subtype( XMLSchemaComplexTypeOP setting ); // used iff type_ is xsctt_group and reference_type_ is the empty string.
	void add_subelement( XMLSchemaElementOP  subelement );
	void add_attribute( XMLSchemaAttribute attribute );
	void min_occurs( int setting );
	void max_occurs( int setting );

	void write_definition( int indentation, std::ostream & os ) const;

private:
	std::string name_;
	XMLSchemaComplexTypeType type_;
	std::string reference_type_; // used if this is simply a reference to an existing type
	XMLSchemaComplexTypeOP subtype_; // used iff type_ is xsctt_group and reference_type_ is the empty string.
	std::list< XMLSchemaElementOP > subelements_;
	std::list< XMLSchemaAttribute > attributes_;
	int min_occurs_;
	int max_occurs_;
	bool suppress_complex_type_output_;
};


/// @brief An XMLSchema element, e.g. <FixbbMover name="fixbb" task_operations="ex1,ex2"/> or
/// <And><Chain id="A"/></And>
class XMLSchemaElement : public utility::pointer::ReferenceCount {
public:
	XMLSchemaElement();
	void name( std::string const & setting );
	void type_name( XMLSchemaType setting );
	void group_name( std::string const & setting );
	// I can't remember why I thought this was important void element_type_def( XMLSchemaAttributeOP setting );
	void element_type_def( XMLSchemaComplexTypeOP setting );
	void restriction_type_def( XMLSchemaRestrictionOP setting );
	void min_occurs( int setting );
	void max_occurs( int setting );

	void write_definition( int indentation, std::ostream & os ) const;

private:

	std::string name_;
	XMLSchemaElementCategory category_;
	XMLSchemaType type_name_;                     // used if category_ is xs_element_is_type_reference
	std::string ref_name_;                        // used if category_ is xs_element_is_group_reference
	XMLSchemaComplexTypeOP complex_type_;         // used if category_ is xs_element_is_complex_type_w_definition
	XMLSchemaRestrictionOP restriction_type_def_; // used if category_ is xs_element_is_restriction_w_definition
	int min_occurs_;
	int max_occurs_;
};

class XMLSchemaDefinition : public utility::pointer::ReferenceCount
{
public:
	XMLSchemaDefinition();
	virtual ~XMLSchemaDefinition();
	void add_top_level_element( std::string const & element_name, std::string const & definition );
	bool has_top_level_element( std::string const & element_name ) const;
	std::string full_definition() const;

private:
	void validate_new_top_level_element( std::string const & element_name, std::string const & definition );

	std::list< std::string > elements_in_order_;
	std::map< std::string, std::string > top_level_elements_;
};


void indent_w_spaces( int indentation, std::ostream & os );
std::string restriction_type_name( XMLSchemaRestrictionType type );
std::ostream & operator << ( std::ostream & os, XMLSchemaRestrictionType type );
std::string xs_complex_type_type_name( XMLSchemaComplexTypeType xsctt );
std::ostream & operator << ( std::ostream & os, XMLSchemaComplexTypeType type );


}
}

#endif
