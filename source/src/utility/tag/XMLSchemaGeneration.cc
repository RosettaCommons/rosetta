// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/utility/tag/XMLSchemaGeneration.cc
/// @brief  Definitions of the classes used to define an XML Schema
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <utility/tag/XMLSchemaGeneration.hh>

// Package headers
#include <utility/tag/Tag.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

namespace utility {
namespace tag {

void indent_w_spaces( int indentation, std::ostream & os ) {
	for ( int ii = 1; ii <= indentation; ++ii ) { os << " "; }
}

/// @details Use "-1" for unspecified, use -2 for "unbounded"
void write_min_occurs_max_occurs_if_necessary( int min_occurs, int max_occurs, std::ostream & os )
{
	if ( min_occurs != xsminmax_unspecified ) {
		os << " minOccurs=\"" << min_occurs << "\"";
	}
	if ( max_occurs == xsminmax_unbounded ) {
		os << " maxOccurs=\"unbounded\"";
	} else if ( max_occurs > 0 ) {
		os << " maxOccurs=\"" << max_occurs << "\"";
	}

}

XMLSchemaType::XMLSchemaType() :
	type_( xs_string )
{}

XMLSchemaType::XMLSchemaType( XMLSchemaDataType setting ) :
	type_( setting )
{}

XMLSchemaType::XMLSchemaType( std::string const & custom_type ) :
	type_( xs_custom ),
	custom_type_name_( custom_type )
{}

XMLSchemaType::XMLSchemaType( char const * custom_type ) :
	type_( xs_custom ),
	custom_type_name_( custom_type )
{}

void XMLSchemaType::type( XMLSchemaDataType setting ) { type_ = setting; }
void XMLSchemaType::custom_type_name( std::string const & setting ) { type_ = xs_custom; custom_type_name_ = setting; }

std::string XMLSchemaType::type_name() const {
	switch ( type_ ) {
	case xs_string : return "xs:string"; break;
	case xs_decimal : return "xs:decimal"; break;
	case xs_integer : return "xs:integer"; break;
	case xs_boolean : return "xs:boolean"; break;
	case xs_date : return "xs:date"; break;
	case xs_time : return "xs:time"; break;
	case xs_custom : return custom_type_name_;
	}
	return "error!"; // appease compiler
}


XMLSchemaAttribute::XMLSchemaAttribute() :
	is_required_( false )
{}

XMLSchemaAttribute::XMLSchemaAttribute(
	std::string const & name,
	XMLSchemaType type,
	bool is_required
) :
	name_( name ),
	type_( type ),
	is_required_( is_required )
{}

XMLSchemaAttribute::XMLSchemaAttribute(
	std::string const & name,
	XMLSchemaType type,
	std::string const & default_value
) :
	name_( name ),
	type_( type ),
	default_value_( default_value ),
	is_required_( false )
{}

XMLSchemaAttribute::XMLSchemaAttribute(
	std::string const & name,
	XMLSchemaType type,
	char const * default_value
) :
	name_( name ),
	type_( type ),
	default_value_( default_value ),
	is_required_( false )
{}

void XMLSchemaAttribute::name( std::string const & setting ) { name_ = setting; }
void XMLSchemaAttribute::type( XMLSchemaType setting ) { type_ = setting; }
void XMLSchemaAttribute::is_required( bool setting ) { is_required_ = setting; }
void XMLSchemaAttribute::default_value( std::string const & setting ) { default_value_ = setting; }

XMLSchemaType const & XMLSchemaAttribute::type() const { return type_; }

void XMLSchemaAttribute::write_definition( int indentation, std::ostream & os ) const {
	indent_w_spaces( indentation, os );
	os << "<xs:attribute name=\"" << name_ << "\" type=\"";
	os << type_.type_name();
	os << "\"";
	if ( default_value_ != "" ) os << " default=\"" << default_value_ << "\"";
	if ( is_required_ ) os << " use=\"required\"";
	os << "/>\n";
}


std::string
restriction_type_name( XMLSchemaRestrictionType type )
{
	switch ( type ) {
	case xsr_enumeration :     return "xs:enumeration";
	case xsr_fractionDigits :  return "xs:fractionDigits";
	case xsr_length :          return "xs:length";
	case xsr_maxExclusive :    return "xs:maxExclusive";
	case xsr_maxInclusive :    return "xs:maxInclusive";
	case xsr_maxLength :       return "xs:maxLength";
	case xsr_minExclusive :    return "xs:minExclusive";
	case xsr_minInclusive :    return "xs:minInclusive";
	case xsr_minLength :       return "xs:minLength";
	case xsr_pattern :         return "xs:pattern";
	case xsr_totalDigits :     return "xs:totalDigits";
	case xsr_whitespace :      return "xs:whitespace";
	}
	return "error!"; // appease compiler
}

std::ostream & operator << ( std::ostream & os, XMLSchemaRestrictionType type ) {
	os << restriction_type_name( type );
	return os;
}

XMLSchemaRestriction::XMLSchemaRestriction() {}

void XMLSchemaRestriction::name( std::string const & setting ) { name_ = setting; }
void XMLSchemaRestriction::base_type( XMLSchemaType setting ) { base_type_ = setting; }
void XMLSchemaRestriction::add_restriction( XMLSchemaRestrictionType type, std::string const & value ) {
	restrictions_.push_back( std::make_pair( type, value ) );
}

void XMLSchemaRestriction::write_definition( int indentation, std::ostream & os ) const {
	indent_w_spaces( indentation, os );
	os << "<xs:simpleType";
	if ( name_ != "" ) {
		os << " name=\"" << name_ << "\"";
	}
	os << ">\n";
	indent_w_spaces( indentation+1, os );
	os << "<xs:restriction base=\"" << base_type_.type_name() << "\">\n";
	for ( std::list< std::pair< XMLSchemaRestrictionType, std::string > >::const_iterator
			iter = restrictions_.begin(),
			iter_end = restrictions_.end();
			iter != iter_end; ++iter ) {
		indent_w_spaces( indentation+2, os );
		os << "<" << iter->first << " value=\"" << iter->second << "\"/>\n";
	}
	indent_w_spaces( indentation+1, os );
	os << "</xs:restriction>\n";
	indent_w_spaces( indentation, os );
	os << "</xs:simpleType>\n";
}

std::string xs_complex_type_type_name( XMLSchemaComplexTypeType xsctt ) {
	switch ( xsctt ) {
	case xsctt_empty :    return "xs:empty"; // this isn't a real type!
	case xsctt_sequence : return "xs:sequence";
	case xsctt_all :      return "xs:all";
	case xsctt_choice :   return "xs:choice";
	case xsctt_group :    return "xs:group";
	}
	return "error!"; // appease compiler
}

std::ostream & operator << ( std::ostream & os, XMLSchemaComplexTypeType type ) {
	os << xs_complex_type_type_name( type );
	return os;
}


XMLSchemaComplexType::XMLSchemaComplexType() :
	type_( xsctt_empty ),
	min_occurs_( xsminmax_unspecified ),
	max_occurs_( xsminmax_unspecified ),
	suppress_complex_type_output_( false )
{}
void XMLSchemaComplexType::name( std::string const & setting ) { name_ = setting; }
void XMLSchemaComplexType::type( XMLSchemaComplexTypeType setting ) { type_ = setting; }
void XMLSchemaComplexType::reference_type( std::string const & setting ) { reference_type_ = setting; }
void XMLSchemaComplexType::subtype( XMLSchemaComplexTypeOP setting ) { type_ = xsctt_group; subtype_ = setting; subtype_->suppress_complex_type_output_ = true; }
void XMLSchemaComplexType::add_subelement( XMLSchemaElementOP  subelement ) { assert( type_ != xsctt_group ); subelements_.push_back( subelement ); }
void XMLSchemaComplexType::add_attribute( XMLSchemaAttribute attribute ) { assert( type_ != xsctt_group ); attributes_.push_back( attribute ); }
void XMLSchemaComplexType::min_occurs( int setting ) { min_occurs_ = setting; }
void XMLSchemaComplexType::max_occurs( int setting ) { max_occurs_ = setting; }

void XMLSchemaComplexType::write_definition( int indentation, std::ostream & os ) const {
	if ( type_ == xsctt_group ) {
		if ( reference_type_ == "" ) {
			indent_w_spaces( indentation, os );
			os << "<xs:group name=\"" << name_ << "\"";
			write_min_occurs_max_occurs_if_necessary( min_occurs_, max_occurs_, os );
			os << ">\n";
			subtype_->write_definition( indentation+1, os );
			indent_w_spaces( indentation, os );
			os << "</xs:group>\n";
		} else {
			indent_w_spaces( indentation, os );
			os << "<xs:group ref=\"" << reference_type_ << "\"";
			os << ">\n";
		}
	} else {
		if ( ! suppress_complex_type_output_ ) {
			indent_w_spaces( indentation, os );
			os << "<xs:complexType";
			if ( name_ != "" ) {
				os << " name=\"" << name_ << "\"";
			}
			os << " mixed=\"true\">\n";
			++indentation;
		}
		if ( ! subelements_.empty() ) {
			indent_w_spaces( indentation, os );
			os << "<" << type_;
			write_min_occurs_max_occurs_if_necessary( min_occurs_, max_occurs_, os );
			os << ">\n";
			++indentation;
			for ( std::list< XMLSchemaElementOP >::const_iterator
					iter = subelements_.begin(),
					iter_end = subelements_.end();
					iter != iter_end; ++iter ) {
				(*iter)->write_definition( indentation, os );
			}
			--indentation;
			indent_w_spaces( indentation, os );
			os << "</" << type_ << ">\n";
		}
		if ( ! attributes_.empty() ) {
			for ( std::list< XMLSchemaAttribute >::const_iterator
					iter = attributes_.begin(),
					iter_end = attributes_.end();
					iter != iter_end; ++iter ) {
				iter->write_definition( indentation, os );
			}
		}
		if ( ! suppress_complex_type_output_ ) {
			--indentation;
			indent_w_spaces( indentation, os );
			os << "</xs:complexType>\n";
		}
	}
}

///////////////////////////////////////////// XMLSchemaElement ///////////////////////////////////////////////////

XMLSchemaElement::XMLSchemaElement() : category_( xs_element_is_type_reference ), min_occurs_( xsminmax_unspecified ), max_occurs_( xsminmax_unspecified ) {}
void XMLSchemaElement::name( std::string const & setting ) { name_ = setting; }
void XMLSchemaElement::type_name( XMLSchemaType setting ) { type_name_ = setting; category_ = xs_element_is_type_reference; }
void XMLSchemaElement::group_name( std::string const & setting ) { ref_name_ = setting; category_ = xs_element_is_group_reference; }
void XMLSchemaElement::element_type_def( XMLSchemaComplexTypeOP setting ) { complex_type_ = setting; category_ = xs_element_is_complex_type_w_definition; }
void XMLSchemaElement::restriction_type_def( XMLSchemaRestrictionOP setting ) { restriction_type_def_ = setting; category_ = xs_element_is_restriction_w_definition; }
void XMLSchemaElement::min_occurs( int setting ) { min_occurs_ = setting; }
void XMLSchemaElement::max_occurs( int setting ) { max_occurs_ = setting; }

void XMLSchemaElement::write_definition( int indentation, std::ostream & os ) const {
	if ( category_ == xs_element_is_group_reference ) {
		indent_w_spaces( indentation, os );
		os << "<xs:group ref=\"" << ref_name_ << "\"" ;
		write_min_occurs_max_occurs_if_necessary( min_occurs_, max_occurs_, os );
		os << "/>\n";
	} else if ( category_ == xs_element_is_type_reference ) {
		indent_w_spaces( indentation, os );
		os << "<xs:element name=\"" << name_ << "\" type=\"" << type_name_.type_name() << "\"" ;
		write_min_occurs_max_occurs_if_necessary( min_occurs_, max_occurs_, os );
		os << "/>\n";
	} else {
		indent_w_spaces( indentation, os );
		os << "<xs:element name=\"" << name_ << "\"";
		write_min_occurs_max_occurs_if_necessary( min_occurs_, max_occurs_, os );
		os << ">\n";
		if ( category_ == xs_element_is_complex_type_w_definition ) {
			complex_type_->write_definition( indentation+1, os );
		} else {
			restriction_type_def_->write_definition( indentation+1, os );
		}

		indent_w_spaces( indentation, os );
		os << "</xs:element>\n";
	}
}


//////////////////////////////// XMLSchemaDefinition /////////////////////////////////////////////

XMLSchemaDefinition::XMLSchemaDefinition() {}

XMLSchemaDefinition::~XMLSchemaDefinition() {}

void XMLSchemaDefinition::add_top_level_element( std::string const & element_name, std::string const & definition )
{
	validate_new_top_level_element( element_name, definition );
	if ( top_level_elements_.find( element_name ) == top_level_elements_.end() ) {
		top_level_elements_[ element_name ] = definition;
		elements_in_order_.push_back( element_name );
	} else if ( top_level_elements_[ element_name ] != definition ) {
		throw utility::excn::EXCN_Msg_Exception( "Name collision in creation of XML Schema definition: top level element with name \"" +
			element_name + "\" already has been defined.\nOld definition:\n" + top_level_elements_[ element_name ] +
			"New definition:\n" + definition );
	}
}

bool XMLSchemaDefinition::has_top_level_element( std::string const & element_name ) const
{
	return top_level_elements_.find( element_name ) != top_level_elements_.end();
}


std::string XMLSchemaDefinition::full_definition() const
{
	std::ostringstream oss;
	oss << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">\n\n";
	for ( std::list< std::string >::const_iterator iter = elements_in_order_.begin(),
			iter_end = elements_in_order_.end(); iter != iter_end; ++iter ) {
		oss << top_level_elements_.find( *iter )->second << "\n";
	}
	oss << "</xs:schema>\n";
	return oss.str();
}

void XMLSchemaDefinition::validate_new_top_level_element( std::string const & element_name, std::string const & definition )
{
	// parse as TagOP, make sure that the "name" option is the same as the element_name
	std::istringstream def_stream( definition );
	utility::tag::TagCOP element_tag = utility::tag::Tag::create( def_stream );
	if ( ! element_tag->hasOption( "name" ) ) {
		throw utility::excn::EXCN_Msg_Exception( "top level tag \"" + element_name + "\" does not have a name attribute (a.k.a. option)." );
	}
	if ( element_tag->getOption< std::string >( "name" ) != element_name ) {
		throw utility::excn::EXCN_Msg_Exception( "top level tag with presumed name \"" + element_name + "\" has an actual name of \""
			+ element_tag->getOption< std::string >( "name" ) + "\"" );
	}
}

}
}
