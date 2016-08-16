// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/utility/tag/XMLSchemaGeneration.cc
/// @brief  Definitions of the classes used to define an XML Schema
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <utility/tag/XMLSchemaGeneration.hh>

// Package headers
#include <utility/tag/Tag.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

namespace utility {
namespace tag {

std::string name_for_common_type( XMLSchemaCommonType common_type )
{
	switch ( common_type ) {
	case xsct_int_cslist : return "int_cslist";
	case xsct_int_wsslist : return "int_wsslist";
	case xsct_real_cslist : return "real_cslist";
	case xsct_real_wsslist : return "real_wsslist";
	case xsct_non_negative_integer : return "non_negative_integer";
	case xsct_rosetta_bool : return "rosetta_bool";
	case xsct_bool_wsslist : return "bool_wsslist";
	case xsct_none :
		throw utility::excn::EXCN_Msg_Exception( "Error in requesting name for xsct_none;" );
		break;
	}
	return "Error in requesting name for common type";
}

std::ostream & operator << ( std::ostream & os, XMLSchemaCommonType common_type )
{
	os << name_for_common_type( common_type );
	return os;
}


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
	type_( xs_string ),
	common_type_( xsct_none )
{}

XMLSchemaType::XMLSchemaType( XMLSchemaDataType setting ) :
	type_( setting ),
	common_type_( xsct_none )
{}

XMLSchemaType::XMLSchemaType( XMLSchemaCommonType setting ) :
	type_( xs_common ),
	common_type_( setting )
{}

XMLSchemaType::XMLSchemaType( std::string const & custom_type ) :
	type_( xs_custom ),
	common_type_( xsct_none ),
	custom_type_name_( custom_type )
{}

XMLSchemaType::XMLSchemaType( char const * custom_type ) :
	type_( xs_custom ),
	common_type_( xsct_none ),
	custom_type_name_( custom_type )
{}

void XMLSchemaType::type( XMLSchemaDataType setting ) { type_ = setting; common_type_ = xsct_none; }
void XMLSchemaType::common_type( XMLSchemaCommonType setting ) { type_ = xs_common; common_type_ = setting; }
void XMLSchemaType::custom_type_name( std::string const & setting ) { type_ = xs_custom; common_type_ = xsct_none; custom_type_name_ = setting; }

std::string XMLSchemaType::type_name() const {
	switch ( type_ ) {
	case xs_string : return "xs:string"; break;
	case xs_decimal : return "xs:decimal"; break;
	case xs_integer : return "xs:integer"; break;
	case xs_boolean : return "xs:boolean"; break;
	case xs_date : return "xs:date"; break;
	case xs_time : return "xs:time"; break;
	case xs_common : return name_for_common_type( common_type_ ); break;
	case xs_custom : return custom_type_name_;
	}
	return "error!"; // appease compiler
}

XMLSchemaDataType XMLSchemaType::type() const { return type_; }
XMLSchemaCommonType XMLSchemaType::common_type() const { return common_type_; }

XMLSchemaAttribute::XMLSchemaAttribute() :
	is_required_( false )
{}

XMLSchemaAttribute::XMLSchemaAttribute(
	std::string const & name,
	XMLSchemaType type
) :
	name_( name ),
	type_( type ),
	is_required_( false )
{}

// XMLSchemaAttribute::XMLSchemaAttribute(
//  std::string const & name,
//  XMLSchemaType type,
//  std::string const & default_value
// ) :
//  name_( name ),
//  type_( type ),
//  default_value_( default_value ),
//  is_required_( false )
// {}

XMLSchemaAttribute XMLSchemaAttribute::required_attribute(
	std::string const & name,
	XMLSchemaType type
)
{
	XMLSchemaAttribute attribute( name, type );
	attribute.is_required_ = true;
	return attribute;
}

XMLSchemaAttribute
XMLSchemaAttribute::required_attribute(
	std::string const & name,
	XMLSchemaType type,
	std::string const & description
)
{
	XMLSchemaAttribute attribute( name, type );
	attribute.is_required( true ).description( description );
	return attribute;
}


XMLSchemaAttribute
XMLSchemaAttribute::attribute_w_default(
	std::string const & name,
	XMLSchemaType type,
	std::string const & default_value
)
{
	XMLSchemaAttribute attribute( name, type );
	attribute.default_value( default_value );
	return attribute;
}


XMLSchemaAttribute
XMLSchemaAttribute::attribute_w_default(
	std::string const & name,
	XMLSchemaType type,
	std::string const & default_value,
	std::string const & description
)
{
	XMLSchemaAttribute attribute( name, type );
	attribute.default_value( default_value ).description( description );
	return attribute;
}

XMLSchemaAttribute & XMLSchemaAttribute::name( std::string const & setting ) { name_ = setting; return *this; }
XMLSchemaAttribute & XMLSchemaAttribute::type( XMLSchemaType setting ) { type_ = setting; return *this; }
XMLSchemaAttribute & XMLSchemaAttribute::is_required( bool setting ) { is_required_ = setting; return *this; }
XMLSchemaAttribute & XMLSchemaAttribute::default_value( std::string const & setting ) { default_value_ = setting; return *this; }
XMLSchemaAttribute & XMLSchemaAttribute::description( std::string const & setting ) { description_ = setting; return *this; }

XMLSchemaType const &
XMLSchemaAttribute::type() const { return type_; }

std::string const & XMLSchemaAttribute::element_name() const {
	return name_;
}

void XMLSchemaAttribute::write_definition( int indentation, std::ostream & os ) const {
	indent_w_spaces( indentation, os );
	os << "<xs:attribute name=\"" << name_ << "\" type=\"";
	os << type_.type_name();
	os << "\"";
	if ( default_value_ != "" ) os << " default=\"" << default_value_ << "\"";
	if ( is_required_ ) os << " use=\"required\"";
	if ( description_ != "" ) os << " desc=\"" << description_ << "\"";
	os << "/>\n";
}

std::string real_regex_pattern() {
	return "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?";
}

/// @brief the set of all strings recognized as booleans according to Rosetta (in utility/string_util.cc)
std::string rosetta_bool_string() {
	return "true|True|TRUE|t|1|on|On|ON|y|Y|yes|Yes|YES|false|False|FALSE|f|F|0|off|Off|OFF|n|N|no|No|NO";
}

void
activate_common_simple_type(
	utility::tag::XMLSchemaDefinition & xsd,
	XMLSchemaCommonType common_type
)
{
	if ( common_type == xsct_int_cslist ) {
		XMLSchemaRestriction int_cslist;
		int_cslist.name( name_for_common_type( common_type ));
		int_cslist.base_type( xs_string );
		int_cslist.add_restriction( xsr_pattern, "[-+]?[0-9]+(,[-+]?[0-9]+)*" );

		xsd.add_top_level_element( int_cslist );
	} else  if ( common_type == xsct_int_wsslist ) {
		XMLSchemaRestriction int_wsslist;
		int_wsslist.name( name_for_common_type( common_type ));
		int_wsslist.base_type( xs_string );
		int_wsslist.add_restriction( xsr_pattern, "[-+]?[0-9]+(\\s+[-+]?[0-9]+)*" );

		xsd.add_top_level_element( int_wsslist );
	} else if ( common_type == xsct_real_cslist ) {
		XMLSchemaRestriction real_cslist;
		real_cslist.name( name_for_common_type( common_type ));
		real_cslist.base_type( xs_string );
		//real_cslist.add_restriction( xsr_pattern, "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?(,[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)" );
		real_cslist.add_restriction( xsr_pattern, real_regex_pattern() + "(," + real_regex_pattern() + ")" );
		xsd.add_top_level_element( real_cslist );
	} else if ( common_type == xsct_real_wsslist ) {
		XMLSchemaRestriction real_wsslist;
		real_wsslist.name( name_for_common_type( common_type ));
		real_wsslist.base_type( xs_string );
		real_wsslist.add_restriction( xsr_pattern, real_regex_pattern() + "(\\s+" + real_regex_pattern() + ")" );
		xsd.add_top_level_element( real_wsslist );
	} else if ( common_type == xsct_bool_wsslist ) {
		XMLSchemaRestriction bool_wsslist;
		bool_wsslist.name( name_for_common_type( common_type ));
		bool_wsslist.base_type( xs_string );
		bool_wsslist.add_restriction( xsr_pattern, rosetta_bool_string() + "(\\s+" + rosetta_bool_string() + ")" );
		xsd.add_top_level_element( bool_wsslist );
	} else if ( common_type == xsct_non_negative_integer ) {
		XMLSchemaRestriction nonneg_int;
		nonneg_int.name( name_for_common_type( common_type ));
		nonneg_int.base_type( xs_integer );
		nonneg_int.add_restriction( xsr_minInclusive, "0" );
		xsd.add_top_level_element( nonneg_int );
	} else if ( common_type == xsct_rosetta_bool ) {
		XMLSchemaRestriction boolean;
		boolean.name( name_for_common_type( common_type ));
		boolean.base_type( xs_string );
		boolean.add_restriction( xsr_pattern, rosetta_bool_string() );
		xsd.add_top_level_element( boolean );
	}

	/* else if ( common_type == xsct_zero_or_one ) {
	XMLSchemaRestriction zero_or_one;
	zero_or_one.name( name_for_common_type( common_type ));
	zero_or_one.base_type( xs_integer );
	zero_or_one.add_restriction( xsr_minInclusive, "0" );
	zero_or_one.add_restriction( xsr_maxInclusive, "1" );
	xsd.add_top_level_element( zero_or_one );
	}*/

}


void XMLSchemaAttribute::prepare_for_output( XMLSchemaDefinition & xsd ) const
{
	if ( type_.common_type() != xsct_none ) {
		activate_common_simple_type( xsd, type_.common_type() );
	}
}

AttributeList &
operator + ( AttributeList & attributes, XMLSchemaAttribute const & attribute_to_append )
{
	attributes.push_back( attribute_to_append );
	return attributes;
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

std::string const & XMLSchemaRestriction::element_name() const {
	return name_;
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

void XMLSchemaRestriction::prepare_for_output( XMLSchemaDefinition & xsd ) const
{
	if ( base_type_.common_type() != xsct_none ) {
		activate_common_simple_type( xsd, base_type_.common_type() );
	}
}

XMLSchemaParticle::XMLSchemaParticle() : min_occurs_( xsminmax_unspecified ), max_occurs_( xsminmax_unspecified ) {}
XMLSchemaParticle::~XMLSchemaParticle() {}

XMLSchemaParticle & XMLSchemaParticle::min_occurs( int setting )
{
	min_occurs_ = setting;
	return *this;
}

XMLSchemaParticle & XMLSchemaParticle::max_occurs( int setting )
{
	max_occurs_ = setting;
	return *this;
}

int XMLSchemaParticle::min_occurs() const { return min_occurs_; }
int XMLSchemaParticle::max_occurs() const { return max_occurs_; }

XMLSchemaModelGroup::XMLSchemaModelGroup() :
	type_( xsmgt_sequence )
{}

XMLSchemaModelGroup::~XMLSchemaModelGroup() {}

XMLSchemaModelGroup::XMLSchemaModelGroup( XMLSchemaModelGroupType type ) :
	type_( type )
{}

XMLSchemaModelGroup::XMLSchemaModelGroup( XMLSchemaModelGroupType type, std::list< XMLSchemaParticleCOP > const & particles ) :
	type_( type ),
	particles_( particles )
{
	validate_content();
}

XMLSchemaModelGroup::XMLSchemaModelGroup( std::string const & group_name ) :
	type_( xsmgt_group ),
	group_name_( group_name )
{}



XMLSchemaModelGroup & XMLSchemaModelGroup::type( XMLSchemaModelGroupType type )
{
	type_ = type;
	validate_content();
	return *this;
}

XMLSchemaModelGroup & XMLSchemaModelGroup::append_particle( XMLSchemaParticleCOP particle )
{
	particles_.push_back( particle );
	validate_content();
	return *this;
}

XMLSchemaModelGroup &
XMLSchemaModelGroup::append_particles( std::list< XMLSchemaParticleCOP > const & particles )
{
	for ( std::list< XMLSchemaParticleCOP >::const_iterator iter = particles.begin(); iter != particles.end(); ++iter ) {
		particles_.push_back( *iter );
	}
	validate_content();
	return *this;
}

XMLSchemaModelGroup & XMLSchemaModelGroup::group_name( std::string const & name )
{
	type_ = xsmgt_group;
	group_name_ = name;
	return *this;
}

std::string const &
XMLSchemaModelGroup::element_name() const {
	return group_name_;
}

void XMLSchemaModelGroup::write_definition( int indentation, std::ostream & os ) const
{
	if ( type_ != xsmgt_group && particles_.empty() ) return;

	indent_w_spaces( indentation, os );
	switch ( type_ ) {
	case xsmgt_sequence :
	case xsmgt_choice :
	case xsmgt_all :
		os << "<" << type_;
		break;
	case xsmgt_group :
		os << "<" << type_;
		if ( group_name_ != "" ) {
			if ( particles_.empty() ) {
				os << " ref=\"";
			} else {
				os << " name=\"";
			}
			os << group_name_ << "\"";
		}
		break;
	}

	if ( ! is_group_holding_all() ) { write_min_occurs_max_occurs_if_necessary( min_occurs(), max_occurs(), os ); }

	if ( particles_.empty() ) {
		os << "/>\n";
		return;
	} else {
		os << ">\n";
	}

	for ( std::list< XMLSchemaParticleCOP >::const_iterator
			iter = particles_.begin(); iter != particles_.end(); ++iter ) {
		(*iter)->write_definition( indentation+1, os );
	}

	indent_w_spaces( indentation, os );
	os << "</" << type_ << ">\n";
}

void XMLSchemaModelGroup::prepare_for_output( XMLSchemaDefinition & xsd ) const
{
	for ( std::list< XMLSchemaParticleCOP >::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
		(*iter)->prepare_for_output( xsd );
	}
}


void
XMLSchemaModelGroup::validate_content() const
{
	bool found_all = false;
	for ( std::list< XMLSchemaParticleCOP >::const_iterator
			iter = particles_.begin(); iter != particles_.end(); ++iter ) {
		XMLSchemaParticleCOP particle(*iter);
		XMLSchemaModelGroupCOP model_group = utility::pointer::dynamic_pointer_cast< XMLSchemaModelGroup const >( particle );
		if ( model_group ) {
			if ( model_group->type_ == xsmgt_all ) {
				if ( type_ != xsmgt_all && type_ != xsmgt_group ) {
					throw utility::excn::EXCN_Msg_Exception( "In XML Schema, xs:all model groups may not be nested beneath anything besides an xs:all" );
				}
				found_all = true;
			} else if ( model_group->is_group_holding_all() ) {
				throw utility::excn::EXCN_Msg_Exception( "In XML Schema, xs:groups holding an xs:all group may only be nested directly beneath a complexType" );
			} else {
				if ( type_ == xsmgt_all ) {
					throw utility::excn::EXCN_Msg_Exception( "In XML Schema, only xs:all model groups may be nested beneath an xs:all group" );
				}
			}
		} // else -- we're looking at an XMLSchemaElement, and it's always ok to hold an XMLSchemaElement.
	}
	if ( found_all && type_ == xsmgt_group ) {
		if ( particles_.size() > 1 ) {
			throw utility::excn::EXCN_Msg_Exception( "In XML Schema, an xs:all may be the only child of an xs:group if it is a child at all" );
		}
	}
}

bool XMLSchemaModelGroup::is_group_holding_all() const
{
	if ( type_ != xsmgt_group ) return false;
	for ( std::list< XMLSchemaParticleCOP >::const_iterator
			iter = particles_.begin(); iter != particles_.end(); ++iter ) {
		XMLSchemaParticleCOP particle(*iter);
		XMLSchemaModelGroupCOP model_group = utility::pointer::dynamic_pointer_cast< XMLSchemaModelGroup const >( particle );
		if ( model_group ) {
			if ( model_group->type_ == xsmgt_all ) return true;
		}
	}
	return false;
}


std::string xs_model_group_name( XMLSchemaModelGroupType xsmgt ) {
	switch ( xsmgt ) {
	case xsmgt_sequence : return "xs:sequence";
	case xsmgt_all :      return "xs:all";
	case xsmgt_choice :   return "xs:choice";
	case xsmgt_group :    return "xs:group";
	}
	return "error!"; // appease compiler
}

std::ostream & operator << ( std::ostream & os, XMLSchemaModelGroupType type ) {
	os << xs_model_group_name( type );
	return os;
}


XMLSchemaComplexType::XMLSchemaComplexType() {}

XMLSchemaComplexType & XMLSchemaComplexType::name( std::string const & setting ) { name_ = setting; return *this; }
XMLSchemaComplexType & XMLSchemaComplexType::set_model_group( XMLSchemaModelGroupCOP model_group ) { model_group_ = model_group; return *this; }
XMLSchemaComplexType & XMLSchemaComplexType::add_attribute( XMLSchemaAttribute attribute ) { attributes_.push_back( attribute ); return *this; }
XMLSchemaComplexType & XMLSchemaComplexType::add_attributes( AttributeList const & attributes ) {
	for ( AttributeList::const_iterator iter = attributes.begin(); iter != attributes.end(); ++iter ) {
		add_attribute( *iter );
	}
	return *this;
}

std::string const & XMLSchemaComplexType::element_name() const {
	return name_;
}

void XMLSchemaComplexType::write_definition( int indentation, std::ostream & os ) const
{

	// write out header
	indent_w_spaces( indentation, os );
	os << "<xs:complexType";
	if ( name_ != "" ) {
		os << " name=\"" << name_ << "\"";
	}
	os << " mixed=\"true\">\n";

	// write out model group
	if ( model_group_ ) {
		model_group_->write_definition( indentation+1, os );
	}

	// write out attributes
	for ( std::list< XMLSchemaAttribute >::const_iterator
			iter = attributes_.begin(),
			iter_end = attributes_.end();
			iter != iter_end; ++iter ) {
		iter->write_definition( indentation+1, os );
	}

	// write out footer
	indent_w_spaces( indentation, os );
	os << "</xs:complexType>\n";
}

void XMLSchemaComplexType::prepare_for_output( XMLSchemaDefinition & xsd ) const
{
	if ( model_group_ ) {
		model_group_->prepare_for_output( xsd );
	}
	for ( std::list< XMLSchemaAttribute >::const_iterator iter = attributes_.begin(); iter != attributes_.end(); ++iter ) {
		iter->prepare_for_output( xsd );
	}

}

///////////////////////////////////////////// XMLSchemaElement ///////////////////////////////////////////////////

XMLSchemaElement::XMLSchemaElement() : category_( xs_element_is_type_reference ) {}
XMLSchemaElement & XMLSchemaElement::name( std::string const & setting ) { name_ = setting; return *this; }
XMLSchemaElement & XMLSchemaElement::set_abstract() { category_ = xs_element_is_abstract; return *this; }
XMLSchemaElement & XMLSchemaElement::substitution_group( std::string const & setting ) { substitution_group_ = setting; return *this; }
XMLSchemaElement & XMLSchemaElement::type_name( XMLSchemaType setting ) { type_name_ = setting; category_ = xs_element_is_type_reference; return *this; }
XMLSchemaElement & XMLSchemaElement::reference_name( std::string const & setting ) { category_ = xs_element_is_element_reference; name_ = setting; return *this; }

XMLSchemaElement & XMLSchemaElement::element_type_def( XMLSchemaComplexTypeOP setting ) { complex_type_ = setting; category_ = xs_element_is_complex_type_w_definition; return *this; }

std::string const & XMLSchemaElement::element_name() const {
	return name_;
}

void XMLSchemaElement::write_definition( int indentation, std::ostream & os ) const {
	if ( category_ == xs_element_is_abstract ) {
		indent_w_spaces( indentation, os );
		os << "<xs:element name=\"" << name_ << "\" abstract=\"true\"/>\n";
	} else if ( category_ == xs_element_is_type_reference ) {
		indent_w_spaces( indentation, os );
		os << "<xs:element name=\"" << name_ << "\" type=\"" << type_name_.type_name() << "\"" ;
		write_min_occurs_max_occurs_if_necessary( min_occurs(), max_occurs(), os );
		if ( substitution_group_ != "" ) {
			os << " substitutionGroup=\"" << substitution_group_ << "\"";
		}
		os << "/>\n";
	} else if ( category_ == xs_element_is_element_reference ) {
		indent_w_spaces( indentation, os );
		os << "<xs:element ref=\"" << name_ << "\"/>\n";
	} else {
		indent_w_spaces( indentation, os );
		os << "<xs:element name=\"" << name_ << "\"";
		write_min_occurs_max_occurs_if_necessary( min_occurs(), max_occurs(), os );
		if ( substitution_group_ != "" ) {
			os << " substitutionGroup=\"" << substitution_group_ << "\"";
		}
		os << ">\n";
		complex_type_->write_definition( indentation+1, os );

		indent_w_spaces( indentation, os );
		os << "</xs:element>\n";
	}
}

void XMLSchemaElement::prepare_for_output( XMLSchemaDefinition & xsd ) const
{
	if ( type_name_.common_type() != xsct_none ) { activate_common_simple_type( xsd, type_name_.common_type() ); }
	if ( complex_type_ ) complex_type_->prepare_for_output( xsd );
}


//////////////////////////////// XMLSchemaDefinition /////////////////////////////////////////////

XMLSchemaDefinition::XMLSchemaDefinition() {}

XMLSchemaDefinition::~XMLSchemaDefinition() {}

void XMLSchemaDefinition::add_top_level_element( XMLSchemaTopLevelElement const & element )
{
	// handle the common simple types by recursing through the input element.
	element.prepare_for_output( *this );

	std::ostringstream oss;
	element.write_definition( 0, oss );
	std::string definition = oss.str();

	validate_new_top_level_element( element.element_name(), definition );
	if ( top_level_elements_.find( element.element_name() ) == top_level_elements_.end() ) {
		top_level_elements_[ element.element_name() ] = definition;
		elements_in_order_.push_back( element.element_name() );
	} else if ( top_level_elements_[ element.element_name() ] != definition ) {
		throw utility::excn::EXCN_Msg_Exception( "Name collision in creation of XML Schema definition: top level element with name \"" +
			element.element_name() + "\" already has been defined.\nOld definition:\n" + top_level_elements_[ element.element_name() ] +
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

XMLSchemaSimpleSubelementList::XMLSchemaSimpleSubelementList() {}

XMLSchemaSimpleSubelementList::~XMLSchemaSimpleSubelementList() {}

XMLSchemaSimpleSubelementList::XMLSchemaSimpleSubelementList( XMLSchemaSimpleSubelementList const & src ) :
	ct_naming_func_for_simple_subelements_( src.ct_naming_func_for_simple_subelements_ ),
	elements_( src.elements_ )
{}

XMLSchemaSimpleSubelementList &
XMLSchemaSimpleSubelementList::operator = ( XMLSchemaSimpleSubelementList const & rhs )
{
	if ( this != &rhs ) {
		ct_naming_func_for_simple_subelements_ = rhs.ct_naming_func_for_simple_subelements_;
		elements_ = rhs.elements_;
	}
	return *this;
}

XMLSchemaSimpleSubelementList &
XMLSchemaSimpleSubelementList::complex_type_naming_func( DerivedNameFunction const & naming_function )
{
	ct_naming_func_for_simple_subelements_ = naming_function;
	return *this;
}

XMLSchemaSimpleSubelementList &
XMLSchemaSimpleSubelementList::add_simple_subelement( std::string const & name, AttributeList const & attributes )
{
	ElementSummary summary;
	summary.element_type = ElementSummary::ct_simple;
	summary.element_name = name;
	summary.ct_name = "";
	summary.attributes = attributes;
	summary.min_or_max_occurs_set = false;
	summary.min_occurs = xsminmax_unspecified;
	summary.max_occurs = xsminmax_unspecified;
	elements_.push_back( summary );
	return *this;
}

XMLSchemaSimpleSubelementList &
XMLSchemaSimpleSubelementList::add_simple_subelement(
	std::string const & name,
	AttributeList const & attributes,
	int min_occurs,
	int max_occurs
)
{
	ElementSummary summary;
	summary.element_type = ElementSummary::ct_simple;
	summary.element_name = name;
	summary.ct_name = "";
	summary.attributes = attributes;
	summary.min_or_max_occurs_set = true;
	summary.min_occurs = min_occurs;
	summary.max_occurs = max_occurs;
	elements_.push_back( summary );
	return *this;
}

XMLSchemaSimpleSubelementList &
XMLSchemaSimpleSubelementList::add_already_defined_subelement(
	std::string const & name,
	DerivedNameFunction const & ct_naming_function
)
{
	ElementSummary summary;
	summary.element_type = ElementSummary::ct_ref;
	summary.element_name = name;
	summary.ct_name = ct_naming_function( name );
	summary.min_or_max_occurs_set = false;
	summary.min_occurs = xsminmax_unspecified;
	summary.max_occurs = xsminmax_unspecified;
	elements_.push_back( summary );
	return *this;
}

XMLSchemaSimpleSubelementList &
XMLSchemaSimpleSubelementList::add_already_defined_subelement(
	std::string const & name,
	DerivedNameFunction const & ct_naming_function,
	int min_occurs,
	int max_occurs
)
{
	ElementSummary summary;
	summary.element_type = ElementSummary::ct_ref;
	summary.element_name = name;
	summary.ct_name = ct_naming_function( name );
	summary.min_or_max_occurs_set = true;
	summary.min_occurs = min_occurs;
	summary.max_occurs = max_occurs;
	elements_.push_back( summary );
	return *this;
}

XMLSchemaSimpleSubelementList &
XMLSchemaSimpleSubelementList::add_group_subelement(
	NameFunction const & group_name_function
)
{
	ElementSummary summary;
	summary.element_type = ElementSummary::ct_group;
	summary.element_name = "";
	summary.ct_name = group_name_function();
	summary.min_or_max_occurs_set = false;
	summary.min_occurs = xsminmax_unspecified;
	summary.max_occurs = xsminmax_unspecified;
	elements_.push_back( summary );
	return *this;
}

XMLSchemaSimpleSubelementList &
XMLSchemaSimpleSubelementList::add_group_subelement(
	NameFunction const & group_name_function,
	int min_occurs,
	int max_occurs
)
{
	ElementSummary summary;
	summary.element_type = ElementSummary::ct_group;
	summary.element_name = "";
	summary.ct_name = group_name_function();
	summary.min_or_max_occurs_set = true;
	summary.min_occurs = min_occurs;
	summary.max_occurs = max_occurs;
	elements_.push_back( summary );
	return *this;
}

bool
XMLSchemaSimpleSubelementList::simple_element_naming_func_has_been_set() const
{
	return ct_naming_func_for_simple_subelements_;
}

XMLSchemaSimpleSubelementList::DerivedNameFunction
XMLSchemaSimpleSubelementList::naming_func() const
{
	return ct_naming_func_for_simple_subelements_;
}

std::string
XMLSchemaSimpleSubelementList::complex_typename_for_element( std::string const & element_name ) const
{
	return ct_naming_func_for_simple_subelements_( element_name );
}

std::list< XMLSchemaSimpleSubelementList::ElementSummary > const &
XMLSchemaSimpleSubelementList::element_list() const
{
	return elements_;
}


class XMLSchemaComplexTypeGeneratorImpl : public utility::pointer::ReferenceCount
{
public:
	typedef boost::function< std::string ( std::string const & ) >     DerivedNameFunction;
	typedef boost::function< std::string () >                          NameFunction;
	typedef std::list< XMLSchemaSimpleSubelementList::ElementSummary > ElementSummaries;

	enum SetOfSubelementsBehavior { ss_repeatable, ss_optional, ss_required, ss_pick_one_opt, ss_pick_one_req };

	typedef std::list< std::pair< XMLSchemaSimpleSubelementList, SetOfSubelementsBehavior > > SubelementSets;

public:
	XMLSchemaComplexTypeGeneratorImpl();
	~XMLSchemaComplexTypeGeneratorImpl();

	void element_name( std::string const & );
	void complex_type_naming_func( DerivedNameFunction const & naming_function );
	void add_attribute( XMLSchemaAttribute const & attribute );
	void add_attributes( AttributeList const & attributes );

	/// @brief set subelements as repeatable (and optional); setting the sub-elements replaces any sub-elements that were previously set.
	/// These repeatable subelements are allowed to appear in any order from the point of view of the XML Schema, which is
	/// not to say that the order in which they appear cannot matter to the code reading these subelements.  The group_name
	/// string must be unique to Rosetta's XSD; it is needed in order to define an xs:group and then refer to that group within
	/// this xs:complexType.
	void set_subelements_repeatable(
		XMLSchemaSimpleSubelementList const & subelements,
		int min_occurs,
		int max_occurs
	);

	void set_subelements_pick_one_required( XMLSchemaSimpleSubelementList const & subelements );
	void set_subelements_pick_one_optional( XMLSchemaSimpleSubelementList const & subelements );

	/// @brief set subelements as single-appearance (and required); setting the sub-elements replaces any sub-elements that were previously set.
	/// These single-appearance subelements are allowed to appear in any order from the point of view of the XML Schema, which is
	/// not to say that the order in which they appear cannot matter to the code reading these subelements.
	void set_subelements_single_appearance_required( XMLSchemaSimpleSubelementList const & subelements );

	/// @brief set subelements as single-appearance (and optional); setting the sub-elements replaces any sub-elements that were previously set.
	/// These single-appearance subelements are allowed to appear in any order from the point of view of the XML Schema, which is
	/// not to say that the order in which they appear cannot matter to the code reading these subelements.
	void set_subelements_single_appearance_optional( XMLSchemaSimpleSubelementList const & subelements );


	void set_subelements_single_appearance_required_and_ordered(
		XMLSchemaSimpleSubelementList const & subelements
	);

	/// @brief Add a subelement list to a growing set of ordered subelements, where elements in
	/// this set are labeled as being repeatable.
	/// This function corresponds to case 8 in the list of behaviors above.
	void add_ordered_subelement_set_as_repeatable(
		XMLSchemaSimpleSubelementList const & subelements
	);

	/// @brief Add a subelement list to a growing set of ordered subelements, where elements in
	/// this set are labeled as being optional.  There should be only a single element in the
	/// input subelement list.
	/// This function corresponds to case 8 in the list of behaviors above.
	void add_ordered_subelement_set_as_optional(
		XMLSchemaSimpleSubelementList const & subelements
	);

	/// @brief Add a subelement list to a growing set of ordered subelements, where elements in
	/// this set are labeled as being requried.  There should be only a single element in the
	/// input subelement list.
	/// This function corresponds to case 8 in the list of behaviors above.
	void add_ordered_subelement_set_as_required(
		XMLSchemaSimpleSubelementList const & subelements
	);

	/// @brief Add a subelement list to a growing set of ordered subelements, where elements in
	/// this set are labeled as "pick one or none".  There should be more than one element in the
	/// input subelement list.
	/// This function corresponds to case 8 in the list of behaviors above.
	void add_ordered_subelement_set_as_pick_one_optional(
		XMLSchemaSimpleSubelementList const & subelements
	);

	/// @brief Add a subelement list to a growing set of ordered subelements, where elements in
	/// this set are labeled as "pick exactly one".  There should be more than one element in the
	/// input subelement list.
	/// This function corresponds to case 8 in the list of behaviors above.
	void add_ordered_subelement_set_as_pick_one_required(
		XMLSchemaSimpleSubelementList const & subelements
	);

	CTGenSubelementBehavior subelement_behavior() const;

	void write_complex_type_to_schema( XMLSchemaDefinition & xsd );

private:

	bool only_one_subelement_and_that_subelement_is_a_ct_group() const;
	bool only_one_subelement_and_that_subelement_is_a_ct_group( XMLSchemaSimpleSubelementList const & subelements ) const;
	void prepare_subelement_repeatable( XMLSchemaDefinition & xsd, XMLSchemaComplexType & complex_type );
	void prepare_subelement_choice_req( XMLSchemaDefinition & xsd, XMLSchemaComplexType & complex_type );
	void prepare_subelement_choice_opt( XMLSchemaDefinition & xsd, XMLSchemaComplexType & complex_type );
	void prepare_subelement_single_instance_required( XMLSchemaDefinition & xsd, XMLSchemaComplexType & complex_type );
	void prepare_subelement_single_instance_optional( XMLSchemaDefinition & xsd, XMLSchemaComplexType & complex_type );
	void prepare_subelement_single_instance_required_ordered( XMLSchemaDefinition & xsd, XMLSchemaComplexType & complex_type );
	void prepare_sequence_of_subelement_sets( XMLSchemaDefinition & xsd, XMLSchemaComplexType & complex_type );
	XMLSchemaParticleOP create_subelement( XMLSchemaSimpleSubelementList::ElementSummary const & summary, XMLSchemaDefinition & xsd );

private:

	CTGenSubelementBehavior subelement_behavior_;
	std::string element_name_;
	DerivedNameFunction complex_type_naming_function_;

	XMLSchemaSimpleSubelementList subelements_;
	SubelementSets subelement_sets_;

	AttributeList attributes_;
	int repeatable_min_occurs_;
	int repeatable_max_occurs_;
};

XMLSchemaComplexTypeGeneratorImpl::XMLSchemaComplexTypeGeneratorImpl() :
	subelement_behavior_( se_none ),
	repeatable_min_occurs_( 0 ),
	repeatable_max_occurs_( xsminmax_unbounded )
{}


XMLSchemaComplexTypeGeneratorImpl::~XMLSchemaComplexTypeGeneratorImpl() {}

void XMLSchemaComplexTypeGeneratorImpl::element_name( std::string const & element_name ) {
	element_name_ = element_name;
}

void XMLSchemaComplexTypeGeneratorImpl::complex_type_naming_func(
	DerivedNameFunction const & naming_function
)
{
	complex_type_naming_function_ = naming_function;
}

void XMLSchemaComplexTypeGeneratorImpl::add_attribute( XMLSchemaAttribute const & attribute )
{
	attributes_.push_back( attribute );
}

void XMLSchemaComplexTypeGeneratorImpl::add_attributes( AttributeList const & attributes )
{
	for ( AttributeList::const_iterator iter = attributes.begin(); iter != attributes.end(); ++iter ) {
		attributes_.push_back( *iter );
	}
}

void
XMLSchemaComplexTypeGeneratorImpl::set_subelements_repeatable(
	XMLSchemaSimpleSubelementList const & subelements,
	int min_occurs,
	int max_occurs
)
{
	subelement_behavior_ = se_repeatable;
	subelements_ = subelements;
	repeatable_min_occurs_ = min_occurs;
	repeatable_max_occurs_ = max_occurs;
	subelement_sets_.clear();
}

void
XMLSchemaComplexTypeGeneratorImpl::set_subelements_pick_one_required( XMLSchemaSimpleSubelementList const & subelements )
{
	subelement_behavior_ = se_choice_req;
	subelements_ = subelements;
	subelement_sets_.clear();
}

void
XMLSchemaComplexTypeGeneratorImpl::set_subelements_pick_one_optional( XMLSchemaSimpleSubelementList const & subelements )
{
	subelement_behavior_ = se_choice_opt;
	subelements_ = subelements;
	subelement_sets_.clear();
}


void
XMLSchemaComplexTypeGeneratorImpl::set_subelements_single_appearance_required( XMLSchemaSimpleSubelementList const & subelements )
{
	if ( subelements.element_list().size() > 1 ) {
		ElementSummaries const & elements( subelements_.element_list() );

		for ( ElementSummaries::const_iterator iter = elements.begin(); iter != elements.end(); ++iter ) {
			if ( iter->element_type == XMLSchemaSimpleSubelementList::ElementSummary::ct_group ) {
				throw utility::excn::EXCN_Msg_Exception( "set_subelement_single_appearance_required cannot be used with"
					" more than one element if any of the subelements are \"group\" references\n" );
			}
		}
	}

	subelement_behavior_ = se_single_req;
	subelements_ = subelements;
	subelement_sets_.clear();
}

void
XMLSchemaComplexTypeGeneratorImpl::set_subelements_single_appearance_optional( XMLSchemaSimpleSubelementList const & subelements )
{
	if ( subelements.element_list().size() > 1 ) {
		ElementSummaries const & elements( subelements_.element_list() );

		for ( ElementSummaries::const_iterator iter = elements.begin(); iter != elements.end(); ++iter ) {
			if ( iter->element_type == XMLSchemaSimpleSubelementList::ElementSummary::ct_group ) {
				throw utility::excn::EXCN_Msg_Exception( "set_subelement_single_appearance_optional cannot be used with"
					" more than one element if any of the subelements are \"group\" references\n" );
			}
		}
	}

	subelement_behavior_ = se_single_opt;
	subelements_ = subelements;
	subelement_sets_.clear();
}

void
XMLSchemaComplexTypeGeneratorImpl::set_subelements_single_appearance_required_and_ordered(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	subelement_behavior_ = se_single_req_ordered;
	subelements_ = subelements;
	subelement_sets_.clear();
}

void
XMLSchemaComplexTypeGeneratorImpl::add_ordered_subelement_set_as_repeatable(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	subelement_behavior_ = se_ordered_sets;
	subelement_sets_.push_back( std::make_pair( subelements, ss_repeatable ));
}


void
XMLSchemaComplexTypeGeneratorImpl::add_ordered_subelement_set_as_optional(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	subelement_behavior_ = se_ordered_sets;
	if ( subelements.element_list().size() != 1 ) {
		throw utility::excn::EXCN_Msg_Exception(
			"When calling XMLSchemaComplexTypeGeneratorImpl::add_ordered_subelement_set_as_optional, "
			" the input subelement list should contain only a single subelement." );
	}
	subelement_sets_.push_back( std::make_pair( subelements, ss_optional ));
}

void
XMLSchemaComplexTypeGeneratorImpl::add_ordered_subelement_set_as_required(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	subelement_behavior_ = se_ordered_sets;
	if ( subelements.element_list().size() != 1 ) {
		throw utility::excn::EXCN_Msg_Exception(
			"When calling XMLSchemaComplexTypeGeneratorImpl::add_ordered_subelement_set_as_required, "
			" the input subelement list should contain only a single subelement." );
	}
	subelement_sets_.push_back( std::make_pair( subelements, ss_required ));
}

void
XMLSchemaComplexTypeGeneratorImpl::add_ordered_subelement_set_as_pick_one_optional(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	subelement_behavior_ = se_ordered_sets;
	subelement_sets_.push_back( std::make_pair( subelements, ss_pick_one_opt ));
}

void
XMLSchemaComplexTypeGeneratorImpl::add_ordered_subelement_set_as_pick_one_required(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	subelement_behavior_ = se_ordered_sets;
	subelement_sets_.push_back( std::make_pair( subelements, ss_pick_one_req ));
}

CTGenSubelementBehavior XMLSchemaComplexTypeGeneratorImpl::subelement_behavior() const
{
	return subelement_behavior_;
}

void
XMLSchemaComplexTypeGeneratorImpl::write_complex_type_to_schema( XMLSchemaDefinition & xsd )
{
	debug_assert( element_name_ != "" );
	debug_assert( complex_type_naming_function_ );

	XMLSchemaComplexType complex_type;
	complex_type.name( complex_type_naming_function_( element_name_ ));

	switch ( subelement_behavior_ ) {
	case se_none :
		break;
	case se_repeatable :
		prepare_subelement_repeatable( xsd, complex_type );
		break;
	case se_choice_req :
		prepare_subelement_choice_req( xsd, complex_type );
		break;
	case se_choice_opt :
		prepare_subelement_choice_opt( xsd, complex_type );
		break;
	case se_single_req :
		prepare_subelement_single_instance_required( xsd, complex_type );
		break;
	case se_single_opt :
		prepare_subelement_single_instance_optional( xsd, complex_type );
		break;
	case se_single_req_ordered :
		prepare_subelement_single_instance_required_ordered( xsd, complex_type );
		break;
	case se_ordered_sets :
		prepare_sequence_of_subelement_sets( xsd, complex_type );
		break;
	}

	complex_type.add_attributes( attributes_ );
	xsd.add_top_level_element( complex_type );
}

bool
XMLSchemaComplexTypeGeneratorImpl::only_one_subelement_and_that_subelement_is_a_ct_group() const
{
	return only_one_subelement_and_that_subelement_is_a_ct_group( subelements_ );
}

bool
XMLSchemaComplexTypeGeneratorImpl::only_one_subelement_and_that_subelement_is_a_ct_group(
	XMLSchemaSimpleSubelementList const & subelements
) const
{
	ElementSummaries const & elements( subelements.element_list() );

	bool one_subelement_that_is_ct_group = false;
	if ( elements.size() == 1 ) {
		XMLSchemaSimpleSubelementList::ElementSummary const & first = *elements.begin();
		if ( first.element_type == XMLSchemaSimpleSubelementList::ElementSummary::ct_group ) {
			one_subelement_that_is_ct_group = true;
		}
	}
	return one_subelement_that_is_ct_group;
}

void XMLSchemaComplexTypeGeneratorImpl::prepare_subelement_repeatable( XMLSchemaDefinition & xsd, XMLSchemaComplexType & complex_type )
{
	ElementSummaries const & elements( subelements_.element_list() );

	if ( only_one_subelement_and_that_subelement_is_a_ct_group() ) {
		XMLSchemaModelGroupOP group( new XMLSchemaModelGroup( elements.begin()->ct_name ));
		group->min_occurs( repeatable_min_occurs_ ).max_occurs( repeatable_max_occurs_ );
		complex_type.set_model_group( group );
	} else {

		XMLSchemaModelGroupOP xs_choice_type( new XMLSchemaModelGroup( xsmgt_choice ) );

		for ( ElementSummaries::const_iterator iter = elements.begin(); iter != elements.end(); ++iter ) {
			if ( iter->min_or_max_occurs_set ) {
				throw utility::excn::EXCN_Msg_Exception( "Subelement named " + iter->element_name + " was initilized with "
					"either min_occurs or max_occurs set, but then handed to the XMLSchemaComplexTypeGenerator through the "
					"set_sublements_repeatable function, which will override the min/max occurence settings." );
			}
			XMLSchemaParticleOP elem = create_subelement( *iter, xsd );
			xs_choice_type->append_particle( elem );
		}
		xs_choice_type->min_occurs( repeatable_min_occurs_ ).max_occurs( repeatable_max_occurs_ );

		// finalize the xs:complexType
		complex_type.set_model_group( xs_choice_type );
	}
}

void
XMLSchemaComplexTypeGeneratorImpl::prepare_subelement_choice_req(
	XMLSchemaDefinition & xsd,
	XMLSchemaComplexType & complex_type
)
{
	ElementSummaries const & elements( subelements_.element_list() );

	XMLSchemaModelGroupOP model_group( new XMLSchemaModelGroup );
	model_group->type( xsmgt_choice );
	for ( ElementSummaries::const_iterator iter = elements.begin(); iter != elements.end(); ++iter ) {
		if ( iter->min_or_max_occurs_set ) {
			throw utility::excn::EXCN_Msg_Exception( "Subelement named " + iter->element_name + " was initilized with "
				"either min_occurs or max_occurs set, but then handed to the XMLSchemaComplexTypeGenerator through the "
				"set_sublements_pick_one function, which will override the min/max occurence settings." );
		}
		XMLSchemaParticleOP elem = create_subelement( *iter, xsd );
		model_group->append_particle( elem );
	}
	complex_type.set_model_group( model_group );
}

void
XMLSchemaComplexTypeGeneratorImpl::prepare_subelement_choice_opt(
	XMLSchemaDefinition & xsd,
	XMLSchemaComplexType & complex_type
)
{
	ElementSummaries const & elements( subelements_.element_list() );

	XMLSchemaModelGroupOP model_group( new XMLSchemaModelGroup );
	model_group->type( xsmgt_choice );
	for ( ElementSummaries::const_iterator iter = elements.begin(); iter != elements.end(); ++iter ) {
		if ( iter->min_or_max_occurs_set ) {
			throw utility::excn::EXCN_Msg_Exception( "Subelement named " + iter->element_name + " was initilized with "
				"either min_occurs or max_occurs set, but then handed to the XMLSchemaComplexTypeGenerator through the "
				"set_sublements_pick_one_optional function, which will override the min/max occurence settings." );
		}
		XMLSchemaParticleOP elem = create_subelement( *iter, xsd );
		model_group->append_particle( elem );
	}
	complex_type.set_model_group( model_group );
}


void XMLSchemaComplexTypeGeneratorImpl::prepare_subelement_single_instance_required(
	XMLSchemaDefinition & xsd,
	XMLSchemaComplexType & complex_type
)
{
	ElementSummaries const & elements( subelements_.element_list() );

	XMLSchemaModelGroupOP model_group( new XMLSchemaModelGroup );
	if ( only_one_subelement_and_that_subelement_is_a_ct_group() ) {
		model_group->type( xsmgt_choice );
	} else {
		model_group->type( xsmgt_all );
	}

	for ( ElementSummaries::const_iterator iter = elements.begin(); iter != elements.end(); ++iter ) {
		XMLSchemaParticleOP elem = create_subelement( *iter, xsd );
		model_group->append_particle( elem );
	}
	complex_type.set_model_group( model_group );
}


void XMLSchemaComplexTypeGeneratorImpl::prepare_subelement_single_instance_optional(
	XMLSchemaDefinition & xsd,
	XMLSchemaComplexType & complex_type
)
{
	ElementSummaries const & elements( subelements_.element_list() );

	XMLSchemaModelGroupOP model_group( new XMLSchemaModelGroup );
	if ( only_one_subelement_and_that_subelement_is_a_ct_group() ) {
		model_group->type( xsmgt_choice );
	} else {
		model_group->type( xsmgt_all );
	}

	for ( ElementSummaries::const_iterator iter = elements.begin(); iter != elements.end(); ++iter ) {
		XMLSchemaParticleOP elem = create_subelement( *iter, xsd );
		if ( iter->min_or_max_occurs_set ) {
			// Error checking. You are allowed to specify a min_occurs of 1 to say that a particular element is required
			// among a list of otherwise optional elements.
			if ( iter->max_occurs != 1 && iter->max_occurs != xsminmax_unspecified ) {
				throw utility::excn::EXCN_Msg_Exception( "Subelement named " + iter->element_name + " was initilized with "
					"max_occurs set to " + utility::to_string( iter->max_occurs )  + ", but then handed to the "
					"XMLSchemaComplexTypeGenerator through the set_sublements_single_instance_optional function, "
					"which only allows a max occurence setting of 1." );
			}
			if ( iter->min_occurs > 1 || iter->min_occurs == xsminmax_unbounded ) {
				throw utility::excn::EXCN_Msg_Exception( "Subelement named " + iter->element_name + " was initilized with "
					"min_occurs set to " + utility::to_string( iter->max_occurs )  + ", but then handed to the "
					"XMLSchemaComplexTypeGenerator through the set_sublements_single_instance_optional function, "
					"which only allows at most a single occurrence." );
			}
		} else {
			elem->min_occurs( 0 ).max_occurs( 1 );
		}
		model_group->append_particle( elem );
	}
	complex_type.set_model_group( model_group );
}


void
XMLSchemaComplexTypeGeneratorImpl::prepare_subelement_single_instance_required_ordered(
	XMLSchemaDefinition & xsd,
	XMLSchemaComplexType & complex_type
)
{
	ElementSummaries const & elements( subelements_.element_list() );

	XMLSchemaModelGroupOP seq( new XMLSchemaModelGroup( xsmgt_sequence ));
	for ( ElementSummaries::const_iterator iter = elements.begin(); iter != elements.end(); ++iter ) {
		XMLSchemaParticleOP elem = create_subelement( *iter, xsd );
		seq->append_particle( elem );
	}
	complex_type.set_model_group( seq );
}

void
XMLSchemaComplexTypeGeneratorImpl::prepare_sequence_of_subelement_sets(
	XMLSchemaDefinition & xsd,
	XMLSchemaComplexType & complex_type
)
{
	XMLSchemaModelGroupOP seq( new XMLSchemaModelGroup( xsmgt_sequence ));
	std::set< std::string > all_element_names;
	for ( SubelementSets::const_iterator iter = subelement_sets_.begin(); iter != subelement_sets_.end(); ++iter ) {
		ElementSummaries const & elements( iter->first.element_list() );
		switch ( iter->second ) {
		case ss_repeatable :
			if ( only_one_subelement_and_that_subelement_is_a_ct_group( iter->first ) ) {
				XMLSchemaModelGroupOP group( new XMLSchemaModelGroup( elements.begin()->ct_name ));
				group->min_occurs( 0 ).max_occurs( xsminmax_unbounded );
				seq->append_particle( group );
			} else {
				XMLSchemaModelGroupOP choice( new XMLSchemaModelGroup( xsmgt_choice ));
				for ( ElementSummaries::const_iterator elem_iter = elements.begin(); elem_iter != elements.end(); ++elem_iter ) {
					if ( elem_iter->min_or_max_occurs_set ) {
						throw utility::excn::EXCN_Msg_Exception( "Subelement named " + elem_iter->element_name + " was initilized with "
							"either min_occurs or max_occurs set, but then handed to the XMLSchemaComplexTypeGenerator through the "
							"add_ordered_sublement_set_as_repeatable function, which will override the min/max occurence settings." );
					}
					if ( all_element_names.count( elem_iter->element_name ) != 0 ) {
						throw utility::excn::EXCN_Msg_Exception( "Subelement named " + elem_iter->element_name + " appears more than"
							" once in subelement lists for complex type for an element named \"" + element_name_ + "\" which"
							" was added through the add_ordered_sublement_set_as_repeatable function" );
					}
					if ( elem_iter->element_type != XMLSchemaSimpleSubelementList::ElementSummary::ct_group ) all_element_names.insert( elem_iter->element_name );
					XMLSchemaParticleOP elem = create_subelement( *elem_iter, xsd );
					choice->append_particle( elem );
				}
				choice->min_occurs( 0 ).max_occurs( xsminmax_unbounded );
				seq->append_particle( choice );
			}
			break;
		case ss_optional :
			{

			assert( elements.size() == 1 ); // this was already checked for in add_ordered_sublement_set_as_optional
			XMLSchemaSimpleSubelementList::ElementSummary elem1 = *elements.begin();
			if ( elem1.min_or_max_occurs_set ) {
				throw utility::excn::EXCN_Msg_Exception( "Subelement named " + elem1.element_name + " was initilized with "
					"either min_occurs or max_occurs set, but then handed to the XMLSchemaComplexTypeGenerator through the "
					"add_ordered_sublement_set_as_optional function, which will override the min/max occurence settings." );
			}
			if ( all_element_names.count( elem1.element_name ) != 0 ) {
				throw utility::excn::EXCN_Msg_Exception( "Subelement named " + elem1.element_name + " appears more than"
					" once in subelement lists for complex type for an element named \"" + element_name_ + "\" which"
					" was added through the add_ordered_sublement_set_as_optional function" );
			}
			if ( elem1.element_type != XMLSchemaSimpleSubelementList::ElementSummary::ct_group ) all_element_names.insert( elem1.element_name );
			XMLSchemaParticleOP elem = create_subelement( elem1, xsd );
			elem->min_occurs( 0 ).max_occurs( 1 );
			seq->append_particle( elem );
		}

			break;
		case ss_required :
			{
			assert( elements.size() == 1 ); // this was already checked for in add_ordered_sublement_set_as_required
			XMLSchemaSimpleSubelementList::ElementSummary elem1 = *elements.begin();
			if ( elem1.min_or_max_occurs_set ) {
				throw utility::excn::EXCN_Msg_Exception( "Subelement named " + elem1.element_name + " was initilized with "
					"either min_occurs or max_occurs set, but then handed to the XMLSchemaComplexTypeGenerator through the "
					"add_ordered_sublement_set_as_required function, which will override the min/max occurence settings." );
			}
			if ( all_element_names.count( elem1.element_name ) != 0 ) {
				throw utility::excn::EXCN_Msg_Exception( "Subelement named " + elem1.element_name + " appears more than"
					" once in subelement lists for complex type for an element named \"" + element_name_ + "\" which"
					" was added through the add_ordered_sublement_set_as_required function" );
			}
			if ( elem1.element_type != XMLSchemaSimpleSubelementList::ElementSummary::ct_group ) all_element_names.insert( elem1.element_name );
			XMLSchemaParticleOP elem = create_subelement( elem1, xsd );
			elem->min_occurs( 1 ).max_occurs( 1 );
			seq->append_particle( elem );
		}
			break;
		case ss_pick_one_opt :
			if ( only_one_subelement_and_that_subelement_is_a_ct_group( iter->first ) ) {
				XMLSchemaModelGroupOP group( new XMLSchemaModelGroup( elements.begin()->ct_name ));
				group->min_occurs( 0 ).max_occurs( 1 );
				seq->append_particle( group );
			} else {
				XMLSchemaModelGroupOP choice( new XMLSchemaModelGroup( xsmgt_choice ));
				for ( ElementSummaries::const_iterator elem_iter = elements.begin(); elem_iter != elements.end(); ++elem_iter ) {
					if ( elem_iter->min_or_max_occurs_set ) {
						throw utility::excn::EXCN_Msg_Exception( "Subelement named " + elem_iter->element_name + " was initilized with "
							"either min_occurs or max_occurs set, but then handed to the XMLSchemaComplexTypeGenerator through the "
							"add_ordered_sublement_set_as_pick_one_or_none function, which will override the min/max occurence settings." );
					}
					if ( all_element_names.count( elem_iter->element_name ) != 0 ) {
						throw utility::excn::EXCN_Msg_Exception( "Subelement named " + elem_iter->element_name + " appears more than"
							" once in subelement lists for complex type for an element named \"" + element_name_ + "\" which"
							" was added through the add_ordered_sublement_set_as_pick_one_or_none function" );
					}
					if ( elem_iter->element_type != XMLSchemaSimpleSubelementList::ElementSummary::ct_group ) all_element_names.insert( elem_iter->element_name );
					XMLSchemaParticleOP elem = create_subelement( *elem_iter, xsd );
					choice->append_particle( elem );
				}
				choice->min_occurs( 0 ).max_occurs( 1 );
				seq->append_particle( choice );
			}
			break;
		case ss_pick_one_req :
			if ( only_one_subelement_and_that_subelement_is_a_ct_group( iter->first ) ) {
				XMLSchemaModelGroupOP group( new XMLSchemaModelGroup( elements.begin()->ct_name ));
				group->min_occurs( 1 ).max_occurs( 1 );
				seq->append_particle( group );
			} else {
				XMLSchemaModelGroupOP choice( new XMLSchemaModelGroup( xsmgt_choice ));
				for ( ElementSummaries::const_iterator elem_iter = elements.begin(); elem_iter != elements.end(); ++elem_iter ) {
					if ( elem_iter->min_or_max_occurs_set ) {
						throw utility::excn::EXCN_Msg_Exception( "Subelement named " + elem_iter->element_name + " was initilized with "
							"either min_occurs or max_occurs set, but then handed to the XMLSchemaComplexTypeGenerator through the "
							"add_ordered_sublement_set_as_pick_one function, which will override the min/max occurence settings." );
					}
					if ( all_element_names.count( elem_iter->element_name ) != 0 ) {
						throw utility::excn::EXCN_Msg_Exception( "Subelement named " + elem_iter->element_name + " appears more than"
							" once in subelement lists for complex type for an element named \"" + element_name_ + "\" which"
							" was added through the add_ordered_sublement_set_as_pick_one function" );
					}
					if ( elem_iter->element_type != XMLSchemaSimpleSubelementList::ElementSummary::ct_group ) all_element_names.insert( elem_iter->element_name );
					XMLSchemaParticleOP elem = create_subelement( *elem_iter, xsd );
					choice->append_particle( elem );
				}
				choice->min_occurs( 1 ).max_occurs( 1 );
				seq->append_particle( choice );
			}
			break;
		}
	}
	complex_type.set_model_group( seq );
}

XMLSchemaParticleOP XMLSchemaComplexTypeGeneratorImpl::create_subelement(
	XMLSchemaSimpleSubelementList::ElementSummary const & summary,
	XMLSchemaDefinition & xsd
)
{
	switch ( summary.element_type ) {
	case XMLSchemaSimpleSubelementList::ElementSummary::ct_simple :
		{
		XMLSchemaElementOP element( new XMLSchemaElement );
		if ( subelements_.simple_element_naming_func_has_been_set() ) {

			// write the complex type to the schema
			XMLSchemaComplexTypeGenerator ctgen;
			ctgen
				.element_name( summary.element_name )
				.complex_type_naming_func( subelements_.naming_func() )
				.add_attributes( summary.attributes );
			ctgen.write_complex_type_to_schema( xsd );

			// initialize the element to refer to the new type
			element->name( summary.element_name );
			element->type_name( subelements_.naming_func()( summary.element_name ));
		} else {
			// if no naming function has been provided, then the type for this element must be provided inline
			XMLSchemaComplexTypeOP ct( new XMLSchemaComplexType );
			ct->add_attributes( summary.attributes );
			element->element_type_def( ct );
			element->name( summary.element_name );
		}
		if ( summary.min_or_max_occurs_set ) {
			if ( summary.min_occurs != xsminmax_unspecified ) { element->min_occurs( summary.min_occurs ); }
		}
		return element;
	}
		break;
	case XMLSchemaSimpleSubelementList::ElementSummary::ct_ref :
		{
		XMLSchemaElementOP element( new XMLSchemaElement );

		element->name( summary.element_name );
		element->type_name( summary.ct_name );
		if ( summary.min_or_max_occurs_set ) {
			if ( summary.min_occurs != xsminmax_unspecified ) { element->min_occurs( summary.min_occurs ); }
			if ( summary.max_occurs != xsminmax_unspecified ) { element->max_occurs( summary.max_occurs ); }
		}
		return element;
	}
		break;
	case XMLSchemaSimpleSubelementList::ElementSummary::ct_group :
		{
		XMLSchemaModelGroupOP group( new XMLSchemaModelGroup );
		group->group_name( summary.ct_name );
		if ( summary.min_or_max_occurs_set ) {
			if ( summary.min_occurs != xsminmax_unspecified ) { group->min_occurs( summary.min_occurs ); }
			if ( summary.max_occurs != xsminmax_unspecified ) { group->max_occurs( summary.max_occurs ); }
		}
		return group;
	}
		break;
	}
	return XMLSchemaModelGroupOP(); // appease compiler
}

XMLSchemaComplexTypeGenerator::XMLSchemaComplexTypeGenerator() : pimpl_( new XMLSchemaComplexTypeGeneratorImpl ) {}
XMLSchemaComplexTypeGenerator::~XMLSchemaComplexTypeGenerator() { delete pimpl_; }
XMLSchemaComplexTypeGenerator::XMLSchemaComplexTypeGenerator( XMLSchemaComplexTypeGenerator const & src ) :
	pimpl_( new XMLSchemaComplexTypeGeneratorImpl( *src.pimpl_ ))
{
}

XMLSchemaComplexTypeGenerator & XMLSchemaComplexTypeGenerator::operator= ( XMLSchemaComplexTypeGenerator const & rhs )
{
	if ( this != & rhs ) {
		*pimpl_ = *rhs.pimpl_;
	}
	return *this;
}

XMLSchemaComplexTypeGenerator & XMLSchemaComplexTypeGenerator::element_name( std::string const & name )
{
	pimpl_->element_name( name );
	return *this;
}

XMLSchemaComplexTypeGenerator & XMLSchemaComplexTypeGenerator::complex_type_naming_func(
	DerivedNameFunction const & naming_function
)
{
	pimpl_->complex_type_naming_func( naming_function );
	return *this;
}

XMLSchemaComplexTypeGenerator & XMLSchemaComplexTypeGenerator::add_attribute( XMLSchemaAttribute const & attribute )
{
	pimpl_->add_attribute( attribute );
	return *this;
}

XMLSchemaComplexTypeGenerator & XMLSchemaComplexTypeGenerator::add_required_name_attribute() {
	pimpl_->add_attribute( required_name_attribute() );
	return *this;
}

XMLSchemaComplexTypeGenerator & XMLSchemaComplexTypeGenerator::add_optional_name_attribute() {
	pimpl_->add_attribute( optional_name_attribute() );
	return *this;
}

XMLSchemaComplexTypeGenerator & XMLSchemaComplexTypeGenerator::add_attributes( AttributeList const & attributes )
{
	pimpl_->add_attributes( attributes );
	return *this;
}

XMLSchemaComplexTypeGenerator &
XMLSchemaComplexTypeGenerator::set_subelements_repeatable(
	XMLSchemaSimpleSubelementList const & subelements,
	int min_occurs,
	int max_occurs
)
{
	pimpl_->set_subelements_repeatable( subelements, min_occurs, max_occurs );
	return *this;
}

XMLSchemaComplexTypeGenerator &
XMLSchemaComplexTypeGenerator::set_subelements_pick_one( XMLSchemaSimpleSubelementList const & subelements )
{
	pimpl_->set_subelements_pick_one_required( subelements );
	return *this;
}

XMLSchemaComplexTypeGenerator &
XMLSchemaComplexTypeGenerator::set_subelements_pick_one_or_none( XMLSchemaSimpleSubelementList const & subelements )
{
	pimpl_->set_subelements_pick_one_optional( subelements );
	return *this;
}


/// @details Why can you not have two or more "group" subelements that are required and whose order is unspecified?
/// Because the folks writing the XML Schema specification decided it would be too difficult for parsers to
/// handle such cases.  See the section on "Limitations of unordered content models" on this page:
/// http://docstore.mik.ua/orelly/xml/schema/ch07_04.htm
/// There are two solutions:
/// 1) drop the flexibility that comes with allowing elements to appear in any order and instead specify a rigid order
/// 2) change your format so that the group elements are beneath regular elements, e.g.
///
/// <OperateOnResidueSubset>
///  < *some residue selector* />
///  < *some residue-level-task-operation* />
/// </OperateOnResidueSubset>
///
/// would instead become
///
/// <OperateOnResidueSubset>
///  <Selector>
///   < *some residue selector* />
///  </Selector>
///  <RLTO>
///   < *some residue-level-task-operation* />
///  </RLTO>
/// </OperateOnResidueSubset>
///
/// which is not as pretty, but would allow order flexibility.
///
/// The reason you can get away with having exactly one group subelement in a "single appearance required"
/// call, is because instead of generating an <xs:all> block, the code instead generates an <xs:choice> block.
/// Of course order does not matter when there is only a single subelement.
XMLSchemaComplexTypeGenerator &
XMLSchemaComplexTypeGenerator::set_subelements_single_appearance_required( XMLSchemaSimpleSubelementList const & subelements )
{
	pimpl_->set_subelements_single_appearance_required( subelements );
	return *this;
}

XMLSchemaComplexTypeGenerator &
XMLSchemaComplexTypeGenerator::set_subelements_single_appearance_optional( XMLSchemaSimpleSubelementList const & subelements )
{
	pimpl_->set_subelements_single_appearance_optional( subelements );
	return *this;
}

XMLSchemaComplexTypeGenerator &
XMLSchemaComplexTypeGenerator::set_subelements_single_appearance_required_and_ordered(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	pimpl_->set_subelements_single_appearance_required_and_ordered( subelements );
	return *this;
}


XMLSchemaComplexTypeGenerator &
XMLSchemaComplexTypeGenerator::add_ordered_subelement_set_as_repeatable(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	pimpl_->add_ordered_subelement_set_as_repeatable( subelements );
	return *this;
}

XMLSchemaComplexTypeGenerator &
XMLSchemaComplexTypeGenerator::add_ordered_subelement_set_as_optional(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	pimpl_->add_ordered_subelement_set_as_optional( subelements );
	return *this;
}


XMLSchemaComplexTypeGenerator &
XMLSchemaComplexTypeGenerator::add_ordered_subelement_set_as_required(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	pimpl_->add_ordered_subelement_set_as_required( subelements );
	return *this;
}


XMLSchemaComplexTypeGenerator &
XMLSchemaComplexTypeGenerator::add_ordered_subelement_set_as_pick_one_or_none(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	pimpl_->add_ordered_subelement_set_as_pick_one_optional( subelements );
	return *this;
}


XMLSchemaComplexTypeGenerator &
XMLSchemaComplexTypeGenerator::add_ordered_subelement_set_as_pick_one(
	XMLSchemaSimpleSubelementList const & subelements
)
{
	pimpl_->add_ordered_subelement_set_as_pick_one_required( subelements );
	return *this;
}


void
XMLSchemaComplexTypeGenerator::write_complex_type_to_schema( XMLSchemaDefinition & xsd )
{
	pimpl_->write_complex_type_to_schema( xsd );
}

CTGenSubelementBehavior XMLSchemaComplexTypeGenerator::subelement_behavior() const
{
	return pimpl_->subelement_behavior();
}


XMLSchemaRestriction
integer_range_restriction( std::string const & name, int lower, int upper )
{
	XMLSchemaRestriction restriction;
	restriction.name( name );
	restriction.base_type( xs_integer );
	restriction.add_restriction( xsr_minInclusive, utility::to_string( lower ) );
	restriction.add_restriction( xsr_maxInclusive, utility::to_string( upper ) );
	return restriction;
}


XMLSchemaAttribute
required_name_attribute()
{
	return XMLSchemaAttribute::required_attribute( "name", xs_string );
}

XMLSchemaAttribute
optional_name_attribute()
{
	return XMLSchemaAttribute( "name", xs_string );
}


/// @details Most XML tags have a "name" attribute; this function does not require that
/// the name be provided.
void
append_name_and_attributes_to_complex_type(
	AttributeList const & attributes,
	XMLSchemaComplexType & type_definition
)
{
	type_definition.add_attribute( optional_name_attribute() );
	type_definition.add_attributes( attributes );
}

/// @details Most XML tags in Rosetta have a "name" attribute; this function appends an attribute "name" and
/// states that the attribute is required
void
append_required_name_and_attributes_to_complex_type(
	AttributeList const & attributes,
	XMLSchemaComplexType & type_definition
)
{
	type_definition.add_attribute( required_name_attribute() );
	type_definition.add_attributes( attributes );
}

}
}
