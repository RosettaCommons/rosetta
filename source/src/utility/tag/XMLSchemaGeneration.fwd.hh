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

#ifndef INCLUDED_utility_tag_XMLSchemaGeneration_FWD_HH
#define INCLUDED_utility_tag_XMLSchemaGeneration_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <list>

namespace utility {
namespace tag {

class XMLSchemaType;
typedef utility::pointer::shared_ptr< XMLSchemaType > XMLSchemaTypeOP;
typedef utility::pointer::shared_ptr< XMLSchemaType const > XMLSchemaTypeCOP;

class XMLSchemaAttribute;
typedef utility::pointer::shared_ptr< XMLSchemaAttribute > XMLSchemaAttributeOP;
typedef utility::pointer::shared_ptr< XMLSchemaAttribute const > XMLSchemaAttributeCOP;

typedef std::list< utility::tag::XMLSchemaAttribute > AttributeList;

class XMLSchemaElement;
typedef utility::pointer::shared_ptr< XMLSchemaElement > XMLSchemaElementOP;
typedef utility::pointer::shared_ptr< XMLSchemaElement const > XMLSchemaElementCOP;

class XMLSchemaRestriction;
typedef utility::pointer::shared_ptr< XMLSchemaRestriction > XMLSchemaRestrictionOP;
typedef utility::pointer::shared_ptr< XMLSchemaRestriction const > XMLSchemaRestrictionCOP;

class XMLSchemaParticle;
typedef utility::pointer::shared_ptr< XMLSchemaParticle > XMLSchemaParticleOP;
typedef utility::pointer::shared_ptr< XMLSchemaParticle const > XMLSchemaParticleCOP;

class XMLSchemaModelGroup;
typedef utility::pointer::shared_ptr< XMLSchemaModelGroup > XMLSchemaModelGroupOP;
typedef utility::pointer::shared_ptr< XMLSchemaModelGroup const > XMLSchemaModelGroupCOP;

class XMLSchemaComplexType;
typedef utility::pointer::shared_ptr< XMLSchemaComplexType > XMLSchemaComplexTypeOP;
typedef utility::pointer::shared_ptr< XMLSchemaComplexType const > XMLSchemaComplexTypeCOP;

class XMLSchemaElementType;
typedef utility::pointer::shared_ptr< XMLSchemaElementType > XMLSchemaElementTypeOP;
typedef utility::pointer::shared_ptr< XMLSchemaElementType const > XMLSchemaElementTypeCOP;

class XMLSchemaDefinition;
typedef utility::pointer::shared_ptr< XMLSchemaDefinition > XMLSchemaDefinitionOP;
typedef utility::pointer::shared_ptr< XMLSchemaDefinition const > XMLSchemaDefinitionCOP;

class XMLSchemaSimpleSubelementList;
typedef utility::pointer::shared_ptr< XMLSchemaSimpleSubelementList > XMLSchemaSimpleSubelementListOP;
typedef utility::pointer::shared_ptr< XMLSchemaSimpleSubelementList const > XMLSchemaSimpleSubelementListCOP;

class XMLComplexTypeSchemaGenerator;
typedef utility::pointer::shared_ptr< XMLComplexTypeSchemaGenerator > XMLComplexTypeSchemaGeneratorOP;
typedef utility::pointer::shared_ptr< XMLComplexTypeSchemaGenerator const > XMLComplexTypeSchemaGeneratorCOP;


enum XMLSchemaDataType {
	xs_string,
	xs_decimal,
	xs_integer,
	xs_boolean,
	xs_date,
	xs_time,
	xs_custom
};

enum XMLSchemaMinOccursMaxOccurs {
	xsminmax_unbounded = -2,
	xsminmax_unspecified = -1
};

enum XMLSchemaRestrictionType {
	xsr_enumeration,
	xsr_fractionDigits,
	xsr_length,
	xsr_maxExclusive,
	xsr_maxInclusive,
	xsr_maxLength,
	xsr_minExclusive,
	xsr_minInclusive,
	xsr_minLength,
	xsr_pattern,
	xsr_totalDigits,
	xsr_whitespace
};

enum XMLSchemaModelGroupType {
	xsmgt_sequence,
	xsmgt_all,
	xsmgt_choice,
	xsmgt_group
};

enum XMLSchemaElementCategory {
	xs_element_is_type_reference,
	xs_element_is_abstract,
	xs_element_is_complex_type_w_definition,
	xs_element_is_element_reference
};

}
}

#endif
