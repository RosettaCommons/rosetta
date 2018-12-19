// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/FalseResidueSelector.cc
/// @brief  The FalseResidueSelector creates an appropriate all-false vector
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <core/select/residue_selector/FalseResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/assert.hh>

// C++ headers

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {


FalseResidueSelector::FalseResidueSelector() :
	ResidueSelector()
{}

FalseResidueSelector::~FalseResidueSelector() = default;

/// @brief Copy constructor
///
FalseResidueSelector::FalseResidueSelector( FalseResidueSelector const &) :
	ResidueSelector()
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP FalseResidueSelector::clone() const { return ResidueSelectorOP( new FalseResidueSelector(*this) ); }

ResidueSubset
FalseResidueSelector::apply( core::pose::Pose const & pose ) const
{
	return ResidueSubset( pose.size(), false );
}

void FalseResidueSelector::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap&
)
{}

std::string FalseResidueSelector::get_name() const {
	return FalseResidueSelector::class_name();
}

std::string FalseResidueSelector::class_name() {
	return "False";
}

void
FalseResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	xsd_type_definition_w_attributes( xsd, class_name(), "A residue selector which doesn't select any residues.", attributes );
}

ResidueSelectorOP
FalseResidueSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new FalseResidueSelector );
}

std::string
FalseResidueSelectorCreator::keyname() const {
	return FalseResidueSelector::class_name();
}

void
FalseResidueSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	FalseResidueSelector::provide_xml_schema( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::FalseResidueSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::FalseResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::FalseResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::FalseResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_FalseResidueSelector )
#endif // SERIALIZATION
