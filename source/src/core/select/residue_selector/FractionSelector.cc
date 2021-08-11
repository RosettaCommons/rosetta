// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/FractionSelector.hh
/// @brief  Selects only either the first or the last residue selected by a provided ResidueSelector
/// @author frankdt (frankdt@email.unc.edu)

// Unit headers
#include <core/select/residue_selector/FractionSelector.hh>
#include <core/select/residue_selector/FractionSelectorCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions

// Project headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "protocols.pose_sewing.residue_selectors.FractionSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
FractionSelector::FractionSelector():
	core::select::residue_selector::ResidueSelector()
{
}

/// @brief Destructor.
///
FractionSelector::~FractionSelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
//FractionSelector::FractionSelector(FractionSelector const & src):
// core::select::residue_selector::ResidueSelector( src )
//{
//}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
FractionSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP( new FractionSelector(*this) );
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
FractionSelector::ResidueSubset
FractionSelector::apply(
	core::pose::Pose const & pose//pose
) const {
	// TODO: remove the line below and use this function to create a residue subset
	//       based on the input pose
	core::select::residue_selector::ResidueSubset selection( pose.total_residue(), false );
	if ( selector_ != nullptr ) {
		selection = selector_->apply( pose );
	}
	ResidueSubset matching_set( pose.size(), false );
	core::Size remaining = N_;
	core::Size selected = 0;
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( selection[i] ) {
			++selected;
		}
	}
	core::Size remaining_fraction = core::Size(fraction_ * selected);
	if ( remaining_fraction < remaining ) {
		remaining = remaining_fraction;
	}
	if ( first_ ) {
		for ( core::Size i=1; i<=pose.size(); ++i ) {
			if ( selection[i] && remaining > 0 ) {
				matching_set[ i ] = true;
				--remaining;
			}
		}
	} else {
		for ( core::Size i=pose.size(); i>=1; --i ) {
			if ( selection[i] && remaining > 0 ) {
				matching_set[ i ] = true;
				--remaining;
			}
		}
	}
	return matching_set;
}
void
FractionSelector::set_first(bool first){
	first_ = first;
}

void
FractionSelector::set_N(core::Size N){
	N_ = N;
}

void
FractionSelector::set_fraction(core::Real fraction){
	fraction_ = fraction;
}
void
FractionSelector::set_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
FractionSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{
	if ( tag->hasOption("first") ) {
		first_ = tag->getOption<bool>("first");
	}
	if ( tag->hasOption("N") ) {
		N_ = tag->getOption<core::Size>("N");
	}
	if ( tag->hasOption("fraction") ) {
		fraction_ = tag->getOption<core::Real>("fraction");
	}
	if ( tag->hasOption( "selector" ) ) {
		selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "selector" ), datamap );
	}
}

std::string FractionSelector::get_name() const
{
	return FractionSelector::class_name();
}

std::string FractionSelector::class_name()
{
	return "FractionSelector";
}

void FractionSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "first", xsct_rosetta_bool, "Selects residues starting from the N terminus if true, and the C terminus if false/", "true" )
		+ XMLSchemaAttribute::attribute_w_default(  "N", xsct_positive_integer, "Hard cap on the number of residues selected", "1" )
		+ XMLSchemaAttribute::attribute_w_default(  "fraction", xsct_real, "Decimal fraction of the residues selected by the provided selector that are to be selected by this selector.", "1" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attributes, "selector",  "residue selector that defines the residues from which the specified fraction is to be selected" );
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "Selects a specified fraction of a provided selector, either in absolute residue count or in terms of a total of the entire selection, ", attributes );

}

core::select::residue_selector::ResidueSelectorOP
FractionSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP( new FractionSelector );
}

std::string
FractionSelectorCreator::keyname() const {
	return FractionSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
FractionSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	FractionSelector::provide_xml_schema( xsd );
}

} //protocols
} //pose_sewing
} //residue_selectors

