// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/VirtualResidueSelector.hh
/// @brief  a residue selector that wraps is_virtual_residue
/// @author frankdt (frankdt@email.unc.edu)

// Unit headers
#include <core/select/residue_selector/VirtualResidueSelector.hh>
#include <core/select/residue_selector/VirtualResidueSelectorCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

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


static basic::Tracer TR( "protocols.pose_sewing.residue_selectors.VirtualResidueSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
VirtualResidueSelector::VirtualResidueSelector():
	core::select::residue_selector::ResidueSelector()
{
}

/// @brief Destructor.
///
VirtualResidueSelector::~VirtualResidueSelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
//VirtualResidueSelector::VirtualResidueSelector(VirtualResidueSelector const & src):
// core::select::residue_selector::ResidueSelector( src )
//{
//}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
VirtualResidueSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP( new VirtualResidueSelector(*this) );
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
VirtualResidueSelector::ResidueSubset
VirtualResidueSelector::apply(core::pose::Pose const & pose) const {

	ResidueSubset matching_set( pose.size(), false );
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		matching_set[ i ] = pose.residue(i).is_virtual_residue();
	}
	return matching_set;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
VirtualResidueSelector::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap & )
{
}

std::string VirtualResidueSelector::get_name() const
{
	return VirtualResidueSelector::class_name();
}

std::string VirtualResidueSelector::class_name()
{
	return "VirtualResidueSelector";
}

void VirtualResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Syntax Example:
	using namespace utility::tag;
	AttributeList attributes;
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "Selects virtual residues.", attributes );

}

core::select::residue_selector::ResidueSelectorOP
VirtualResidueSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP( new VirtualResidueSelector );
}

std::string
VirtualResidueSelectorCreator::keyname() const {
	return VirtualResidueSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
VirtualResidueSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	VirtualResidueSelector::provide_xml_schema( xsd );
}

} //protocols
} //pose_sewing
} //residue_selectors

