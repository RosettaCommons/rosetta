// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/AsymmetricUnitSelector.hh
/// @brief  A residue selector for selecting the master subunit in a symmetrical pose.  If not symmetrical, will select the whole pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/select/residue_selector/AsymmetricUnitSelector.hh>
#include <core/select/residue_selector/AsymmetricUnitSelectorCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

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


static basic::Tracer TR( "core.select.residue_selector.AsymmetricUnitSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
AsymmetricUnitSelector::AsymmetricUnitSelector():
	core::select::residue_selector::ResidueSelector()
{
}

/// @brief Destructor.
///
AsymmetricUnitSelector::~AsymmetricUnitSelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
//AsymmetricUnitSelector::AsymmetricUnitSelector(AsymmetricUnitSelector const & src):
// core::select::residue_selector::ResidueSelector( src )
//{
//}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
AsymmetricUnitSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP( new AsymmetricUnitSelector(*this) );
}


/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
AsymmetricUnitSelector::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap & )
{
}

std::string AsymmetricUnitSelector::get_name() const
{
	return AsymmetricUnitSelector::class_name();
}

std::string AsymmetricUnitSelector::class_name()
{
	return "AsymmetricUnitSelector";
}

void AsymmetricUnitSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Syntax Example:
	using namespace utility::tag;
	AttributeList attributes;
	//attributes
	// + XMLSchemaAttribute::attribute_w_default(  "select_positive_phi",      xsct_rosetta_bool, "Description of first option here!", "true" )
	// + XMLSchemaAttribute::attribute_w_default(  "ignore_unconnected_upper", xsct_rosetta_bool, "Description of second option here!", "true" );
	//core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "Description of residue selector here!", attributes );

	xsd_type_definition_w_attributes( xsd, class_name(),
		"A ResidueSelector that selects only the master subunit of a symmetric pose or all of the residues in a non-symmetric pose", attributes );

}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
AsymmetricUnitSelector::ResidueSubset
AsymmetricUnitSelector::apply(
	core::pose::Pose const & pose
) const {
	using namespace core::conformation::symmetry;

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		utility::vector1< bool > subset(pose.size(), false);

		auto const & symm_conf (
			dynamic_cast<SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( symm_info->bb_is_independent(i) ) {
				subset[i] = true;
			}
		}
		return subset;
	}

	utility::vector1< bool > subset( pose.size(), true);
	return subset;

}

core::select::residue_selector::ResidueSelectorOP
AsymmetricUnitSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP( new AsymmetricUnitSelector );
}

std::string
AsymmetricUnitSelectorCreator::keyname() const {
	return AsymmetricUnitSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
AsymmetricUnitSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	AsymmetricUnitSelector::provide_xml_schema( xsd );
}


} //core
} //select
} //residue_selector

#ifdef    SERIALIZATION

// See the serialization documentation here.  There is a script you can run.
//  https://wiki.rosettacommons.org/index.php/SerializationFAQ

template< class Archive >
void
core::select::residue_selector::AsymmetricUnitSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	// You need to add "arc( CEREAL_NVP( datamember_name ) );" calls here for each of your data members.
}

template< class Archive >
void
core::select::residue_selector::AsymmetricUnitSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	// You need to add "arc( datamember_name );" calls here for each of your data members.
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::AsymmetricUnitSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::AsymmetricUnitSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_AsymmetricUnitSelector )
#endif // SERIALIZATION
