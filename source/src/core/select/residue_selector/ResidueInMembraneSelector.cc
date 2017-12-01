// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/ResidueInMembraneSelector.hh
/// @brief  The ResidueInMembraneSelector selects residues either in or out of the membrane
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

// Unit headers
#include <core/select/residue_selector/ResidueInMembraneSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// package headers
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <utility/assert.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

ResidueInMembraneSelector::ResidueInMembraneSelector():
	select_in_membrane_(true)
{
}

/// @brief Copy constructor
///
ResidueInMembraneSelector::ResidueInMembraneSelector( ResidueInMembraneSelector const &src) :
	select_in_membrane_( src.select_in_membrane_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP ResidueInMembraneSelector::clone() const { return ResidueSelectorOP( new ResidueInMembraneSelector(*this) ); }


ResidueInMembraneSelector::~ResidueInMembraneSelector() {}

ResidueSubset
ResidueInMembraneSelector::apply( core::pose::Pose const & pose ) const
{
	ResidueSubset subset( pose.total_residue(), false );
	// quit here if there are no residues in the pose
	if ( pose.total_residue() == 0 ) {
		return subset;
	}

	core::Real thickness = pose.conformation().membrane_info()->membrane_thickness();
	for ( core::Size res_i=1; res_i <= pose.total_residue(); res_i++ ) {
		core::conformation::Residue const & rsd = pose.residue(res_i);
		if ( ! rsd.is_protein() ) continue;

		Real const z_position( pose.conformation().membrane_info()->residue_z_position( pose.conformation(), rsd.seqpos() ) );


		bool in_membrane;
		if ( -thickness <= z_position && z_position <= thickness ) {
			in_membrane = true;
		} else {
			in_membrane = false;
		}

		if ( select_in_membrane_ ) {
			subset[res_i] = in_membrane;
		} else {
			subset[res_i] = !in_membrane;
		}

	}
	return subset;
}

void
ResidueInMembraneSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &)
{
	if ( tag->hasOption( "select_in_membrane" ) ) {
		set_select_in_membrane( tag->getOption< bool >( "select_in_membrane" ) );
	} else {
		set_select_in_membrane( true );
	}
}


void
ResidueInMembraneSelector::set_select_in_membrane( bool const & select_in_membrane ) {
	select_in_membrane_ = select_in_membrane;
}

std::string ResidueInMembraneSelector::get_name() const {
	return ResidueInMembraneSelector::class_name();
}

std::string ResidueInMembraneSelector::class_name() {
	return "ResidueInMembrane";
}

void
ResidueInMembraneSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "select_in_membrane", xsct_rosetta_bool, "whether choose resiudes in membrane or out" );
	xsd_type_definition_w_attributes( xsd, class_name(), "choose membrane or extramembrane residues", attributes );
}

ResidueSelectorOP
ResidueInMembraneSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new ResidueInMembraneSelector );
}

std::string
ResidueInMembraneSelectorCreator::keyname() const {
	return ResidueInMembraneSelector::class_name();
}

void ResidueInMembraneSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueInMembraneSelector::provide_xml_schema( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::ResidueInMembraneSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( select_in_membrane_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::ResidueInMembraneSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( select_in_membrane_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::ResidueInMembraneSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::ResidueInMembraneSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_ResidueInMembraneSelector )
#endif // SERIALIZATION
