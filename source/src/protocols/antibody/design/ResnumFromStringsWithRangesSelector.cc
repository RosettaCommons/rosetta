// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/ResnumFromStringsWithRangesSelector.hh
/// @brief Select residues with the protocols::antibody::design::get_resnums_from_strings_with_ranges() approach.
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <protocols/antibody/design/ResnumFromStringsWithRangesSelector.hh>
#include <protocols/antibody/design/util.hh>

#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.antibody.design.ResnumFromStringsWithRangesSelector" );

namespace protocols {
namespace antibody {
namespace design {

/// @brief Constructor.
///
ResnumFromStringsWithRangesSelector::ResnumFromStringsWithRangesSelector() = default;

ResnumFromStringsWithRangesSelector::ResnumFromStringsWithRangesSelector( utility::vector1<std::string> const & pdb_residues ):
	ResidueSelector(),
	pdb_residues_(pdb_residues)
{}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
ResnumFromStringsWithRangesSelector::ResidueSelectorOP
ResnumFromStringsWithRangesSelector::clone() const {
	return utility::pointer::make_shared< ResnumFromStringsWithRangesSelector >(*this);
}

std::string ResnumFromStringsWithRangesSelector::get_name() const
{
	return ResnumFromStringsWithRangesSelector::class_name();
}

std::string ResnumFromStringsWithRangesSelector::class_name()
{
	return "ResnumFromStringsWithRangesSelector";
}


/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
ResnumFromStringsWithRangesSelector::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap & )
{
	utility_exit_with_message("ResnumFromStringsWithRangesSelector not currently set up for XML usage.");
}

void ResnumFromStringsWithRangesSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & )
{
	//Syntax Example:
	//using namespace utility::tag;
	//AttributeList attributes;
	//attributes
	// + XMLSchemaAttribute::attribute_w_default(  "select_positive_phi",      xsct_rosetta_bool, "true" )
	// + XMLSchemaAttribute::attribute_w_default(  "ignore_unconnected_upper", xsct_rosetta_bool, "true" );
	//xsd_type_definition_w_attributes( xsd, class_name(), attributes );

}


/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ResnumFromStringsWithRangesSelector::ResidueSubset
ResnumFromStringsWithRangesSelector::apply(
	core::pose::Pose const & pose
) const {

	return get_resnums_from_strings_with_ranges( pose, pdb_residues_ );
}


} //design
} //antibody
} //protocols

#ifdef    SERIALIZATION

/// @brief serialization method
template< class Archive >
void
protocols::antibody::design::ResnumFromStringsWithRangesSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( pdb_residues_ ) ); // utility::vector1<std::string
}

/// @brief deserialization method
template< class Archive >
void
protocols::antibody::design::ResnumFromStringsWithRangesSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( pdb_residues_ ); // utility::vector1<std::string
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::antibody::design::ResnumFromStringsWithRangesSelector );
CEREAL_REGISTER_TYPE( protocols::antibody::design::ResnumFromStringsWithRangesSelector )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_antibody_design_ResnumFromStringsWithRangesSelector )
#endif // SERIALIZATION


