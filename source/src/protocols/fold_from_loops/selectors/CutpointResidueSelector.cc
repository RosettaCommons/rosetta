// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fold_from_loops/CutpointResidueSelector.cc
/// @brief  Selects those residues that are amino acids
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

// Unit headers
#include <protocols/fold_from_loops/selectors/CutpointResidueSelector.hh>
#include <protocols/fold_from_loops/selectors/CutpointResidueSelectorCreator.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace fold_from_loops {
namespace selectors {


CutpointResidueSelector::CutpointResidueSelector() :
	ResidueSelector(),
	use_foldtree_( true )
{}

CutpointResidueSelector::CutpointResidueSelector( bool use_foldtree ) :
	ResidueSelector(),
	use_foldtree_( use_foldtree )
{}

CutpointResidueSelector::~CutpointResidueSelector() {}

/// @brief Copy constructor
///
CutpointResidueSelector::CutpointResidueSelector( CutpointResidueSelector const &) :
	ResidueSelector()
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
core::select::residue_selector::ResidueSelectorOP CutpointResidueSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP( new CutpointResidueSelector(*this) );
}

core::select::residue_selector::ResidueSubset
CutpointResidueSelector::apply( core::pose::Pose const & pose ) const
{
	core::select::residue_selector::ResidueSubset aminoacids;
	for ( core::Size i = 1; i <= pose.size(); i++ ) {
		if ( use_foldtree_ ) {
			aminoacids.push_back(
				pose.fold_tree().is_cutpoint(i) &&
				!pose.residue(i).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT )
			);
		} else {
			aminoacids.push_back(
				pose.residue(i).has_variant_type( core::chemical::CUTPOINT_LOWER ) &&
				!pose.residue(i).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) &&
				!pose.residue(i).has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ));
		}
	}
	return aminoacids;
}

void CutpointResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap&
)
{
	use_foldtree( tag->getOption< bool >( "use_foldtree", true ) );
}

std::string CutpointResidueSelector::get_name() const {
	return CutpointResidueSelector::class_name();
}

std::string CutpointResidueSelector::class_name() {
	return "CutpointResidueSelector";
}

void
CutpointResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute::attribute_w_default( "use_foldtree", xsct_rosetta_bool,
		"Selector requires the residue to be cutpoint in the FoldTree and as variant type. Otherwise, just based on variant type " , std::to_string( true ) );
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "Selects all amino acid residues that are cutpoints of the FoldTree", attributes );
}

core::select::residue_selector::ResidueSelectorOP
CutpointResidueSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP( new CutpointResidueSelector );
}

std::string
CutpointResidueSelectorCreator::keyname() const {
	return CutpointResidueSelector::class_name();
}

void
CutpointResidueSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	CutpointResidueSelector::provide_xml_schema( xsd );
}

} //namespace selectors
} //namespace fold_from_loops
} //namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::fold_from_loops::selectors::CutpointResidueSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( use_foldtree_ ) ); // bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::fold_from_loops::selectors::CutpointResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( use_foldtree_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::fold_from_loops::selectors::CutpointResidueSelector );
CEREAL_REGISTER_TYPE( protocols::fold_from_loops::selectors::CutpointResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_fold_from_loops_CutpointResidueSelector )
#endif // SERIALIZATION
