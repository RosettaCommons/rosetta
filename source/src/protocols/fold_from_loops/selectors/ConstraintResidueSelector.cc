// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fold_from_loops/ConstraintResidueSelector.cc
/// @brief  Selects those residues that have distance/angle/dihedral constraints
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

// Unit headers
#include <protocols/fold_from_loops/selectors/ConstraintResidueSelector.hh>
#include <protocols/fold_from_loops/selectors/ConstraintResidueSelectorCreator.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
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


ConstraintResidueSelector::ConstraintResidueSelector() :
	ResidueSelector(),
	distance_( default_distance() ),
	angle_( default_angle() ),
	dihedral_( default_dihedral() )
{}

ConstraintResidueSelector::~ConstraintResidueSelector() {}

/// @brief Copy constructor
///
ConstraintResidueSelector::ConstraintResidueSelector( ConstraintResidueSelector const &) :
	ResidueSelector()
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
core::select::residue_selector::ResidueSelectorOP ConstraintResidueSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP( new ConstraintResidueSelector(*this) );
}

core::select::residue_selector::ResidueSubset
ConstraintResidueSelector::apply( core::pose::Pose const & pose ) const
{
	core::select::residue_selector::ResidueSubset aminoacids( pose.size(), false);
	core::scoring::constraints::ConstraintSetCOP constraints = pose.constraint_set();
	if ( !constraints->is_empty() ) {
		for ( auto cst : constraints->get_all_constraints() ) {
			std::string csttype=cst->type();
			if ( distance_ ) {
				if ( csttype == "AtomPair" ) {
					aminoacids[cst->atom(1).rsd()] = true;
					aminoacids[cst->atom(2).rsd()] = true;
				}
			}
			if ( angle_ ) {
				if ( csttype == "Angle" ) {
					aminoacids[cst->atom(1).rsd()] = true;
					aminoacids[cst->atom(2).rsd()] = true;
					aminoacids[cst->atom(3).rsd()] = true;
				}
			}
			if ( dihedral_ ) {
				if ( csttype == "Dihedral" ) {
					aminoacids[cst->atom(1).rsd()] = true;
					aminoacids[cst->atom(2).rsd()] = true;
					aminoacids[cst->atom(3).rsd()] = true;
					aminoacids[cst->atom(4).rsd()] = true;
				}
			}
		}
	}
	return aminoacids;
}

void ConstraintResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap&
)
{
	distance( tag->getOption< bool >( "distance", default_distance() ) );
	angle( tag->getOption< bool >( "abgle", default_angle() ) );
	dihedral( tag->getOption< bool >( "dihedral", default_dihedral() ) );

}

std::string ConstraintResidueSelector::get_name() const {
	return ConstraintResidueSelector::class_name();
}

std::string ConstraintResidueSelector::class_name() {
	return "ConstraintResidueSelector";
}

void
ConstraintResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "distance", xsct_rosetta_bool, "Select distance (csttype=AtomPair) constraints.", std::to_string( default_distance() ) )
		+ XMLSchemaAttribute::attribute_w_default( "angle", xsct_rosetta_bool, "Select angle (csttype=Angle) constraints.", std::to_string( default_angle() ) )
		+ XMLSchemaAttribute::attribute_w_default( "dihedral", xsct_rosetta_bool, "Select dihedral (csttype=Dihedral) constraints.", std::to_string( default_dihedral() ) );
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "Selects those residues that have distance/angle/dihedral constraints", attributes );
}

core::select::residue_selector::ResidueSelectorOP
ConstraintResidueSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP( new ConstraintResidueSelector );
}

std::string
ConstraintResidueSelectorCreator::keyname() const {
	return ConstraintResidueSelector::class_name();
}

void
ConstraintResidueSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ConstraintResidueSelector::provide_xml_schema( xsd );
}

} //namespace selectors
} //namespace fold_from_loops
} //namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::fold_from_loops::selectors::ConstraintResidueSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( distance_ ) ); // bool
	arc( CEREAL_NVP( angle_ ) ); // bool
	arc( CEREAL_NVP( dihedral_ ) ); // bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::fold_from_loops::selectors::ConstraintResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( distance_ );
	arc( angle_ );
	arc( dihedral_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::fold_from_loops::selectors::ConstraintResidueSelector );
CEREAL_REGISTER_TYPE( protocols::fold_from_loops::selectors::ConstraintResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_fold_from_loops_ConstraintResidueSelector )
#endif // SERIALIZATION
