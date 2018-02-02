// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fold_from_loops/ProteinResidueSelector.cc
/// @brief  Selects those residues that are amino acids
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

// Unit headers
#include <protocols/fold_from_loops/selectors/ProteinResidueSelector.hh>
#include <protocols/fold_from_loops/selectors/ProteinResidueSelectorCreator.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
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


ProteinResidueSelector::ProteinResidueSelector() :
	ResidueSelector()
{}

ProteinResidueSelector::~ProteinResidueSelector() {}

/// @brief Copy constructor
///
ProteinResidueSelector::ProteinResidueSelector( ProteinResidueSelector const &) :
	ResidueSelector()
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
core::select::residue_selector::ResidueSelectorOP ProteinResidueSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP( new ProteinResidueSelector(*this) );
}

core::select::residue_selector::ResidueSubset
ProteinResidueSelector::apply( core::pose::Pose const & pose ) const
{
	core::select::residue_selector::ResidueSubset aminoacids;
	for ( core::Size i = 1; i <= pose.size(); i++ ) {
		aminoacids.push_back( pose.residue(i).is_protein() );
	}
	return aminoacids;
}

void ProteinResidueSelector::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap&
)
{}

std::string ProteinResidueSelector::get_name() const {
	return ProteinResidueSelector::class_name();
}

std::string ProteinResidueSelector::class_name() {
	return "ProteinResidueSelector";
}

void
ProteinResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "Selects all amino acid residues from the Pose", attributes );
}

core::select::residue_selector::ResidueSelectorOP
ProteinResidueSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP( new ProteinResidueSelector );
}

std::string
ProteinResidueSelectorCreator::keyname() const {
	return ProteinResidueSelector::class_name();
}

void
ProteinResidueSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ProteinResidueSelector::provide_xml_schema( xsd );
}

} //namespace selectors
} //namespace fold_from_loops
} //namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::fold_from_loops::selectors::ProteinResidueSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::fold_from_loops::selectors::ProteinResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::fold_from_loops::selectors::ProteinResidueSelector );
CEREAL_REGISTER_TYPE( protocols::fold_from_loops::selectors::ProteinResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_fold_from_loops_ProteinResidueSelector )
#endif // SERIALIZATION
