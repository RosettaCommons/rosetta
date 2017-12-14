// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueSelector.cc
/// @brief  The ResidueSelector class identifies a subset of residues from a Pose
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/select/residue_selector/ResidueSelector.hh>

// Utility headers
#include <utility/tag/Tag.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

ResidueSelector::ResidueSelector() = default;

ResidueSelector::~ResidueSelector() = default;

/// @details Noop implementation in the base class in the case that a derived
/// class has no need to read data from an input tag
void ResidueSelector::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap &
)
{}

} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::ResidueSelector::save( Archive & ) const {}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::ResidueSelector::load( Archive & ) {}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::ResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::ResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_ResidueSelector )
#endif // SERIALIZATION
