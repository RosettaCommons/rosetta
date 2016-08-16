// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>

#include <core/pack/rotamer_set/WaterAnchorInfo.hh>

// utility headers
#include <string>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {

/// @details Auto-generated virtual destructor
WaterAnchorInfo::~WaterAnchorInfo() {}

Size
WaterAnchorInfo::anchor_residue() const {
	return anchor_residue_;
}

void
WaterAnchorInfo::anchor_residue( Size const rsd ) {
	anchor_residue_ = rsd;
}

bool
WaterAnchorInfo::attaches_to_residue_type( ResidueType
	const & type ) const {
	return type.aa() == aa_ && type.has( anchor_atom_name_ );
}

Size
WaterAnchorInfo::anchor_atom( ResidueType const & type
) const {
	debug_assert( attaches_to_residue_type( type ) );
	return type.atom_index( anchor_atom_name_ );
}

void
WaterAnchorInfo::anchor_atom( std::string const & name
) {
	anchor_atom_name_ = name;
}

void
WaterAnchorInfo::aa( AA const & aa_in ) {
	aa_ = aa_in;
}

Size
WaterAnchorInfo::nstep() const {
	return nstep_;
}

void
WaterAnchorInfo::nstep( Size const nstep_in ) {
	nstep_ = nstep_in;
}

} // namespace rotamer_set
} // namespace pack
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::rotamer_set::WaterAnchorInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( anchor_residue_ ) ); // Size
	arc( CEREAL_NVP( anchor_atom_name_ ) ); // std::string
	arc( CEREAL_NVP( aa_ ) ); // AA
	arc( CEREAL_NVP( nstep_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::rotamer_set::WaterAnchorInfo::load( Archive & arc ) {
	arc( anchor_residue_ ); // Size
	arc( anchor_atom_name_ ); // std::string
	arc( aa_ ); // AA
	arc( nstep_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::rotamer_set::WaterAnchorInfo );
CEREAL_REGISTER_TYPE( core::pack::rotamer_set::WaterAnchorInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_rotamer_set_WaterAnchorInfo )
#endif // SERIALIZATION
