// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief A Chi record object for a MutableResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <core/chemical/MutableChiRecord.hh>
#include <core/chemical/MutableResidueType.hh>

// Project headers

// Numeric headers

// Utility headers

//#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>

// C++ headers

#ifdef    SERIALIZATION
#include <core/chemical/ResidueGraphTypes.srlz.hh>

// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp> // For std::pair
#endif // SERIALIZATION

namespace core {
namespace chemical {

static basic::Tracer TR( "core.chemical.MutableChiRecord" );

//////////////////////////////////////////////////////////////////////////////////////////////////

MutableChiRecord::MutableChiRecord(VD atm1, VD atm2, VD atm3, VD atm4):
	chi_atoms_( { atm1, atm2, atm3, atm4 } )
{}

MutableChiRecord::MutableChiRecord(utility::vector1< VD > const & atm_vec):
	chi_atoms_( atm_vec )
{
	debug_assert( chi_atoms_.size() == 4 );
}

void
MutableChiRecord::set_proton_chi(
	utility::vector1< Real > const & dihedral_samples,
	utility::vector1< Real > const & extra_samples
) {
	is_proton_chi_ = true;
	proton_chi_samples_ = dihedral_samples;
	proton_chi_extra_samples_ = extra_samples;
}

void
MutableChiRecord::add_chi_rotamer( Real const mean, Real const sdev ) {
	chi_rotamers_.push_back( std::make_pair( mean, sdev ) );
}

void
MutableChiRecord::clear_chi_rotamers() {
	chi_rotamers_.clear();
}

void
MutableChiRecord::remap_atom_vds( std::map< VD, VD > const & old_to_new ) {
	for ( core::Size ii(1); ii <= chi_atoms_.size(); ++ii ) {
		if ( old_to_new.count( chi_atoms_[ ii ] ) ) {
			chi_atoms_[ ii ] = old_to_new.at( chi_atoms_[ ii ] );
		}
	}
}

} // chemical
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::MutableChiRecord::save( Archive & arc ) const {
	SERIALIZE_VD_VECTOR( arc, chi_atoms_ ); // EXEMPT chi_atoms_
	arc( CEREAL_NVP( is_proton_chi_ ) );
	arc( CEREAL_NVP( proton_chi_samples_ ) );
	arc( CEREAL_NVP( proton_chi_extra_samples_ ) );
	arc( CEREAL_NVP( chi_rotamers_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::MutableChiRecord::load( Archive & arc ) {
	DESERIALIZE_VD_VECTOR( arc, chi_atoms_ ); // EXEMPT chi_atoms_
	arc( is_proton_chi_ );
	arc( proton_chi_samples_ );
	arc( proton_chi_extra_samples_ );
	arc( chi_rotamers_ );

	// IMPORTANT - You probably need to call remap_atom_vds() on this class afterward
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::MutableChiRecord );
#endif // SERIALIZATION
