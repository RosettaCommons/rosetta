// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/CacheableResidueTypeSets.cc
/// @brief A (Pose-cacheable) container for ResidueTypeSets
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/CacheableResidueTypeSets.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/PoseResidueTypeSet.hh>

#include <basic/Tracer.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "core.chemical.CacheableResidueTypeSets" );


namespace core {
namespace chemical {

CacheableResidueTypeSets::CacheableResidueTypeSets():
	basic::datacache::CacheableData(),
	res_type_sets_(TYPE_SET_MODES_LENGTH, nullptr)
{}

CacheableResidueTypeSets::~CacheableResidueTypeSets()= default;

CacheableResidueTypeSets::CacheableResidueTypeSets( CacheableResidueTypeSets const & other ) :
	basic::datacache::CacheableData(*this),
	// Shallow copy of the PoseResidueTypeSets
	res_type_sets_( other.res_type_sets_ )
{}

// Default is a shallow copy of the data members.
CacheableResidueTypeSets & CacheableResidueTypeSets::operator= (CacheableResidueTypeSets const & ) = default;

basic::datacache::CacheableDataOP
CacheableResidueTypeSets::clone() const {
	return basic::datacache::CacheableDataOP( new CacheableResidueTypeSets( *this ) );
}

void
CacheableResidueTypeSets::clear() {
	res_type_sets_.clear();
	res_type_sets_.resize( TYPE_SET_MODES_LENGTH, nullptr );
}

/// @brief Do we have a 'mode' ResidueTypeSet already instantiated?
bool
CacheableResidueTypeSets::has_res_type_set( TypeSetMode mode ) const {
	debug_assert( core::Size(mode) >= 1 );
	debug_assert( mode <= TYPE_SET_MODES_LENGTH );
	return res_type_sets_[ mode ] != nullptr;
}

/// @brief Return a ResidueTypeSet of the appropriate type,
/// If one doesn't already exist, return a null pointer.
PoseResidueTypeSetCOP
CacheableResidueTypeSets::get_res_type_set( TypeSetMode mode /*= FULL_ATOM_t*/ ) const {
	debug_assert( core::Size(mode) >= 1 );
	debug_assert( mode <= TYPE_SET_MODES_LENGTH );
	return res_type_sets_[ mode ];
}

/// @brief Replace the current ResidueTypeSet of the appropriate type
/// with the given RTS.
void
CacheableResidueTypeSets::set_res_type_set( PoseResidueTypeSetOP rts, TypeSetMode mode /*= INVALID_t*/ ) {
	if ( mode == INVALID_t ) {
		debug_assert( rts != nullptr );
		mode = rts->mode();
	}
	debug_assert( core::Size(mode) >= 1 );
	debug_assert( mode <= TYPE_SET_MODES_LENGTH );
	res_type_sets_[ mode ] = rts;
}

} //core
} //chemical

#ifdef    SERIALIZATION

template< class Archive >
void
core::chemical::CacheableResidueTypeSets::save( Archive & arc ) const {
	arc( CEREAL_NVP( res_type_sets_ ) ); // utility::vector1< PoseResidueTypeSetOP >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::CacheableResidueTypeSets::load( Archive & arc ) {
	arc( res_type_sets_ ); // utility::vector1< PoseResidueTypeSetOP >
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::CacheableResidueTypeSets );
CEREAL_REGISTER_TYPE( core::chemical::CacheableResidueTypeSets )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_CacheableResidueTypeSets )
#endif // SERIALIZATION




