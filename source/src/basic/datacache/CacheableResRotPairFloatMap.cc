// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableResRotPairFloatMap.cc
/// @brief  A CacheableData that stores a std::map<ResRotPair,float>
/// @author Brian Coventry
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Added caching of floats indexed by rotamer memory address,
/// used _only_ for self-interactions in the symmetric case.

// unit headers
#include <basic/datacache/CacheableResRotPairFloatMap.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

#include <numeric/MathMatrix.srlz.hh>

// Cereal headers
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/map.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

CacheableResRotPairFloatMap::CacheableResRotPairFloatMap( )
:
	CacheableData(),
	shallow_copy_( false )
{
	mapOP_ = utility::pointer::make_shared< std::unordered_map< ResRotPair, float, ResRotPairHasher > >();
	four_int_indexed_mapOP_ = utility::pointer::make_shared< std::map< std::tuple< platform::Size, platform::Size, platform::Size, platform::Size>, float > >();
}

CacheableResRotPairFloatMap::CacheableResRotPairFloatMap( CacheableResRotPairFloatMap const & ot ) : CacheableData( ot ) {
	mapOP_ = utility::pointer::make_shared< std::unordered_map< ResRotPair, float, ResRotPairHasher > >();
	four_int_indexed_mapOP_ = utility::pointer::make_shared< std::map< std::tuple< platform::Size, platform::Size, platform::Size, platform::Size>, float > >();
	*this = ot;
}

CacheableResRotPairFloatMap &
CacheableResRotPairFloatMap::operator=( CacheableResRotPairFloatMap const & ot ) {
	debug_assert( mapOP_ );
	debug_assert( four_int_indexed_mapOP_ );
	shallow_copy_ = ot.shallow_copy_;
	if ( shallow_copy_ ) {
		mapOP_ = ot.mapOP_;
		four_int_indexed_mapOP_ = ot.four_int_indexed_mapOP_;
	} else {
		// Copy constructor doesn't work for unordered_map. But assignment does.
		*mapOP_ = *ot.mapOP_;
		*four_int_indexed_mapOP_ = *ot.four_int_indexed_mapOP_;
	}

	return *this;
}

CacheableResRotPairFloatMap::~CacheableResRotPairFloatMap() = default;

CacheableDataOP
CacheableResRotPairFloatMap::clone() const {
	return utility::pointer::make_shared< CacheableResRotPairFloatMap >(*this);
}

CacheableResRotPairFloatMapOP
CacheableResRotPairFloatMap::shared_from_this() {
	return utility::pointer::static_pointer_cast<CacheableResRotPairFloatMap>( CacheableData::shared_from_this() );
}

void
CacheableResRotPairFloatMap::set_shallow_copy( bool shallow ) {
	shallow_copy_ = shallow;
}

std::unordered_map< ResRotPair, float, ResRotPairHasher > &
CacheableResRotPairFloatMap::map(){
	debug_assert( mapOP_ );
	return *mapOP_;
}

std::unordered_map< ResRotPair, float, ResRotPairHasher > const &
CacheableResRotPairFloatMap::map() const {
	debug_assert( mapOP_ );
	return *mapOP_;
}

/// @brief Nonconst access to map of four ints --> float value.
/// @details Used only to store energies of a rotamer interacting with its own symmetric copies in the symmetric case.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
std::map< std::tuple< platform::Size, platform::Size, platform::Size, platform::Size>, float > &
CacheableResRotPairFloatMap::four_int_indexed_map() {
	debug_assert( mapOP_ );
	return *four_int_indexed_mapOP_;
}

/// @brief Const access to map of four ints --> float value.
/// @details Used only to store energies of a rotamer interacting with its own symmetric copies in the symmetric case.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
std::map< std::tuple< platform::Size, platform::Size, platform::Size, platform::Size>, float > const &
CacheableResRotPairFloatMap::four_int_indexed_map() const {
	debug_assert( mapOP_ );
	return *four_int_indexed_mapOP_;
}

} // namespace datacache
} // namespace basic


#ifdef    SERIALIZATION
/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::ResRotPair::save( Archive & arc ) const {
	arc( CEREAL_NVP( first_res ) );
	arc( CEREAL_NVP( first_rot ) );
	arc( CEREAL_NVP( second_res ) );
	arc( CEREAL_NVP( second_rot ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::ResRotPair::load( Archive & arc ) {
	arc( first_res );
	arc( first_rot );
	arc( second_res );
	arc( second_rot );
}

/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::CacheableResRotPairFloatMap::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( mapOP_ ) ); // std::unordered_map<uint64_t, MathMatrix<float>>
	arc( CEREAL_NVP( four_int_indexed_mapOP_ ) );
	arc( CEREAL_NVP( shallow_copy_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::CacheableResRotPairFloatMap::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( mapOP_ ); // std::unordered_map<uint64_t, MathMatrix<float>>
	arc( four_int_indexed_mapOP_ );
	arc( shallow_copy_ );
}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::CacheableResRotPairFloatMap );
CEREAL_REGISTER_TYPE( basic::datacache::CacheableResRotPairFloatMap )

CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_CacheableResRotPairFloatMap )
#endif // SERIALIZATION
