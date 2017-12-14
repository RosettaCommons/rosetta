// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/trie/TrieCollection.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/trie/TrieCollection.hh>

// Package Headers
#include <core/scoring/trie/RotamerTrieBase.hh>

// STL Headers

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace trie {

RotamerTrieBaseCOP
TrieCollection::trie( Size index ) const
{
	return tries_[ index ];
}

void
TrieCollection::total_residue( Size total_residue_in )
{
	tries_.resize( total_residue_in );
	std::fill( tries_.begin(), tries_.end(), RotamerTrieBaseOP( nullptr ) );
}

Size
TrieCollection::total_residue() const
{
	return tries_.size();
}

void
TrieCollection::trie( Size index, RotamerTrieBaseOP new_trie )
{
	tries_[ index ] = new_trie;
}

basic::datacache::CacheableDataOP
TrieCollection::clone() const
{
	return basic::datacache::CacheableDataOP( new TrieCollection( *this ) );
}


} // namespace trie
} // namespace scoring
} // namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::trie::TrieCollection::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( tries_ ) ); // utility::vector1<RotamerTrieBaseOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::trie::TrieCollection::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( tries_ ); // utility::vector1<RotamerTrieBaseOP>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::trie::TrieCollection );
CEREAL_REGISTER_TYPE( core::scoring::trie::TrieCollection )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_trie_TrieCollection )
#endif // SERIALIZATION
