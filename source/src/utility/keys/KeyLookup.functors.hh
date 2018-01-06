// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/KeyLookup.functors.hh
/// @brief  utility::keys::KeyLookup functor declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note Convenience header for inclusion/injection into key collection namespaces
/// @note Must typedef the collection key type to KeyType above point of inclusion

#include <utility/keys/KeyLookup.hh>

/// @brief Key collection iterators
typedef  utility::keys::KeyLookup< KeyType >::const_iterator  const_iterator;
typedef  utility::keys::KeyLookup< KeyType >::const_iterator  ConstIterator;

/// @brief Lookup functors
extern utility::keys::lookup::has  < KeyType > const has  ; // Provides has( id ) and has( key )
extern utility::keys::lookup::key  < KeyType > const key  ; // Provides key( id ) and key( key )
extern utility::keys::lookup::gen  < KeyType > const gen  ; // Provides gen( id, identifier, code ) and gen( key )
extern utility::keys::lookup::n_key< KeyType > const n_key; // Provides n_key()
extern utility::keys::lookup::begin< KeyType > const begin; // Provides begin()
extern utility::keys::lookup::end  < KeyType > const end  ; // Provides end()
