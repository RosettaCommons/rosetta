// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/trie/TrieCollection.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/trie/TrieCollection.hh>

// Package Headers
#include <core/scoring/trie/RotamerTrieBase.hh>

// STL Headers

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
	std::fill( tries_.begin(), tries_.end(), RotamerTrieBaseOP( 0 ) );
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

