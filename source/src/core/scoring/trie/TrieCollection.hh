// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/trie/TrieCollection.hh
/// @brief Container class for storing rotamer tries in an Energies object.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_trie_TrieCollection_hh
#define INCLUDED_core_scoring_trie_TrieCollection_hh

// Unit Headers
#include <core/scoring/trie/TrieCollection.fwd.hh>

// Package Headers
#include <core/scoring/trie/RotamerTrieBase.fwd.hh>

// Project Headers
#include <basic/datacache/CacheableData.hh>

// Utility Headers

#include <utility/vector1_bool.hh>


namespace core {
namespace scoring {
namespace trie {

class TrieCollection : public basic::datacache::CacheableData
{
public:
	RotamerTrieBaseCOP
	trie( Size index ) const;

	void total_residue( Size );
	Size total_residue() const;

	void trie( Size index, RotamerTrieBaseOP new_trie );

	virtual basic::datacache::CacheableDataOP clone() const;

private:
	utility::vector1< RotamerTrieBaseOP > tries_;

};


} // namespace trie
} // namespace scoring
} // namespace core

#endif
