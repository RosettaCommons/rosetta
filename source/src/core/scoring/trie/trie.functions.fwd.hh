// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/trie/trie.functions.fwd.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_trie_trie_functions_fwd_hh
#define INCLUDED_core_scoring_trie_trie_functions_fwd_hh

// Package Headers
#include <core/scoring/trie/CPDataCorrespondence.fwd.hh>
#include <core/scoring/trie/RotamerTrieBase.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>


// Project Headers
#include <core/conformation/Residue.fwd.hh>
//#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/types.hh>

#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {
namespace trie {

/// @brief Create a trie where cpdata is initialized to reflect
/// path distances to connection-point atoms.
template < class AT, class CPDAT >
RotamerTrieBaseOP
create_trie(
	conformation::RotamerSetBase const & rotset,
	AT    const & /* dummy variable for type identification */,
	CPDAT const & /* dummy variable for type identification */,
	CPDataCorrespondence const & cpdata_map,
	Distance atomic_interaction_cutoff

);


template < class CPDAT >
void
initialize_cpdata_for_atom(
	CPDAT & cpdata,
	Size atom_index,
	conformation::Residue const & res,
	CPDataCorrespondence const & cpdata_map
);


} // namespace trie
} // namespace scoring
} // namespace core

#endif
