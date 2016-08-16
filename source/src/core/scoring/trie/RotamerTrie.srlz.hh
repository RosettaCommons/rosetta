// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/trie/RotamerTrie.srlz.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_trie_RotamerTrie_SRLZ_HH
#define INCLUDED_core_scoring_trie_RotamerTrie_SRLZ_HH

#ifdef SERIALIZATION

// Unit headers
#include <core/scoring/trie/RotamerTrie.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>

namespace core {
namespace scoring {
namespace trie {

template< class AT, class CPDATA >
template< class Archive >
void
TrieNode< AT, CPDATA >::save( Archive & arc ) const
{
	arc( atom_ );
	arc( cp_data_ );
	arc( subtree_intxn_sphere_radius_sqr_ );
	arc( first_atom_in_branch_ );
	arc( is_hydrogen_ );
	arc( is_term_ );
	arc( sibling_ );
	arc( rotamers_in_subtree_ );

}

template < class AT, class CPDATA >
template < class Archive >
void
TrieNode< AT, CPDATA >::load( Archive & arc )
{
	arc( atom_ );
	arc( cp_data_ );
	arc( subtree_intxn_sphere_radius_sqr_ );
	arc( first_atom_in_branch_ );
	arc( is_hydrogen_ );
	arc( is_term_ );
	arc( sibling_ );
	arc( rotamers_in_subtree_ );

}


template < class AT, class CPDATA >
template < class Archive >
void
RotamerTrie< AT, CPDATA >::save( Archive & arc ) const
{
	arc( trie_ );
	arc( num_total_atoms_ );
	arc( num_heavyatoms_ );
	arc( max_atoms_per_rotamer_ );
	arc( num_unique_rotamers_ );
	arc( num_total_rotamers_ );
	arc( total_rotamers_2_unique_rotamers_ );
	arc( max_branch_depth_ );
	arc( max_heavyatom_depth_ );
	arc( max_atom_depth_ );
}

template < class AT, class CPDATA >
template < class Archive >
void
RotamerTrie< AT, CPDATA >::load( Archive & arc )
{
	arc( trie_ );
	arc( num_total_atoms_ );
	arc( num_heavyatoms_ );
	arc( max_atoms_per_rotamer_ );
	arc( num_unique_rotamers_ );
	arc( num_total_rotamers_ );
	arc( total_rotamers_2_unique_rotamers_ );
	arc( max_branch_depth_ );
	arc( max_heavyatom_depth_ );
	arc( max_atom_depth_ );
}

}
}
}

#endif // SERIALIZATION
#endif
