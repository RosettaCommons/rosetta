// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/trie/trie_vs_trie.fwd.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_trie_trie_vs_trie_fwd_hh
#define INCLUDED_core_scoring_trie_trie_vs_trie_fwd_hh

#include <core/scoring/trie/RotamerTrie.fwd.hh>

#include <core/types.hh>


// Utility Headers
#include <utility/vector1.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>

namespace core {
namespace scoring {
namespace trie {


template < class AT, class CPDAT1, class CPDAT2, class CPFXN, class SFXN >
void
trie_vs_trie(
	RotamerTrie< AT, CPDAT1 > & trie1,
	RotamerTrie< AT, CPDAT2 > & trie2,
	CPFXN & count_pair,
	SFXN & score_function,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
);


void
convert_inorder_table_to_original_order(
	utility::vector1< Size > const & total_rotamers_2_unique_rotamers_1,
	utility::vector1< Size > const & total_rotamers_2_unique_rotamers_2,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2A< core::PackerEnergy > const & inorder_table
);

} // namespace trie
} // namespace scoring
} // namespace core


#endif
