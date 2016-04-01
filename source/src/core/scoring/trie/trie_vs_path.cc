// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/trie/trie_vs_trie.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_packer_trie_trie_vs_trie_HH
#define INCLUDED_core_pack_packer_trie_trie_vs_trie_HH

// Unit Headers
#include <core/scoring/trie/trie_vs_path.hh>

// Project Headers
#include <core/types.hh>

// STL Headers

// Utility Headers

namespace core {
namespace scoring {
namespace trie {

void
convert_inorder_vector_to_original_order(
	utility::vector1< Size > const & total_rotamers_2_unique_rotamers,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< Energy > const & temp_vector
) {
	Size const nrots1 = total_rotamers_2_unique_rotamers.size();

	for ( Size ii = 1; ii <= nrots1; ++ii ) {
		Size const ii_inorder = total_rotamers_2_unique_rotamers[ ii ];
		pair_energy_vector[ ii ] = temp_vector[  ii_inorder ];
	}
	//std::cout << "finished converting inorder vector" << std::endl;
}

} // namespace trie
} // namespace scoring
} // namespace core


#endif
