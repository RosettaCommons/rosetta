// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/trie/trie_vs_trie.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/trie/trie_vs_trie.hh>

// Project Headers
#include <core/types.hh>

namespace core {
namespace scoring {
namespace trie {

void
convert_inorder_table_to_original_order(
	utility::vector1< Size > const & total_rotamers_2_unique_rotamers_1,
	utility::vector1< Size > const & total_rotamers_2_unique_rotamers_2,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2A< core::PackerEnergy > const & rot_rot_table
) {
	using namespace ObjexxFCL;

	Size const nrots1 = total_rotamers_2_unique_rotamers_1.size();
	Size const nrots2 = total_rotamers_2_unique_rotamers_2.size();

	FArray2A< core::PackerEnergy > etable_proxy( pair_energy_table, nrots2, nrots1 );
	for ( Size ii = 1; ii <= nrots1; ++ii ) {
		Size const ii_inorder = total_rotamers_2_unique_rotamers_1[ ii ];
		for ( Size jj = 1; jj <= nrots2; ++jj ) {
			Size const jj_inorder = total_rotamers_2_unique_rotamers_2[ jj ];
			etable_proxy( jj, ii ) = rot_rot_table( jj_inorder, ii_inorder );
			//std::cout << "Writing inorder cell (" << jj_inorder << ", " << ii_inorder << ") to orig_order cell (";
			//std::cout << jj << ", " << ii << ") with value: " << rot_rot_table( jj_inorder, ii_inorder ) << std::endl;
		}
	}
	//std::cout << "finished converting inorder table" << std::endl;
}


} // namespace trie
} // namespace scoring
} // namespace core

