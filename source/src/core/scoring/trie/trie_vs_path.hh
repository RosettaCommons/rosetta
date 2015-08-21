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

#ifndef INCLUDED_core_scoring_trie_trie_vs_path_hh
#define INCLUDED_core_scoring_trie_trie_vs_path_hh

// Package Headers

// Project Headers
#include <core/types.hh>

// STL Headers

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>

// Utility Headers

#include <core/scoring/trie/RotamerTrie.fwd.hh>
#include <utility/vector1_bool.hh>


namespace core {
namespace scoring {
namespace trie {

void
convert_inorder_vector_to_original_order(
	utility::vector1< Size > const & total_rotamers_2_unique_rotamers,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< Energy > const & temp_vector
);


/// @brief trie vs path algorithm, templated on the Atom type that each TrieAtom is templated on,
/// along with the Count Pair data (two tries must share the same Atom type to be used for the
/// trie-vs-trie algorithm, but may contain different peices of count-pair data), on the kind of
/// count pair function used, and finally, on the score function itself.
///
template < class AT, class CPDAT1, class CPDAT2, class CPFXN, class SFXN >
void
trie_vs_path(
	RotamerTrie< AT, CPDAT1 > const & trie1,
	RotamerTrie< AT, CPDAT2 > const & trie2,
	CPFXN & count_pair,
	SFXN & score_function,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	using namespace ObjexxFCL;

	DistanceSquared const hydrogen_interaction_cutoff = score_function.hydrogen_interaction_cutoff2();

	Size trie2_num_atoms = trie2.atoms().size();
	Size const trie2_num_heavyatoms = trie2.num_heavy_atoms();

	Size const trie1_num_atoms = trie1.atoms().size();

	//trie_node* this_trie = get_trie_for_curr_bb();
	//trie_node* bg_trie = bg.get_trie_for_curr_bb();

	typename utility::vector1< TrieNode < AT, CPDAT1 > > const & trie1_atoms = trie1.atoms();
	typename utility::vector1 < TrieNode < AT, CPDAT2 > > const & trie2_atoms = trie2.atoms();

	Size const trie1_max_branch_depth = trie1.max_branch_depth();

	utility::vector1< Energy > energy_sum( trie1.max_atom_depth() + 1 );
	energy_sum[1] = 0;
	//FArray2D_bool parent_heavy_wi_hydrogen_cutoff( trie2_num_heavyatoms, trie1.max_heavyatom_depth(), true);
	//FArray2D_bool r_heavy_skip_s_subtree(trie2_num_heavyatoms, trie1.max_heavyatom_depth(), false);
	FArray2D_int parent_heavy_wi_hydrogen_cutoff( trie2_num_heavyatoms, trie1.max_heavyatom_depth(), true);
	FArray2D_int r_heavy_skip_s_subtree(trie2_num_heavyatoms, trie1.max_heavyatom_depth(), false);


	utility::vector1< Size > r_heavy_depth_stack( trie1_max_branch_depth );
	r_heavy_depth_stack[1] = 0;

	Size r_num_rotamers_seen = 0;
	Size r_curr_stack_top = 2;
	Size s_heavy_depth;

	for ( Size ii = 1; ii <= trie1_num_atoms; ++ii ) {
		TrieNode< AT, CPDAT1 > const & r = trie1_atoms[ ii ];

		if ( r.first_atom_in_branch() ) --r_curr_stack_top;
		if ( r.has_sibling() ) { //push - copy stack downwards
			++r_curr_stack_top;
			energy_sum[ r_curr_stack_top ] = energy_sum[ r_curr_stack_top - 1];
			r_heavy_depth_stack[ r_curr_stack_top ] = r_heavy_depth_stack[ r_curr_stack_top-1 ];
		}

		//++r_tree_depth_stack[ r_curr_stack_top ];
		if ( ! r.is_hydrogen() ) ++r_heavy_depth_stack[ r_curr_stack_top ];

		//FArray1A_bool parent_wi_h_dist
		// ( parent_heavy_wi_hydrogen_cutoff( 1,r_heavy_depth_stack[ r_curr_stack_top ] ) );
		//parent_wi_h_dist.dimension(trie2_num_heavyatoms);

		//FArray1A_int r_heavy_skip_s_subtree_proxy
		// ( r_heavy_skip_s_subtree(1, r_heavy_depth_stack[ r_curr_stack_top ] ));
		//r_heavy_skip_s_subtree_proxy.dimension(trie2_num_heavyatoms);

		s_heavy_depth = 0;

		if ( r.is_hydrogen() ) {
			for ( Size jj = 1; jj <= trie2_num_atoms; ++jj ) {
				//trie_node & s = bg_trie[jj];
				TrieNode< AT, CPDAT2 > const & s  = trie2_atoms[ jj ];
				Real weight( 1.0 ); Energy e( 0.0 ); Size path_dist(0);
				if ( s.is_hydrogen() ) {
					if ( parent_heavy_wi_hydrogen_cutoff( s_heavy_depth, r_heavy_depth_stack[ r_curr_stack_top ]  ) &&
							count_pair( r.cp_data(), s.cp_data(), weight, path_dist ) ) {
						e = score_function.hydrogenatom_hydrogenatom_energy( r.atom(), s.atom(), path_dist );
						energy_sum[ r_curr_stack_top ] += weight * e;
					}
				} else {
					++s_heavy_depth;

					if ( r_heavy_skip_s_subtree(s_heavy_depth,r_heavy_depth_stack[ r_curr_stack_top ]) ) break;

					if ( parent_heavy_wi_hydrogen_cutoff(s_heavy_depth, r_heavy_depth_stack[ r_curr_stack_top ] ) &&
							count_pair( r.cp_data(), s.cp_data(), weight, path_dist ) ) {
						e = score_function.hydrogenatom_heavyatom_energy( r.atom(), s.atom(), path_dist );
						energy_sum[ r_curr_stack_top ] += e * weight;
					}

				}
			}
		} else { //else r is a heavy atom
			for ( Size jj = 1; jj <= trie2_num_heavyatoms; ++jj ) {
				r_heavy_skip_s_subtree( jj, r_heavy_depth_stack[ r_curr_stack_top ] ) = false;
			}
			for ( Size jj = 1; jj <= trie2_num_atoms; ++jj ) {
				TrieNode< AT, CPDAT2 > const & s  = trie2_atoms[ jj ];
				Real weight( 1.0 ); Energy e( 0.0 ); Size path_dist(0);
				if ( s.is_hydrogen() ) {
					if ( parent_heavy_wi_hydrogen_cutoff( s_heavy_depth, r_heavy_depth_stack[ r_curr_stack_top ] ) &&
							count_pair( r.cp_data(), s.cp_data(), weight, path_dist ) ) {
						e = score_function.heavyatom_hydrogenatom_energy( r.atom(), s.atom(), path_dist );
						energy_sum[ r_curr_stack_top ] += weight * e;
					}
				} else {
					++s_heavy_depth;
					DistanceSquared d2;

					if ( count_pair( r.cp_data(), s.cp_data(), weight, path_dist ) ) {
						e = score_function.heavyatom_heavyatom_energy(r.atom(), s.atom(), d2, path_dist);
						//std::cout << "hv/hv atom pair energy: " << ii << " & " << jj << " = " << weight * e << "( unweighted: " << e << ") estack: " << energy_stack[ s_curr_stack_top ] << std::endl;
						energy_sum[ r_curr_stack_top ] += weight * e;
					} else {
						/// compute d2
						d2 = r.atom().xyz().distance_squared( s.atom().xyz() );
					}

					parent_heavy_wi_hydrogen_cutoff( s_heavy_depth, r_heavy_depth_stack[ r_curr_stack_top ] ) = (d2 < hydrogen_interaction_cutoff);
					if ( d2 > s.subtree_interaction_sphere_square_radius() ) {
						r_heavy_skip_s_subtree(s_heavy_depth, r_heavy_depth_stack[ r_curr_stack_top ]) = true;
						break;
					}
				}
			}
		}

		if ( r.is_rotamer_terminal() ) {
			++r_num_rotamers_seen;
			temp_vector[r_num_rotamers_seen] = energy_sum[ r_curr_stack_top ];
		}
	}

	convert_inorder_vector_to_original_order(
		trie1.total_rotamers_2_unique_rotamers(),
		pair_energy_vector,
		temp_vector );

}

/* accidental path-vs-path alg in mini format -- use this later
{


Energy esum = 0;
Size s_heavy_depth = 0;

Size const trie1_num_atoms = trie1.atoms().size();
Size const trie2_num_atoms = trie2.atoms().size();
Size const trie2_num_heavyatoms = trie2.num_heavy_atoms();
Size const trie2_num_heavyatoms_p_1 = trie2_num_heavyatoms + 1;

//trie_node* rt_trie = rt.get_trie_for_curr_bb();
//trie_node* this_trie = get_trie_for_curr_bb();

typename utility::vector1< TrieNode < AT, CPDAT1 > > const & trie1_atoms = trie1.atoms();
typename utility::vector1 < TrieNode < AT, CPDAT2 > > const & trie2_atoms = trie2.atoms();

//FArray1D_bool parent_wi_h_dist( rt_num_heavyatoms, false );
bool  parent_wi_h_dist[ 40 ]; //<--- replace this magic number ASAP
//FArray1D_bool r_heavy_skip_s_subtree( rt_num_heavyatoms, false );
Size r_heavy_skip_s_subtree_depth = rt_num_heavyatoms_p_1;

for (int Size = 1; ii <= trie1_num_atoms; ++ii )
{
//trie_node & r = this_trie[ii];
TrieNode< AT, CPDAT1 > const & r = trie1_atoms[ ii ];
if ( r.is_hydrogen() && r_heavy_skip_s_subtree_depth == 1)
{
continue;
}

//unsigned short * r_class_masks = NULL;
//for (int eq_class = 0; eq_class < 6; ++eq_class)
//{ if (this_equivalence_class_bitmasks[eq_class] & r.flags2_)
// {  r_class_masks = rt_equivalence_class_bitmasks[eq_class];
//  break;
// }
//}
//debug_assert( r_class_masks );

s_heavy_depth = 0;
if ( r.is_hydrogen() ) {
for ( Size jj = 1; jj <= trie2_num_atoms; ++jj ) {
//trie_node & s = rt_trie[jj];
TrieNode< AT, CPDAT2 > const & s  = trie2_atoms[ jj ];
Real weight( 1.0 ), e( 0.0 );
if ( s.is_hydrogen() ) {
if (parent_wi_h_dist[ s_heavy_depth ] &&
count_pair( r.cp_data(), s.cp_data(), weight)) {
e = score_function.hydrogenatom_hydrogenatom_energy(r.atom(), s.atom() );
esum += weight * e;
}
} else {
++s_heavy_depth;

if (r_heavy_skip_s_subtree_depth == s_heavy_depth) break;

if ( parent_wi_h_dist[s_heavy_depth] &&
count_pair( r.cp_data(), s.cp_data(), weight)) {
e = score_function.hydrogenatom_heavyatom_energy( r.atom(), s.atom() );
esum += e * weight;
}
}
}
} else {
r_heavy_skip_s_subtree_depth = rt_num_heavyatoms_p_1;
for ( Size jj = 1; jj <= rt_num_atoms; ++jj) {
trie_node & s = rt_trie[jj];
Real weight( 1.0 );
if ( s.is_hydrogen()) {
if (parent_wi_h_dist[ s_heavy_depth ]  &&
count_pair( r.cp_data(), s.cp_data(), weight) ) {
Energy e = score_function.heavyatom_hydrogenatom_energy( r.atom(), s.atom() );
esum += e * weight;
}
}
else
{
++s_heavy_depth;
DistanceSquared d2;
if ( count_pair( r.cp_data(), s.cp_data(), weight) ) {
Energy e = score_function.heavyatom_heavyatom_energy(r.atom(), s.atom(), d2);
//std::cout << "hv/hv atom pair energy: " << ii << " & " << jj << " = " << weight * e << "( unweighted: " << e << ") estack: " << energy_stack[ s_curr_stack_top ] << std::endl;
esum += weight * e;
} else {
/// compute d2
d2 = r.atom().xyz().distance_squared( s.atom().xyz() );
}

parent_wi_h_dist[s_heavy_depth] = (d2 < hydrogen_interaction_cutoff);
if (d2 > s.subtree_interaction_sphere_square_radius())
{
r_heavy_skip_s_subtree_depth = s_heavy_depth;
break;
}
}
}
}
}

return esum;
}
*/

} // namespace trie
} // namespace scoring
} // namespace core


#endif
