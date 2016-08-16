// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/trie/trie_vs_trie.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_trie_trie_vs_trie_hh
#define INCLUDED_core_scoring_trie_trie_vs_trie_hh

// Unit Headers
#include <core/scoring/trie/trie_vs_trie.fwd.hh>

// Package Headers

// Project Headers
#include <core/types.hh>


// STL Headers

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2A.hh>

// Utility Headers

#include <utility/vector1_bool.hh>


namespace core {
namespace scoring {
namespace trie {

//typedef core::PackerEnergy core::PackerEnergy;

/// @brief trie vs trie algorithm, templated on the Atom type that each TrieAtom is templated on,
/// along with the Count Pair data (two tries must share the same Atom type to be used for the
/// trie-vs-trie algorithm, but may contain different peices of count-pair data), on the kind of
/// count pair function used, and finally, on the score function itself.
///
template < class AT, class CPDAT1, class CPDAT2, class CPFXN, class SFXN >
void
trie_vs_trie(
	RotamerTrie< AT, CPDAT1 > const & trie1,
	RotamerTrie< AT, CPDAT2 > const & trie2,
	CPFXN & count_pair,
	SFXN & score_function,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	/*
	std::cout << "I made it" << std::endl;
	score_function.print();
	trie1.print();
	trie2.print();
	CPFXN::print();
	*/

	using namespace ObjexxFCL;

	ObjexxFCL::FArray2A< core::PackerEnergy > rot_rot_table(
		temp_table,
		trie2.num_unique_rotamers(),
		trie1.num_unique_rotamers() );

	DistanceSquared const hydrogen_interaction_cutoff = score_function.hydrogen_interaction_cutoff2();

	typename utility::vector1< TrieNode < AT, CPDAT1 > > const & trie1_atoms = trie1.atoms();
	Size const trie1_natoms = trie1.atoms().size();
	//Size const trie1_num_unique_rotamers = trie1.num_unique_rotamers();

	typename utility::vector1 < TrieNode < AT, CPDAT2 > > const & trie2_atoms = trie2.atoms();
	Size const trie2_natoms = trie2.atoms().size();

	Size const trie2_num_heavyatoms = trie2.num_heavy_atoms();
	Size const trie2_num_unique_rotamers = trie2.num_unique_rotamers();

	Size const trie1_max_branch_depth = trie1.max_branch_depth();
	Size const trie1_max_heavyatom_depth = trie1.max_heavyatom_depth();

	FArray2D< core::PackerEnergy > at_v_rot_stack( trie2_num_unique_rotamers, trie1_max_branch_depth, 0.0);

	FArray2D_int parent_heavy_wi_hydrogen_cutoff(trie2_num_heavyatoms, trie1_max_heavyatom_depth, true);
	FArray2D_int r_heavy_skip_s_subtree(trie2_num_heavyatoms, trie1_max_heavyatom_depth, false);

	utility::vector1< core::PackerEnergy > energy_stack( std::max( trie1.max_atom_depth() + 1, trie2.max_atom_depth() + 1 ));
	energy_stack[1] = 0;
	//energy_stack[0] = 0;

	utility::vector1< Size >  r_heavy_depth_stack( trie1.max_branch_depth() + 1 );
	utility::vector1< Size >  s_heavy_depth_stack( trie2.max_branch_depth() + 1 );

	//utility::vector1< bool > parent_heavy_wi_hcut_stack( trie2.max_heavyatom_depth() );
	int * parent_heavy_wi_hcut_stack = new int[ trie2.max_heavyatom_depth() + 2 ]; // 17% of time spent looking up utility::vector1<bool>

	utility::vector1< Size > s_sibling_stack( trie2.max_atom_depth() + 1 );
	s_sibling_stack[1] = trie2_natoms + 1; // out-of-range sibling
	for ( Size ii = 2; ii <= trie2.max_atom_depth(); ++ii ) {
		s_sibling_stack[ii] = 0;
	}

	Size r_curr_stack_top = 2; //immediately decremented to 1
	Size s_curr_stack_top;
	Size r_rotamers_seen = 0;
	Size s_rotamers_seen;
	Size s_heavyatoms_seen;

	r_heavy_depth_stack[1] = 0;

	for ( Size ii = 1; ii <= trie1_natoms; ++ii ) {
		//std::cout << "tvt with ii = " << ii << std::endl;
		/*typename*/ TrieNode< AT, CPDAT1 > const & r = trie1_atoms[ ii ];
		energy_stack[1] = 0;

		s_rotamers_seen = 0;
		s_heavyatoms_seen = 0;

		if ( r.first_atom_in_branch() ) --r_curr_stack_top;

		if ( r.has_sibling() ) { //push - copy stack downwards

			++r_curr_stack_top;
			//std::cout << "r.has_sibling(): new stack top: " << r_curr_stack_top << std::endl;

			// when I previously dimensioned by trie1_num_unique_rotamers (a bug) I didn't get an assertion failure here!
			// wtf?!
			//FArray1A< core::PackerEnergy > at_rot_array_d_proxy(at_v_rot_stack(1, r_curr_stack_top), trie2_num_unique_rotamers);
			//FArray1A< core::PackerEnergy > at_rot_array_dminus1_proxy(at_v_rot_stack(1, r_curr_stack_top-1), trie2_num_unique_rotamers);
			//at_rot_array_d_proxy = at_rot_array_dminus1_proxy;
			for ( Size jj = 1, litop = at_v_rot_stack.index( 1, r_curr_stack_top ),
					liprev = at_v_rot_stack.index( 1, r_curr_stack_top - 1 );
					jj <= trie2_num_unique_rotamers; ++jj, ++litop, ++liprev ) {
				at_v_rot_stack[ litop ] = at_v_rot_stack[ liprev ];
			}

			r_heavy_depth_stack[ r_curr_stack_top ] = r_heavy_depth_stack[ r_curr_stack_top-1 ];
			//std::cout << "r_heavy_depth_stack[" << r_curr_stack_top << "] = " << r_heavy_depth_stack[ r_curr_stack_top ] << std::endl;
		}

		//++r_tree_depth_stack[ r_curr_stack_top ];

		if ( ! r.is_hydrogen() ) ++r_heavy_depth_stack[ r_curr_stack_top ];


		//FArray1A< core::PackerEnergy > at_rot_array_proxy(at_v_rot_stack(1, r_curr_stack_top), trie2_num_unique_rotamers);

		//FArray1A_int parent_wi_h_dist( parent_heavy_wi_hydrogen_cutoff( 1,
		// r_heavy_depth_stack[ r_curr_stack_top ] ), trie2_num_heavyatoms );

		//FArray1A_int r_heavy_skip_s_subtree_proxy
		// ( r_heavy_skip_s_subtree(1, r_heavy_depth_stack[ r_curr_stack_top ] ), trie2_num_heavyatoms);

		if ( r.is_hydrogen() ) {
			s_curr_stack_top = 2;
			s_heavy_depth_stack[1] = 0;
			//s_tree_depth_stack[1] = 0;

			for ( Size jj = 1; jj <= trie2_natoms; ++jj ) {
				//std::cerr << "  ii: " << ii << " is H, jj : " << jj << " " << s_curr_stack_top << " ";
				//std::cerr << "s_sibling_stack [";
				//for ( Size kk = 1; kk <= s_curr_stack_top; ++kk )
				//  std::cerr << s_sibling_stack[kk] << " ";
				//std::cerr << " " <<  std::endl;


				TrieNode< AT, CPDAT2 > const & s  = trie2_atoms[ jj ];

				if ( s.first_atom_in_branch() ) {
					--s_curr_stack_top;
				}

				if ( s.has_sibling() ) { //push - copy stack downwards
					++s_curr_stack_top;
					energy_stack[s_curr_stack_top] = energy_stack[s_curr_stack_top - 1];
					s_heavy_depth_stack[ s_curr_stack_top] = s_heavy_depth_stack[ s_curr_stack_top - 1];
					//s_tree_depth_stack[ s_curr_stack_top]
					// = s_tree_depth_stack[ s_curr_stack_top - 1];
					s_sibling_stack[s_curr_stack_top] = s.sibling();
				}
				//++s_tree_depth_stack[ s_curr_stack_top];

				if ( !s.is_hydrogen() ) {
					++s_heavyatoms_seen;
					++s_heavy_depth_stack[s_curr_stack_top];
					parent_heavy_wi_hcut_stack[s_heavy_depth_stack[s_curr_stack_top]]
						= parent_heavy_wi_hydrogen_cutoff(s_heavyatoms_seen, r_heavy_depth_stack[ r_curr_stack_top ] );

					if ( r_heavy_skip_s_subtree(s_heavyatoms_seen, r_heavy_depth_stack[ r_curr_stack_top ]) ) {
						if ( energy_stack[s_curr_stack_top] != 0.0f ) {
							for ( Size kk = s_rotamers_seen + 1;
									kk <= s_rotamers_seen + s.num_rotamers_in_subtree();
									++kk ) {
								at_v_rot_stack(kk, r_curr_stack_top) += energy_stack[s_curr_stack_top];
							}
						}
						s_rotamers_seen += s.num_rotamers_in_subtree();
						jj = s_sibling_stack[s_curr_stack_top] - 1;
						continue;
					}
				}

				Real weight(1.0); Size path_dist(0);
				if ( parent_heavy_wi_hcut_stack[s_heavy_depth_stack[s_curr_stack_top] ] &&
						count_pair( r.cp_data(), s.cp_data(), weight, path_dist) ) {
					if ( s.is_hydrogen() ) {
						core::PackerEnergy e = score_function.hydrogenatom_hydrogenatom_energy(r.atom(), s.atom(), path_dist );
						energy_stack[ s_curr_stack_top ] += weight * e;
						//std::cout << "h/h atom pair energy: " << ii << " & " << jj << " = " << weight * e << "( unweighted: " << e <<  ") estack: " << energy_stack[ s_curr_stack_top ] << std::endl;

					} else {
						core::PackerEnergy e = score_function.hydrogenatom_heavyatom_energy( r.atom(), s.atom(), path_dist );
						energy_stack[ s_curr_stack_top ] += weight * e;
						//std::cout << "h/hv atom pair energy: " << ii << " & " << jj << " = " << weight * e << "( unweighted: " << e << ") estack: " << energy_stack[ s_curr_stack_top ] << std::endl;
					}
				}

				if ( s.is_rotamer_terminal() ) {
					++s_rotamers_seen;
					at_v_rot_stack(s_rotamers_seen, r_curr_stack_top) += energy_stack[ s_curr_stack_top ];
				}
			}
		} else { // r is a heavy atom
			s_curr_stack_top = 2;
			s_heavy_depth_stack[1] = 0;
			//s_tree_depth_stack[1] = 0;

			for ( Size jj = 1,
					liskip = r_heavy_skip_s_subtree.index( 1, r_heavy_depth_stack[ r_curr_stack_top ] );
					jj <= trie2_num_heavyatoms; ++jj, ++liskip ) {
				r_heavy_skip_s_subtree[ liskip ] = false;
			}

			for ( Size jj = 1; jj <= trie2_natoms; ++jj ) {
				//std::cerr << "! ii: " << ii << " is H, jj : " << jj << " " << s_curr_stack_top <<
				// " " << s_rotamers_seen << " " << s_heavyatoms_seen << std::endl;
				//trie_node & s = rt_trie[jj];
				TrieNode< AT, CPDAT2 > const & s  = trie2_atoms[ jj ];

				//std::cerr << "s_sibling_stack [";
				//for ( Size kk = 1; kk <= s_curr_stack_top; ++kk )
				//  std::cerr << s_sibling_stack[kk] << " ";
				//std::cerr << " " <<  std::endl;

				if ( s.first_atom_in_branch() ) --s_curr_stack_top;

				if ( s.has_sibling() ) { //push - copy stack downwards
					++s_curr_stack_top;
					energy_stack[s_curr_stack_top] = energy_stack[s_curr_stack_top - 1];
					s_heavy_depth_stack[ s_curr_stack_top] = s_heavy_depth_stack[ s_curr_stack_top - 1];
					//s_tree_depth_stack[ s_curr_stack_top]
					// = s_tree_depth_stack[ s_curr_stack_top - 1];
					s_sibling_stack[s_curr_stack_top] = s.sibling();
				}
				//++s_tree_depth_stack[ s_curr_stack_top];

				if ( s.is_hydrogen() ) {
					Real weight(1.0); Size path_dist(0);
					if ( parent_heavy_wi_hcut_stack[s_heavy_depth_stack[s_curr_stack_top]] &&
							count_pair( r.cp_data(), s.cp_data(), weight, path_dist ) ) {
						core::PackerEnergy e = score_function.heavyatom_hydrogenatom_energy( r.atom(), s.atom(), path_dist );
						energy_stack[ s_curr_stack_top ] += weight * e;
						//std::cout << "hv/h atom pair energy: " << ii << " & " << jj << " = " << weight * e  << "( unweighted: " << e << ") estack: " << energy_stack[ s_curr_stack_top ] << std::endl;
						/*,trie1.res_id(), trie2.res_id(),
						tri1.num_neighbors(), trie2.num_neighbors(),
						Whbond);*/
					}

				} else { // s is a heavy atom
					++s_heavyatoms_seen;
					++s_heavy_depth_stack[s_curr_stack_top];

					DistanceSquared d2(0.0);
					Real weight = 1;
					Size path_dist(0);

					if ( count_pair( r.cp_data(), s.cp_data(), weight, path_dist) ) {
						core::PackerEnergy e = score_function.heavyatom_heavyatom_energy(r.atom(), s.atom(), d2, path_dist);
						//std::cout << "hv/hv atom pair energy: " << ii << " & " << jj << " = " << weight * e << "( unweighted: " << e << ") estack: " << energy_stack[ s_curr_stack_top ] << std::endl;
						energy_stack[s_curr_stack_top] += weight * e;
					} else {
						/// compute d2
						d2 = r.atom().xyz().distance_squared( s.atom().xyz() );
					}

					parent_heavy_wi_hcut_stack[s_heavy_depth_stack
						[s_curr_stack_top]] =
						parent_heavy_wi_hydrogen_cutoff(s_heavyatoms_seen, r_heavy_depth_stack[ r_curr_stack_top ]) =
						(d2 < hydrogen_interaction_cutoff);

					if ( d2 > s.subtree_interaction_sphere_square_radius() ) {
						r_heavy_skip_s_subtree(s_heavyatoms_seen, r_heavy_depth_stack[ r_curr_stack_top ]) = true;
						if ( energy_stack[s_curr_stack_top] != 0. ) {
							for ( Size kk = s_rotamers_seen + 1,
									li_avrstack = at_v_rot_stack.index( kk, r_curr_stack_top );
									kk <= s_rotamers_seen + s.num_rotamers_in_subtree();
									++kk, ++li_avrstack ) {
								at_v_rot_stack[ li_avrstack ] += energy_stack[s_curr_stack_top];
							}
						}
						//std::cout << "heavy atom / subtree prune: " << ii << " " << jj << " : jumping to sibling " << s_sibling_stack[s_curr_stack_top] << std::endl;
						jj = s_sibling_stack[s_curr_stack_top] - 1;
						s_rotamers_seen += s.num_rotamers_in_subtree();
						continue;
					}
				}
				if ( s.is_rotamer_terminal() ) {
					++s_rotamers_seen;
					//std::cout << "s.is_rotamer_terminal " << s_rotamers_seen << " with energy: " << energy_stack[ s_curr_stack_top ] << std::endl;
					at_v_rot_stack(s_rotamers_seen, r_curr_stack_top) += energy_stack[s_curr_stack_top];
				}
			}
		}
		if ( r.is_rotamer_terminal() ) {
			++r_rotamers_seen;
			//std::cout << "terminal #" << r_rotamers_seen << ": writing at_rot_array_proxy: ";
			//for ( Size jj = 1; jj <= trie2_num_unique_rotamers; ++jj ) {
			// std::cout << " " << at_rot_array_proxy( jj );
			//}
			//std::cout << std::endl;
			//FArray1A< core::PackerEnergy > rot_rot_table_row( rot_rot_table(1, r_rotamers_seen), trie2_num_unique_rotamers);
			//rot_rot_table_row.dimension(trie2_num_unique_rotamers);
			//rot_rot_table_row = at_rot_array_proxy;
			for ( Size jj = 1; jj <= trie2_num_unique_rotamers; ++jj ) {
				rot_rot_table( jj, r_rotamers_seen ) = at_v_rot_stack( jj, r_curr_stack_top );
			}


			//std::cout << "rot rot table, column " << r_rotamers_seen << std::endl;
			//for ( Size jj = 1; jj <= trie2_num_unique_rotamers; ++jj ) {
			// std::cout << " " << rot_rot_table( jj, r_rotamers_seen );
			//}
			//std::cout << std::endl;

		}


	}
	convert_inorder_table_to_original_order(
		trie1.total_rotamers_2_unique_rotamers(),
		trie2.total_rotamers_2_unique_rotamers(),
		pair_energy_table,
		rot_rot_table );
	//std::cout << "complete trie-vs-trie" << std::endl;

	delete [] parent_heavy_wi_hcut_stack;
}


} // namespace trie
} // namespace scoring
} // namespace core


#endif
