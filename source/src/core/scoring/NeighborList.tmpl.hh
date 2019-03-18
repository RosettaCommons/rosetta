// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_NeighborList_tmpl_hh
#define INCLUDED_core_scoring_NeighborList_tmpl_hh

// Unit headers
#include <core/scoring/NeighborList.hh>

// Package headers
//#include <core/scoring/EnergyGraph.hh> // necessary?
//#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/prof.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/PointGraph.hh>

namespace core {
namespace scoring {

///////////////////////////////////////////////////////////////////////////////
/// @details const so that it may be called within setup_for_scoring
/// T_Etable class must implement the following functions:
/// bool defines_intrares_energy( EnergyMap const & ) const;
/// CountPairFunctionCOP get_intrares_countpair(
///   conformation::Residue const &,
///   pose::Pose const &,
///   ScoreFunction const & ) const;
/// CountPairFunctionCOP get_count_pair_function(
///   Size &,
///   Size &,
///   pose::Pose const &,
///   ScoreFunction const & ) const;

template < class T_Etable >
void
NeighborList::setup(
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	T_Etable const & etable_method
) const
{
	using namespace etable::count_pair;
	using namespace utility::graph;

	PROF_START( basic::SETUP_NBLIST );

	Size const nres( pose.size() );

	if ( auto_update_  ) {

		if ( reference_coords_.size() == 0 ) {
			reference_coords_.resize( nres );
			wide_reference_coords_.resize( nres );
			atom_needs_update_from_wide_.resize( nres );
			//atom_has_been_updated_from_wide_.resize( nres );
			for ( Size ii = 1; ii <= nres; ++ii ) {
				Size const ii_natoms = pose.residue( ii ).natoms();
				reference_coords_[ ii ].resize( ii_natoms );
				wide_reference_coords_[ ii ].resize( ii_natoms );
				atom_needs_update_from_wide_[ ii ].resize( ii_natoms, 0 );
				//atom_has_been_updated_from_wide_[ ii ].resize( ii_natoms, 0 );
			}
		}
		debug_assert( reference_coords_.size() == nres );
		for ( Size ii = 1; ii <= nres; ++ii ) {
			debug_assert( reference_coords_[ ii ].size() == pose.residue( ii ).natoms() );
			debug_assert( wide_reference_coords_[ ii ].size() == pose.residue( ii ).natoms() );
			for ( Size jj = 1; jj <= reference_coords_[ ii ].size(); ++jj ) {
				reference_coords_[ ii ][ jj ] = pose.residue( ii ).xyz( jj );
				wide_reference_coords_[ ii ][ jj ] = pose.residue( ii ).xyz( jj );
				/// if these aren't zero, then the logic for updating atom nblists from the wide-nblists has failed
				debug_assert( atom_needs_update_from_wide_[ ii ][ jj ] == 0 );
				//debug_assert( atom_has_been_updated_from_wide_[ ii ][ jj ] == 0 );
			}
		}
	}

	/////////////////////////
	// dimension the nblists, or on a subsquenct setup, remove the stale data from
	// the nblists
	bool const first_time_setup( nblist_.size() == 0 );
	debug_assert( first_time_setup || nblist_.size() == pose.size() );
	if ( first_time_setup ) { nblist_.resize( nres ); upper_nblist_.resize( nres ); intrares_upper_nblist_.resize( nres ); }
	if ( auto_update_ && first_time_setup ) wide_nblist_.resize( nres );
	for ( Size i=1; i<= nres; ++i ) {
		Size const natoms( pose.residue(i).natoms() );
		if ( first_time_setup ) { nblist_[i].resize( natoms ); upper_nblist_[i].resize( natoms ); intrares_upper_nblist_[i].resize( natoms ); }
		if ( auto_update_ && first_time_setup ) wide_nblist_[i].resize( natoms );
		for ( Size j=1; j<= natoms; ++j ) {
			nblist_[i][j].clear();
			upper_nblist_[i][j].clear();
			intrares_upper_nblist_[i][j].clear();
			if ( auto_update_ ) wide_nblist_[i][j].clear();
		}
	}

	/// Detect residue neighbors
	utility::vector1< bool > residue_mask = pose.conformation().get_residue_mask();
	core::conformation::PointGraphOP residue_point_graph( new core::conformation::PointGraph );
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *residue_point_graph );
	core::conformation::find_neighbors_restricted<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>(
		residue_point_graph,
		2 * pose::pose_max_nbr_radius( pose ) +
		XX_cutoff_ +
		( auto_update_ ? 2 * wide_nblist_extension_ : 0 ) ,
		residue_mask
	);
	//  find_neighbors(
	//   residue_point_graph,
	//   2 * pose::pose_max_nbr_radius( pose ) +
	//   XX_cutoff_ +
	//   ( auto_update_ ? 2 * wide_nblist_extension_ : 0 )
	//  );


	////////////////////
	// fill the nblist
	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & ires( pose.residue( i ) );
		int const imap( domain_map_(i) );

		if ( !residue_mask[i] ) continue;
		core::Real weight_func = pose.conformation().get_residue_weight(i,i);

		// are we defining intraresidue interactions? if so, add them
		// Do we only include intraresidue pair energies in the neighborlist if imap == 0?
		if ( etable_method.defines_intrares_energy( sfxn.weights() ) && imap == 0 ) {

			CountPairFunctionCOP count_pair( etable_method.get_intrares_countpair( ires, pose, sfxn) );

			for ( int ii=1, ii_end = ires.natoms(); ii<= ii_end; ++ii ) {
				conformation::Atom const & iatom( ires.atom(ii) );
				bool const iatom_is_hydrogen( ires.atom_is_hydrogen( ii ) );
				for ( int jj=ii+1, jj_end = ires.natoms(); jj<= jj_end; ++jj ) {

					// this code is duplicated -- a candidate for refactoring...
					conformation::Atom const & jatom( ires.atom(jj) );
					bool const jatom_is_hydrogen( ires.atom_is_hydrogen( jj ) );

					Real weight(1.0);
					Size path_dist( 0 );
					if ( count_pair->count( ii, jj, weight, path_dist ) ) {
						Real const dist_sq( iatom.xyz().distance_squared( jatom.xyz() ));
						Real const cutoff = atom_pair_cutoff( iatom_is_hydrogen, jatom_is_hydrogen );
						if ( dist_sq <= cutoff ) {
							declare_atoms_neighbors( id::AtomID( ii, i ), id::AtomID( jj, i ), path_dist, weight, weight_func );
						} // distance check

						if ( auto_update_ ) {
							Real const wide_cutoff
								( ( iatom_is_hydrogen && jatom_is_hydrogen ) ?
								HH_cutoff_wide_ : ( ( iatom_is_hydrogen || jatom_is_hydrogen ) ?
								XH_cutoff_wide_ : XX_cutoff_wide_ ) );
							if ( dist_sq <= wide_cutoff ) {
								wide_nblist_[i][ii].push_back( AtomNeighbor( i, jj, path_dist, weight ) );
								wide_nblist_[i][jj].push_back( AtomNeighbor( i, ii, path_dist, weight ) );
							}
						}

					} // count_pair check

				} // jj  = ii=1,ires.natoms()
			} // ii  = 1,ires.natoms()
		}

		Real const ireach( pose.residue( i ).nbr_radius() + sqrt_XX_cutoff_ + ( auto_update_ ? 2 * wide_nblist_extension_ : 0 ) );

		// Iterate across the neighbors of residue i
		//for ( graph::Graph::EdgeListConstIter
		//  iru  = pose.energies().energy_graph().get_node(i)->const_upper_edge_list_begin(),
		//  irue = pose.energies().energy_graph().get_node(i)->const_upper_edge_list_end();
		//  iru != irue; ++iru ) {
		for ( core::conformation::PointGraph::UpperEdgeListConstIter
				iru = residue_point_graph->get_vertex( i ).upper_edge_list_begin(),
				irue = residue_point_graph->get_vertex( i ).upper_edge_list_end();
				iru != irue; ++iru ) {

			//EnergyEdge const * edge( static_cast< EnergyEdge const *>(*iru));
			//EnergyEdge const * edge( utility::down_cast< EnergyEdge const * > (*iru) );

			Size const j( iru->upper_vertex() );
			Distance const ijsqrdist( iru->data().dsq() );
			Distance const ijreach( ireach + pose.residue( j ).nbr_radius() );
			/// ignore residues with either negative nbr_radii, or residue pairs that are sufficiently separated.
			if ( ijreach < 0 || ijsqrdist > ijreach * ijreach ) continue;

			//count_pair::CountPairFunctionCOP count_pair( edge->count_pair_function() );
			if ( imap == domain_map_(j) && imap != 0 ) continue;

			CountPairFunctionCOP count_pair( etable_method.get_count_pair_function( i, j, pose, sfxn) );

			conformation::Residue const & jres( pose.residue( j ) );
			core::Real weight_func2 = pose.conformation().get_residue_weight(i,j);

			for ( int ii=1, ii_end = ires.natoms(); ii<= ii_end; ++ii ) {
				conformation::Atom const & iatom( ires.atom(ii) );
				bool const iatom_is_hydrogen( ires.atom_is_hydrogen( ii ) );
				for ( int jj=1, jj_end = jres.natoms(); jj<= jj_end; ++jj ) {
					conformation::Atom const & jatom( jres.atom(jj) );
					bool const jatom_is_hydrogen( jres.atom_is_hydrogen( jj ) );

					Real weight(1.0);
					Size path_dist( 0 );
					if ( count_pair->count( ii, jj, weight, path_dist ) ) {
						Real const dist_sq( iatom.xyz().distance_squared( jatom.xyz() ));
						Real const cutoff( atom_pair_cutoff( iatom_is_hydrogen, jatom_is_hydrogen ));
						if ( dist_sq <= cutoff ) {
							declare_atoms_neighbors( id::AtomID( ii, i ), id::AtomID( jj, j ), path_dist, weight, weight_func2 );
						} // distance check

						if ( auto_update_ ) {
							Real const wide_cutoff
								( ( iatom_is_hydrogen && jatom_is_hydrogen ) ?
								HH_cutoff_wide_ : ( ( iatom_is_hydrogen || jatom_is_hydrogen ) ?
								XH_cutoff_wide_ : XX_cutoff_wide_ ) );
							if ( dist_sq <= wide_cutoff ) {
								wide_nblist_[i][ii].push_back( AtomNeighbor( j, jj, path_dist, weight, weight_func2 ) );
								wide_nblist_[j][jj].push_back( AtomNeighbor( i, ii, path_dist, weight, weight_func2 ) );
							}
						}

					}   // count_pair check
				}     // jj  = 1,jres.natoms()
			}       // ii  = 1,ires.natoms()
		}         // iru = { rsd nbrs of ires }
	}           // i   = 1,nres
	PROF_STOP ( basic::SETUP_NBLIST );
}

/// @brief If auto_update_, ensure that no atom in the pose has not moved too much
/// since the last time the neighborlist was updated.  The neighborlist
/// tracks the starting coords for all atoms, and then updates
template < class T_Etable >
void
NeighborList::prepare_for_scoring(
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	T_Etable const & etable_method
) const
{
	//std::cout << "prepare_for_scoring: " << auto_update_ << std::endl;
	if ( !auto_update_ ) return;

	++n_prepare_for_scorings_;

	atoms_to_update_.clear();
	Size const nres = pose.size();

	bool update_narrow = false; // true if any atom has moved > sqrt( move_tolerance_sqr_ ) from reference_coords_
	bool update_wide = false; // true if any atom has moved > sqrt( wide_move_tolerance_sqr_ ) from wide_reference_coords_

	debug_assert( atom_needs_update_from_wide_.size() == nres );
	for ( Size ii = 1; ii <= nres; ++ii ) {
		debug_assert( reference_coords_[ ii ].size() == pose.residue( ii ).natoms() );
		debug_assert( atom_needs_update_from_wide_[ ii ].size() ==  pose.residue( ii ).natoms() );
		for ( Size jj = 1; jj <= reference_coords_[ ii ].size(); ++jj ) {
			DistanceSquared dsqr_from_ref = reference_coords_[ ii ][ jj ].distance_squared( pose.residue( ii ).xyz( jj ));
			if ( dsqr_from_ref > move_tolerance_sqr_ ) {
				DistanceSquared dsqr_from_wide_ref = wide_reference_coords_[ ii ][ jj ].distance_squared( pose.residue( ii ).xyz( jj ));
				if ( dsqr_from_wide_ref > wide_move_tolerance_sqr_ ) {
					//std::cout << "Atom " << pose.residue( ii ).atom_name( jj ) << " on residue " << pose.residue( ii ).name() << " ";
					//std::cout << ii << " moved " << dsqr_from_wide_ref;
					//std::cout << " which is greater than wide move tolerance " << sqrt( wide_move_tolerance_sqr_ ) << std::endl;

					update_wide = true;
					break;
				}
				//std::cout << "Atom " << pose.residue( ii ).atom_name( jj ) << " on residue " << pose.residue( ii ).name() << " ";
				//std::cout << ii << " moved " << reference_coords_[ ii ][ jj ].distance( pose.residue( ii ).xyz( jj ));
				//std::cout << " which is greater than " << sqrt( move_tolerance_sqr_ ) << std::endl;
				update_narrow = true;
				atom_needs_update_from_wide_[ ii ][ jj ] = 1;
				atoms_to_update_.push_back( id::AtomID( jj, ii ) );
			}
		}
		if ( update_wide ) break;
	}
	if ( update_wide ) {
		/// Any atoms that we thought needed narrow-from-wide updates need to be reset.
		for ( Size ii = 1; ii <= atoms_to_update_.size(); ++ii ) {
			id::AtomID ii_atom( atoms_to_update_[ ii ] );
			atom_needs_update_from_wide_[ ii_atom.rsd() ][ ii_atom.atomno() ] = 0;
		}

		//std::cout << "Updating entire neighborlist!" << std::endl;
		++n_full_updates_;
		setup( pose, sfxn, etable_method );
	} else if ( update_narrow ) {
		//std::cout << "Updating neighborlist from wide neighborlist" << std::endl;
		++n_update_from_wide_;
		update_from_wide_nblist( pose );
	}
}

}
}

#endif
