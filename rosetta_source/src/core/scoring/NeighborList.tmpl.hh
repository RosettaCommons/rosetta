// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
//XRW_B_T1
//#include <core/coarse/Translator.fwd.hh>
//XRW_E_T1
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
//#include <core/kinematics/MinimizerMapBase.fwd.hh> // Is this necessary?

#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
//#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
//XRW_B_T1
//#include <core/scoring/etable/CoarseEtableEnergy.fwd.hh>
//XRW_E_T1
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <basic/MetricValue.fwd.hh>
// AUTO-REMOVED #include <basic/OStream.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>


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
	using namespace graph;

	PROF_START( basic::SETUP_NBLIST );

	Size const nres( pose.total_residue() );

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
		assert( reference_coords_.size() == nres );
		for ( Size ii = 1; ii <= nres; ++ii ) {
			assert( reference_coords_[ ii ].size() == pose.residue( ii ).natoms() );
			assert( wide_reference_coords_[ ii ].size() == pose.residue( ii ).natoms() );
			for ( Size jj = 1; jj <= reference_coords_[ ii ].size(); ++jj ) {
				reference_coords_[ ii ][ jj ] = pose.residue( ii ).xyz( jj );
				wide_reference_coords_[ ii ][ jj ] = pose.residue( ii ).xyz( jj );
				/// if these aren't zero, then the logic for updating atom nblists from the wide-nblists has failed
				assert( atom_needs_update_from_wide_[ ii ][ jj ] == 0 );
				//assert( atom_has_been_updated_from_wide_[ ii ][ jj ] == 0 );
			}
		}
	}

	/////////////////////////
	// dimension the nblists, or on a subsquenct setup, remove the stale data from
	// the nblists
	bool const first_time_setup( nblist_.size() == 0 );
	assert( first_time_setup || nblist_.size() == pose.total_residue() );
	if ( first_time_setup ) { nblist_.resize( nres ); upper_nblist_.resize( nres ); intrares_upper_nblist_.resize( nres ); }
	if ( auto_update_ && first_time_setup ) wide_nblist_.resize( nres );
	for ( Size i=1; i<= nres; ++i ) {
		Size const natoms( pose.residue(i).natoms() );
		if ( first_time_setup ) { nblist_[i].resize( natoms ); upper_nblist_[i].resize( natoms ); intrares_upper_nblist_[i].resize( natoms ); }
		if ( auto_update_ && first_time_setup ) wide_nblist_[i].resize( natoms );;
		for ( Size j=1; j<= natoms; ++j ) {
			nblist_[i][j].clear();
			upper_nblist_[i][j].clear();
			intrares_upper_nblist_[i][j].clear();
			if ( auto_update_ ) wide_nblist_[i][j].clear();
		}
	}

	/// Detect residue neighbors
	utility::vector1< bool > residue_mask = pose.conformation().get_residue_mask();
	core::conformation::PointGraphOP residue_point_graph = new core::conformation::PointGraph;
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *residue_point_graph );
	core::conformation::find_neighbors_restricted<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>(
		residue_point_graph,
		2 * pose::pose_max_nbr_radius( pose ) +
		XX_cutoff_ +
		( auto_update_ ? 2 * wide_nblist_extension_ : 0 ) ,
		residue_mask
	);
// 	find_neighbors(
// 		residue_point_graph,
// 		2 * pose::pose_max_nbr_radius( pose ) +
// 		XX_cutoff_ +
// 		( auto_update_ ? 2 * wide_nblist_extension_ : 0 )
// 	);


	////////////////////
	// fill the nblist
	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & ires( pose.residue( i ) );
		int const imap( domain_map_(i) );

		if (!residue_mask[i]) continue;
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
		//		iru  = pose.energies().energy_graph().get_node(i)->const_upper_edge_list_begin(),
		//		irue = pose.energies().energy_graph().get_node(i)->const_upper_edge_list_end();
		//		iru != irue; ++iru ) {
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
			core::Real weight_func = pose.conformation().get_residue_weight(i,j);

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
							declare_atoms_neighbors( id::AtomID( ii, i ), id::AtomID( jj, j ), path_dist, weight, weight_func );
						} // distance check

						if ( auto_update_ ) {
							Real const wide_cutoff
								( ( iatom_is_hydrogen && jatom_is_hydrogen ) ?
									HH_cutoff_wide_ : ( ( iatom_is_hydrogen || jatom_is_hydrogen ) ?
									XH_cutoff_wide_ : XX_cutoff_wide_ ) );
							if ( dist_sq <= wide_cutoff ) {
								wide_nblist_[i][ii].push_back( AtomNeighbor( j, jj, path_dist, weight, weight_func ) );
								wide_nblist_[j][jj].push_back( AtomNeighbor( i, ii, path_dist, weight, weight_func ) );
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
	Size const nres = pose.total_residue();

	bool update_narrow = false; // true if any atom has moved > sqrt( move_tolerance_sqr_ ) from reference_coords_
	bool update_wide = false; // true if any atom has moved > sqrt( wide_move_tolerance_sqr_ ) from wide_reference_coords_

	assert( atom_needs_update_from_wide_.size() == nres );
	for ( Size ii = 1; ii <= nres; ++ii ) {
		assert( reference_coords_[ ii ].size() == pose.residue( ii ).natoms() );
		assert( atom_needs_update_from_wide_[ ii ].size() ==  pose.residue( ii ).natoms() );
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
