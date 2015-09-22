// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/SymmMinimalistInteractionGraph.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <core/pack/interaction_graph/SymmMinimalistInteractionGraph.hh>

/// Debugging headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

#include <iostream>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1A.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>

namespace core {
namespace pack {
namespace interaction_graph {

static THREAD_LOCAL basic::Tracer T( "core.pack.interaction_graph.symm_symmin_ig", basic::t_error );


/// @brief main constructor, no default or copy constructors
SymmMinimalistNode::SymmMinimalistNode(
	InteractionGraphBase * owner,
	int node_id,
	int num_states
) :
	SymmOnTheFlyNode( owner, node_id, num_states ),
	current_state_( 0 ),
	curr_state_one_body_energy_( 0.0f ),
	curr_state_total_energy_( 0.0f ),
	alternate_state_( 0 ),
	alternate_state_one_body_energy_( 0 ),
	alternate_state_total_energy_( 0 ),
	alternate_state_is_being_considered_( false ),
	already_prepped_for_simA_( false )
{
}

SymmMinimalistNode::~SymmMinimalistNode()
{}

void
SymmMinimalistNode::prepare_for_simulated_annealing()
{
	if ( ! get_edge_vector_up_to_date() ) update_internal_vectors();
	already_prepped_for_simA_ = true;
	return;
}

void
SymmMinimalistNode::print() const
{

	T << "SymmMinimalistNode " << get_node_index() << " with " << get_num_states() << " states" << std::endl;
	T << "curr_state " << current_state_ << " ";
	T << "Curr One Body Energy: " << curr_state_one_body_energy_ << std::endl;
	T << "Curr Two Body Energies:";
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		T << " " << get_index_of_adjacent_node(ii) << ":" << curr_state_two_body_energies_[ ii ];
	}
	T << std::endl;

	if ( ! alternate_state_is_being_considered_ ) return;
	T << "Alt One Body Energy: " << alternate_state_one_body_energy_ << std::endl;
	T << "Alt Two Body Energies:";
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		T << " " << get_index_of_adjacent_node(ii) << ":" << alternate_state_two_body_energies_[ ii ];
	}
	T << std::endl  << "-----------------" << std::endl;


}

unsigned int
SymmMinimalistNode::count_static_memory() const
{
	return sizeof( SymmMinimalistNode );
}

unsigned int
SymmMinimalistNode::count_dynamic_memory() const
{
	unsigned int total_memory = SymmOnTheFlyNode::count_dynamic_memory();

	total_memory += neighbors_curr_state_.size() * sizeof( int );
	total_memory += curr_state_two_body_energies_.size() * sizeof( core::PackerEnergy );
	total_memory += alternate_state_two_body_energies_.size() * sizeof( core::PackerEnergy );

	return total_memory;
}


/// @brief puts the symminNode in the unassigned state
void
SymmMinimalistNode::assign_zero_state()
{
	current_state_ = 0;
	alternate_state_ = 0;
	alternate_state_is_being_considered_ = false;

	curr_state_one_body_energy_ = 0.0f;
	std::fill(
		curr_state_two_body_energies_.begin(),
		curr_state_two_body_energies_.end(),
		0.0f);
	curr_state_total_energy_ = 0.0f;

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_symmin_edge(ii)->
			acknowledge_state_zeroed( get_node_index() );
	}

	return;
}


//// @brief assigns a new state to the Node
void
SymmMinimalistNode::assign_state(int new_state)
{
	debug_assert( new_state >= 0 && new_state <= get_num_states());

	if ( new_state == 0 ) assign_zero_state();
	else {
		//T << "assign_state: node -  " << get_node_index() << " new state " << new_state << "...";
		current_state_ = new_state;
		curr_state_one_body_energy_ = get_one_body_energy( current_state_ );
		curr_state_total_energy_ = curr_state_one_body_energy_;
		alternate_state_is_being_considered_ = false;

		for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
			get_incident_symmin_edge(ii)->acknowledge_state_change(
				get_node_index(),
				current_state_,
				curr_state_two_body_energies_[ii]
			);

			curr_state_total_energy_ += curr_state_two_body_energies_[ ii ];
		}
		//T<< "..done" << std::endl;
	}

	//if ( debug ) {
	// get_on_the_fly_owner()->non_const_pose().replace_residue( get_rotamer(current_state_).seqpos(), get_rotamer( current_state_ ), false );
	// get_on_the_fly_owner()->score_function()( get_on_the_fly_owner()->non_const_pose() );
	//}
}


void
SymmMinimalistNode::partial_assign_state( int new_state )
{
	if ( new_state == 0 ) {
		assign_zero_state();
		return;
	}

	current_state_ = new_state;

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_symmin_edge(ii)->acknowledge_partial_state_change(
			get_node_index(),
			current_state_
		);
	}
	alternate_state_is_being_considered_ = false;
}


void SymmMinimalistNode::complete_state_assignment()
{
	if ( current_state_ == 0 ) return;

	curr_state_total_energy_ = curr_state_one_body_energy_ =
		get_one_body_energy( current_state_ );
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		curr_state_two_body_energies_[ ii ] =
			get_incident_symmin_edge( ii )->
			get_energy_following_partial_state_assignment();
		curr_state_total_energy_ += curr_state_two_body_energies_[ ii ];
	}
}


core::PackerEnergy
SymmMinimalistNode::project_deltaE_for_substitution
(
	int alternate_state,
	core::PackerEnergy & prev_node_energy
)
{
	alternate_state_is_being_considered_ = true;

	alternate_state_ = alternate_state;

	alternate_state_one_body_energy_ = get_one_body_energy( alternate_state );
	alternate_state_total_energy_ = alternate_state_one_body_energy_;
	prev_node_energy = curr_state_total_energy_;


	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		if ( neighbors_curr_state_[ ii ] != 0 ) {
			alternate_state_two_body_energies_[ ii ] = get_incident_symmin_edge( ii )->get_energy_for_alt_state( get_node_index() );
		} else {
			alternate_state_two_body_energies_[ ii ] = 0;
		}

		alternate_state_total_energy_ += alternate_state_two_body_energies_[ ii ];
	}

	// if ( debug && ! get_owner()->any_vertex_state_unassigned() ) {
	// debug_assert( get_rotamer(alternate_state_).seqpos() == get_rotamer(current_state_).seqpos() );
	//  get_on_the_fly_owner()->non_const_pose().replace_residue( get_rotamer(alternate_state_).seqpos(), get_rotamer( alternate_state_ ), false);
	//  Real score_after = get_on_the_fly_owner()->score_function()( get_on_the_fly_owner()->non_const_pose() );
	//  /// Now handled automatically.  get_on_the_fly_owner()->score_function().accumulate_residue_total_energies( get_on_the_fly_owner()->non_const_pose() );
	//  Real rep_after = get_on_the_fly_owner()->pose().energies().residue_total_energies( get_rotamer(alternate_state_).seqpos() )[ scoring::fa_rep ];
	//
	//  get_on_the_fly_owner()->non_const_pose().replace_residue( get_rotamer(current_state_).seqpos(), get_rotamer( current_state_ ), false);
	//  Real score_before = get_on_the_fly_owner()->score_function()( get_on_the_fly_owner()->non_const_pose() );
	//  /// Now handled automatically.  get_on_the_fly_owner()->score_function().accumulate_residue_total_energies( get_on_the_fly_owner()->non_const_pose() );
	//  Real rep_before = get_on_the_fly_owner()->pose().energies().residue_total_energies( get_rotamer(alternate_state_).seqpos() )[ scoring::fa_rep ];
	//
	//  Real actual_score_delta = score_after - score_before;
	//  Real projected_score_delta = alternate_state_total_energy_ - curr_state_total_energy_;
	//  Real delta_delta = actual_score_delta - projected_score_delta;
	//
	//  if ( (std::abs( delta_delta ) > 0.001 && std::abs( delta_delta / score_after ) > 10E-5) &&
	//   (rep_after < 4 && rep_before < 4) ) {
	//
	//   T << "Score before: " << score_before << " Score after " << score_after << " delta: " << actual_score_delta;
	//   T << " projected delta: " << projected_score_delta << " delta delta: " << delta_delta << " rep: " << rep_before << " " << rep_after <<  std::endl;
	//
	//
	//   /// LOOK AT CURRENT ENERGIES
	//   Size const seqpos( get_rotamer(alternate_state_).seqpos() );
	//   T << "Problem rotamer substitution at " << seqpos << ": from " << get_rotamer( current_state_).name() << " to " << get_rotamer(alternate_state_).name() << std::endl;
	//   T << "CURR One body energies: ";
	//   T << get_on_the_fly_owner()->score_function().weights().dot( get_on_the_fly_owner()->pose().energies().onebody_energies( get_rotamer(alternate_state_).seqpos() ) ) << std::endl;
	//   T << "internal one body energies: " << curr_state_one_body_energy_ << std::endl;
	//   T << "location: curr_state_one_body_energy_  " << & curr_state_one_body_energy_ << std::endl;
	//
	//   { //scope
	//   scoring::EnergyGraph const & energygraph = get_on_the_fly_owner()->pose().energies().energy_graph();
	//   for ( core::graph::Graph::EdgeListConstIter
	//     iter = energygraph.get_node( seqpos )->const_edge_list_begin(),
	//     iter_end = energygraph.get_node( seqpos)->const_edge_list_end();
	//     iter != iter_end; ++iter ) {
	//    bool corresponding_edge_found_in_ig( false );
	//    scoring::EnergyMap const tbemap( (static_cast< scoring::EnergyEdge const * > (*iter))->fill_energy_map() );
	//    Size const other_node_index = (*iter)->get_other_ind( seqpos );
	//    Real const real_energy = get_on_the_fly_owner()->score_function().weights().dot( tbemap );
	//    for ( Size ii = 1; ii <= (Size) get_num_incident_edges(); ++ii ) {
	//     if ( (Size) get_index_of_adjacent_node( ii ) != other_node_index ) continue;
	//     corresponding_edge_found_in_ig = true;
	//     if ( std::abs( real_energy - curr_state_two_body_energies_[ ii ]) > 0.001 ) {
	//      T << "Other residue: " << get_adjacent_symmin_node( ii )->get_rotamer( neighbors_curr_state_[ ii ]).name() << std::endl;
	//      T << "CURR Real score: edge to " << other_node_index << " energy: " << get_on_the_fly_owner()->score_function().weights().dot( tbemap ) << std::endl;
	//      T << "CURR Predicted score: edge to " << get_index_of_adjacent_node( ii ) << " energy: " << curr_state_two_body_energies_[ ii ] << std::endl;
	//      T << "CURR Real - Predicted: " << real_energy - curr_state_two_body_energies_[ ii ] << std::endl;
	//
	//      tbemap.show_nonzero( T );
	//      T << std::endl;
	//
	//      int const this_aa( curr_state_sparse_mat_info_.get_aa_type());
	//      int const other_aa( neighbors_curr_state_sparse_info_[ii].get_aa_type() );
	//      T << "Sparse matrix info: (this,other): " ;
	//      T << (int)  get_incident_symmin_edge( ii )->get_sparse_aa_neighbor_info()( this_aa, other_aa );
	//      T << " (other,this): ";
	//      T << (int) get_incident_symmin_edge( ii )->get_sparse_aa_neighbor_info()( other_aa, this_aa ) << std::endl;
	//
	//      core::PackerEnergy recomputed = compute_rotamer_pair_energy( ii, current_state_, neighbors_curr_state_[ ii ] );
	//      T << "Recomputed energy: " << recomputed << std::endl;
	//
	//      scoring::EnergyMap tbemap;
	//      get_on_the_fly_owner()->score_function().eval_ci_2b(
	//       get_on_the_fly_owner()->pose().residue( get_node_index() ),
	//       get_on_the_fly_owner()->pose().residue( other_node_index ),
	//       get_on_the_fly_owner()->pose(),
	//       tbemap );
	//      get_on_the_fly_owner()->score_function().eval_cd_2b(
	//       get_on_the_fly_owner()->pose().residue( get_node_index() ),
	//       get_on_the_fly_owner()->pose().residue( other_node_index ),
	//       get_on_the_fly_owner()->pose(),
	//       tbemap );
	//      T << "Rescored from pose: " << get_on_the_fly_owner()->score_function().weights().dot( tbemap ) << std::endl;
	//
	//      tbemap.zero();
	//      get_on_the_fly_owner()->score_function().eval_ci_2b(
	//       get_rotamer( current_state_ ),
	//       get_on_the_fly_owner()->pose().residue( other_node_index ),
	//       get_on_the_fly_owner()->pose(),
	//       tbemap );
	//      get_on_the_fly_owner()->score_function().eval_cd_2b(
	//       get_rotamer( current_state_ ),
	//       get_on_the_fly_owner()->pose().residue( other_node_index ),
	//       get_on_the_fly_owner()->pose(),
	//       tbemap );
	//      T << "Rescored combo: " << get_on_the_fly_owner()->score_function().weights().dot( tbemap ) << std::endl;
	//
	//
	//      tbemap.zero();
	//      get_on_the_fly_owner()->score_function().eval_ci_2b(
	//       get_adjacent_symmin_node( ii )->get_rotamer( neighbors_curr_state_[ ii ]),
	//       get_on_the_fly_owner()->pose().residue( get_node_index() ),
	//       get_on_the_fly_owner()->pose(),
	//       tbemap );
	//      get_on_the_fly_owner()->score_function().eval_cd_2b(
	//       get_adjacent_symmin_node( ii )->get_rotamer( neighbors_curr_state_[ ii ]),
	//       get_on_the_fly_owner()->pose().residue( get_node_index() ),
	//       get_on_the_fly_owner()->pose(),
	//       tbemap );
	//      T << "Rescored combo 2: " << get_on_the_fly_owner()->score_function().weights().dot( tbemap ) << std::endl;
	//
	//
	//      /// Check if order dependence is causing a bug -- res1 and res2 should not have to be ordered in the
	//      /// residue pair energy calls
	//      tbemap.zero();
	//      get_on_the_fly_owner()->score_function().eval_ci_2b(
	//       get_on_the_fly_owner()->pose().residue( other_node_index ),
	//       get_rotamer( current_state_ ),
	//       get_on_the_fly_owner()->pose(),
	//       tbemap );
	//      get_on_the_fly_owner()->score_function().eval_cd_2b(
	//       get_on_the_fly_owner()->pose().residue( other_node_index ),
	//       get_rotamer( current_state_ ),
	//       get_on_the_fly_owner()->pose(),
	//       tbemap );
	//      T << "Rescored combo swapped: " << get_on_the_fly_owner()->score_function().weights().dot( tbemap ) << std::endl;
	//
	//      /// Check if order dependence is causing a bug -- res1 and res2 should not have to be ordered in the
	//      /// residue pair energy calls
	//      tbemap.zero();
	//      get_on_the_fly_owner()->score_function().eval_ci_2b(
	//       get_on_the_fly_owner()->pose().residue( get_node_index() ),
	//       get_adjacent_symmin_node( ii )->get_rotamer( neighbors_curr_state_[ ii ]),
	//       get_on_the_fly_owner()->pose(),
	//       tbemap );
	//      get_on_the_fly_owner()->score_function().eval_cd_2b(
	//       get_on_the_fly_owner()->pose().residue( get_node_index() ),
	//       get_adjacent_symmin_node( ii )->get_rotamer( neighbors_curr_state_[ ii ]),
	//       get_on_the_fly_owner()->pose(),
	//       tbemap );
	//      T << "Rescored combo 2 swapped: " << get_on_the_fly_owner()->score_function().weights().dot( tbemap ) << std::endl;
	//
	//      //These references are useful in GDB if you need to debug.
	//      //
	//      //conformation::Residue const & res_in_pose = get_on_the_fly_owner()->pose().residue( get_node_index() );
	//      //conformation::Residue const & res_on_node = get_rotamer( current_state_ );
	//      //conformation::Residue const & other_res_in_pose = get_on_the_fly_owner()->pose().residue( other_node_index );
	//      //conformation::Residue const & other_res_on_node = get_adjacent_symmin_node( ii )->get_rotamer( neighbors_curr_state_[ ii ]);
	//      //SymmMinimalistNode * neighbor = get_adjacent_symmin_node( ii );
	//
	//      break;
	//     }
	//     if ( !corresponding_edge_found_in_ig ) {
	//      T << "Did not find edge in energy map to " << other_node_index << " with energy " << real_energy << " in the interaction graph!" << std::endl;
	//     }
	//    }
	//   }
	//   }// end scope
	//
	//   /// Look at interaction graph edges that are absent from the energy graph
	//   { //scope
	//   scoring::EnergyGraph const & energygraph = get_on_the_fly_owner()->pose().energies().energy_graph();
	//   for ( Size ii = 1; ii <= (Size) get_num_incident_edges(); ++ii ) {
	//    Size const other_node_index = (Size) get_index_of_adjacent_node( ii );
	//    if ( curr_state_two_body_energies_[ ii ] == 0 ) continue;
	//    bool found_similar_edge( false );
	//
	//    for ( core::graph::Graph::EdgeListConstIter
	//      iter = energygraph.get_node( seqpos )->const_edge_list_begin(),
	//      iter_end = energygraph.get_node( seqpos)->const_edge_list_end();
	//      iter != iter_end; ++iter ) {
	//     if ( other_node_index != (Size) (*iter)->get_other_ind( seqpos ) ) continue;
	//     found_similar_edge = true;
	//
	//     //scoring::EnergyMap const & tbemap( (static_cast< scoring::EnergyEdge const * > (*iter))->energy_map() );
	//     //Real const real_energy = get_on_the_fly_owner()->score_function().weights().dot( tbemap );
	//
	//    }
	//    if ( ! found_similar_edge ) {
	//     T << "Edge in lmig CUR to node " << other_node_index << " with energy: " << curr_state_two_body_energies_[ ii ] << " absent from energy graph!" << std::endl;
	//    }
	//   }
	//   } // end scope
	//
	//   /// Place the alternate rotamer on the pose and rescore.
	//   get_on_the_fly_owner()->non_const_pose().replace_residue( get_rotamer(alternate_state_).seqpos(), get_rotamer( alternate_state_ ), false);
	//   get_on_the_fly_owner()->score_function()( get_on_the_fly_owner()->non_const_pose() );
	//
	//   T << "ALT One body energies: ";
	//   T << get_on_the_fly_owner()->score_function().weights().dot( get_on_the_fly_owner()->pose().energies().onebody_energies( get_rotamer(alternate_state_).seqpos() ) ) << std::endl;
	//   T << "internal one body energies: " << alternate_state_one_body_energy_ << std::endl;
	//   T << "location: alternate_state_one_body_energy_  " << & alternate_state_one_body_energy_ << std::endl;
	//
	//
	//   { //scope
	//   scoring::EnergyGraph const & energygraph = get_on_the_fly_owner()->pose().energies().energy_graph();
	//   for ( core::graph::Graph::EdgeListConstIter
	//     iter = energygraph.get_node( seqpos )->const_edge_list_begin(),
	//     iter_end = energygraph.get_node( seqpos)->const_edge_list_end();
	//     iter != iter_end; ++iter ) {
	//    bool corresponding_edge_found_in_ig( false );
	//    scoring::EnergyMap const tbemap( (static_cast< scoring::EnergyEdge const * > (*iter))->fill_energy_map() );
	//    Size const other_node_index = (*iter)->get_other_ind( seqpos );
	//    Real const real_energy = get_on_the_fly_owner()->score_function().weights().dot( tbemap );
	//    for ( Size ii = 1; ii <= (Size) get_num_incident_edges(); ++ii ) {
	//     if ( (Size) get_index_of_adjacent_node( ii ) != other_node_index ) continue;
	//     corresponding_edge_found_in_ig = true;
	//     if ( std::abs( real_energy - alternate_state_two_body_energies_[ ii ]) > 0.001 ) {
	//      T << "ALT Real score: edge to " << other_node_index << " energy: " << get_on_the_fly_owner()->score_function().weights().dot( tbemap ) << std::endl;
	//      T << "ALT Predicted score: edge to " << get_index_of_adjacent_node( ii ) << " energy: " << alternate_state_two_body_energies_[ ii ] << std::endl;
	//      T << "ALT Real - Predicted: " << real_energy - alternate_state_two_body_energies_[ ii ] << std::endl;
	//      tbemap.show_nonzero( T );
	//      T << std::endl;
	//
	//      int const this_aa( alt_state_sparse_mat_info_.get_aa_type());
	//      int const other_aa( neighbors_curr_state_sparse_info_[ii].get_aa_type() );
	//      T << "Sparse matrix info: (this,other): " ;
	//      T << (int) get_incident_symmin_edge( ii )->get_sparse_aa_neighbor_info()( this_aa, other_aa );
	//      T << " (other,this): ";
	//      T << (int) get_incident_symmin_edge( ii )->get_sparse_aa_neighbor_info()( other_aa, this_aa ) << std::endl;
	//
	//      core::PackerEnergy recomputed = compute_rotamer_pair_energy( ii, alternate_state_, neighbors_curr_state_[ ii ] );
	//      T << "Recomputed energy: " << recomputed << std::endl;
	//
	//      break;
	//     }
	//     if ( !corresponding_edge_found_in_ig ) {
	//      T << "Did not find edge in energy map to " << other_node_index << " with energy " << real_energy << " in the interaction graph!" << std::endl;
	//     }
	//
	//    }
	//   }
	//   }// end scope
	//   /// Look at interaction graph edges that are absent from the energy graph
	//   { //scope
	//   scoring::EnergyGraph const & energygraph = get_on_the_fly_owner()->pose().energies().energy_graph();
	//   for ( Size ii = 1; ii <= (Size) get_num_incident_edges(); ++ii ) {
	//    Size const other_node_index = (Size) get_index_of_adjacent_node( ii );
	//    if ( alternate_state_two_body_energies_[ ii ] == 0 ) continue;
	//    bool found_similar_edge( false );
	//
	//    for ( core::graph::Graph::EdgeListConstIter
	//      iter = energygraph.get_node( seqpos )->const_edge_list_begin(),
	//      iter_end = energygraph.get_node( seqpos)->const_edge_list_end();
	//      iter != iter_end; ++iter ) {
	//     if ( other_node_index != (Size) (*iter)->get_other_ind( seqpos ) ) continue;
	//     found_similar_edge = true;
	//
	//     //scoring::EnergyMap const & tbemap( (static_cast< scoring::EnergyEdge const * > (*iter))->energy_map() );
	//     //Real const real_energy = get_on_the_fly_owner()->score_function().weights().dot( tbemap );
	//
	//    }
	//    if ( ! found_similar_edge ) {
	//     T << "Edge in lmig ALT to node " << other_node_index << " with energy: " << alternate_state_two_body_energies_[ ii ] << " absent from energy graph!" << std::endl;
	//    }
	//   }
	//   } // end scope
	//
	//
	//   get_on_the_fly_owner()->non_const_pose().replace_residue( get_rotamer(alternate_state_).seqpos(), get_rotamer( current_state_ ), false);
	//   get_on_the_fly_owner()->score_function()( get_on_the_fly_owner()->non_const_pose() );
	//  }
	// } // end debug

	return alternate_state_total_energy_ - curr_state_total_energy_;

}


/// @brief commits the last substitution that was considered by this Node
void
SymmMinimalistNode::commit_considered_substitution()
{
	debug_assert( alternate_state_is_being_considered_ );

	current_state_ = alternate_state_;
	curr_state_one_body_energy_ = alternate_state_one_body_energy_;
	curr_state_total_energy_ = alternate_state_total_energy_;

	//copies from [1] to end
	//utility::vector1< core::PackerEnergy >::iterator alt_position1 = alternate_state_two_body_energies_.begin();
	//utility::vector1< core::PackerEnergy >::iterator curr_position1 = curr_state_two_body_energies_.begin();

	std::copy( alternate_state_two_body_energies_.begin(),
		alternate_state_two_body_energies_.end(),
		curr_state_two_body_energies_.begin() );

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_symmin_edge(ii)->acknowledge_substitution(
			get_node_index(),
			curr_state_two_body_energies_[ ii ],
			current_state_
		);
	}

	alternate_state_is_being_considered_ = false;

	// if ( debug ) {
	//  get_on_the_fly_owner()->non_const_pose().replace_residue( get_rotamer(current_state_).seqpos(), get_rotamer( current_state_ ), false );
	//  get_on_the_fly_owner()->score_function()( get_on_the_fly_owner()->non_const_pose() );
	// }

	return;
}

void
SymmMinimalistNode::acknowledge_last_substititon_not_committed()
{
	alternate_state_is_being_considered_ = false;
}

core::PackerEnergy
SymmMinimalistNode::compute_pair_energy_for_current_state(
	int edge_making_energy_request
)
{
	return compute_rotamer_pair_energy(
		edge_making_energy_request,
		current_state_,
		neighbors_curr_state_[ edge_making_energy_request ]
	);
}

core::PackerEnergy
SymmMinimalistNode::compute_pair_energy_for_alternate_state(
	int edge_making_energy_request
)
{
	return compute_rotamer_pair_energy(
		edge_making_energy_request,
		alternate_state_,
		neighbors_curr_state_[ edge_making_energy_request ] );
}


void
SymmMinimalistNode::acknowledge_neighbors_partial_state_substitution(
	int edge_to_altered_neighbor,
	int other_node_new_state
)
{
	curr_state_total_energy_ = 0;
	curr_state_two_body_energies_[ edge_to_altered_neighbor ] = 0;
	neighbors_curr_state_[ edge_to_altered_neighbor ] = other_node_new_state;
}

void
SymmMinimalistNode::print_internal_energies() const
{
	T << "curr_state " << current_state_ << " ";
	T << "curr_state_one_body_energy_ ";
	T << curr_state_one_body_energy_ << " ";
	T << "curr_state_total_energy_" << curr_state_total_energy_ << " ";
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		T << "(" << get_index_of_adjacent_node(ii) << ":" << curr_state_two_body_energies_[ ii ] << ") ";
	}
	T << std::endl;
}


void
SymmMinimalistNode::update_internal_energy_sums()
{
	debug_assert( get_edge_vector_up_to_date() );
	curr_state_total_energy_ = 0;
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		curr_state_total_energy_ += get_incident_symmin_edge(ii)->get_current_two_body_energy();
	}
	curr_state_total_energy_ += curr_state_one_body_energy_;
	return;
}

void SymmMinimalistNode::update_internal_vectors()
{
	NodeBase::update_edge_vector();
	neighbors_curr_state_.resize( get_num_incident_edges());

	curr_state_two_body_energies_.resize( get_num_incident_edges());
	alternate_state_two_body_energies_.resize( get_num_incident_edges());
	return;
}


//-----------------------------------------------------------------//


SymmMinimalistEdge::SymmMinimalistEdge(
	InteractionGraphBase* owner,
	int first_node_ind,
	int second_node_ind
):
	SymmOnTheFlyEdge( owner, first_node_ind, second_node_ind),
	curr_state_energy_( 0.0f ),
	partial_state_assignment_( false ),
	preped_for_sim_annealing_( false )
{
}

SymmMinimalistEdge::~SymmMinimalistEdge() {}

core::PackerEnergy SymmMinimalistEdge::get_two_body_energy( int const , int const ) const
{
	throw utility::excn::EXCN_Msg_Exception( "Method unimplemented: SymmMinimalistEdge::get_two_body_energy" );
	return 0.0;
}

void
SymmMinimalistEdge::declare_energies_final()
{}

void
SymmMinimalistEdge::prepare_for_simulated_annealing()
{
	if ( preped_for_sim_annealing_ ) return;
	preped_for_sim_annealing_ = true;
}

unsigned int
SymmMinimalistEdge::count_static_memory() const
{
	return sizeof( SymmMinimalistEdge );
}


unsigned int
SymmMinimalistEdge::count_dynamic_memory() const
{
	unsigned int total_memory = SymmOnTheFlyEdge::count_dynamic_memory();
	return total_memory;
}

/// @details  DANGER: this will not update the cached energies on the nodes this edge is incident upon.
void
SymmMinimalistEdge::set_edge_weight( Real weight )
{
	edge_weight( weight );
}

core::PackerEnergy
SymmMinimalistEdge::get_current_two_body_energy() const
{
	return curr_state_energy_;
}


void
SymmMinimalistEdge::acknowledge_state_change(
	int node_ind,
	int new_state,
	core::PackerEnergy & new_energy
)
{
	int node_substituted =  ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = get_symmin_node( 0 )->
		compute_pair_energy_for_current_state(
		get_edges_position_in_nodes_edge_vector( 0 ) );

	new_energy = curr_state_energy_;

	get_symmin_node( node_not_substituted )->acknowledge_neighbors_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		new_state
	);
}


void
SymmMinimalistEdge::acknowledge_state_zeroed( int node_ind )
{
	int node_substituted = ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = 0;

	get_symmin_node( node_not_substituted )->acknowledge_neighbors_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		0
	);
	return;
}


void SymmMinimalistEdge::acknowledge_partial_state_change(
	int node_ind,
	int new_state
)
{
	int node_substituted =  ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = 0;

	get_symmin_node( node_not_substituted )->acknowledge_neighbors_partial_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		new_state );
	partial_state_assignment_ = true;
}


core::PackerEnergy
SymmMinimalistEdge::get_energy_following_partial_state_assignment()
{
	if ( partial_state_assignment_
			&& get_symmin_node(0)->get_current_state() != 0
			&& get_symmin_node(1)->get_current_state() != 0 ) {

		curr_state_energy_ = get_symmin_node( 0 )->
			compute_pair_energy_for_current_state(
			get_edges_position_in_nodes_edge_vector( 0 ) );
		partial_state_assignment_ = false;
	}
	return curr_state_energy_;
}

core::PackerEnergy
SymmMinimalistEdge::get_energy_for_alt_state(
	int changing_node_index
)
{
	int const node_changing = changing_node_index == get_node_index( 0 ) ? 0 : 1;
	alt_state_energy_ = get_symmin_node( node_changing )->
		compute_pair_energy_for_alternate_state(
		get_edges_position_in_nodes_edge_vector( node_changing ));

	return alt_state_energy_;
}

int SymmMinimalistEdge::get_two_body_table_size() const
{
	return 1;
}

void
SymmMinimalistEdge::print_current_energy() const
{
	T << "SymmMinimalistEdge: " << get_node_index( 0 ) << "/" << get_node_index( 1 );
	T << " energy= " << curr_state_energy_ << std::endl;
}


//-------------------------------------------------------------------//

SymmMinimalistInteractionGraph::SymmMinimalistInteractionGraph(
	int numNodes
) :
	SymmOnTheFlyInteractionGraph( numNodes ),
	first_time_prepping_for_simA_( true ),
	num_commits_since_last_update_( 0 ),
	total_energy_current_state_assignment_( 0.0 ),
	total_energy_alternate_state_assignment_( 0.0 ),
	node_considering_alt_state_( 0 ),
	have_not_committed_last_substitution_( false )
{
}


SymmMinimalistInteractionGraph::~SymmMinimalistInteractionGraph() {}

void
SymmMinimalistInteractionGraph::blanket_assign_state_0()
{
	have_not_committed_last_substitution_ = false;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_symmin_node( ii )->assign_zero_state();
	}
	total_energy_current_state_assignment_ = 0;
}


core::PackerEnergy
SymmMinimalistInteractionGraph::set_state_for_node(int node_ind, int new_state)
{
	have_not_committed_last_substitution_ = false;
	get_symmin_node( node_ind )->assign_state( new_state );
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}


core::PackerEnergy
SymmMinimalistInteractionGraph::set_network_state(
	ObjexxFCL::FArray1_int & node_states
)
{
	have_not_committed_last_substitution_ = false;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_symmin_node( ii )->partial_assign_state( node_states( ii ) );
	}
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_symmin_node( ii )->complete_state_assignment();
	}
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}


void
SymmMinimalistInteractionGraph::consider_substitution(
	int node_ind,
	int new_state,
	core::PackerEnergy & delta_energy,
	core::PackerEnergy & prev_energy_for_node
)
{
	if ( have_not_committed_last_substitution_ ) {
		get_symmin_node( node_considering_alt_state_ )->acknowledge_last_substititon_not_committed();
	}

	node_considering_alt_state_ = node_ind;

	delta_energy = get_symmin_node( node_ind )->project_deltaE_for_substitution( new_state, prev_energy_for_node );

	total_energy_alternate_state_assignment_ = total_energy_current_state_assignment_ + delta_energy;
	have_not_committed_last_substitution_ = true;
}

core::PackerEnergy
SymmMinimalistInteractionGraph::commit_considered_substitution()
{
	have_not_committed_last_substitution_ = false;
	get_symmin_node( node_considering_alt_state_ )->commit_considered_substitution();

	total_energy_current_state_assignment_ = total_energy_alternate_state_assignment_;

	++num_commits_since_last_update_;
	if ( num_commits_since_last_update_ == COMMIT_LIMIT_BETWEEN_UPDATES ) {
		update_internal_energy_totals();
	}

	return total_energy_alternate_state_assignment_;
}


core::PackerEnergy
SymmMinimalistInteractionGraph::get_energy_current_state_assignment()
{
	//T << "Num rotamer pair energy calculations performed: " << SymmMinimalistNode::num_rpe_calcs << std::endl;
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}

/// @brief O(1) total energy report.  Protected read access for derived classes.
core::PackerEnergy
SymmMinimalistInteractionGraph::get_energy_PD_current_state_assignment()
{
	return total_energy_current_state_assignment_;
}

int
SymmMinimalistInteractionGraph::get_edge_memory_usage() const
{
	int sum = 0;
	for ( std::list< EdgeBase* >::const_iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		sum += ((SymmMinimalistEdge*) *iter)->get_two_body_table_size();
	}
	return sum;
}

void
SymmMinimalistInteractionGraph::print_current_state_assignment() const
{
	T << "State Assignment: " << std::endl;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		T << "Node " << ii << " state " << get_symmin_node(ii)->get_current_state() << std::endl;
		get_symmin_node(ii)->print();
	}

	for ( std::list< EdgeBase* >::const_iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		((SymmMinimalistEdge*) (*iter))->print_current_energy();
	}
	T << "Energy: " << total_energy_current_state_assignment_ << std::endl;
}


void
SymmMinimalistInteractionGraph::set_errorfull_deltaE_threshold( core::PackerEnergy ) {}

core::PackerEnergy
SymmMinimalistInteractionGraph::get_energy_sum_for_vertex_group( int group_id )
{
	core::PackerEnergy esum = 0;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		if ( get_vertex_member_of_energy_sum_group( ii, group_id ) ) {
			esum += get_symmin_node( ii )->get_one_body_energy_current_state();
		}
	}

	for ( std::list< EdgeBase* >::iterator edge_iter = get_edge_list_begin();
			edge_iter != get_edge_list_end(); ++edge_iter ) {
		int first_node_ind = (*edge_iter)->get_first_node_ind();
		int second_node_ind = (*edge_iter)->get_second_node_ind();

		if ( get_vertex_member_of_energy_sum_group( first_node_ind, group_id )
				&& get_vertex_member_of_energy_sum_group( second_node_ind, group_id ) ) {
			esum += ((SymmMinimalistEdge*) (*edge_iter))->get_current_two_body_energy();
		}
	}

	return esum;
}

void
SymmMinimalistInteractionGraph::prepare_for_simulated_annealing()
{
	if ( first_time_prepping_for_simA_ ) {
		first_time_prepping_for_simA_ = false;
	}
	InteractionGraphBase::prepare_for_simulated_annealing();
}

unsigned int
SymmMinimalistInteractionGraph::count_static_memory() const
{
	return sizeof( SymmMinimalistInteractionGraph );
}

unsigned int
SymmMinimalistInteractionGraph::count_dynamic_memory() const
{
	unsigned int total_memory = SymmOnTheFlyInteractionGraph::count_dynamic_memory();
	return total_memory;
}


NodeBase*
SymmMinimalistInteractionGraph::create_new_node( int node_index, int num_states )
{
	return new SymmMinimalistNode( this, node_index, num_states );
}


EdgeBase*
SymmMinimalistInteractionGraph::create_new_edge( int index1, int index2 )
{
	return new SymmMinimalistEdge( this, index1, index2 );
}

void
SymmMinimalistInteractionGraph::update_internal_energy_totals()
{
	total_energy_current_state_assignment_ = 0;

	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		total_energy_current_state_assignment_ += get_symmin_node( ii )->
			get_one_body_energy_current_state();
	}

	for ( std::list<EdgeBase*>::iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		total_energy_current_state_assignment_ +=
			((SymmMinimalistEdge*) *iter)->get_current_two_body_energy();
	}

	num_commits_since_last_update_ = 0;
	return;
}

} // namespace interaction_graph
} // namespace pack
} // namespace core
