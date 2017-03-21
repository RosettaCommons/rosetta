// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/interaction_graph/FlexbbIteractionGraph.cc
/// @brief  Declaration for flexible-backbone-packing interaction graph interface & base classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Unit headers
#include <protocols/flexpack/interaction_graph/FlexbbInteractionGraph.hh>

/// Package headers
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSets.hh>

/// Project headers
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

/// C++ headers
#include <iostream>

#include <protocols/flexpack/rotamer_set/FlexbbRotamerSet.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1A.hh>


namespace protocols {
namespace flexpack {
namespace interaction_graph {

using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer TR( "protocols.flexpack.interaction_graph" );

FlexbbNode::FlexbbNode(
	FlexbbInteractionGraph * owner,
	int node_id,
	int num_states
) :
	parent( owner, node_id, num_states ),
	num_aa_types_( get_flexbbig_owner()->get_num_aa_types() ),
	num_bb_( 1 ),
	num_states_for_bb_( 1, num_states ),
	state_offsets_for_bb_( 1, 0 ),
	state_info_( num_states + 1 ),
	one_body_energies_( num_states, 0.0 ),
	current_state_( 0 ),
	curr_state_one_body_energy_( 0.0 ),
	curr_state_total_energy_( 0.0 ),
	alternate_state_( 0 ),
	alternate_state_one_body_energy_( 0.0 ),
	alternate_state_total_energy_( 0.0 ),
	alternate_state_is_being_considered_( false )
{
	/// Default: for nodes without backbone moves, they will never have "set_num_states_per_backbone" called,
	/// so set their state_info_ to the detaul backbone #1.
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		state_info_[ ii ].set_bb( 1 );
	}
	//std::cout << "FlexbbNode ctor: num_aa_types_: " << num_aa_types_ << std::endl;
}

FlexbbNode::~FlexbbNode() {}

void
FlexbbNode::print() const {
	std::cout << "FlexbbNode " << get_node_index() << " #bb " << num_bb_
		<< " curr: " << current_state_ << " " << state_info_[ current_state_ ].get_aa_type()
		<< " " << state_info_[ current_state_ ].get_bb() << " energies: 1b: "
		<< curr_state_one_body_energy_ << " 2b: ";
	for ( Size ii = 1; ii <= curr_state_two_body_energies_.size(); ++ii ) {
		std::cout << curr_state_two_body_energies_[ ii ] << " ";
	}
	std::cout << std::endl;
}

/// @details Requires that num_aa_types_ is set first; num_aa_types_ is currently set in the node's ctor.
void
FlexbbNode::set_num_distinct_backbones( int nbbconfs ) {
	debug_assert( nbbconfs > 0 ); // cannot have fewer than 1 backbone conformation.
	debug_assert( num_aa_types_ != 0 );
	num_bb_ = nbbconfs;
	num_states_for_bb_.resize( nbbconfs, 0 );
	state_offsets_for_bb_.resize( nbbconfs, 0 );
	num_states_for_aa_for_bb_.dimension( num_aa_types_, nbbconfs );
	num_states_for_aa_for_bb_ = 0;
	state_offsets_for_aa_for_bb_.dimension( num_aa_types_, nbbconfs );
	state_offsets_for_aa_for_bb_ = 0;
	closest_state_on_alt_bb_.dimension( num_aa_types_, nbbconfs );
	closest_state_on_alt_bb_ = 0;
}

/// @details precondition: num_bb_ must have been previously set and must be equal to the number
/// of backbones implicitly identified in num_states_for_bb by its size.
void
FlexbbNode::set_num_states_per_backbone( utility::vector1< int > const & num_states_for_bb )
{
	debug_assert( num_states_for_bb_.size() == num_states_for_bb.size() );
	debug_assert( num_states_for_bb_.size() == state_offsets_for_bb_.size() );

	num_states_for_bb_ = num_states_for_bb;
	state_offsets_for_bb_[ 1 ] = 0;
	for ( int ii = 2; ii <= num_bb_; ++ii ) {
		state_offsets_for_bb_[ ii ] = state_offsets_for_bb_[ ii  - 1 ] + num_states_for_bb_[ ii - 1 ];
	}
	int bbindex = 1;
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		state_info_[ ii ].set_bb( bbindex );
		if ( bbindex < num_bb_ && state_offsets_for_bb_[ bbindex + 1 ] == ii ) ++bbindex;
	}

}

//int
//FlexbbNode::get_bb_for_state( int  state ) const {
// return state_info_[ state ].get_bb();
//}

int
FlexbbNode::get_num_states_for_bb( int  bbconf ) const {
	return num_states_for_bb_[ bbconf ];
}

int
FlexbbNode::get_state_offset_for_bb( int  bbconf ) const {
	return state_offsets_for_bb_[ bbconf ];
}

void
FlexbbNode::get_states_on_curr_bb(
	utility::vector1< Size > & rotlist,
	int rotamer_number_offset
) const
{
	Size const initial_size = rotlist.size();
	int currbb = state_info_[ current_state_ ].get_bb();
	if ( currbb == 0 ) currbb = 1;
	rotlist.reserve( initial_size + num_states_for_bb_[ currbb ] );
	for ( int ii = 1; ii <= num_states_for_bb_[ currbb ]; ++ii ) {
		rotlist.push_back( rotamer_number_offset + state_offsets_for_bb_[ currbb ] + ii );
	}
}

void
FlexbbNode::get_all_states(
	utility::vector1< Size > & state_list,
	int offset
) const
{
	Size const initial_size = state_list.size();
	state_list.reserve( initial_size + get_num_states() );
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		state_list.push_back( offset + ii );
	}
}

void
FlexbbNode::set_closest_states_on_other_bbs( ObjexxFCL::FArray2D_int const & closest_state_on_alt_bb )
{
	closest_state_on_alt_bb_ = closest_state_on_alt_bb;
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		for ( int jj = 1; jj <= num_bb_; ++jj ) {
			if ( state_info_[ ii ].get_bb() != jj  && closest_state_on_alt_bb_( jj, ii ) == 0 ) {
				TR << "Note state " << ii << " on FlexbbNode " << get_node_index() <<
					" has no close alternate state on backbone " << jj << std::endl;
				debug_assert( state_info_[ closest_state_on_alt_bb_( jj, ii ) ].get_bb() == jj );
			}
			if ( closest_state_on_alt_bb_( jj, ii ) > get_num_states() ) {
				std::cerr << "ALT STATE ON CLOSEST BB OUT OF RANGE: " << ii << " " << jj << " " << closest_state_on_alt_bb_( jj, ii ) << std::endl;
				utility_exit();
			}
		}
	}
}

/// @details Precondition: num bb have already been set.
void
FlexbbNode::set_amino_acid_types( utility::vector1< int > const & aatypes )
{
	int state_index_for_aa = 1;
	int last_aa = -1;
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		debug_assert( aatypes[ ii ] > 0 );
		if ( aatypes[ ii ] != last_aa ) {
			int temp_last_aa = last_aa;
			last_aa = aatypes[ ii ];
			if ( ii != 1 ) {
				state_offsets_for_aa_for_bb_( last_aa, state_info_[ ii ].get_bb() ) =
					state_offsets_for_aa_for_bb_( temp_last_aa, state_info_[ ii - 1 ].get_bb() )
					+ state_index_for_aa;
			}
			state_index_for_aa = 1;
		}
		state_info_[ ii ].set_aa_type( last_aa );
		state_info_[ ii ].set_state_ind_for_this_aa_type( state_index_for_aa );
		++num_states_for_aa_for_bb_( last_aa, state_info_[ ii ].get_bb() );
	}
	//for ( Size ii = 1; ii <= num_bb_; ++ii ) {
	// for ( Size jj = 1; jj <= num_aa_types_; ++jj ) {
	//  std::cout << get_node_index() << " " << ii << " " << jj << " " << num_states_for_aa_for_bb_( jj, ii ) << std::endl;
	// }
	//}
}

FArray1A_int
FlexbbNode::getNumStatesPerAAPerBB( int bb )
{
	return num_states_for_aa_for_bb_( 1, bb );
}


void
FlexbbNode::add_to_one_body_energies( FArray1< PackerEnergy > & energies )
{
	debug_assert( energies.size() == one_body_energies_.size() );
	for ( Size ii = 1; ii <= one_body_energies_.size(); ++ii ) {
		// if ( get_node_index() == 17 ) { std::cout << "FlexbbNode::add_to_one_body_energies " << ii << " " << energies( ii ) << " " << one_body_energies_[ ii ] << " " << one_body_energies_[ ii ] + energies( ii )  << std::endl; }
		one_body_energies_[ ii ] += energies( ii );
	}
}

void
FlexbbNode::add_to_one_body_energy( int state, PackerEnergy energy )
{
	//if ( get_node_index() == 17 ) { std::cout << "FlexbbNode::add_to_one_body_energy " << state << " " << energy  << " " << one_body_energies_[ state ] << " " << one_body_energies_[ state ] + energy << std::endl; }
	one_body_energies_[ state ] += energy;
}

/// @details Setter.  This function is not called anywhere and should be removed...
void
FlexbbNode::update_one_body_energy( int state, PackerEnergy energy )
{
	one_body_energies_[ state ] = energy;
}

void FlexbbNode::zero_one_body_energies()
{
	std::fill( one_body_energies_.begin(), one_body_energies_.end(), 0.0f );
}

bool FlexbbNode::state_unassigned() const {
	return current_state_ == 0;
}


FlexbbNode::PackerEnergy
FlexbbNode::get_one_body_energy( int state ) const
{
	return one_body_energies_[ state ];
}

void
FlexbbNode::prepare_for_simulated_annealing()
{
	update_internal_vectors();
}

int
FlexbbNode::get_current_state() const
{
	return current_state_;
}

int
FlexbbNode::get_backbone_for_current_state() const
{
	return state_info_[ current_state_ ].get_bb();
}

void
FlexbbNode::write_current_state_to_state_array( FArray1_int & nodes_states)
{
	nodes_states[ get_node_index() ] = current_state_;
}

/// @details Edges must be informed of the alternate state for both
/// of their nodes; this should happen before any edge is asked for
/// an alternate-state energy.
bool
FlexbbNode::inform_edges_of_alt_state_before_bbjump()
{
	//std::cout << "Node " << get_node_index() << " inform_edges_of_alt_state_before_bbjump() " << std::endl;
	alternate_state_is_being_considered_ = true;

	if ( current_state_ == 0 || alternate_state_ == 0 ) return false;

	get_flexbbig_owner()->increment_count_nodes_in_flexseg();

	for ( int ii = 1; ii  <= get_num_incident_edges(); ++ii ) {
		get_incident_flexbb_edge( ii )->set_alt_state(
			get_node_index(), alternate_state_, alternate_state_info_ );
	}

	return true;
}

/// @details In a bbmove, did the backbone actually move?  Should be called only
/// by the node who was initially contacted
//bool
//FlexbbNode::bb_move_actually_kept_original_bb() const {
// debug_assert( node_contacted_by_graph_about_bb_move_ );
// return resolved_considered_bb_move_;
//}

int
FlexbbNode::get_num_states_for_aa_type_for_bb(int aa_type, int bb) const
{
	return num_states_for_aa_for_bb_( aa_type, bb );
}

void
FlexbbNode::print_internal_energies()
{
	// stubbed out
}

/// @details Numerical noise creeps in to the total_energy as energies are added and subtracted
/// due to neighbors accepting rotamer substitutions.
void
FlexbbNode::update_internal_energy_sums()
{
	curr_state_total_energy_ = curr_state_one_body_energy_;
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		curr_state_total_energy_ += curr_state_two_body_energies_[ ii ];
	}
}

unsigned int
FlexbbNode::count_dynamic_memory() const
{
	unsigned int total = parent::count_dynamic_memory();
	total += num_states_for_bb_.size() * sizeof( int );
	total += state_offsets_for_bb_.size() * sizeof( int );
	total += num_states_for_aa_for_bb_.size() * sizeof( int );
	total += closest_state_on_alt_bb_.size() * sizeof( int );
	return total;
}


void FlexbbNode::update_internal_vectors()
{
	parent::update_edge_vector();

	curr_state_two_body_energies_.resize(      get_num_incident_edges() );
	alternate_state_two_body_energies_.resize( get_num_incident_edges() );
	edge_connects_flexsegmate_.resize(         get_num_incident_edges() );

	std::fill( curr_state_two_body_energies_.begin(), curr_state_two_body_energies_.end(), 0.0 );
	std::fill( alternate_state_two_body_energies_.begin(), alternate_state_two_body_energies_.end(), 0.0 );
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		edge_connects_flexsegmate_[ ii ] = get_incident_flexbb_edge( ii )->get_nodes_from_same_flexseg();
	}

	return;

}

void
FlexbbNode::inform_edges_considered_fixedbb_substition_uncommitted()
{
	alternate_state_is_being_considered_ = false;
	for ( int ii = 1; ii  <= get_num_incident_edges(); ++ii ) {
		get_incident_flexbb_edge( ii )->reset_alternate_states_for_uncommited_substitution();
	}

}


//bool
//FlexbbNode::energies_already_projected() const {
//std::cout << "energies already projected? " << get_node_index() << " " << projected_energies_for_bb_move_ << std::endl;
//bool temp = projected_energies_for_bb_move_;
//projected_energies_for_bb_move_ = true;
//return temp;
//}

/*
void
FlexbbNode::reset_bbmove_bookkeeping_data() {
resolved_considered_bb_move_ = true;
told_edges_alt_state_for_bb_move_ = false;
projected_energies_for_bb_move_ = false;
node_contacted_by_graph_about_bb_move_ = false;
}
*/

void
FlexbbNode::reset_all_rotamer_substitution_bookkeeping_data()
{
	alternate_state_is_being_considered_ = false;
}

void
FlexbbNode::partial_state_assignment( int new_state )
{
	alternate_state_is_being_considered_ = true;
	current_state_ = new_state;
	current_state_info_ = state_info_[ new_state ];
	alternate_state_ = new_state;
	alternate_state_info_ = state_info_[ new_state ];
	if ( new_state != 0 ) {
		alternate_state_one_body_energy_ = one_body_energies_[ new_state ];
	} else {
		alternate_state_one_body_energy_  = 0;
	}
	alternate_state_total_energy_ = alternate_state_one_body_energy_;
}

void
FlexbbNode::inform_incident_edges_about_partial_state_assignment()
{
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_flexbb_edge( ii )->acknowledge_partial_state_assignment(
			get_node_index(), current_state_, current_state_info_ );
	}
}


void FlexbbNode::copy_alternate_to_current()
{
	debug_assert( alternate_state_is_being_considered_ );

	current_state_ = alternate_state_;
	current_state_info_ = alternate_state_info_;
	curr_state_one_body_energy_ = alternate_state_one_body_energy_;
	curr_state_total_energy_ = alternate_state_total_energy_;

	std::copy( alternate_state_two_body_energies_.begin(),
		alternate_state_two_body_energies_.end(),
		curr_state_two_body_energies_.begin() );

	alternate_state_is_being_considered_ = false;

}

void FlexbbNode::have_edges_copy_alternate_to_current()
{
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_flexbb_edge( ii )->note_state_substitution_accepted();
	}
}

void FlexbbNode::have_edges_copy_alternate_to_current_following_flexbb_accept()
{
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		if ( ! edge_connects_flexsegmate()[ii] || get_index_of_adjacent_node( ii ) > get_node_index() ) {
			get_incident_flexbb_edge( ii )->note_state_substitution_accepted();
		}
	}
}

//// EDGE

FlexbbEdge::FlexbbEdge(
	FlexbbInteractionGraph * owner,
	int first_node_ind,
	int second_node_ind
) :
	parent( owner, first_node_ind, second_node_ind ),
	nodes_part_of_same_flexseg_( get_flexbbig_owner()->nodes_from_same_flexseg( first_node_ind, second_node_ind ) ),
	cur_energy_( 0.0 ),
	alt_energy_( 0.0 ),
	alt_e_up_to_date_( true ),
	nodes_considering_bb_move_( false )
{
	nodes_num_bb_[ 0 ] = get_flexbb_node( 0 )->get_num_distinct_backbones();
	nodes_num_bb_[ 1 ] = get_flexbb_node( 1 )->get_num_distinct_backbones();
	nodes_cur_state_[ 0 ] = 0;
	nodes_cur_state_[ 1 ] = 0;
	nodes_alt_state_[ 0 ] = 0;
	nodes_alt_state_[ 1 ] = 0;
}

FlexbbEdge::~FlexbbEdge() {}

/// @details during backbone-moving rotamer substitutions, the nodes on
/// flexible segments must inform their edges of their new states.
void
FlexbbEdge::set_alt_state(
	int node_index,
	int alternate_state,
	FlexbbSparseMatrixIndex const & alt_state_info
)
{
	alt_e_up_to_date_ = false;
	int node_changed = which_node( node_index );

	debug_assert( nodes_part_of_same_flexseg_ || ( nodes_cur_state_[ !node_changed ] == nodes_alt_state_[ !node_changed ] ));

	nodes_alt_state_[ node_changed ] = alternate_state;
	nodes_alt_info_[  node_changed ] = alt_state_info;

	nodes_considering_bb_move_ = nodes_part_of_same_flexseg_ &&
		nodes_alt_info_[ node_changed ].get_bb() != nodes_cur_info_[ node_changed ].get_bb();
}

void
FlexbbEdge::acknowledge_partial_state_assignment(
	int node_index,
	int new_state,
	FlexbbSparseMatrixIndex const & state_info
)
{
	alt_e_up_to_date_ = false;
	nodes_considering_bb_move_ = false;

	int node_changed = which_node( node_index );
	nodes_cur_state_[ node_changed ] = new_state;
	nodes_cur_info_[  node_changed ] = state_info;
	nodes_alt_state_[ node_changed ] = new_state;
	nodes_alt_info_[  node_changed ] = state_info;

	cur_energy_ = -1234;
	alt_energy_ = -1234;
}

void
FlexbbEdge::note_state_substitution_accepted()
{
	nodes_cur_state_[ 0 ] = nodes_alt_state_[ 0 ];
	nodes_cur_info_[  0 ] = nodes_alt_info_[  0 ];
	nodes_cur_state_[ 1 ] = nodes_alt_state_[ 1 ];
	nodes_cur_info_[  1 ] = nodes_alt_info_[  1 ];
	cur_energy_ = alt_energy_;
	nodes_considering_bb_move_ = false;
}

void
FlexbbEdge::reset_alternate_states_for_uncommited_substitution()
{
	nodes_alt_state_[ 0 ] = nodes_cur_state_[ 0 ];
	nodes_alt_info_[  0 ] = nodes_cur_info_[  0 ];
	nodes_alt_state_[ 1 ] = nodes_cur_state_[ 1 ];
	nodes_alt_info_[  1 ] = nodes_cur_info_[  1 ];
	nodes_considering_bb_move_ = false;
	alt_energy_ = cur_energy_;
	alt_e_up_to_date_ = true;
}


unsigned int
FlexbbEdge::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}

void
FlexbbEdge::copy_alternate_to_current()
{
	debug_assert( alt_e_up_to_date_ );

	nodes_cur_state_[ 0 ] = nodes_alt_state_[ 0 ];
	nodes_cur_state_[ 1 ] = nodes_alt_state_[ 1 ];
	nodes_cur_info_[ 0 ] = nodes_alt_info_[ 0 ];
	nodes_cur_info_[ 1 ] = nodes_alt_info_[ 1 ];
	cur_energy_ = alt_energy_;

	nodes_considering_bb_move_ = false;
}

void
FlexbbEdge::set_node_state_to_zero( int which_node )
{
	nodes_cur_state_[ which_node ] = nodes_alt_state_[ which_node ] = 0;
	nodes_cur_info_[ which_node ] = nodes_alt_info_[ which_node ] = FlexbbSparseMatrixIndex();
	cur_energy_ = alt_energy_ = 0;
	alt_e_up_to_date_ = true;
	nodes_considering_bb_move_ = false;

}

/// GRAPH

FlexbbInteractionGraph::~FlexbbInteractionGraph() {}
FlexbbInteractionGraph::FlexbbInteractionGraph(int num_nodes) :
	parent( num_nodes ),
	num_aa_types_( 0 ),
	total_energy_current_state_assignment_( 0.0 ),
	total_energy_alternate_state_assignment_( 0.0 ),
	node_considering_alt_state_( 0 ),
	flexseg_considering_alt_bb_( 0 ),
	num_flexible_segments_( 0 ),
	num_total_bb_( 0 ),
	flexseg_for_moltenres_( num_nodes, 0 ),
	enforce_bb_contiguity_( true ),
	last_sub_attempted_backbone_move_(false),
	last_considered_backbone_sub_valid_( true ),
	last_considered_backbone_sub_unresolved_( false ),
	last_considered_fixedbb_sub_unresolved_( false ),
	num_nodes_changing_state_( 0 ),
	num_commits_since_last_update_( 0 )
{}

void
FlexbbInteractionGraph::initialize( core::pack::rotamer_set::RotamerSetsBase const & rot_sets )
{
	using namespace rotamer_set;
	using namespace core::conformation;
	using namespace core;

	debug_assert( dynamic_cast< FlexbbRotamerSets const * > ( & rot_sets ) );
	rotamer_set::FlexbbRotamerSets const & flex_sets( static_cast< rotamer_set::FlexbbRotamerSets const & > ( rot_sets ) );

	/// Set Graph data:
	/// 1. flexseg and bb id mappings.
	set_num_flexsegs( flex_sets.nflexible_segments() );
	set_total_num_backbones( flex_sets.nbackbone_conformations());

	if ( num_flexible_segments_ != 0 ) { flexseg_bb_offset_[ 1 ] = 0; }
	for ( int ii = 1; ii <= num_flexible_segments_; ++ii ) {
		//flexseg_representative_[ ii ] = flex_sets.flexsegment_start_moltenresid( ii );
		/// List all the moltenresidues in this flexible segment.
		flexseg_members_[ ii ].resize( flex_sets.flexsegment_stop_moltenresid( ii ) - flex_sets.flexsegment_start_moltenresid( ii ) + 1 );
		for ( Size jj = 1; jj <= flexseg_members_[ ii ].size(); ++jj ) {
			flexseg_members_[ ii ][ jj ] = flex_sets.flexsegment_start_moltenresid( ii ) - 1 + jj;
		}
		Size const ii_nbb = flex_sets.nbbconfs_for_moltenres( flexseg_members_[ ii ][ 1 ] );
		num_bb_alternatives_for_flexseg_[ ii ] = ii_nbb;
		if ( ii > 1 ) {
			flexseg_bb_offset_[ ii ] = flexseg_bb_offset_[ ii - 1 ] + num_bb_alternatives_for_flexseg_[ ii - 1];
		}
	}
	int flexseg_id( 1 );
	for ( int ii = 1; ii <= num_total_bb_; ++ii ) {
		if ( flexseg_id < num_flexible_segments_ && ii > flexseg_bb_offset_[ flexseg_id + 1 ] ) {
			++flexseg_id;
		}
		flexseg_for_bb_[ ii ] = flexseg_id;
	}
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		flexseg_for_moltenres_[ ii ] = flex_sets.flexsegid_for_moltenres( ii );
	}

	// 2. determine max # of residue types, asking the RotamerSets how many rotamer type groups they have
	Size max_nresgroups = 0;
	for ( Size ii = 1; ii <= flex_sets.nmoltenres(); ++ii ) {
		for ( Size jj = 1; jj <= flex_sets.nbbconfs_for_moltenres( ii ); ++jj ) {
			Size jj_nresgroups =  flex_sets.rotset_for_moltenres( ii, jj )->get_n_residue_groups();
			//std::cout << ii << " " << jj << " jj_nrestypes: " << jj_nrestypes << " " << flex_sets.rotset_for_moltenres( ii, jj ).get() << std::endl;
			if ( jj_nresgroups > max_nresgroups ) max_nresgroups = jj_nresgroups;
		}
	}
	//std::cout << "MAX_NRESTYPES: " << max_nrestypes << std::endl;
	num_aa_types_ = max_nresgroups;

	/// 3. create nodes, inform them of their number of backbone conformations,
	/// and break down their rotamers by amino acid type.
	for ( Size ii = 1; ii <= flex_sets.nmoltenres(); ++ii ) {
		Size const ii_num_states = flex_sets.nrotamers_for_moltenres( ii );
		Size const ii_nbb =        flex_sets.nbbconfs_for_moltenres( ii );
		set_num_states_for_node( ii, ii_num_states ); // allocate the node.
		get_flexbb_node( ii )->set_num_distinct_backbones( ii_nbb );
		get_flexbb_node( ii )->set_num_states_per_backbone( flex_sets.num_states_per_backbone_for_moltenres( ii ) );

		// figure out which residue-type group each rotamer is a member of
		utility::vector1< int > aatype_for_state( ii_num_states, 0 );
		Size curr_resgroup = 1;
		Size count_for_resgroup = 1;
		Size const ii_nresgroups = flex_sets.rotset_for_moltenres( ii )->get_n_residue_groups();
		Size which_bb = 1;
		for ( Size jj = 1; jj <= ii_num_states; ++jj ) {
			if ( which_bb + 1 <= ii_nbb && jj == flex_sets.local_rotid_start_for_moltenres_in_bbconf( ii, which_bb + 1 ) + 1 ) {
				++which_bb;
				curr_resgroup = 1;
				count_for_resgroup = 1;
			}
			aatype_for_state[ jj ] = curr_resgroup;
			++count_for_resgroup;
			while ( count_for_resgroup > flex_sets.rotset_for_moltenres( ii, which_bb )->get_n_rotamers_for_residue_group( curr_resgroup ) ) {
				// increment curr_restype and skip over restypes with 0 rotamers
				++curr_resgroup;
				count_for_resgroup = 1;
				if ( curr_resgroup > ii_nresgroups ) break;
			}
		}
		set_aatypes_for_node( ii, aatype_for_state );
	}

	/// 4. Figure out for each rotamer on each residue with multiple backbone conformations
	/// what the closest rotamer on would be for each alternate backbone conformation.
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		Size const ii_nbb( get_flexbb_node( ii )->get_num_distinct_backbones()  );
		if ( ii_nbb == 1 ) continue;

		FArray2D_int best_match( ii_nbb, flex_sets.nrotamers_for_moltenres( ii ), 0 );
		for ( Size jj = 1; jj <= ii_nbb; ++jj ) {
			Size jjoffset = flex_sets.local_rotid_start_for_moltenres_in_bbconf( ii, jj );
			FlexbbRotamerSetCOP jj_flexset( flex_sets.rotset_for_moltenres( ii, jj ));
			for ( Size kk = 1; kk <= ii_nbb; ++kk ) {
				FlexbbRotamerSetCOP kk_flexset( flex_sets.rotset_for_moltenres(ii, kk ));
				if ( jj == kk ) continue; // skip diagonal.

				utility::vector1< Size > kk_rotamers_matching_ll;
				for ( int ll = 1; ll <= get_flexbb_node( ii )->get_num_states_for_bb( jj ); ++ll ) {
					if ( ll == 1 || jj_flexset->rotamer( ll )->name() != jj_flexset->rotamer( ll - 1 )->name() ) {
						kk_rotamers_matching_ll.clear();
						for ( Size mm = 1; mm <= kk_flexset->num_rotamers(); ++mm ) {
							if ( jj_flexset->rotamer( ll )->name() == kk_flexset->rotamer( mm )->name() ) {
								kk_rotamers_matching_ll.push_back( mm );
							}
						}
					}
					ResidueCOP llrotop = jj_flexset->rotamer( ll );
					Residue const & llrot( *llrotop );

					/// Now iterate across the rotamers with the same name and compute rms.
					Size best_match_index( 0 );
					Real best_rms( 12345 );
					for ( Size mm = 1; mm <= kk_rotamers_matching_ll.size(); ++mm ) {
						Size const mm_rotid = kk_rotamers_matching_ll[ mm ];
						ResidueCOP mmrotop = kk_flexset->rotamer( mm_rotid );
						Residue const & mmrot( *mmrotop );
						debug_assert( llrot.name() == mmrot.name() );

						Real mm_rms( 0.0 );
						for ( Size nn = 1; nn <= llrot.natoms(); ++nn ) {
							mm_rms += llrot.xyz( nn ).distance_squared( mmrot.xyz( nn ));
						}
						if ( mm == 1 || best_rms > mm_rms ) {
							best_match_index = mm_rotid;
							best_rms = mm_rms;
						}
					}
					best_match( kk, ll + jjoffset ) =
						flex_sets.local_rotid_start_for_moltenres_in_bbconf( ii, kk ) + best_match_index;

				}
			}
		}
		get_flexbb_node( ii )->set_closest_states_on_other_bbs( best_match );

	}
}

void
FlexbbInteractionGraph::set_num_flexsegs(int num_flexsegs)
{
	num_flexible_segments_ = num_flexsegs;
	//flexseg_representative_.resize( num_flexible_segments_ );
	flexseg_members_.resize( num_flexible_segments_ );
	//std::fill( flexseg_representative_.begin(), flexseg_representative_.end(), 0 );
	flexseg_bb_offset_.resize( num_flexible_segments_ );
	std::fill( flexseg_bb_offset_.begin(), flexseg_bb_offset_.end(), 0 );
	num_bb_alternatives_for_flexseg_.resize( num_flexible_segments_ );
	std::fill( num_bb_alternatives_for_flexseg_.begin(), num_bb_alternatives_for_flexseg_.end(), 0 );

}

void
FlexbbInteractionGraph::set_total_num_backbones( int num_backbones )
{
	num_total_bb_ = num_backbones;
	flexseg_for_bb_.resize( num_backbones );
	std::fill( flexseg_for_bb_.begin(), flexseg_for_bb_.end(), 0 );
	flexseg_bb_offset_.resize( num_backbones );
	std::fill( flexseg_bb_offset_.begin(), flexseg_bb_offset_.end(), 0 );
}


bool
FlexbbInteractionGraph::nodes_from_same_flexseg( int node1, int node2 ) const
{
	return flexseg_for_moltenres_[ node1 ] == flexseg_for_moltenres_[ node2 ];
}

/// @details Only call once. -- depricated!
//void
//FlexbbInteractionGraph::set_representitive_node_for_flexseg( int flexseg, int node_index)
//{
// debug_assert( flexseg_representative_[ flexseg ] == 0 );
// flexseg_representative_[ flexseg ] = node_index;
//}

/// @details Must include initial backbone as well all the alternative backbones.
void
FlexbbInteractionGraph::set_num_bb_for_node( int node, int numbb)
{
	get_flexbb_node( node )->set_num_distinct_backbones( numbb );
}

void
FlexbbInteractionGraph::set_num_states_per_backbone_for_node(
	int node, utility::vector1< int > const & states_per_bb
)
{
	get_flexbb_node( node )->set_num_states_per_backbone( states_per_bb );
}

int
FlexbbInteractionGraph::get_num_states_per_backbone_for_node( int node, int bb ) const
{
	return get_flexbb_node( node )->get_num_states_for_bb( bb );
}

int
FlexbbInteractionGraph::get_bb_for_state( int node, int state ) const
{
	return get_flexbb_node( node )->get_bb_for_state( state );
}

void
FlexbbInteractionGraph::set_aatypes_for_node(
	int node_ind,
	utility::vector1< int > const & aatypes
)
{
	get_flexbb_node( node_ind )->set_amino_acid_types( aatypes );
}

void
FlexbbInteractionGraph::set_closest_states_on_other_bbs(
	int node_ind,
	FArray2D_int const & closest_states
)
{
	get_flexbb_node( node_ind )->set_closest_states_on_other_bbs( closest_states );
}

/*
void
FlexbbInteractionGraph::set_edge_connecting_nodes_on_same_flexseg( int node1, int node2 )
{
FlexbbEdge * edge = find_flexbb_edge( node1, node2 );
if ( ! edge && std::abs( node1 - node2 ) == 1 ) {
add_edge( node1, node2 );
edge = find_flexbb_edge( node1, node2 );
} else if ( !edge ) {
return;
}

edge->set_nodes_from_same_flexseg( true );

}
*/

void
FlexbbInteractionGraph::get_accessible_states(
	Subsitution move_mode,
	utility::vector1< Size > & rotlist
) const
{
	rotlist.clear();
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		if ( move_mode == SC_ONLY ) {
			get_flexbb_node( ii )->get_states_on_curr_bb( rotlist, get_node_state_offset( ii ) );
		} else {
			get_flexbb_node( ii )->get_all_states( rotlist, get_node_state_offset( ii ) );
		}
	}
}

void
FlexbbInteractionGraph::get_backbone_list(
	utility::vector1< Size > & bblist
) const
{
	bblist.resize( num_total_bb_ );
	for ( Size ii = 1; ii <= bblist.size(); ++ii ) {
		bblist[ ii ] = ii;
	}
}

/// @brief Is the backbone conformation (in the global enumertion of backbone conformations) already
/// assigned to the network?  False if any residue on the flexible segment that this bbid corresponds to
/// is assigned state 0.
bool
FlexbbInteractionGraph::get_backbone_currently_assigned( int bbid ) const
{
	int flexsegid = flexseg_for_bb_[ bbid ];
	//int representative_node = flexseg_representative_[ flexsegid ];
	int representative_node = flexseg_members_[ flexsegid ][ 1 ];
	return bbid == get_flexbb_node( representative_node )->get_backbone_for_current_state();
}

/// @brief FlexbbNodes will ask: am I allowed to have a state that breaks the backbone?
/// There are brief periods when the backbone is "broken" as the graph assigns new states to
/// nodes on the same flexible segment.
bool
FlexbbInteractionGraph::get_enforce_bb_contiguity() const
{
	return enforce_bb_contiguity_;
}

/// @brief I don't remember what this is for... prev: increment_count_nodes_in_moving_fragment
void
FlexbbInteractionGraph::increment_count_nodes_in_flexseg()
{
	++num_nodes_changing_state_;
}

unsigned int
FlexbbInteractionGraph::count_dynamic_memory() const
{
	unsigned int total = parent::count_dynamic_memory();
	//total += flexseg_representative_.size() * sizeof( int );
	total += flexseg_members_.size() * sizeof( utility::vector1< Size > );
	for ( Size ii = 1; ii <= flexseg_members_.size(); ++ii ) {
		total += flexseg_members_[ ii ].size() * sizeof( Size );
	}
	total += num_bb_alternatives_for_flexseg_.size() * sizeof( int );
	total += flexseg_for_bb_.size() * sizeof( int );
	total += flexseg_bb_offset_.size() * sizeof( int );
	return total;
}

void
FlexbbInteractionGraph::set_total_energy_current_state_assignment( PackerEnergy setting )
{
	total_energy_current_state_assignment_ = setting;
	++num_commits_since_last_update_;

	if ( num_commits_since_last_update_ == COMMIT_LIMIT_BETWEEN_UPDATES ) {
		update_internal_energy_totals();
	}
}

void
FlexbbInteractionGraph::set_enforce_bb_contiguity( bool setting )
{
	enforce_bb_contiguity_ = setting;
}


void
FlexbbInteractionGraph::note_last_considered_substitution_resolved()
{
	last_considered_backbone_sub_unresolved_ = false;
}

void
FlexbbInteractionGraph::reset_node_in_moving_flexseg_count()
{
	/// The node contacted by the graph about the backbone move
	/// does not invoke increment on this count.   Therefore,
	/// the count for the number of changing nodes starts at one.
	num_nodes_changing_state_ = 1;
}

int
FlexbbInteractionGraph::get_num_nodes_changing_state() const
{
	return num_nodes_changing_state_;
}

void FlexbbInteractionGraph::update_internal_energy_totals()
{
	//PackerEnergy prev = total_energy_current_state_assignment_;
	total_energy_current_state_assignment_ = 0;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		total_energy_current_state_assignment_ += get_flexbb_node( ii )->curr_state_one_body_energy();
		get_flexbb_node( ii )->update_internal_energy_sums();
	}
	for ( std::list<EdgeBase*>::iterator
			iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		total_energy_current_state_assignment_ +=
			cast_flexbb_edge(*iter)->cur_energy();
	}
	num_commits_since_last_update_ = 0;
	//std::cout << "Updated current energy totals: " << total_energy_current_state_assignment_ << " delta: " << total_energy_current_state_assignment_ - prev << std::endl;
}


}
}
}

