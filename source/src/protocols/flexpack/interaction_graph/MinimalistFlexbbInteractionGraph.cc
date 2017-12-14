// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/interaction_graph/MinimalistFlexbbInteractionGraph.cc
/// @brief  Class implementation for minimimalist on-the-fly RPE calculating FlexbbInteractionGraph
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/flexpack/interaction_graph/MinimalistFlexbbInteractionGraph.hh>

/// Project headers

/// C++ headers
#include <iostream>

#include <utility/exit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace flexpack {
namespace interaction_graph {

MinimalistFlexbbNode::MinimalistFlexbbNode(
	MinimalistFlexbbInteractionGraph * owner,
	int node_id,
	int num_states
) :
	parent( owner, node_id, num_states )
{}

MinimalistFlexbbNode::~MinimalistFlexbbNode() = default;

/// Virtual functions from NodeBase
void
MinimalistFlexbbNode::assign_zero_state()
{
	set_current_state( 0 );
	set_curr_state_one_body_energy( 0.0 );
	set_curr_state_total_energy( 0.0 );
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		set_curr_state_two_body_energies( ii , 0.0 );
		get_incident_minimalistflexbb_edge( ii )->acknowledge_state_zeroed( get_node_index() );
	}
	reset_all_rotamer_substitution_bookkeeping_data();
}

void
MinimalistFlexbbNode::prepare_for_simulated_annealing()
{
	parent::prepare_for_simulated_annealing();
}

/*
void
MinimalistFlexbbNode::add_to_one_body_energy( int state, PackerEnergy energy );

void
MinimalistFlexbbNode::update_one_body_energy( int state, PackerEnergy energy);

void
MinimalistFlexbbNode::zero_one_body_energies();
*/

void
MinimalistFlexbbNode::print() const
{
	parent::print();
}

//bool
//MinimalistFlexbbNode::state_unassigned() const
//{
// return current_state() == 0;
//}


MinimalistFlexbbNode::PackerEnergy
MinimalistFlexbbNode::project_deltaE_for_substitution(
	int alternate_state,
	PackerEnergy & prev_energy_for_node
)
{
	debug_assert( get_bb_for_state( current_state() ) == 0 || get_bb_for_state( alternate_state ) == get_bb_for_state( current_state() ));
	set_considering_alternate_state();
	set_alternate_state( alternate_state );
	set_alternate_state_one_body_energy( one_body_energies()[ alternate_state ]);

	set_alternate_state_total_energy( alternate_state_one_body_energy() );
	prev_energy_for_node = curr_state_total_energy();

	/// for debugging -- make sure that curr_state_total_energy
	/// has not accumulated too much numerical noise.
	//Real curr_state_actual_total = curr_state_one_body_energy();
	//PackerEnergy two_body_energy_sum( 0.0 );

	//std::cout << "\npdE: " << get_node_index() << " alt state: " << alternate_state << " amino acid " << rotamer( alternate_state ).aa() ;
	//std::cout << " coord: " << rotamer( alternate_state ).xyz( rotamer( alternate_state ).nheavyatoms() ).x() <<
	//     " " << rotamer( alternate_state ).xyz( rotamer( alternate_state ).nheavyatoms() ).y() <<
	//     " " << rotamer( alternate_state ).xyz( rotamer( alternate_state ).nheavyatoms() ).z() << std::endl;
	//std::cout << "alt state one body energy: " << one_body_energies()[ alternate_state ] << " vs curr: " << one_body_energies()[ current_state() ] << std::endl;

	for ( int ii = 1; ii <= get_num_edges_to_smaller_indexed_nodes(); ++ii ) {
		get_incident_flexbb_edge( ii )->set_alt_state(
			get_node_index(), alternate_state, alternate_state_info() );
		set_alternate_state_two_body_energies(
			ii, get_incident_otfflexbb_edge( ii )->compute_samebbconf_alternate_state_energy_second_node());
		inc_alternate_state_total_energy( alternate_state_two_body_energies( ii ));

		//std::cout << "\n" << get_index_of_adjacent_node( ii ) << " " << alternate_state_two_body_energies( ii );
		/*
		two_body_energy_sum += alternate_state_two_body_energies( ii );
		if ( get_node_index() == 17 && alternate_state == 63 ) {
		std::cout << "Two-Body Energy Predicted: " << get_node_index() << " " << get_index_of_adjacent_node( ii ) << " " << alternate_state_two_body_energies( ii ) << std::endl;
		}
		curr_state_actual_total += curr_state_two_body_energies()[ ii ];
		if ( std::abs( curr_state_two_body_energies()[ ii ] - get_incident_otfflexbb_edge( ii )->cur_energy()) > 1e-6 ) {
		std::cout << "CURRENT STATE ENERGY DISCREPANCY: " << get_node_index() << " " << get_index_of_adjacent_node( ii ) <<
		" local: " << curr_state_two_body_energies()[ ii ] <<
		" on edge: " << get_incident_otfflexbb_edge( ii )->cur_energy() <<
		" diff: " << curr_state_two_body_energies()[ ii ] - get_incident_otfflexbb_edge( ii )->cur_energy() <<
		std::endl;
		}
		*/
	}

	for ( int ii = get_num_edges_to_smaller_indexed_nodes() + 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_flexbb_edge( ii )->set_alt_state(
			get_node_index(), alternate_state, alternate_state_info() );
		set_alternate_state_two_body_energies(
			ii, get_incident_otfflexbb_edge( ii )->compute_samebbconf_alternate_state_energy_first_node() );
		inc_alternate_state_total_energy( alternate_state_two_body_energies( ii ));
		/*
		if ( get_node_index() == 17 && alternate_state == 63 ) {
		std::cout << "Two-Body Energy Predicted: " << get_node_index() << " " << get_index_of_adjacent_node( ii ) << " " << alternate_state_two_body_energies( ii ) << std::endl;
		}
		//std::cout << "\n" << get_index_of_adjacent_node( ii ) << " " << alternate_state_two_body_energies( ii );
		two_body_energy_sum += alternate_state_two_body_energies( ii );
		curr_state_actual_total += curr_state_two_body_energies()[ ii ];
		if ( std::abs( curr_state_two_body_energies()[ ii ] - get_incident_otfflexbb_edge( ii )->cur_energy()) > 1e-6 ) {
		std::cout << "CURRENT STATE ENERGY DISCREPANCY: " << get_node_index() << " " << get_index_of_adjacent_node( ii ) <<
		" local: " << curr_state_two_body_energies()[ ii ] <<
		" on edge: " << get_incident_otfflexbb_edge( ii )->cur_energy() <<
		" diff: " << curr_state_two_body_energies()[ ii ] - get_incident_otfflexbb_edge( ii )->cur_energy() <<
		std::endl;
		}
		*/
	}
	//std::cout << "\n Two body energy total: " << two_body_energy_sum << std::endl;
	//std::cout << "\n Alt state total energy: " << alternate_state_total_energy() << " vs curr: " << curr_state_total_energy() << std::endl;

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	get_otfflexbbig_owner()->debug_note_considered_substitution( rotamer( alternate_state ), alternate_state - state_offsets_for_bb()[ alternate_state_info().get_bb() ] );
	//get_otfflexbbig_owner()->debug_note_considered_substitution( rotamer( alternate_state ), alternate_state );
#endif

	/*if ( std::abs( curr_state_actual_total - curr_state_total_energy()) > 1e-6 ) {
	std::cout << "Inaccurate current energy: " << get_node_index() << " " << alternate_state << " "
	<< curr_state_actual_total << " " << curr_state_total_energy() << std::endl;
	}
	if ( get_node_index() == 17 && alternate_state == 63 ) {
	std::cout << "Alternate_state_total_energy: " << alternate_state_total_energy() << " " << two_body_energy_sum
	<< " " << alternate_state_one_body_energy() + two_body_energy_sum
	<< " " << alternate_state_total_energy() - curr_state_total_energy()
	<< " " << alternate_state_total_energy() - curr_state_actual_total
	<< std::endl;
	}*/

	return alternate_state_total_energy() - curr_state_total_energy();

}

/*
MinimalistFlexbbNode::PackerEnergy
MinimalistFlexbbNode::project_deltaE_with_backbone_move(
int alternate_state,
PackerEnergy & prev_energy_for_flexseg,
bool & valid_motion
)
{
//std::cout << "project_deltaE_with_backbone_move: " << get_node_index() << " " << alternate_state << " " << curr_bb() << std::endl;
prev_energy_for_flexseg = 0;
register_contacted_node_for_bb_jump();

if ( get_bb_for_state( alternate_state ) == get_bb_for_state( current_state() ) ) {
//std::cout << "project_deltaE_with_backbone_move: same backbone " << get_node_index() << " " << alternate_state << " " << curr_bb() << std::endl;
return project_deltaE_for_substitution( alternate_state, prev_energy_for_flexseg );
}

set_alternate_state( alternate_state );

int const alt_bb = alternate_state_info().get_bb();

/// FlexbbNode base class call -- visit all nodes in this flexseg and their adjacent edges
/// in a DFS traversal.  If any node in the flexseg has state 0, quit.
if ( ! prepare_for_bb_jump( alt_bb ) ) {
//std::cout << "INVALID MOTION" << std::endl;
valid_motion = false; return 0.0;
}

PackerEnergy alt_energy_for_flexseg = get_altE_for_bb_move( prev_energy_for_flexseg );
//std::cout << "Finished flexseg energy calc: " << alt_energy_for_flexseg << " " << prev_energy_for_flexseg << std::endl;
return alt_energy_for_flexseg - prev_energy_for_flexseg;

}


MinimalistFlexbbNode::PackerEnergy
MinimalistFlexbbNode::project_deltaE_for_backbone_move(
int alt_bb,
PackerEnergy & prev_energy_for_flexseg,
bool & valid_motion
)
{
//std::cout << "Considering bb move: " << get_node_index() << " " << alt_bb << " " << curr_bb() << std::endl;
if ( current_state() == 0 || alt_bb == curr_bb() ) { valid_motion = false; return 0.0; }
set_alternate_state( closest_state_on_alt_bb()( alt_bb, current_state() ) );
return project_deltaE_with_backbone_move( alternate_state(), prev_energy_for_flexseg, valid_motion );
}
*/

bool
MinimalistFlexbbNode::prepare_for_altbb_move_to_state( int alt_state )
{
	if ( alt_state == 0 ) return false;
	if ( alt_state > get_num_states() ) {
		std::cerr << "CANNOT MOVE TO ALT_STATE " << alt_state << " WHEN THERE ARE ONLY " << get_num_states() << " STATES" << std::endl;
		utility_exit();
		return false;
	}

	set_considering_alternate_state();
	set_alternate_state( alt_state );
	set_alternate_state_one_body_energy( one_body_energies()[ alt_state ]);

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	get_otfflexbbig_owner()->debug_note_considered_substitution( rotamer( alt_state ), alt_state - state_offsets_for_bb()[ alternate_state_info().get_bb() ] );
#endif

	return inform_edges_of_alt_state_before_bbjump();
}

bool
MinimalistFlexbbNode::prepare_for_altbb_move_to_closest_state( int alt_bb )
{
	if ( current_state() == 0 ) return false;
	if ( current_state() > get_num_states() ) { std::cerr << "Current state out-of-range: " << get_num_states() << " " << current_state() << std::endl; return false; }
	if ( alt_bb == 0 || alt_bb > get_num_distinct_backbones() ) { std::cerr << "ALTERNATE BACKBONE OUT-OF-RANGE: " << get_num_distinct_backbones() << " " << alt_bb << std::endl; return false; }
	if ( closest_state_on_alt_bb()( alt_bb, current_state() ) == 0 ) return false;

	if ( closest_state_on_alt_bb()( alt_bb, current_state() ) > get_num_states() ) {
		std::cerr << "CANNOT2 MOVE TO ALT_STATE " << closest_state_on_alt_bb()( alt_bb, current_state() ) << " WHEN THERE ARE ONLY " << get_num_states() << " STATES; " << alt_bb << " " << current_state() << " " << curr_bb() << std::endl;
		utility_exit();
		return false;
	}
	return prepare_for_altbb_move_to_state( closest_state_on_alt_bb()( alt_bb, current_state() ) );
}

/// @details This function call will be parallelized with openMP.
/// The calculations performed by different nodes may be performed in
/// any order and in parallel without any data collision.  After this function
/// completes, a second pass across the nodes in the graph must be
/// performed to compute the node alternate_state_total_energy --
/// but this second traversal is fast.
MinimalistFlexbbNode::PackerEnergy
MinimalistFlexbbNode::get_frag_energy_for_alt_bb_state()
{
	PackerEnergy alt_frag_etotal( alternate_state_one_body_energy() );
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		if ( edge_connects_flexsegmate()[ ii ] ) {
			if ( get_index_of_adjacent_node( ii ) >  get_node_index() ) {
				set_alternate_state_two_body_energies( ii,
					get_incident_minimalistflexbb_edge( ii )->get_alt_stateE());

				alt_frag_etotal  += alternate_state_two_body_energies()[ ii ];
			}
		} else {
			set_alternate_state_two_body_energies( ii, get_incident_minimalistflexbb_edge( ii )->get_alt_stateE() );
			alt_frag_etotal += alternate_state_two_body_energies()[ ii ];
		}
	}
	return alt_frag_etotal;
}

/// @details This method could be parallelized with openMP, but probably would
/// not be worth parallelizing.  It is the second step of computing the deltaE
/// for a backbone move.  The alternate_state_total_energy is computed here, as
/// well as the curr_frag_total_energy which is returned.
MinimalistFlexbbNode::PackerEnergy
MinimalistFlexbbNode::get_frag_energy_for_curr_bb_state_and_finalize_alt_energy_total()
{
	PackerEnergy curr_frag_etotal( curr_state_one_body_energy() );
	set_alternate_state_total_energy( alternate_state_one_body_energy() );

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		if ( edge_connects_flexsegmate()[ ii ] ) {
			if ( get_index_of_adjacent_node( ii ) <  get_node_index() ) {
				/// this energy has already been computed
				set_alternate_state_two_body_energies( ii,
					get_incident_minimalistflexbb_edge( ii )->get_alt_stateE());
			} else {
				curr_frag_etotal += curr_state_two_body_energies()[ ii ];
			}
			inc_alternate_state_total_energy( alternate_state_two_body_energies()[ ii ] );
		} else {
			inc_alternate_state_total_energy( alternate_state_two_body_energies()[ ii ] );
			curr_frag_etotal += curr_state_two_body_energies()[ ii ];
		}
	}

	return curr_frag_etotal;
}


void
MinimalistFlexbbNode::commit_considered_substitution()
{
	debug_assert( alternate_state_is_being_considered() );

	copy_alternate_to_current();
	have_edges_copy_alternate_to_current();
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_minimalistflexbb_edge( ii )->acknowledge_substitution( get_node_index() );
	}
}

void
MinimalistFlexbbNode::acknowledge_neighbors_substitution(
	int which_edge,
	PackerEnergy alternate_twobody_energy
)
{
	inc_curr_state_total_energy( alternate_twobody_energy - curr_state_two_body_energies()[ which_edge ] );
	set_curr_state_two_body_energies( which_edge, alternate_twobody_energy );
	/// DEBUG:
	//PackerEnergy true_curr_state_total_energy = curr_state_one_body_energy();
	//for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
	// true_curr_state_total_energy += curr_state_two_body_energies()[ ii ];
	//}
	//debug_assert( std::abs( (true_curr_state_total_energy - curr_state_total_energy())/std::max( PackerEnergy(1.0), std::abs(curr_state_total_energy())) ) < 1e-5 );
}


void MinimalistFlexbbNode::commit_considered_substitution( ObjexxFCL::FArray1_int & state_on_node )
{
	commit_considered_substitution();
	state_on_node( get_node_index() ) = current_state();
	return;
}

//// @details Cannot be parallelized since it's updating current_total_energy_ on other nodes.
void
MinimalistFlexbbNode::commit_alt_bb_substitution( ObjexxFCL::FArray1_int & state_on_node )
{
	copy_alternate_to_current();
	have_edges_copy_alternate_to_current_following_flexbb_accept();
	state_on_node( get_node_index() ) = alternate_state();
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_minimalistflexbb_edge( ii )->acknowledge_substitution( get_node_index() );
	}
}

void
MinimalistFlexbbNode::resolve_uncommitted_substitution()
{
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_minimalistflexbb_edge( ii )->reset_alternate_states_for_uncommited_substitution();
	}
	set_alternate_state( current_state() );
}

MinimalistFlexbbNode::PackerEnergy
MinimalistFlexbbNode::assign_state( int new_state )
{
	PackerEnergy prev_total = curr_state_total_energy();
	partially_assign_state( new_state );
	complete_partial_state_assignment();
	return curr_state_total_energy() - prev_total;
}

void
MinimalistFlexbbNode::partially_assign_state( int new_state )
{
	//set_considering_alternate_state();
	partial_state_assignment( new_state );
	inform_incident_edges_about_partial_state_assignment();

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	get_otfflexbbig_owner()->debug_note_considered_substitution( rotamer( new_state ), alternate_state() - state_offsets_for_bb()[ alternate_state_info().get_bb() ] );
	//get_otfflexbbig_owner()->debug_note_considered_substitution( rotamer( new_state ), alternate_state() );
#endif

}

void
MinimalistFlexbbNode::complete_partial_state_assignment()
{
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		set_alternate_state_two_body_energies( ii,
			get_incident_minimalistflexbb_edge( ii )->get_alt_stateE());
		inc_alternate_state_total_energy( alternate_state_two_body_energies()[ ii ]);
	}
	have_edges_copy_alternate_to_current();
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_minimalistflexbb_edge( ii )->acknowledge_substitution( get_node_index() );
	}
	copy_alternate_to_current();
}


unsigned int
MinimalistFlexbbNode::count_static_memory() const
{
	return sizeof( MinimalistFlexbbNode );
}

unsigned int
MinimalistFlexbbNode::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}

/*
MinimalistFlexbbNode::PackerEnergy
MinimalistFlexbbNode::get_altE_for_bb_move( PackerEnergy & curr_frag_etotal )
{
if ( energies_already_projected() ) return 0.0;

//std::cout << "get_altE_for_bb_move " << get_node_index() << std::endl;

PackerEnergy alt_frag_etotal( 0.0 );

set_alternate_state_one_body_energy( one_body_energies()[ alternate_state() ]);
//if ( get_node_index() == 17 && alternate_state() == 63 ) {
// std::cout << "SETTING ALT STATE 11 (63) ON NODE 17: " << alternate_state_one_body_energy() << std::endl;
//}

set_alternate_state_total_energy( alternate_state_one_body_energy() );

alt_frag_etotal += alternate_state_one_body_energy();
curr_frag_etotal += curr_state_one_body_energy();

for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
if ( edge_connects_flexsegmate()[ ii ] ) {
alt_frag_etotal +=
get_adjacent_minimalistflexbb_node( ii )->get_altE_for_bb_move( curr_frag_etotal );
//std::cerr << get_node_index() << " alt_frag_etotal( " << ii << ", " << alt_frag_etotal  << ") " << std::endl;
//std::cerr << get_node_index() << " curr_frag_etotal( " << ii << ", " << curr_frag_etotal << ") " << std::endl;

set_alternate_state_two_body_energies( ii,
get_incident_minimalistflexbb_edge( ii )->get_alt_stateE());
inc_alternate_state_total_energy( alternate_state_two_body_energies()[ ii ]);

if ( count_energy_to_node_in_my_fragtotalE( ii ) ) {

alt_frag_etotal  += alternate_state_two_body_energies()[ ii ];
curr_frag_etotal += curr_state_two_body_energies()[ ii ];

}
} else {

//if ( ii < get_num_edges_to_smaller_indexed_nodes() ) {
// set_alternate_state_two_body_energies(
//  ii, get_incident_otfflexbb_edge( ii )->compute_samebbconf_alternate_state_energy_second_node());
//} else {
// set_alternate_state_two_body_energies(
//  ii, get_incident_otfflexbb_edge( ii )->compute_samebbconf_alternate_state_energy_first_node() );
//}

set_alternate_state_two_body_energies( ii, get_incident_minimalistflexbb_edge( ii )->get_alt_stateE() );

alt_frag_etotal += alternate_state_two_body_energies()[ ii ];
inc_alternate_state_total_energy( alternate_state_two_body_energies()[ ii ] );
curr_frag_etotal += curr_state_two_body_energies()[ii];

}
}

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
get_otfflexbbig_owner()->debug_note_considered_substitution( rotamer( alternate_state() ), alternate_state() - state_offsets_for_bb()[ alternate_state_info().get_bb() ] );
//get_otfflexbbig_owner()->debug_note_considered_substitution( rotamer( alternate_state() ), alternate_state() );
#endif

return alt_frag_etotal;
}
*/

/// EDGE

MinimalistFlexbbEdge::MinimalistFlexbbEdge(
	MinimalistFlexbbInteractionGraph * owner,
	int node1,
	int node2
) :
	parent( owner, node1, node2 )
{}

MinimalistFlexbbEdge::~MinimalistFlexbbEdge() = default;

void MinimalistFlexbbEdge::declare_energies_final()
{}

void MinimalistFlexbbEdge::prepare_for_simulated_annealing()
{
	//std::cout << "MinFlexbbEdge::prepare_for_simulated_annealing" << std::endl;
	parent::prepare_for_simulated_annealing();
}

void MinimalistFlexbbEdge::set_edge_weight( Real weight )
{
	/// This should also go and reweight energies already stored on the edges...
	/// assert that SA hasn't yet begun for now.
	debug_assert( nodes_cur_state( 0 ) == 0 && nodes_cur_state( 1 ) == 0 );
	/// EdgeBase class protected setter.
	edge_weight( weight );
}


MinimalistFlexbbEdge::PackerEnergy
MinimalistFlexbbEdge::get_alt_stateE()
{
	//std::cout << "get_alt_stateE: ";
	if ( alt_e_up_to_date() ) {
		//std::cout << "... done" << std::endl;
		return alt_energy();
	}

	compute_altbbconf_alternate_state_energy();
	//std::cout << ".. computed: " << alt_energy() << std::endl;
	return alt_energy();
}

void
MinimalistFlexbbEdge::acknowledge_substitution(
	int node_that_changed
)
{
	int node_not_changing( ! which_node( node_that_changed ));
	/// copy_alternate_to_current(); -- this call is handled by flexbb_node::have_edges_copy_alternate_to_current();
	otfedge_note_substitution_accepted();
	get_minimalistflexbb_node( node_not_changing )->acknowledge_neighbors_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_changing ),
		cur_energy() );
}

void
MinimalistFlexbbEdge::acknowledge_state_zeroed( int node_index )
{
	parent::set_node_state_to_zero( which_node( node_index ));
}

unsigned int
MinimalistFlexbbEdge::count_static_memory() const
{
	return sizeof( MinimalistFlexbbEdge );
}

unsigned int
MinimalistFlexbbEdge::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}


/// GRAPH

MinimalistFlexbbInteractionGraph::MinimalistFlexbbInteractionGraph( int num_nodes ) :
	parent( num_nodes )
{}

MinimalistFlexbbInteractionGraph::~MinimalistFlexbbInteractionGraph() = default;

void
MinimalistFlexbbInteractionGraph::initialize( core::pack::rotamer_set::RotamerSetsBase const & rot_sets )
{
	parent::initialize( rot_sets );
}

MinimalistFlexbbInteractionGraph::PackerEnergy
MinimalistFlexbbInteractionGraph::get_one_body_energy_for_node_state(
	int node,
	int state
)
{
	return get_minimalistflexbb_node( node )->get_one_body_energy( state );
}

void
MinimalistFlexbbInteractionGraph::blanket_assign_state_0()
{
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_minimalistflexbb_node( ii )->assign_zero_state();
	}
	note_last_considered_substitution_resolved();
	set_total_energy_current_state_assignment( 0 );
	set_node_considering_alt_state( 0 );
}

MinimalistFlexbbInteractionGraph::PackerEnergy
MinimalistFlexbbInteractionGraph::set_state_for_node(int node_ind, int new_state)
{
	if ( last_considered_substitution_unresolved() ) {
		resolve_uncommitted_substitution();
	}

	set_node_considering_alt_state( 0 );

	PackerEnergy deltaE = get_minimalistflexbb_node( node_ind )->assign_state( new_state );
	set_total_energy_current_state_assignment( total_energy_current_state_assignment() + deltaE );
	update_internal_energy_totals();
	//std::cout << "Accept: " << total_energy_current_state_assignment() << std::endl;
	return total_energy_current_state_assignment();
}

MinimalistFlexbbInteractionGraph::PackerEnergy
MinimalistFlexbbInteractionGraph::set_network_state( ObjexxFCL::FArray1_int & node_states )
{
	if ( last_considered_substitution_unresolved() ) {
		resolve_uncommitted_substitution();
	}

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	PackerEnergy last_total_energy = total_energy_current_state_assignment();
#endif

	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_minimalistflexbb_node( ii )->partially_assign_state( node_states( ii ) );
	}
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_minimalistflexbb_node( ii )->complete_partial_state_assignment();
	}
	update_internal_energy_totals();

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	PackerEnergy current_total_energy = total_energy_current_state_assignment();
	debug_note_projected_deltaE_of_considered_substitution( current_total_energy - last_total_energy, 1.0, false );
	debug_note_accepted_substitution();
#endif
	//std::cout << "Accept: " << total_energy_current_state_assignment() << std::endl;

	return total_energy_current_state_assignment();
}

void
MinimalistFlexbbInteractionGraph::consider_substitution(
	int node_ind,
	int new_state,
	PackerEnergy & delta_energy,
	PackerEnergy & prev_energy_for_node
)
{
	if ( last_considered_substitution_unresolved() ) {
		resolve_uncommitted_substitution();
	}
	note_fixedbb_substitution();

	set_node_considering_alt_state( node_ind );
	delta_energy = get_minimalistflexbb_node( node_ind )->
		project_deltaE_for_substitution( new_state, prev_energy_for_node );

	set_total_energy_alternate_state_assignment( total_energy_current_state_assignment() + delta_energy );

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	debug_note_projected_deltaE_of_considered_substitution(
		delta_energy, get_minimalistflexbb_node( node_ind )->alternate_state_total_energy()  );
#endif

}

MinimalistFlexbbInteractionGraph::PackerEnergy
MinimalistFlexbbInteractionGraph::commit_considered_substitution()
{
	get_minimalistflexbb_node( node_considering_alt_state() )->commit_considered_substitution();
	note_last_considered_substitution_resolved();
	set_total_energy_current_state_assignment( total_energy_alternate_state_assignment() );

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	debug_note_accepted_substitution();
#endif

	//std::cout << "Accept: " << total_energy_current_state_assignment() << std::endl;

	return total_energy_current_state_assignment();
}

MinimalistFlexbbInteractionGraph::PackerEnergy
MinimalistFlexbbInteractionGraph::get_energy_current_state_assignment()
{
	return total_energy_current_state_assignment();
}

int
MinimalistFlexbbInteractionGraph::get_edge_memory_usage() const
{
	int sum = 0;
	for ( auto iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		sum += cast_minimalist_flexbb_edge(*iter)->count_dynamic_memory(); // close enough...
	}
	return sum;

}

void
MinimalistFlexbbInteractionGraph::print_current_state_assignment() const
{
	std::cerr << "Curr States: ";
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		std::cerr << "(" << ii << ", ";
		std::cerr << get_minimalistflexbb_node(ii)->get_current_state() << ") ";
	}
	std::cerr << std::endl;
}

/// @details noop.  This IG does not allow for inaccuracies in its energy function.
void
MinimalistFlexbbInteractionGraph::set_errorfull_deltaE_threshold( PackerEnergy  )
{}


MinimalistFlexbbInteractionGraph::PackerEnergy
MinimalistFlexbbInteractionGraph::get_energy_sum_for_vertex_group( int group_id )
{
	PackerEnergy esum = 0;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		if ( get_vertex_member_of_energy_sum_group( ii, group_id ) ) {
			esum += get_minimalistflexbb_node( ii )->curr_state_one_body_energy();
		}
	}

	for ( auto edge_iter = get_edge_list_begin();
			edge_iter != get_edge_list_end(); ++edge_iter ) {
		int first_node_ind = (*edge_iter)->get_first_node_ind();
		int second_node_ind = (*edge_iter)->get_second_node_ind();

		if ( get_vertex_member_of_energy_sum_group( first_node_ind, group_id )
				&& get_vertex_member_of_energy_sum_group( second_node_ind, group_id ) ) {
			esum += cast_minimalist_flexbb_edge(*edge_iter)->cur_energy();
		}
	}
	return esum;
}

/// Virtual functions from FlexbbInteractionGraph
void
MinimalistFlexbbInteractionGraph::consider_backbone_move(
	int bb_id,
	core::PackerEnergy & delta_energy,
	core::PackerEnergy & prev_flexseg_energy,
	bool & valid_motion,
	int & num_nodes_changing_state
)
{
	if ( last_considered_substitution_unresolved() ) {
		resolve_uncommitted_substitution();
	}

	note_bbjump_substitution();

	/// brace for an invalid move;
	/// early return statements "return" these values to the calling function.
	valid_motion = false;
	delta_energy = 0.0;
	prev_flexseg_energy = 0.0;

	reset_node_in_moving_flexseg_count();
	int moving_flexseg = get_flexseg_for_bb( bb_id );
	set_flexseg_considering_alt_bb( moving_flexseg );
	int altbb_for_flexseg = bb_id - get_flexseg_bb_offset( moving_flexseg );

	int representative = flexseg_members( moving_flexseg )[ 1 ];

	if ( get_flexbb_node( representative )->get_backbone_for_current_state() == altbb_for_flexseg ) {
		return;
	}

	auto const nflexseg_members = (int) flexseg_members( moving_flexseg ).size();
	//std::cout << "bbmove: " << moving_flexseg << " " << bb_id << " " << nflexseg_members;
	//std::cout << std::flush;
	for ( int ii = 1; ii <= nflexseg_members; ++ii ) {
		//std::cout << "." << std::flush;
		int iinode = flexseg_members( moving_flexseg )[ ii ];
		bool move_valid = get_minimalistflexbb_node( iinode )->
			prepare_for_altbb_move_to_closest_state( altbb_for_flexseg );
		if ( ! move_valid ) {
			//std::cout << std::endl;
			return;
		}
	}

	complete_deltaE_prediction_for_bbmove(
		delta_energy, prev_flexseg_energy,
		valid_motion, num_nodes_changing_state );
	//std::cout << "...done" << std::endl;
}

void
MinimalistFlexbbInteractionGraph::consider_bbmove_w_state_substitution(
	int node_ind,
	int new_state,
	core::PackerEnergy & delta_energy,
	core::PackerEnergy & prev_flexseg_energy,
	bool & valid_motion,
	int & num_nodes_changing_state
)
{
	if ( last_considered_substitution_unresolved() ) {
		resolve_uncommitted_substitution();
	}
	reset_node_in_moving_flexseg_count();

	/// brace for an invalid move;
	/// early return statements "return" these values to the calling function.
	valid_motion = false;
	delta_energy = 0.0;
	prev_flexseg_energy = 0.0;

	if ( get_flexbb_node( node_ind )->current_state() == 0 ) {
		return;
	}


	/// Decide: does this flexible-backbone substitution actually need to move the backbone?
	if (  flexseg_for_moltenres( node_ind ) == 0 ||
			get_flexbb_node( node_ind )->state_has_same_backbone_as_current( new_state ) ) {
		valid_motion = true;
		num_nodes_changing_state = 1;
		consider_substitution( node_ind, new_state, delta_energy, prev_flexseg_energy );
		return;
	}


	reset_node_in_moving_flexseg_count();
	note_bbjump_substitution();
	int moving_flexseg = flexseg_for_moltenres( node_ind );
	set_flexseg_considering_alt_bb( moving_flexseg );
	int altbb_for_flexseg = get_flexbb_node( node_ind )->state_info( new_state ).get_bb();

	//std::cout << "bbmove2: " << moving_flexseg << " " << altbb_for_flexseg << " " << flexseg_members(moving_flexseg).size();
	//std::cout << std::flush;

	for ( Size ii = 1; ii <= flexseg_members( moving_flexseg ).size(); ++ii ) {
		int iinode = flexseg_members( moving_flexseg )[ ii ];
		if ( iinode == node_ind ) {
			bool valid = get_minimalistflexbb_node( iinode )->prepare_for_altbb_move_to_state( new_state );
			if ( ! valid ) { return; }
		} else {
			bool move_valid = get_minimalistflexbb_node( iinode )->
				prepare_for_altbb_move_to_closest_state( altbb_for_flexseg );
			if ( ! move_valid ) {
				//std::cout << std::endl;
				return;
			}
		}
	}
	complete_deltaE_prediction_for_bbmove(
		delta_energy, prev_flexseg_energy,
		valid_motion, num_nodes_changing_state );
	//std::cout << std::endl;

}

void
MinimalistFlexbbInteractionGraph::complete_deltaE_prediction_for_bbmove(
	core::PackerEnergy & delta_energy,
	core::PackerEnergy & prev_flexseg_energy,
	bool & valid_motion,
	int & num_nodes_changing_state
)
{
	int moving_flexseg = flexseg_considering_alt_bb();
	auto const nflexseg_members = (int) flexseg_members( moving_flexseg ).size();
	Real total_frag_energy_alt( 0.0 );
	for ( int ii = 1; ii <= nflexseg_members; ++ii ) {
		int iinode = flexseg_members( moving_flexseg )[ ii ];
		total_frag_energy_alt += get_minimalistflexbb_node( iinode )->
			get_frag_energy_for_alt_bb_state();
	}

	Real total_frag_energy_curr( 0.0 );
	for ( int ii = 1; ii <= nflexseg_members; ++ii ) {
		int iinode = flexseg_members( moving_flexseg )[ ii ];
		total_frag_energy_curr += get_minimalistflexbb_node( iinode )->
			get_frag_energy_for_curr_bb_state_and_finalize_alt_energy_total();
	}


	valid_motion = true;
	delta_energy = total_frag_energy_alt - total_frag_energy_curr;
	prev_flexseg_energy = total_frag_energy_curr;

	set_last_considered_backbone_sub_valid( valid_motion );

	set_total_energy_alternate_state_assignment( total_energy_current_state_assignment() + delta_energy );

	num_nodes_changing_state = get_num_nodes_changing_state();

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	debug_note_projected_deltaE_of_considered_substitution( delta_energy, prev_flexseg_energy );
#endif
}


MinimalistFlexbbInteractionGraph::PackerEnergy
MinimalistFlexbbInteractionGraph::commit_considered_backbone_move(
	ObjexxFCL::FArray1_int & rotamer_on_node
)
{
	debug_assert( last_considered_backbone_sub_valid() );
	//get_minimalistflexbb_node( node_considering_alt_state() )->commit_alt_bb_substitution( rotamer_on_node );


	if ( last_considered_substitution_moved_the_backbone() ) {
		int moving_flexseg = flexseg_considering_alt_bb();
		for ( Size ii = 1; ii <= flexseg_members( moving_flexseg ).size(); ++ii ) {
			int iinode = flexseg_members( moving_flexseg )[ ii ];
			get_minimalistflexbb_node( iinode)->commit_alt_bb_substitution( rotamer_on_node );
		}
	} else {
		get_minimalistflexbb_node( node_considering_alt_state() )->commit_considered_substitution( rotamer_on_node );
	}

	note_last_considered_substitution_resolved();
	set_total_energy_current_state_assignment( total_energy_alternate_state_assignment() );

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	debug_note_accepted_substitution();
#endif
	//std::cout << "Accept BBMove: " << total_energy_current_state_assignment() << std::endl;


	return total_energy_current_state_assignment();
}


unsigned int
MinimalistFlexbbInteractionGraph::count_static_memory() const
{
	return sizeof( MinimalistFlexbbInteractionGraph );
}

unsigned int
MinimalistFlexbbInteractionGraph::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}

core::pack::interaction_graph::NodeBase*
MinimalistFlexbbInteractionGraph::create_new_node( int node_index, int num_states)
{
	return new MinimalistFlexbbNode( this, node_index, num_states );
}

core::pack::interaction_graph::EdgeBase*
MinimalistFlexbbInteractionGraph::create_new_edge( int index1, int index2)
{
	return new MinimalistFlexbbEdge( this, index1, index2 );
}


void
MinimalistFlexbbInteractionGraph::resolve_uncommitted_substitution()
{
	debug_assert( last_considered_substitution_unresolved() );
	if ( last_considered_substitution_kept_backbone_fixed() ) {
		get_minimalistflexbb_node( node_considering_alt_state() )->resolve_uncommitted_substitution();
	} else {
		debug_assert( last_considered_substitution_moved_the_backbone() );
		int moving_flexseg = flexseg_considering_alt_bb();
		for ( Size ii = 1; ii <= flexseg_members( moving_flexseg ).size(); ++ii ) {
			int iinode = flexseg_members( moving_flexseg )[ ii ];
			get_minimalistflexbb_node( iinode)->resolve_uncommitted_substitution();
		}
	}
	note_last_considered_substitution_resolved();

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	debug_note_rejected_substitution();
#endif
}

}
}
}


