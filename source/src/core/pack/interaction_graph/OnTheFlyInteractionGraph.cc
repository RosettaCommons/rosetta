// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/OnTheFlyInteractionGraph.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.hh>

// Package headers
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/util.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace interaction_graph {

bool const debug = { false };
static basic::Tracer T("core.pack.interaction_graph.otf_ig", basic::t_error );

/// @brief main constructor, no default or copy constructors
///
/// @detailed
OnTheFlyNode::OnTheFlyNode(
	InteractionGraphBase * owner,
	int node_id,
	int num_states
) :
	FixedBBNode( owner, node_id, num_states ),
	rotamer_set_( 0 ),
	rotamers_( num_states ),
	sc_bounding_spheres_( num_states, std::make_pair( Vector( 0.0 ), Real( 0.0 ) ) ),
	bb_bounding_sphere_( std::make_pair( Vector( 0.0 ), Real( 0.0 ) )),
	num_aa_types_( get_on_the_fly_owner()->get_num_aatypes() ),
	num_states_for_aatype_( num_aa_types_, 0 ),
	state_offset_for_aatype_( num_aa_types_, 0 ),
	sparse_mat_info_for_state_( num_states ),
	one_body_energies_( num_states, 0.0 ),
	distinguish_backbone_and_sidechain_( false ) // by default, do not
{

}

/// @details
OnTheFlyNode::~OnTheFlyNode()
{}

/// @details avoid polymorphic lookup later by nabbing rotamers now.
void
OnTheFlyNode::set_rotamers(
	rotamer_set::RotamerSetCOP rotamers
)
{
	assert( rotamers->num_rotamers() == (Size) get_num_states() );
	rotamer_set_ = rotamers;
	for ( Size ii = 1; ii <= rotamer_set_->num_rotamers(); ++ii ) {
		rotamers_[ ii ] = rotamer_set_->rotamer( ii );
		if ( ii == 1 ) {
			bb_bounding_sphere_.first  = scoring::compute_bb_centroid( *rotamers_[ ii ] );
			bb_bounding_sphere_.second = scoring::compute_bb_radius(   *rotamers_[ ii ], bb_bounding_sphere_.first );
		}
		sc_bounding_spheres_[ ii ].first  = scoring::compute_sc_centroid( *rotamers_[ ii ] );
		sc_bounding_spheres_[ ii ].second = scoring::compute_sc_radius(   *rotamers_[ ii ], sc_bounding_spheres_[ ii ].first );
	}

	// figure out which residue-type group each rotamer is a member of
	Size curr_restype = 1;
	Size count_for_restype = 1;
	Size const num_aa_types_ = rotamers->get_n_residue_types();
	//num_states_for_aatype_.resize( num_aa_types_ );
	for ( Size ii = 1; ii <= rotamers_.size(); ++ii ) {

		sparse_mat_info_for_state_[ ii ].set_aa_type( curr_restype );
		sparse_mat_info_for_state_[ ii ].set_state_ind_for_this_aa_type( count_for_restype );
		++count_for_restype;
		++num_states_for_aatype_[ curr_restype ];
		while ( count_for_restype > rotamers->get_n_rotamers_for_residue_group( curr_restype )) {
			// increment curr_restype and skip over restypes with 0 rotamers
			++curr_restype;
			count_for_restype = 1;
			if ( curr_restype > num_aa_types_ ) break;
			state_offset_for_aatype_[ curr_restype ] = ii;
		}
	}
}

void
OnTheFlyNode::distinguish_backbone_and_sidechain( bool setting )
{
	if ( get_num_incident_edges() != 0 ) {
		utility_exit_with_message( "ERROR:: Must set distinguish_backbone_and_sidechain before adding edges" );
	}
	distinguish_backbone_and_sidechain_ = setting;
}


void
OnTheFlyNode::zero_one_body_energies()
{
	std::fill( one_body_energies_.begin(), one_body_energies_.end(), 0.0f );
}

void
OnTheFlyNode::add_to_one_body_energies( ObjexxFCL::FArray1_float & energy1b )
{
	for ( Size ii = 1; ii <= one_body_energies_.size(); ++ii ) {
		one_body_energies_[ ii ] += energy1b( ii );
	}
}

void
OnTheFlyNode::update_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] = energy;
}


/// @details sets one body energies for a given state
void
OnTheFlyNode::set_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] = energy;
}

/// @details increments one body energies for a given state
void
OnTheFlyNode::add_to_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] +=  energy;
}



void
OnTheFlyNode::zero_one_body_energy( int state )
{
	one_body_energies_[ state ] = 0;
}

/// @details reports the amount of dynamic memory this node allocates
/// recursing to its base class
unsigned int
OnTheFlyNode::count_dynamic_memory() const
{
	unsigned int total_memory = NodeBase::count_dynamic_memory();

	total_memory += rotamers_.size() * sizeof ( conformation::ResidueCOP );
	total_memory += sc_bounding_spheres_.size() * sizeof( BoundingSphere );
	total_memory += num_states_for_aatype_.size() * sizeof( int );
	total_memory += state_offset_for_aatype_.size() * sizeof( int );
	total_memory += sparse_mat_info_for_state_.size() * sizeof( SparseMatrixIndex );
	total_memory += one_body_energies_.size() * sizeof( core::PackerEnergy );

	return total_memory;
}


core::PackerEnergy
OnTheFlyNode::compute_rotamer_pair_energy(
	int edge_making_energy_request,
	int state_this,
	int state_other
) const
{
	using namespace scoring;
	using namespace scoring::methods;

	core::PackerEnergy esum( 0.0 );

	OnTheFlyEdge const & spanning_edge( * get_incident_otf_edge( edge_making_energy_request ) );
	OnTheFlyNode const & neighbor( * get_adjacent_otf_node( edge_making_energy_request ) );

	get_on_the_fly_owner()->note_rpe_calculated();

	if ( spanning_edge.short_range_interactions_exist() ) {
		EnergyMap tbody_emap;

		switch ( spanning_edge.eval_type( get_node_index() )) {
		case ( sc_sc ) :
			scoring::eval_scsc_sr2b_energies(
				get_rotamer( state_this ),
				neighbor.get_rotamer( state_other ),
				sc_bounding_sphere( state_this ).first,
				neighbor.sc_bounding_sphere( state_other ).first,
				sc_bounding_sphere( state_this ).second,
				neighbor.sc_bounding_sphere( state_other ).second,
				get_on_the_fly_owner()->pose(),
				get_on_the_fly_owner()->score_function(),
				tbody_emap );

			/*get_on_the_fly_owner()->score_function().eval_ci_2b_sc_sc(
				get_rotamer( state_this ),
				neighbor.get_rotamer( state_other ),
				get_on_the_fly_owner()->pose(),
				tbody_emap );

			get_on_the_fly_owner()->score_function().eval_cd_2b_sc_sc(
				get_rotamer( state_this ),
				neighbor.get_rotamer( state_other ),
				get_on_the_fly_owner()->pose(),
				tbody_emap
			);*/

		break;
		case ( sc_whole ) :
			get_on_the_fly_owner()->score_function().eval_ci_2b_sc_sc(
				get_rotamer( state_this ),
				neighbor.get_rotamer( state_other ),
				get_on_the_fly_owner()->pose(),
				tbody_emap );

			get_on_the_fly_owner()->score_function().eval_cd_2b_sc_sc(
				get_rotamer( state_this ),
				neighbor.get_rotamer( state_other ),
				get_on_the_fly_owner()->pose(),
				tbody_emap
			);
			/// Also evaluate the bb/sc energy between other/this
			get_on_the_fly_owner()->score_function().eval_ci_2b_bb_sc(
				neighbor.get_rotamer( state_other ),
				get_rotamer( state_this ),
				get_on_the_fly_owner()->pose(),
				tbody_emap );

			get_on_the_fly_owner()->score_function().eval_cd_2b_bb_sc(
				neighbor.get_rotamer( state_other ),
				get_rotamer( state_this ),
				get_on_the_fly_owner()->pose(),
				tbody_emap
			);

		break;
		case ( whole_sc ) :
			get_on_the_fly_owner()->score_function().eval_ci_2b_sc_sc(
				get_rotamer( state_this ),
				neighbor.get_rotamer( state_other ),
				get_on_the_fly_owner()->pose(),
				tbody_emap );

			get_on_the_fly_owner()->score_function().eval_cd_2b_sc_sc(
				get_rotamer( state_this ),
				neighbor.get_rotamer( state_other ),
				get_on_the_fly_owner()->pose(),
				tbody_emap
			);
			/// Also evaluate the bb/sc energy between this/other
			get_on_the_fly_owner()->score_function().eval_ci_2b_bb_sc(
				get_rotamer( state_this ),
				neighbor.get_rotamer( state_other ),
				get_on_the_fly_owner()->pose(),
				tbody_emap );

			get_on_the_fly_owner()->score_function().eval_cd_2b_bb_sc(
				get_rotamer( state_this ),
				neighbor.get_rotamer( state_other ),
				get_on_the_fly_owner()->pose(),
				tbody_emap
			);


		break;
		case ( whole_whole ) :

			get_on_the_fly_owner()->score_function().eval_ci_2b(
				get_rotamer( state_this ),
				neighbor.get_rotamer( state_other ),
				get_on_the_fly_owner()->pose(),
				tbody_emap );

			get_on_the_fly_owner()->score_function().eval_cd_2b(
				get_rotamer( state_this ),
				neighbor.get_rotamer( state_other ),
				get_on_the_fly_owner()->pose(),
				tbody_emap
			);

		break;
		}
		esum += static_cast< core::PackerEnergy > ( get_on_the_fly_owner()->score_function().weights().dot( tbody_emap ) );

	}

	/// Gen-Borne will currently evaluate too many long range energies...
	/// For now, long range interactions are computed for the whole residue, and do not divide it into pieces
	if ( get_incident_otf_edge( edge_making_energy_request )->long_range_interactions_exist() ) {
		EnergyMap emap;
		for ( ScoreFunction::LR_2B_MethodIterator iter = get_on_the_fly_owner()->score_function().long_range_energies_begin(),
					iter_end = get_on_the_fly_owner()->score_function().long_range_energies_end();
					iter != iter_end; ++iter ) {
			(*iter)->residue_pair_energy(
				get_rotamer( state_this ),
				get_adjacent_otf_node( edge_making_energy_request )->get_rotamer( state_other ),
				get_on_the_fly_owner()->pose(),
				get_on_the_fly_owner()->score_function(),
				emap );
		}
		esum += static_cast< core::PackerEnergy > ( get_on_the_fly_owner()->score_function().weights().dot( emap ) );
	}


	if ( get_adjacent_otf_node( edge_making_energy_request )->get_rotamer( state_other ).aa() == chemical::aa_pro )
	{
		esum += get_incident_otf_edge( edge_making_energy_request )->
			get_proline_correction_for_node( get_node_index(), state_this );
	}
	if ( get_rotamer( state_this ).aa() == chemical::aa_pro )
	{
		esum += get_incident_otf_edge( edge_making_energy_request )->
			get_proline_correction_for_node( get_index_of_adjacent_node( edge_making_energy_request ),
			state_other
		);
	}


	return get_incident_edge( edge_making_energy_request )->edge_weight() * esum;

	//	if ( ! get_on_the_fly_owner()->check_empty_weight_map() ) {
	//esum *= get_on_the_fly_owner()->get_residue_weights(this_rotamer.get_resid(), aa_this,
	// other_rotamer.get_resid(), aa_other);

	//return esum;

}


//-------------------------------------------------------------//

OnTheFlyEdge::OnTheFlyEdge(
	InteractionGraphBase * owner,
	int first_node_ind,
	int second_node_ind
)
:
	FixedBBEdge( owner, first_node_ind, second_node_ind ),
	long_range_interactions_exist_( false ),
	short_range_interactions_exist_( false )
{
	bool distinguish_sc_bb[ 2 ];
	for ( int ii = 0; ii < 2; ++ii) {
		proline_corrections_[ ii ].resize( get_num_states_for_node( ii ) );
		std::fill( proline_corrections_[ ii ].begin(), proline_corrections_[ ii ].end(), 0.0f );
		distinguish_sc_bb[ ii ] = get_otf_node( ii )->distinguish_backbone_and_sidechain();
	}

	if ( distinguish_sc_bb[ 0 ] && distinguish_sc_bb[ 1 ] ) {
		eval_types_[ 0 ] = eval_types_[ 1 ] = sc_sc;
	} else if ( ! distinguish_sc_bb[ 0 ] && distinguish_sc_bb[ 1 ] ) {
		eval_types_[ 0 ] = whole_sc;
		eval_types_[ 1 ] = sc_whole;
	} else if ( distinguish_sc_bb[ 0 ] && ! distinguish_sc_bb[ 1 ] ) {
		eval_types_[ 1 ] = whole_sc;
		eval_types_[ 0 ] = sc_whole;
	} else {
		eval_types_[ 1 ] = eval_types_[ 0 ] = whole_whole;
	}
}


OnTheFlyEdge::~OnTheFlyEdge(){}

void
OnTheFlyEdge::set_ProCorrection_values(
	int node_not_necessarily_proline,
	int state,
	core::PackerEnergy bb_nonprobb_E,
	core::PackerEnergy bb_probb_E,
	core::PackerEnergy sc_nonprobb_E,
	core::PackerEnergy sc_probb_E
)
{
	int const which_node = node_not_necessarily_proline == get_node_index( 0 ) ? 0 : 1;

	proline_corrections_[ which_node ][ state ] =
		sc_probb_E + 0.5 * bb_probb_E -
		(sc_nonprobb_E + 0.5 * bb_nonprobb_E);
}


unsigned int
OnTheFlyEdge::count_dynamic_memory() const
{
	unsigned int total_memory_usage = EdgeBase::count_dynamic_memory();

	total_memory_usage += proline_corrections_[ 0 ].size() * sizeof( core::PackerEnergy );
	total_memory_usage += proline_corrections_[ 1 ].size() * sizeof( core::PackerEnergy );

	return total_memory_usage;
}

//-------------------------------------------------------------//

OnTheFlyInteractionGraph::OnTheFlyInteractionGraph(int num_nodes )
:
	FixedBBInteractionGraph( num_nodes ),
	num_aa_types_( 0 )
{}

OnTheFlyInteractionGraph::~OnTheFlyInteractionGraph()
{}

bool
OnTheFlyInteractionGraph::distinguish_backbone_and_sidechain_for_node( int node_ind ) const
{
	return get_on_the_fly_node( node_ind )->distinguish_backbone_and_sidechain();
}

void
OnTheFlyInteractionGraph::distinguish_backbone_and_sidechain_for_node( int node_ind, bool setting )
{
	get_on_the_fly_node( node_ind )->distinguish_backbone_and_sidechain( setting );
}


void
OnTheFlyInteractionGraph::initialize(
	rotamer_set::RotamerSetsBase const & rot_sets_base
)
{
	rotamer_set::RotamerSets const & rot_sets( static_cast< rotamer_set::RotamerSets const & > (rot_sets_base) );

	// determine max # of residue types
	Size max_nresgroups = 0;
	for ( Size ii = 1; ii <= rot_sets.nmoltenres(); ++ii ) {
		Size ii_nresgroups =  rot_sets.rotamer_set_for_moltenresidue( ii )->get_n_residue_groups();
		if ( ii_nresgroups > max_nresgroups ) max_nresgroups = ii_nresgroups;
	}

	//"aa types" means "distinct groups of rotamers" -- this ig has no idea
	// what an amino acid is or why they might be different from one another
	num_aa_types_ = max_nresgroups;


	for ( Size ii = 1; ii <= rot_sets.nmoltenres(); ++ii ) {
		Size const ii_num_states = rot_sets.rotamer_set_for_moltenresidue( ii )->num_rotamers();
		set_num_states_for_node( ii, ii_num_states );
		get_on_the_fly_node( ii )->set_rotamers( rot_sets.rotamer_set_for_moltenresidue( ii ) );
	}

}

void
OnTheFlyInteractionGraph::set_score_function( ScoreFunction const & sfxn )
{
	score_function_ = sfxn.clone();
	if ( pose_ ) (*score_function_)(*pose_); // rescore the pose with the input score function
}

void
OnTheFlyInteractionGraph::set_pose( pose::Pose const & pose )
{
	pose_ = new pose::Pose( pose );
	if ( score_function_ ) (*score_function_)(*pose_);
}

void
OnTheFlyInteractionGraph::zero_one_body_energy_for_node_state(
	int node_ind,
	int state
)
{
	get_on_the_fly_node( node_ind )->zero_one_body_energy( state );
}

void
OnTheFlyInteractionGraph::add_to_one_body_energy_for_node_state(
	int node_ind,
	int state,
	core::PackerEnergy one_body_energy
)
{
	get_on_the_fly_node( node_ind )->add_to_one_body_energy( state, one_body_energy);
}

void
OnTheFlyInteractionGraph::set_one_body_energy_for_node_state(
	int node_ind,
	int state,
	core::PackerEnergy one_body_energy
)
{
	get_on_the_fly_node( node_ind )->set_one_body_energy( state, one_body_energy);
}


core::PackerEnergy
OnTheFlyInteractionGraph::get_one_body_energy_for_node_state( int node, int state)
{
	return get_on_the_fly_node( node )->get_one_body_energy( state );
}

void
OnTheFlyInteractionGraph::reset_rpe_calculations_count()
{
	num_rpe_calcs_ = 0;
}


Size
OnTheFlyInteractionGraph::get_num_rpe_calculations_count() const
{
	return num_rpe_calcs_;
}


void
OnTheFlyInteractionGraph::set_sparse_aa_info_for_edge(
	int node1,
	int node2,
	FArray2_bool const & sparse_conn_info
)
{
	OnTheFlyEdge* edge = (OnTheFlyEdge*) find_edge( node1, node2 );
	if (edge != 0) {
		edge->set_sparse_aa_info( sparse_conn_info );
	}
}

void
OnTheFlyInteractionGraph::set_ProCorrection_values_for_edge(
	int node1,
	int node2,
	int node_not_neccessarily_proline,
	int state,
	core::PackerEnergy bb_nonprobb_E,
	core::PackerEnergy bb_probb_E,
	core::PackerEnergy sc_nonprobb_E,
	core::PackerEnergy sc_probb_E
)
{
	OnTheFlyEdge* edge = (OnTheFlyEdge*) find_edge( node1, node2 );
	if (edge != 0)
	{
		edge->set_ProCorrection_values(
			node_not_neccessarily_proline, state,
			bb_nonprobb_E, bb_probb_E, sc_nonprobb_E, sc_probb_E );
	}
}

void
OnTheFlyInteractionGraph::note_short_range_interactions_exist_for_edge(
	int node1,
	int node2
)
{
	OnTheFlyEdge* edge = (OnTheFlyEdge*) find_edge( node1, node2 );
	if (edge != 0) {
		edge->note_short_range_interactions_exist();
	}
}

void
OnTheFlyInteractionGraph::note_long_range_interactions_exist_for_edge(
	int node1,
	int node2
)
{
	OnTheFlyEdge* edge = (OnTheFlyEdge*) find_edge( node1, node2 );
	if (edge != 0) {
		edge->note_long_range_interactions_exist( );
	}
}


unsigned int
OnTheFlyInteractionGraph::count_dynamic_memory() const
{
	unsigned int total_memory = InteractionGraphBase::count_dynamic_memory();
	return total_memory;
}

} // namespace interaction_graph
} // namespace pack
} // namespace core
