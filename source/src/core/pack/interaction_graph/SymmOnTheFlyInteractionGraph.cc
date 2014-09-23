// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/SymmOnTheFlyInteractionGraph.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/interaction_graph/SymmOnTheFlyInteractionGraph.hh>

// Package headers
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
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

static thread_local basic::Tracer T( "core.pack.interaction_graph.otf_ig", basic::t_error );

/// @brief main constructor, no default or copy constructors
///
/// @detailed
SymmOnTheFlyNode::SymmOnTheFlyNode(
	InteractionGraphBase * owner,
	int node_id,
	int num_states
) :
	FixedBBNode( owner, node_id, num_states ),
	rotamer_set_( /* 0 */ ),
	rotamers_( num_states ),
	sc_bounding_spheres_( num_states, std::make_pair( Vector( 0.0 ), Real( 0.0 ) ) ),
	bb_bounding_sphere_( std::make_pair( Vector( 0.0 ), Real( 0.0 ) )),
	num_res_types_( 0 ), // wait until rotamers are set
	num_restype_groups_( 0 ), // wait until rotamer are set
	//sparse_mat_info_for_state_( num_states ),
	one_body_energies_( num_states, 0.0 ),
	distinguish_backbone_and_sidechain_( false ) // by default, do not
{

}

/// @details
SymmOnTheFlyNode::~SymmOnTheFlyNode()
{}

/// @details avoid polymorphic lookup later by nabbing rotamers now.
void
SymmOnTheFlyNode::set_rotamers(
	rotamer_set::RotamerSetCOP rotamers
)
{
	num_res_types_ = rotamers->get_n_residue_types();
	num_restype_groups_ = rotamers->get_n_residue_groups();

	num_states_for_restype_group_.resize( num_restype_groups_, 0 ),
	state_offset_for_restype_group_.resize( num_restype_groups_, 0 ),

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

	Size const nsubunits = get_on_the_fly_owner()->symm_info()->subunits();
	Size const nres_asu = get_on_the_fly_owner()->symm_info()->num_independent_residues();
	core::pose::Pose const & pose = get_on_the_fly_owner()->pose();
	rotamer_representatives_.resize( nsubunits );
	for ( Size ii = 1; ii <= nsubunits; ++ii ) {
		state_offset_for_restype_group_[ ii ] = rotamers->get_residue_type_begin( ii ) - 1;

		rotamer_representatives_[ ii ].resize( num_res_types_ );
		for ( int jj = 1; jj <= num_res_types_; ++jj ) {
			rotamer_representatives_[ ii ][ jj ] = rotamers->rotamer( rotamers->get_residue_type_begin( jj ) )->clone();
			rotamer_representatives_[ ii ][ jj ]->seqpos( rotamers->resid() + (ii-1)*nres_asu );
			rotamer_representatives_[ ii ][ jj ]->copy_residue_connections_from( pose.residue( rotamer_representatives_[ ii ][ jj ]->seqpos() ) );
		}
	}

	// figure out which residue-type group each rotamer is a member of
	Size curr_restype_group = 1;
	Size count_for_restype_group = 1;
	Size const num_restype_groups = rotamers->get_n_residue_groups();
	for ( Size ii = 1; ii <= rotamers_.size(); ++ii ) {

		//sparse_mat_info_for_state_[ ii ].set_aa_type( curr_restype_group );
		//sparse_mat_info_for_state_[ ii ].set_state_ind_for_this_aa_type( count_for_restype_group );
		++count_for_restype_group;
		while ( count_for_restype_group > rotamers->get_n_rotamers_for_residue_group( curr_restype_group )) {
			// increment curr_restype_group and skip over restypes with 0 rotamers
			++curr_restype_group;
			count_for_restype_group = 1;
			if ( curr_restype_group > num_restype_groups ) break;
			//state_offset_for_aatype_[ curr_restype_group ] = ii;
		}
	}
}

void
SymmOnTheFlyNode::distinguish_backbone_and_sidechain( bool setting )
{
	if ( get_num_incident_edges() != 0 ) {
		utility_exit_with_message( "ERROR:: Must set distinguish_backbone_and_sidechain before adding edges" );
	}
	distinguish_backbone_and_sidechain_ = setting;
}


void
SymmOnTheFlyNode::zero_one_body_energies()
{
	std::fill( one_body_energies_.begin(), one_body_energies_.end(), 0.0f );
}

void
SymmOnTheFlyNode::add_to_one_body_energies( ObjexxFCL::FArray1_float & energy1b )
{
	for ( Size ii = 1; ii <= one_body_energies_.size(); ++ii ) {
		one_body_energies_[ ii ] += energy1b( ii );
	}
}

void
SymmOnTheFlyNode::update_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] = energy;
}


/// @details sets one body energies for a given state
void
SymmOnTheFlyNode::set_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] = energy;
}

/// @details increments one body energies for a given state
void
SymmOnTheFlyNode::add_to_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] +=  energy;
}



void
SymmOnTheFlyNode::zero_one_body_energy( int state )
{
	one_body_energies_[ state ] = 0;
}

/// @details reports the amount of dynamic memory this node allocates
/// recursing to its base class
unsigned int
SymmOnTheFlyNode::count_dynamic_memory() const
{
	unsigned int total_memory = NodeBase::count_dynamic_memory();

	total_memory += rotamers_.size() * sizeof ( conformation::ResidueCOP );
	//total_memory += num_states_for_aatype_.size() * sizeof( int );
	//total_memory += state_offset_for_aatype_.size() * sizeof( int );
	//total_memory += sparse_mat_info_for_state_.size() * sizeof( SparseMatrixIndex );
	total_memory += one_body_energies_.size() * sizeof( core::PackerEnergy );

	return total_memory;
}


/// @details Compute all the energies that are appropriately two body energies for
/// this edge. In the outter loop, pick one of the two nodes on this edge and consider
/// its rotamer as having originated from the asymmetric unit.  In the inner loop,
/// iterate across all subunits, considering each subunit as the origination of the
/// rotamer for the other node (that is, if ii == 1, then "this" node is treated as
/// if its rotamer is originating from the asymmetric unit, and the "other" node
/// takes its rotamer from subunit jj -- which may very well be the asymmetric unit.
/// If however, ii == 2, then the "other" node is treated as if its rotamer is originating
/// from the asymmetric unit and the "this" node's rotamer is taken from jj).
/// Now, inside this inner loop, with ii and jj known, first, ask the edge whether
/// given ii, jj, and the amino acid types of the two rotamers, the two rotamers
/// implied by ii and jj should have their interaction energies calculated.  If not, continue.
/// If so, then ask each edge for a representative rotamer, telling each node which
/// subunit that rotamer should come from.  This may likely require rotating and translating the
/// rotamers from the asymmetric unit into the appropriate symmetric clone
core::PackerEnergy
SymmOnTheFlyNode::compute_rotamer_pair_energy(
	int edge_making_energy_request,
	int state_this,
	int state_other
) const
{
	using namespace scoring;
	using namespace scoring::methods;

	Size const asu_index = 1; // OK will it always be the case that the asymmetic unit is the first subunit?  I'm assuming it always will be.

	core::PackerEnergy esum( 0.0 );

	SymmOnTheFlyEdge const & spanning_edge( * get_incident_otf_edge( edge_making_energy_request ) );
	SymmOnTheFlyNode const & neighbor( * get_adjacent_otf_node( edge_making_energy_request ) );

	get_on_the_fly_owner()->note_rpe_calculated();

	Size nsubunits = get_on_the_fly_owner()->symm_info()->subunits();
	Size this_restype_group = rotamer_set_->get_residue_group_index_for_rotamer( state_this );
	Size other_restype_group = neighbor.rotamer_set_->get_residue_group_index_for_rotamer( state_other );

	EnergyMap tbody_emap;

	for ( Size ii = 1; ii <= 2; ++ii ) {
		bool which_way_is_up = (ii == 1) == (get_node_index() < neighbor.get_node_index());
		Size ii_restype_group = which_way_is_up ? this_restype_group : other_restype_group;
		Size jj_restype_group = which_way_is_up ? other_restype_group : this_restype_group;

		for ( Size jj = 1; jj <= nsubunits; ++jj ) {

			tbody_emap.zero();

			/// iijj_score_multiply does double duty:
			/// 1 ) it indicates whether a pair of amino acids on certain subunits are adjacent
			/// 2 ) it gives the symmetry multiplication factor for the subunit pair
			/// It does not, however, indicate what the multiplication factor should be for long-range two body
			/// interactions if the amino acids are not close enough to qualify for short-range interactionsm
			unsigned char iijj_score_multiply = spanning_edge.residues_adjacent_for_subunit_pair( ii, jj, ii_restype_group, jj_restype_group );

			if ( iijj_score_multiply == 0 && ! spanning_edge.long_range_interactions_exist() ) continue;

			Size this_subunit = ii == 1 ? asu_index : jj;
			Size other_subunit = ii == 1 ? jj : asu_index;

			core::conformation::Residue const & this_rotamer( get_rotamer( state_this, this_subunit ));
			core::conformation::Residue const & other_rotamer( neighbor.get_rotamer( state_other, other_subunit ));

			BoundingSphere this_sc_bounding_sphere = sc_bounding_sphere( state_this, this_subunit );
			BoundingSphere other_sc_bounding_sphere = neighbor.sc_bounding_sphere( state_other, other_subunit );
			BoundingSphere this_bb_bounding_sphere = bb_bounding_sphere( this_subunit );
			BoundingSphere other_bb_bounding_sphere = neighbor.bb_bounding_sphere( other_subunit );

			// evaluate short-ranged interactions for this pair, if appropriate
			if ( iijj_score_multiply != 0 && spanning_edge.short_range_interactions_exist() ) {

				switch ( spanning_edge.eval_type( get_node_index() )) {
				case ( sc_sc ) :
					scoring::eval_scsc_sr2b_energies(
						this_rotamer, other_rotamer,
						this_sc_bounding_sphere.first, other_sc_bounding_sphere.first,
						this_sc_bounding_sphere.second, other_sc_bounding_sphere.second,
						get_on_the_fly_owner()->pose(), get_on_the_fly_owner()->score_function(),
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
						this_rotamer,
						other_rotamer,
						get_on_the_fly_owner()->pose(),
						tbody_emap );

					get_on_the_fly_owner()->score_function().eval_cd_2b_sc_sc(
						this_rotamer,
						other_rotamer,
						get_on_the_fly_owner()->pose(),
						tbody_emap
					);
					/// Also evaluate the bb/sc energy between other/this
					get_on_the_fly_owner()->score_function().eval_ci_2b_bb_sc(
						other_rotamer,
						this_rotamer,
						get_on_the_fly_owner()->pose(),
						tbody_emap );

					get_on_the_fly_owner()->score_function().eval_cd_2b_bb_sc(
						other_rotamer,
						this_rotamer,
						get_on_the_fly_owner()->pose(),
						tbody_emap
					);

				break;
				case ( whole_sc ) :
					get_on_the_fly_owner()->score_function().eval_ci_2b_sc_sc(
						this_rotamer,
						other_rotamer,
						get_on_the_fly_owner()->pose(),
						tbody_emap );

					get_on_the_fly_owner()->score_function().eval_cd_2b_sc_sc(
						this_rotamer,
						other_rotamer,
						get_on_the_fly_owner()->pose(),
						tbody_emap
					);
					/// Also evaluate the bb/sc energy between this/other
					get_on_the_fly_owner()->score_function().eval_ci_2b_bb_sc(
						this_rotamer,
						other_rotamer,
						get_on_the_fly_owner()->pose(),
						tbody_emap );

					get_on_the_fly_owner()->score_function().eval_cd_2b_bb_sc(
						this_rotamer,
						other_rotamer,
						get_on_the_fly_owner()->pose(),
						tbody_emap
					);


				break;
				case ( whole_whole ) :

					get_on_the_fly_owner()->score_function().eval_ci_2b(
						this_rotamer,
						other_rotamer,
						get_on_the_fly_owner()->pose(),
						tbody_emap );

					get_on_the_fly_owner()->score_function().eval_cd_2b(
						this_rotamer,
						other_rotamer,
						get_on_the_fly_owner()->pose(),
						tbody_emap
					);

				break;
				}

				if ( false ) {
					std::cout << get_node_index() << " " << neighbor.get_node_index() << " ii: " << ii << " jj: " << jj << " " << (int) iijj_score_multiply << " * " << static_cast< core::PackerEnergy > ( get_on_the_fly_owner()->score_function().weights().dot( tbody_emap ) ) << " = " <<  iijj_score_multiply * static_cast< core::PackerEnergy > ( get_on_the_fly_owner()->score_function().weights().dot( tbody_emap ) )<< std::endl;
					std::cout << "this cbeta: " << this_rotamer.xyz( "CB" ).x() << " " << this_rotamer.xyz( "CB" ).y() << " " << this_rotamer.xyz( "CB" ).z();
					std::cout << "; other cbeta: " << other_rotamer.xyz( "CB" ).x() << " " << other_rotamer.xyz( "CB" ).y() << " " << other_rotamer.xyz( "CB" ).z() << std::endl;;
				}
				esum += iijj_score_multiply *	static_cast< core::PackerEnergy >
					( get_on_the_fly_owner()->score_function().weights().dot( tbody_emap ) );
			}

			/// Long range interaction energies for this rotamer pair
			if ( spanning_edge.long_range_interactions_exist() ) {
				EnergyMap emap;
				Size lr_iijj_score_multiply( iijj_score_multiply );
				if ( lr_iijj_score_multiply == 0 ) {
					lr_iijj_score_multiply = Size(get_on_the_fly_owner()->symm_info()->score_multiply(
						this_rotamer.seqpos(), other_rotamer.seqpos() ) );
				}
				if ( lr_iijj_score_multiply != 0 ) {
					for ( ScoreFunction::LR_2B_MethodIterator
							iter = get_on_the_fly_owner()->score_function().long_range_energies_begin(),
							iter_end = get_on_the_fly_owner()->score_function().long_range_energies_end();
							iter != iter_end; ++iter ) {
						(*iter)->residue_pair_energy(
							this_rotamer, other_rotamer,
							get_on_the_fly_owner()->pose(),
							get_on_the_fly_owner()->score_function(),
							emap );
					}
					esum += lr_iijj_score_multiply * static_cast< core::PackerEnergy >
						( get_on_the_fly_owner()->score_function().weights().dot( emap ) );
				}
			}
		}
	}


	if ( get_adjacent_otf_node( edge_making_energy_request )->get_rotamer( state_other, 1 ).aa() == chemical::aa_pro ) {
		esum += get_incident_otf_edge( edge_making_energy_request )->
			get_proline_correction_for_node( get_node_index(), state_this );
		//std::cout << "adding proline correction 1: " << get_incident_otf_edge( edge_making_energy_request )->
		//	get_proline_correction_for_node( get_node_index(), state_this ) << std::endl;
	}
	if ( get_rotamer( state_this, 1 ).aa() == chemical::aa_pro ) {
		esum += get_incident_otf_edge( edge_making_energy_request )->
			get_proline_correction_for_node( get_index_of_adjacent_node( edge_making_energy_request ),
			state_other
		);
		//std::cout << "adding proline correction 2: " <<  get_incident_otf_edge( edge_making_energy_request )->
    //  get_proline_correction_for_node( get_index_of_adjacent_node( edge_making_energy_request ), state_other ) << std::endl;
	}

	return get_incident_edge( edge_making_energy_request )->edge_weight() * esum;
}


/// @brief Returns a reference to the rotamer object in the requested subunit.  This reference
/// is valid only until the next call to get_rotamer, and which point, the coordinates inside
/// the requested rotamer may change.
conformation::Residue const &
SymmOnTheFlyNode::get_rotamer( int state, int subunit ) const
{
	/// ASSUMPTION: the asymmetric unit is subunit #1
	if ( subunit == 1 ) {
		return *rotamers_[ state ];
	}

	int state_aa = rotamer_set_->get_residue_type_index_for_rotamer( state );
	//if ( rotamers_currently_reprsented_[ subunit ][ state_aa ] == state ) {
	//	return *rotamer_representatives_[ subunit ][ state_aa ];
	//}
	//rotamers_currently_represented_[ subunit ][ state_aa ] = state;

	conformation::Residue & rotamer = *rotamer_representatives_[ subunit ][ state_aa ];
	conformation::Residue const & asymm_rotamer = *rotamers_[ state ];
	assert( &rotamer.type() == & asymm_rotamer.type() );

	// copy torsions
	rotamer.mainchain_torsions() = asymm_rotamer.mainchain_torsions();
	rotamer.chi() = asymm_rotamer.chi();

	// rotate coordinates
	SymmOnTheFlyInteractionGraph::HTReal transform = get_on_the_fly_owner()->symmetric_transform( subunit );
	for ( Size ii = 1, iiend = rotamer.natoms(); ii <= iiend; ++ii ) {
		rotamer.set_xyz( ii, transform * asymm_rotamer.xyz( ii ) );
	}

	// rotate the act_coord
	rotamer.actcoord() = transform * asymm_rotamer.actcoord();

	// orbital rotation would happen here, too, I think

	return rotamer;
}

/// @brief Returns a bounding sphere for the sidechain of a given state on a particular subunit.
SymmOnTheFlyNode::BoundingSphere
SymmOnTheFlyNode::sc_bounding_sphere( int state, int subunit ) const
{
	BoundingSphere transformed_bs = sc_bounding_spheres_[ state ];
	transformed_bs.first = get_on_the_fly_owner()->symmetric_transform( subunit ) * transformed_bs.first;
	return transformed_bs;
}

/// @brief Returns a bounding sphere for the backbone on a particular subunit.
SymmOnTheFlyNode::BoundingSphere
SymmOnTheFlyNode::bb_bounding_sphere( int subunit ) const
{
	BoundingSphere transformed_bs = bb_bounding_sphere_;
	transformed_bs.first = get_on_the_fly_owner()->symmetric_transform( subunit ) * transformed_bs.first;
	return transformed_bs;
}


//-------------------------------------------------------------//

SymmOnTheFlyEdge::SymmOnTheFlyEdge(
	InteractionGraphBase * owner,
	int first_node_ind,
	int second_node_ind
)
:
	FixedBBEdge( owner, first_node_ind, second_node_ind ),
	restypegroup_adjacency_(
		get_otf_owner()->get_num_restype_groups(),
		get_otf_owner()->get_num_restype_groups(),
		get_otf_owner()->symm_info()->subunits(),
		2,
		(unsigned char) 0  // initial value -- set everything to "false"
	),
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


SymmOnTheFlyEdge::~SymmOnTheFlyEdge(){}

void
SymmOnTheFlyEdge::add_ProCorrection_values(
	int node_not_necessarily_proline,
	int state,
	core::PackerEnergy bb_nonprobb_E,
	core::PackerEnergy bb_probb_E,
	core::PackerEnergy sc_nonprobb_E,
	core::PackerEnergy sc_probb_E
)
{
	int const which_node = node_not_necessarily_proline == get_node_index( 0 ) ? 0 : 1;

	//std::cout << "  SymmOnTheFlyEdge::add_ProCorrection_values" << get_first_node_ind() << " " << get_second_node_ind() << " " << node_not_necessarily_proline << " " << state << " : " << sc_probb_E + 0.5 * bb_probb_E - (sc_nonprobb_E + 0.5 * bb_nonprobb_E) << std::endl;

	proline_corrections_[ which_node ][ state ] +=
		sc_probb_E + 0.5 * bb_probb_E -
		(sc_nonprobb_E + 0.5 * bb_nonprobb_E);
}


unsigned int
SymmOnTheFlyEdge::count_dynamic_memory() const
{
	unsigned int total_memory_usage = EdgeBase::count_dynamic_memory();

	total_memory_usage += proline_corrections_[ 0 ].size() * sizeof( core::PackerEnergy );
	total_memory_usage += proline_corrections_[ 1 ].size() * sizeof( core::PackerEnergy );

	return total_memory_usage;
}

void
SymmOnTheFlyEdge::set_residues_adjacent_for_subunit_pair(
	int asu_node_index,
	int other_node_subunit
)
{
	if ( asu_node_index == 2 && other_node_subunit == 1 ) return;

	//std::cout << "SymmOnTheFlyEdge " << get_first_node_ind() << " " << get_second_node_ind() << " set_residues_adjacent_for_subunit_pair " << asu_node_index << " " << other_node_subunit << std::endl;

	int const asu_index = 1; // ASSUMPTION, the asymmetric unit is unit 1!
	int const nres_asu = get_otf_owner()->symm_info()->num_independent_residues();
	Size const score_multiply =
			Size(get_otf_owner()->symm_info()->score_multiply( 1, (other_node_subunit-1)*nres_asu + 1 ));
	if ( score_multiply == 0 ) return;

	// OK: now iterate across the rotamer representatives for this pair and figure out their cbeta distances
	// and ask whether they should be considered adjacent

	int other_node_index = ( asu_node_index == 1 ) ? 2 : 1;
	SymmOnTheFlyNode const * asu_node = get_otf_node( asu_node_index - 1 );
	SymmOnTheFlyNode const * other_node = get_otf_node( other_node_index - 1 );
	Real const sfxn_reach = get_otf_owner()->score_function().max_atomic_interaction_cutoff();
	for ( int ii = 1; ii <= asu_node->get_num_restype_groups(); ++ii ) {
		int iistate = asu_node->get_state_offset_for_restype_group( ii ) + 1;
		core::conformation::Residue const & iires = asu_node->get_rotamer( iistate, asu_index );
		for ( int jj = 1; jj <= other_node->get_num_restype_groups(); ++jj ) {
			int jjstate = other_node->get_state_offset_for_restype_group( jj ) + 1;
			core::conformation::Residue const & jjres = other_node->get_rotamer( jjstate, other_node_subunit );
			Real ii_jj_radii_sum = iires.nbr_radius() + jjres.nbr_radius();
			Real d2 = iires.xyz( iires.nbr_atom() ).distance_squared( jjres.xyz( jjres.nbr_atom() ));
			//std::cout << "ii: " << ii << " " << iires.name() << " jj: " << jj << " " << jjres.name() << " d2: " << d2 << " sfxn_reach: " << sfxn_reach << " ii_jj_radii_sum: " << ii_jj_radii_sum << ", squared: " << (sfxn_reach + ii_jj_radii_sum )*(sfxn_reach + ii_jj_radii_sum ) << std::endl;
			if ( d2 < (sfxn_reach + ii_jj_radii_sum )*(sfxn_reach + ii_jj_radii_sum ) ) {
				//std::cout << "declared neighbors: " << get_otf_owner()->symm_info()->score_multiply( 1, (other_node_subunit-1)*nres_asu + 1 ) << std::endl;
				restypegroup_adjacency_( jj, ii, other_node_subunit, asu_node_index ) = score_multiply;
			}
		}
	}

}

unsigned char
SymmOnTheFlyEdge::residues_adjacent_for_subunit_pair(
	int which_node, // 1 or 2
	int other_node_subunit, // subunit for othernode
	int whichnode_restypegroup, // restype group for first node
	int othernode_restypegroup // restype group for the other node
) const
{
	return restypegroup_adjacency_( othernode_restypegroup, whichnode_restypegroup, other_node_subunit, which_node );
}

//-------------------------------------------------------------//

SymmOnTheFlyInteractionGraph::SymmOnTheFlyInteractionGraph(int num_nodes )
:
	FixedBBInteractionGraph( num_nodes ),
	num_restype_groups_( 0 )
{}

SymmOnTheFlyInteractionGraph::~SymmOnTheFlyInteractionGraph()
{}

bool
SymmOnTheFlyInteractionGraph::distinguish_backbone_and_sidechain_for_node( int node_ind ) const
{
	return get_on_the_fly_node( node_ind )->distinguish_backbone_and_sidechain();
}

void
SymmOnTheFlyInteractionGraph::distinguish_backbone_and_sidechain_for_node( int node_ind, bool setting )
{
	get_on_the_fly_node( node_ind )->distinguish_backbone_and_sidechain( setting );
}


void
SymmOnTheFlyInteractionGraph::initialize(
	rotamer_set::RotamerSetsBase const & rot_sets_base
)
{
	rotamer_set::RotamerSets const & rot_sets( static_cast< rotamer_set::RotamerSets const & > (rot_sets_base) );

	// determine max # of residue types
	Size max_nrestypes = 0;
	for ( Size ii = 1; ii <= rot_sets.nmoltenres(); ++ii ) {
		Size ii_nrestypes =  rot_sets.rotamer_set_for_moltenresidue( ii )->get_n_residue_types();
		if ( ii_nrestypes > max_nrestypes ) max_nrestypes = ii_nrestypes;
	}

	//"aa types" means "distinct groups of rotamers" -- this ig has no idea
	// what an amino acid is or why they might be different from one another
	num_restype_groups_ = max_nrestypes;

	for ( Size ii = 1; ii <= rot_sets.nmoltenres(); ++ii ) {
		Size const ii_num_states = rot_sets.rotamer_set_for_moltenresidue( ii )->num_rotamers();
		set_num_states_for_node( ii, ii_num_states );
		get_on_the_fly_node( ii )->set_rotamers( rot_sets.rotamer_set_for_moltenresidue( ii ) );
	}

}

void
SymmOnTheFlyInteractionGraph::set_score_function( ScoreFunction const & sfxn )
{
	score_function_ = sfxn.clone();
	if ( pose_ ) (*score_function_)(*pose_); // rescore the pose with the input score function
}

void
SymmOnTheFlyInteractionGraph::set_pose( pose::Pose const & pose )
{
	pose_ = pose::PoseOP( new pose::Pose( pose ) );
	if ( score_function_ ) (*score_function_)(*pose_);

	using conformation::symmetry::SymmetricConformation;
	using id::AtomID;

	SymmetricConformation & symm_conf ( dynamic_cast<SymmetricConformation &> ( pose_->conformation()) );
	symm_info_ = symm_conf.Symmetry_Info();

	// ok, now figure out the HTs to translate from the reference frame to all of
	// symmetric clone frames
	Size const nsubunits = symm_info_->subunits();
	Size const asu_size = symm_info_->num_independent_residues();
	assert( nsubunits > 0 );
	symmetric_transforms_.resize( nsubunits );
	symmetric_transforms_[1].set_identity();

	// 1. find a residue with at least three atoms
	Size orientation_residue = 0;
	for ( Size ii = 1; ii <= asu_size; ++ii ) {
		if ( pose.residue( ii ).natoms() >= 3 ) {
			orientation_residue = ii;
			break;
		}
	}
	if ( orientation_residue == 0 ) {
		/// there are other residues that could be used to generate these
		/// orientation, but that requires a better understanding of the fold tree
		utility_exit_with_message( "Failed to find a residue in the pose with at least three atoms" );
	}

	HTReal frame1(
		pose.xyz( AtomID( 1, orientation_residue )),
		pose.xyz( AtomID( 2, orientation_residue )),
		pose.xyz( AtomID( 3, orientation_residue )) );
	HTReal invframe1 = frame1.inverse();
	for ( Size ii = 2; ii <= nsubunits; ++ii ) {
		Size iiresidue = orientation_residue + (ii-1)*asu_size;
		HTReal iiframe(
			pose.xyz( AtomID( 1, iiresidue )),
			pose.xyz( AtomID( 2, iiresidue )),
			pose.xyz( AtomID( 3, iiresidue )) );
		symmetric_transforms_[ ii ] = iiframe * invframe1;
	}

}

conformation::symmetry::SymmetryInfoCOP
SymmOnTheFlyInteractionGraph::symm_info() const {
	return symm_info_;
}


void
SymmOnTheFlyInteractionGraph::zero_one_body_energy_for_node_state(
	int node_ind,
	int state
)
{
	get_on_the_fly_node( node_ind )->zero_one_body_energy( state );
}

void
SymmOnTheFlyInteractionGraph::add_to_one_body_energy_for_node_state(
	int node_ind,
	int state,
	core::PackerEnergy one_body_energy
)
{
	get_on_the_fly_node( node_ind )->add_to_one_body_energy( state, one_body_energy);
}

void
SymmOnTheFlyInteractionGraph::set_one_body_energy_for_node_state(
	int node_ind,
	int state,
	core::PackerEnergy one_body_energy
)
{
	get_on_the_fly_node( node_ind )->set_one_body_energy( state, one_body_energy);
}


core::PackerEnergy
SymmOnTheFlyInteractionGraph::get_one_body_energy_for_node_state( int node, int state)
{
	return get_on_the_fly_node( node )->get_one_body_energy( state );
}

void
SymmOnTheFlyInteractionGraph::set_residues_adjacent_for_subunit_pair_for_edge(
	int node1,
	int node2,
	int asu_node_index,
	int other_node_subunit
)
{
	SymmOnTheFlyEdge * edge = static_cast< SymmOnTheFlyEdge * > ( find_edge( node1, node2 ) );
	if ( edge ) {
		edge->set_residues_adjacent_for_subunit_pair( asu_node_index, other_node_subunit );
	}
}

void
SymmOnTheFlyInteractionGraph::reset_rpe_calculations_count()
{
	num_rpe_calcs_ = 0;
}


Size
SymmOnTheFlyInteractionGraph::get_num_rpe_calculations_count() const
{
	return num_rpe_calcs_;
}


//void
//SymmOnTheFlyInteractionGraph::set_sparse_aa_info_for_edge(
//	int node1,
//	int node2,
//	FArray2_bool const & sparse_conn_info
//)
//{
//	SymmOnTheFlyEdge* edge = (SymmOnTheFlyEdge*) find_edge( node1, node2 );
//	if (edge != 0) {
//		edge->set_sparse_aa_info( sparse_conn_info );
//	}
//}

void
SymmOnTheFlyInteractionGraph::add_ProCorrection_values_for_edge(
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
	SymmOnTheFlyEdge* edge = (SymmOnTheFlyEdge*) find_edge( node1, node2 );
	if (edge != 0)
	{
		edge->add_ProCorrection_values(
			node_not_neccessarily_proline, state,
			bb_nonprobb_E, bb_probb_E, sc_nonprobb_E, sc_probb_E );
	}
}

void
SymmOnTheFlyInteractionGraph::note_short_range_interactions_exist_for_edge(
	int node1,
	int node2
)
{
	SymmOnTheFlyEdge* edge = (SymmOnTheFlyEdge*) find_edge( node1, node2 );
	if (edge != 0) {
		edge->note_short_range_interactions_exist();
	}
}

void
SymmOnTheFlyInteractionGraph::note_long_range_interactions_exist_for_edge(
	int node1,
	int node2
)
{
	SymmOnTheFlyEdge* edge = (SymmOnTheFlyEdge*) find_edge( node1, node2 );
	if (edge != 0) {
		edge->note_long_range_interactions_exist( );
	}
}


unsigned int
SymmOnTheFlyInteractionGraph::count_dynamic_memory() const
{
	unsigned int total_memory = InteractionGraphBase::count_dynamic_memory();
	return total_memory;
}

} // namespace interaction_graph
} // namespace pack
} // namespace core
