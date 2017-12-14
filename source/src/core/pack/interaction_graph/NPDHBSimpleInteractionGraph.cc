// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/NPDHBSimpleInteractionGraph.cc
/// @brief
/// @author Andrew Leaver-Fay (aleavefay@gmail.com)

// Unit headers
#include <core/pack/interaction_graph/NPDHBSimpleInteractionGraph.hh>

// Package headers
#include <core/pack/interaction_graph/NPDHBondInteractionGraph.hh>

// Project headers
#include <utility/graph/Graph.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/util.hh>

// Basic headers
#include <basic/Tracer.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/assert.hh>

namespace core {
namespace pack {
namespace interaction_graph {

using namespace utility::graph;
using namespace scoring;

void
create_hbonds_one_way(
	scoring::hbonds::HBondDatabase const & database,
	scoring::hbonds::HBondOptions const & hbondoptions,
	scoring::hbonds::HBondSet const & hbset,
	NPDHBSimpleInteractionGraph & ig,
	utility::vector1< char > & hbonding_to_res,
	conformation::Residue const & acc_res,
	utility::vector1< NPDHBondOP > & acc_hbonds,
	utility::vector1< utility::vector1< NPDHBondOP > > & acc_atom_hbonds,
	conformation::Residue const & don_res,
	utility::vector1< NPDHBondOP > & don_hbonds,
	utility::vector1< utility::vector1< NPDHBondOP > > & don_atom_hbonds
)
{
	for ( Size ii_at : acc_res.accpt_pos() ) {
		//Size ii_at = acc_res.accpt_pos()[ ii ];
		Vector const & iixyz = acc_res.xyz( ii_at );
		for ( Size jj_at : don_res.Hpos_polar() ) {
			//Size jj_at = don_res.Hpos_polar()[ jj ];
			Vector const & jjxyz = don_res.xyz( jj_at );
			if ( iixyz.distance_squared( jjxyz ) > core::scoring::hbonds::MAX_R2 ) continue;
			Real e = core::scoring::hbonds::hb_energy( database, hbondoptions, hbset, acc_res, ii_at, don_res, jj_at );
			if ( e >= core::scoring::hbonds::MAX_HB_ENERGY ) continue;
			NPDHBondOP hb = ig.unused_hbond();
			hb->energy_ = e;
			hb->don_wt_ = -1234;
			hb->acc_wt_ = -1234;
			hb->don_wt_alt_ = -1234;
			hb->acc_wt_alt_ = -1234;
			hb->don_rsd_ = don_res.seqpos();
			hb->acc_rsd_ = acc_res.seqpos();
			hb->don_atm_ = jj_at;
			hb->acc_atm_ = ii_at;

			// score function weight should be included in the donor and acceptor weight calculation
			scoring::hbonds::HBEvalTuple hb_eval_tup( don_res.atom_base( jj_at ), don_res, ii_at, acc_res );
			hb->sfxn_wt_ = ig.npd_hb_weight( hb_eval_tup.eval_type(), don_res.seqpos() == acc_res.seqpos() );

			acc_hbonds.push_back( hb );
			if ( acc_res.seqpos() != don_res.seqpos() ) don_hbonds.push_back( hb );
			acc_atom_hbonds[ ii_at ].push_back( hb );
			don_atom_hbonds[ jj_at ].push_back( hb );
			hbonding_to_res[ don_res.seqpos() ] = 1;
			hbonding_to_res[ acc_res.seqpos() ] = 1;
		}
	}
}



//static basic::Tracer TR( "core.pack.interaction_graph.NPDHBSimpleInteractionGraph" );

NPDHBSimpleNode::NPDHBSimpleNode( Graph* owner, Size node_id ):
	parent( owner, node_id  ),
	seqpos_( node_id )
{
	initialize();
}

NPDHBSimpleNode::~NPDHBSimpleNode() = default;


void
NPDHBSimpleNode::commit_change_npd()
{
	parent::commit_change();

	curr_hbonds_.swap( alt_hbonds_ );
	curr_atom_hbonds_.swap( alt_atom_hbonds_ );
	for ( auto const & hb : curr_hbonds_ ) {
		hb->don_wt_ = hb->don_wt_alt_;
		hb->acc_wt_ = hb->acc_wt_alt_;
	}

	// return the (now) alternate hbonds to the owner
	for ( auto const & hb : alt_hbonds_ ) {
		npdhb_owner()->return_hbond_to_queue( hb );
	}
	alt_hbonds_.resize( 0 );

	for ( auto iter = edge_list_begin(), iter_end = edge_list_end(); iter != iter_end; ++iter ) {
		NPDHBSimpleNode * nbr( npdhb_cast( (*iter)->get_other_node( get_node_index() ) ) );
		nbr->acknowledge_neighbors_substitution();
	}
}

/// @brief Copy the alternate energies to the current energies for this node
/// and its incident edges.
void
NPDHBSimpleNode::commit_change_no_res_pointer_update()
{
	parent::commit_change_no_res_pointer_update();

	curr_hbonds_.swap( alt_hbonds_ );
	curr_atom_hbonds_.swap( alt_atom_hbonds_ );
	for ( auto const & hb : curr_hbonds_ ) {
		hb->don_wt_ = hb->don_wt_alt_;
		hb->acc_wt_ = hb->acc_wt_alt_;
	}

	// return the (now) alternate hbonds to the owner
	for ( auto const & hb : alt_hbonds_ ) {
		npdhb_owner()->return_hbond_to_queue( hb );
	}
	alt_hbonds_.resize( 0 );

	for ( auto iter = edge_list_begin(), iter_end = edge_list_end(); iter != iter_end; ++iter ) {
		NPDHBSimpleNode * nbr( npdhb_cast( (*iter)->get_other_node( get_node_index() ) ) );
		nbr->acknowledge_neighbors_substitution();
	}

}

void
NPDHBSimpleNode::reject_change_npd( conformation::ResidueCOP res, basic::datacache::BasicDataCache & residue_data_cache )
{
	parent::reject_change( res, residue_data_cache );
	// return the alternate hbonds to the owner
	for ( auto const & hb : alt_hbonds_ ) {
		npdhb_owner()->return_hbond_to_queue( hb );
	}
	alt_hbonds_.resize( 0 );
}


/// @brief Compute the total energy for the current state assignment
Real NPDHBSimpleNode::get_curr_upper_hbond_energies()
{
	Real totalE( 0 );
	for ( auto const & hb : curr_hbonds_ ) {
		debug_assert( hb->don_rsd_ == seqpos_ || hb->acc_rsd_ == seqpos_ );
		if ( (hb->don_rsd_ == seqpos_ && seqpos_ < hb->acc_rsd_) || (hb->acc_rsd_ == seqpos_ && seqpos_ < hb->don_rsd_ ) ) {
			totalE += hb->energy_ * hb->sfxn_wt_ * hb->acc_wt_ * hb->don_wt_;
		}
	}
	return  totalE;
}

/// @brief Compute the change in NPDHBond energy for the given (previously supplied!)
/// rotamer substitution
Real NPDHBSimpleNode::get_npdhb_deltaE_for_substitution()
{
	//npdhb_owner()->reset_hbondind_to_res();
	utility::vector1< char > & hbonding_to_res( npdhb_owner()->hbonding_to_res() );
	std::fill( hbonding_to_res.begin(), hbonding_to_res.end(), 0 );
	for ( auto const & hb : curr_hbonds_ ) {
		debug_assert( hb->don_rsd_ == seqpos_ || hb->acc_rsd_ == seqpos_ );
		hbonding_to_res[ hb->don_rsd_ ] = 1;
		hbonding_to_res[ hb->acc_rsd_ ] = 1;
	}
	alt_hbonds_.clear();
	if ( alt_atom_hbonds_.size() < get_alternate_ref().natoms() ) {
		alt_atom_hbonds_.resize( get_alternate_ref().natoms() );
		for ( Size ii = 1; ii <= alt_atom_hbonds_.size(); ++ii ) {
			alt_atom_hbonds_[ ii ].reserve( 5 );
		}
	}
	for ( Size ii = 1; ii <= alt_atom_hbonds_.size(); ++ii ) {
		alt_atom_hbonds_[ ii ].clear();
	}

	for ( auto iter = edge_list_begin(), iter_end = edge_list_end(); iter != iter_end; ++iter ) {
		NPDHBSimpleNode * nbr( npdhb_cast( (*iter)->get_other_node( get_node_index() ) ) );
		nbr->prepare_for_neighbors_substitution( this );
	}

	for ( auto iter = edge_list_begin(), iter_end = edge_list_end(); iter != iter_end; ++iter ) {
		NPDHBSimpleNode * nbr( npdhb_cast( (*iter)->get_other_node( get_node_index() ) ) );
		nbr->find_hbs_for_nbrs_alt_state_step1( this );
	}

	// Intra residue hbonds, if requested
	// TEMP! Assume that the residue is a protein residue
	// TO DO: figure out how to get this consistent w/ rna & DNA & protein
	if ( ! npdhb_owner()->hbond_options().exclude_intra_res_protein() ) {
		NPDHBSimpleInteractionGraph & ig( *npdhb_owner() );
		scoring::hbonds::HBondDatabase const & database( ig.hbond_database() );
		scoring::hbonds::HBondOptions const & hbondoptions( ig.hbond_options() );
		scoring::hbonds::NPDHBondSet const & hbset( ig.npd_hbond_set() );

		create_hbonds_one_way( database, hbondoptions, hbset, ig, hbonding_to_res,
			get_alternate_ref(), alt_hbonds_, alt_atom_hbonds_,
			get_alternate_ref(), alt_hbonds_, alt_atom_hbonds_
		);
	}

	compute_alt_weights_for_hbonds( false );

	Real deltaE = 0;
	for ( auto iter = edge_list_begin(), iter_end = edge_list_end(); iter != iter_end; ++iter ) {
		NPDHBSimpleNode * nbr( npdhb_cast( (*iter)->get_other_node( get_node_index() ) ) );
		if ( hbonding_to_res[ nbr->seqpos_ ] ) {
			deltaE += nbr->find_hbs_for_nbrs_alt_state_step2( this );
		}
	}

	// iterate across hbs for current and alternate and sum them up
	Real curr_hb( 0 ), alt_hb( 0 );
	for ( auto const & hb : curr_hbonds_ ) {
		curr_hb += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_ * hb->acc_wt_;
	}
	for ( auto const & hb : alt_hbonds_ ) {
		alt_hb += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_alt_ * hb->acc_wt_alt_;
	}
	// Return the -deltaE
	// std::cout << "    NPD HB deltaE: " << ( deltaE + alt_hb - curr_hb ) << std::endl;
	//return -1 * npdhb_owner()->npd_hbond_weight() * ( deltaE + alt_hb - curr_hb );
	return -1 * ( deltaE + alt_hb - curr_hb );
}

void NPDHBSimpleNode::prepare_for_neighbors_substitution( NPDHBSimpleNode * nbr )
{
	alt_hbonds_.clear();
	if ( alt_atom_hbonds_.size() < get_current_ref().natoms() ) {
		alt_atom_hbonds_.resize( get_current_ref().natoms() );
		for ( Size ii = 1; ii <= alt_atom_hbonds_.size(); ++ii ) {
			alt_atom_hbonds_[ ii ].reserve( 5 );
		}
	}
	for ( Size ii = 1; ii <= alt_atom_hbonds_.size(); ++ii ) {
		alt_atom_hbonds_[ ii ].clear();
	}
	for ( auto const & hb : curr_hbonds_ ) {
		debug_assert( hb->don_rsd_ == seqpos_ || hb->acc_rsd_ == seqpos_ );
		if ( hb->don_rsd_ != nbr->seqpos_ && hb->acc_rsd_ != nbr->seqpos_ ) {
			hb->don_wt_alt_ = hb->don_wt_;
			hb->acc_wt_alt_ = hb->acc_wt_;

			alt_hbonds_.push_back( hb );
			if ( hb->don_rsd_ == seqpos_ ) {
				alt_atom_hbonds_[ hb->don_atm_ ].push_back( hb );
			} else {
				alt_atom_hbonds_[ hb->acc_atm_ ].push_back( hb );
			}
		}
	}
}


void NPDHBSimpleNode::find_hbs_for_nbrs_alt_state_step1( NPDHBSimpleNode * nbr )
{
	NPDHBSimpleInteractionGraph & ig( *npdhb_owner() );
	utility::vector1< char > & hbonding_to_res( ig.hbonding_to_res() );
	scoring::hbonds::HBondDatabase const & database( ig.hbond_database() );
	scoring::hbonds::HBondOptions const & hbondoptions( ig.hbond_options() );
	scoring::hbonds::NPDHBondSet const & hbset( ig.npd_hbond_set() );

	create_hbonds_one_way( database, hbondoptions, hbset, ig, hbonding_to_res,
		get_current_ref(), alt_hbonds_, alt_atom_hbonds_,
		nbr->get_alternate_ref(), nbr->alt_hbonds_, nbr->alt_atom_hbonds_
	);

	create_hbonds_one_way( database, hbondoptions, hbset, ig, hbonding_to_res,
		nbr->get_alternate_ref(), nbr->alt_hbonds_, nbr->alt_atom_hbonds_,
		get_current_ref(), alt_hbonds_, alt_atom_hbonds_
	);

	compute_alt_weights_for_hbonds( true );

}

Real NPDHBSimpleNode::find_hbs_for_nbrs_alt_state_step2( NPDHBSimpleNode * nbr )
{
	NPDHBSimpleInteractionGraph & ig( *npdhb_owner() );
	utility::vector1< char > & hbonding_to_res( ig.hbonding_to_res() );

	Real curr_hbE( 0 ), alt_hbE( 0 );
	for ( auto const & hb : curr_hbonds_ ) {
		debug_assert( hb->don_rsd_ == seqpos_ || hb->acc_rsd_ == seqpos_ );
		Size other = hb->don_rsd_ == seqpos_ ? hb->acc_rsd_ : hb->don_rsd_;
		if ( other == nbr->seqpos_ ) continue;
		if ( ! hbonding_to_res[ other ] || seqpos_ <= other ) {
			curr_hbE += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_ * hb->acc_wt_;
		}
	}

	for ( auto const & hb : alt_hbonds_ ) {
		debug_assert( hb->don_rsd_ == seqpos_ || hb->acc_rsd_ == seqpos_ );
		Size other = hb->don_rsd_ == seqpos_ ? hb->acc_rsd_ : hb->don_rsd_;
		if ( other == nbr->seqpos_ ) continue;
		if ( ! hbonding_to_res[ other ] || seqpos_ <= other ) {
			alt_hbE += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_alt_ * hb->acc_wt_alt_;
		}
	}
	return alt_hbE - curr_hbE;
}

void NPDHBSimpleNode::acknowledge_neighbors_substitution()
{
	curr_hbonds_.swap( alt_hbonds_ );
	curr_atom_hbonds_.swap( alt_atom_hbonds_ );
	for ( auto const & hb : curr_hbonds_ ) {
		hb->don_wt_ = hb->don_wt_alt_;
		hb->acc_wt_ = hb->acc_wt_alt_;
	}
}

void NPDHBSimpleNode::compute_alt_weights_for_hbonds(
	bool curr_state
)
{
	if ( ( curr_state && get_current() == nullptr ) || ( ! curr_state && get_alternate() == nullptr ) ) return;

	conformation::Residue const & res( curr_state ? get_current_ref() : get_alternate_ref() );
	compute_alt_weights_for_npd_hbonds( res, alt_atom_hbonds_, tmp_energies_, tmp_weights_ );
}

void
NPDHBSimpleNode::reset_hbs()
{
	curr_hbonds_.clear();
	if ( curr_atom_hbonds_.size() < get_current_ref().natoms() ) {
		curr_atom_hbonds_.resize( get_current_ref().natoms() );
	}
	for ( Size ii = 1; ii <= curr_atom_hbonds_.size(); ++ii ) {
		curr_atom_hbonds_[ ii ].clear();
		curr_atom_hbonds_[ ii ].reserve( 5 );
	}

	alt_hbonds_.clear();
	for ( Size ii = 1; ii <= alt_atom_hbonds_.size(); ++ii ) {
		alt_atom_hbonds_[ ii ].clear();
		alt_atom_hbonds_[ ii ].reserve( 5 );
	}
}


void
NPDHBSimpleNode::find_curr_hbonds_to_upper_neighbors()
{
	NPDHBSimpleInteractionGraph & ig( *npdhb_owner() );
	utility::vector1< char > & hbonding_to_res( ig.hbonding_to_res() );
	scoring::hbonds::HBondDatabase const & database( ig.hbond_database() );
	scoring::hbonds::HBondOptions const & hbondoptions( ig.hbond_options() );
	scoring::hbonds::NPDHBondSet const & hbset( ig.npd_hbond_set() );

	// this should only be called once per node at the very beginning of initialization
	for ( auto iter = upper_edge_list_begin(), iter_end = upper_edge_list_end();
			iter != iter_end; ++iter ) {
		NPDHBSimpleNode * nbr( npdhb_cast( (*iter)->get_other_node( get_node_index() ) ) );
		create_hbonds_one_way( database, hbondoptions, hbset, ig, hbonding_to_res,
			nbr->get_current_ref(), nbr->curr_hbonds_, nbr->curr_atom_hbonds_,
			get_current_ref(), curr_hbonds_, curr_atom_hbonds_
		);

		create_hbonds_one_way( database, hbondoptions, hbset, ig, hbonding_to_res,
			get_current_ref(), curr_hbonds_, curr_atom_hbonds_,
			nbr->get_current_ref(), nbr->curr_hbonds_, nbr->curr_atom_hbonds_
		);
	}

	// Intra residue hbonds, if requested
	// TEMP! Assume that the residue is a protein residue
	// TO DO: figure out how to get this consistent w/ rna & DNA & protein
	if ( ! npdhb_owner()->hbond_options().exclude_intra_res_protein() ) {
		create_hbonds_one_way( database, hbondoptions, hbset, ig, hbonding_to_res,
			get_current_ref(), curr_hbonds_, curr_atom_hbonds_,
			get_current_ref(), curr_hbonds_, curr_atom_hbonds_
		);
	}

}

void
NPDHBSimpleNode::calc_curr_hbond_weights()
{
	compute_alt_weights_for_npd_hbonds( get_current_ref(), curr_atom_hbonds_, tmp_energies_, tmp_weights_ );
	for ( auto const & hb : curr_hbonds_ ) {
		hb->don_wt_ = hb->don_wt_alt_;
		hb->acc_wt_ = hb->acc_wt_alt_;
	}
}


NPDHBSimpleNode * NPDHBSimpleNode::npdhb_cast( utility::graph::Node * node ) const
{
	return static_cast< NPDHBSimpleNode * > ( node );
}

NPDHBSimpleNode const * NPDHBSimpleNode::npdhb_cast( utility::graph::Node const * node ) const
{
	return static_cast< NPDHBSimpleNode const * > ( node );
}

NPDHBSimpleEdge * NPDHBSimpleNode::npdhb_cast( utility::graph::Edge * edge ) const
{
	return static_cast< NPDHBSimpleEdge * > ( edge );
}

NPDHBSimpleEdge const * NPDHBSimpleNode::npdhb_cast( utility::graph::Edge const * edge ) const
{
	return static_cast< NPDHBSimpleEdge const * > ( edge );
}

NPDHBSimpleInteractionGraph const * NPDHBSimpleNode::npdhb_owner() const
{
	return static_cast< NPDHBSimpleInteractionGraph const * > ( get_owner() );
}

NPDHBSimpleInteractionGraph * NPDHBSimpleNode::npdhb_owner()
{
	return static_cast< NPDHBSimpleInteractionGraph * > ( get_owner() );
}


void
NPDHBSimpleNode::initialize()
{
	parent::initialize();
}

NPDHBSimpleEdge::NPDHBSimpleEdge( Graph* owner, Size res1, Size res2 ):
	parent( owner, res1, res2 )
{}

NPDHBSimpleEdge::~NPDHBSimpleEdge() = default;

NPDHBSimpleInteractionGraph *
NPDHBSimpleEdge::get_npdhb_simple_ig_owner()
{
	return static_cast< NPDHBSimpleInteractionGraph * > ( get_owner() );
}

NPDHBSimpleInteractionGraph const *
NPDHBSimpleEdge::get_npdhb_simple_ig_owner() const
{
	return static_cast< NPDHBSimpleInteractionGraph const * > ( get_owner() );
}


NPDHBSimpleInteractionGraph::NPDHBSimpleInteractionGraph() :
	parent(),
	npd_hbond_weight_( 0 )
{}

NPDHBSimpleInteractionGraph::~NPDHBSimpleInteractionGraph() = default;

void NPDHBSimpleInteractionGraph::set_scorefunction( scoring::ScoreFunction const & sfxn )
{
	parent::set_scorefunction( sfxn );
	npd_hbond_weight_ = std::max( {
		sfxn.weights()[ scoring::npd_hbond_sr_bb ],
		sfxn.weights()[ scoring::npd_hbond_lr_bb ],
		sfxn.weights()[ scoring::npd_hbond_bb_sc ],
		sfxn.weights()[ scoring::npd_hbond_sc ],
		sfxn.weights()[ scoring::npd_hbond ] } );
}

void
NPDHBSimpleInteractionGraph::initialize( pose::Pose const & pose )
{
	parent::initialize( pose ); // calls set_pose_no_initialize internally
	hbonding_to_res_.resize( pose.total_residue(), 0 );
	setup_after_edge_addition();
}

void
NPDHBSimpleInteractionGraph::set_pose_no_initialize( pose::Pose const & pose )
{
	using scoring::EnergiesCacheableDataType::NPD_HBOND_SET;
	using namespace scoring::hbonds;

	npd_hbond_set_ = utility::pointer::static_pointer_cast< NPDHBondSet const > ( pose.energies().data().get_const_ptr( NPD_HBOND_SET ) );
	if ( ! npd_hbond_set_ ) {
		npd_hbond_set_ = NPDHBondSetOP( new NPDHBondSet );
	}

	hbond_options_ = HBondOptionsOP( new HBondOptions( npd_hbond_set_->hbond_options() ));
	hbond_database_ = HBondDatabase::get_database(hbond_options_->params_database_tag());

	parent::set_pose_no_initialize( pose );
	hbonding_to_res_.resize( pose.total_residue(), 0 );
}

void
NPDHBSimpleInteractionGraph::setup_after_edge_addition()
{
	for ( Size ii = 1; ii <= num_nodes(); ++ii ) npdhb_cast( get_node( ii ) )->reset_hbs();
	for ( Size ii = 1; ii <= num_nodes(); ++ii ) npdhb_cast( get_node( ii ) )->find_curr_hbonds_to_upper_neighbors();
	for ( Size ii = 1; ii <= num_nodes(); ++ii ) npdhb_cast( get_node( ii ) )->calc_curr_hbond_weights();

}


void
NPDHBSimpleInteractionGraph::commit_change( Size node_id )
{
	static_cast< NPDHBSimpleNode * >(get_node( node_id ))->commit_change_npd();
}

void
NPDHBSimpleInteractionGraph::reject_change( Size node_id, conformation::ResidueCOP res, basic::datacache::BasicDataCache & residue_data_cache )
{
	static_cast< NPDHBSimpleNode * >( get_node( node_id ) )->reject_change_npd( res, residue_data_cache );
}

scoring::hbonds::HBondDatabase const &
NPDHBSimpleInteractionGraph::hbond_database() const
{
	return *hbond_database_;
}

scoring::hbonds::HBondOptions const &
NPDHBSimpleInteractionGraph::hbond_options() const
{
	return *hbond_options_;
}

scoring::hbonds::NPDHBondSet const &
NPDHBSimpleInteractionGraph::npd_hbond_set() const
{
	return *npd_hbond_set_;
}

/// @details Note, this function returns (currE - altE) which represents
/// the negative of the change in energy for the substition
Real
NPDHBSimpleInteractionGraph::consider_substitution(
	Size node_id,
	conformation::ResidueCOP new_state,
	basic::datacache::BasicDataCache & residue_data_cache
) {
	Real deltaE = parent::consider_substitution( node_id, new_state, residue_data_cache );

	auto* node = static_cast< NPDHBSimpleNode *>(get_node( node_id ));
	deltaE += node->get_npdhb_deltaE_for_substitution();
	return deltaE;
}


Real
NPDHBSimpleInteractionGraph::total_energy() {
	Real total_energy = parent::total_energy(); // get PD total from parent
	for ( Size ii = 1; ii <= pose().total_residue(); ++ii ) {
		total_energy += get_npdhb_simple_node( ii )->get_curr_upper_hbond_energies();
	}
	return total_energy;

}

utility::vector1< char > & NPDHBSimpleInteractionGraph::hbonding_to_res()
{
	return hbonding_to_res_;
}


NPDHBondOP NPDHBSimpleInteractionGraph::unused_hbond()
{
	if ( hbonds_queue_.empty() ) {
		return NPDHBondOP( new NPDHBond );
	}
	NPDHBondOP next = hbonds_queue_.front();
	hbonds_queue_.pop_front();
	return next;
}

void NPDHBSimpleInteractionGraph::return_hbond_to_queue( NPDHBondOP const & hbond )
{
	// mark the hbond as not part of any residue; this is mostly for debugging purposes:
	// if an hbond has been recycled but is still being used in any calculation, then something
	// within the NPDHBondInteractionGraph has gone wrong. There are assertions in the code to
	// make sure that the hbond that's about to be used belongs to the residue that's using it
	// so marking the residue as 0 will trip one of those assertions
	hbond->don_rsd_ = hbond->acc_rsd_ = 0;

	hbonds_queue_.push_back( hbond );
}


Real NPDHBSimpleInteractionGraph::npd_hb_weight(
	scoring::hbonds::HBEvalType eval_type,
	bool intra_res
)
{
	return scoring::hbonds::npd_hb_eval_type_weight( eval_type, scorefunction().weights(), intra_res );
}

void NPDHBSimpleInteractionGraph::delete_edge( Edge * edge )
{
	delete edge;
}

Node* NPDHBSimpleInteractionGraph::create_new_node( platform::Size node_index ){
	return new NPDHBSimpleNode( this, node_index );
}

Edge* NPDHBSimpleInteractionGraph::create_new_edge( platform::Size index1, platform::Size index2 ){
	return new NPDHBSimpleEdge( this, index1, index2 );
}

NPDHBSimpleNode * NPDHBSimpleInteractionGraph::npdhb_cast( utility::graph::Node * node )
{
	return static_cast< NPDHBSimpleNode * > ( node );
}

NPDHBSimpleNode const * NPDHBSimpleInteractionGraph::npdhb_cast( utility::graph::Node const * node ) const
{
	return static_cast< NPDHBSimpleNode const * > ( node );
}



} //namespace interaction_graph
} //namespace pack
} //namespace core
