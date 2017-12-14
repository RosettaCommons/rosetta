// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/SimpleInteractionGraph.cc
/// @brief
/// @author Liz Kellogg (ekellogg@u.washington.edu)

// Unit headers
#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>

// Package headers
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/packer_neighbors.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/util.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/Energies.hh>


// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/graph/Graph.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace interaction_graph {

using namespace utility::graph;
using namespace scoring;

static basic::Tracer TR( "core.pack.interaction_graph.SimpleInteractionGraph" );

SimpleNode::SimpleNode( Graph* owner, Size node_id ):
	Node( owner, node_id  ),
	moved_(false),
	// resnum_( node_id ),
	current_one_body_energy_(0),
	alternate_one_body_energy_(0)
{
	initialize();
}

SimpleNode::~SimpleNode() = default;

void
SimpleNode::initialize()
{
	/**
	SimpleInteractionGraphCOP ig( static_cast< SimpleInteractionGraph const * >(Node::get_owner()));
	runtime_assert( ig );
	set_current( ig->pose().residue( resnum_ ) );
	**/ //better to set up in graph initialization??
}

Real
SimpleNode::one_body_energy() const
{
	if ( moved_ ) {
		return alternate_one_body_energy_;
	} else {
		return current_one_body_energy_;
	}
}

Real
SimpleNode::proposed_one_body_energy() const
{
	return alternate_one_body_energy_;
}

Real
SimpleNode::current_one_body_energy() const
{
	return current_one_body_energy_;
}

void
SimpleNode::commit_change()
{
	current_residue_ = alternate_residue_;
	commit_change_no_res_pointer_update();
}

void
SimpleNode::commit_change_no_res_pointer_update()
{
	moved_ = false;
	current_one_body_energy_ = alternate_one_body_energy_;
	for ( EdgeListIterator edge_itr = edge_list_begin(); edge_itr != edge_list_end(); ++edge_itr ) {
		//update energies
		auto* current_edge(static_cast< SimpleEdge *>(*edge_itr));
		current_edge->commit_change();
	}
}

void
SimpleNode::reject_change( conformation::ResidueCOP res, basic::datacache::BasicDataCache & residue_data_cache )
{
	moved_ = false;
	setup_for_scoring_for_residue( res, residue_data_cache );
}

bool
SimpleNode::moved() const {
	return moved_;
}

void
SimpleNode::set_current( conformation::ResidueCOP res, basic::datacache::BasicDataCache & residue_data_cache )
{
	setup_for_scoring_for_residue( res, residue_data_cache );

	//TR.Debug << "setting res " << res->seqpos() << " to current: " << std::endl;
	if ( ! current_residue_ ) {
		/// This is the first time current_ is being set; calculate backbone centroid and radius
		/*bb_centroid_.zero();
		Size count_n_bb( 0 );
		for ( Size ii = 1; ii <= res->type().first_sidechain_atom() - 1; ++ii ) {
		++count_n_bb;
		bb_centroid_ += res->xyz( ii );
		}
		if ( count_n_bb != 0 ) { bb_centroid_ /= count_n_bb; }
		bb_radius_ = 0.0;
		for ( Size ii = 1; ii <= res->type().first_sidechain_atom() - 1; ++ii ) {
		Real d2 = res->xyz( ii ).distance_squared( bb_centroid_ );
		if ( bb_radius_ < d2 ) bb_radius_ = d2;
		}*/

		bb_centroid_ = scoring::compute_bb_centroid( *res );
		bb_radius_ = scoring::compute_bb_radius( *res, bb_centroid_ );
	}

	current_residue_ = res;
	curr_sc_centroid_ = scoring::compute_sc_centroid( *current_residue_ );
	curr_sc_radius_ =   scoring::compute_sc_radius(   *current_residue_, curr_sc_centroid_ );

	moved_ = false;
	//current_residue_->update_actcoord();
	update_current_one_body_energy();
	//also update edge energies
	for ( EdgeListIterator edge_itr = edge_list_begin(); edge_itr != edge_list_end(); ++edge_itr ) {
		//update energies
		auto* current_edge(static_cast< SimpleEdge *>(*edge_itr));
		current_edge->update_current_energy();
	}
}


void
SimpleNode::set_alternate( conformation::ResidueCOP res, basic::datacache::BasicDataCache & residue_data_cache )
{
	setup_for_scoring_for_residue( res, residue_data_cache );

	//TR.Debug << "setting res " << res->seqpos() << " to alternate: " << std::endl;
	if ( alternate_residue_ != nullptr && TR.Debug.visible() ) {
		// output statements like this are very expensive, even when muted!
		TR.Debug << "setting res " << alternate_residue_->type().name() << " to new-res: " << res->type().name() << std::endl;
	}
	alternate_residue_ = res;

	alt_sc_centroid_ = scoring::compute_sc_centroid( *alternate_residue_ );
	alt_sc_radius_ =   scoring::compute_sc_radius(   *alternate_residue_, alt_sc_centroid_ );

	moved_ = true;
	//alternate_residue_->update_actcoord();
	update_alternate_one_body_energy();
	//also update edge energies
	for ( EdgeListIterator edge_itr = this->edge_list_begin();
			edge_itr != this->edge_list_end();
			++edge_itr ) {
		auto* current_edge(static_cast< SimpleEdge *>(*edge_itr));
		current_edge->update_proposed_energy();
	}
}

void SimpleNode::set_current_no_E_update( conformation::ResidueCOP res, basic::datacache::BasicDataCache & residue_data_cache )
{
	setup_for_scoring_for_residue( res, residue_data_cache );
	current_residue_ = res;
}

void SimpleNode::set_alternate_no_E_update( conformation::ResidueCOP res, basic::datacache::BasicDataCache & residue_data_cache )
{
	setup_for_scoring_for_residue( res, residue_data_cache );
	alternate_residue_ = res;
	moved_ = true;
}

void SimpleNode::update_energies_after_passive_change()
{
	update_alternate_one_body_energy();
	//also update edge energies
	for ( EdgeListIterator edge_itr = this->edge_list_begin();
			edge_itr != this->edge_list_end(); ++edge_itr ) {
		auto* current_edge(static_cast< SimpleEdge *>(*edge_itr));
		current_edge->update_proposed_energy();
	}
}


conformation::ResidueCOP
SimpleNode::get_current() const
{
	return current_residue_;
}

conformation::ResidueCOP
SimpleNode::get_alternate() const
{
	return alternate_residue_;
}

/// @brief return the current state -- must not be 0
conformation::Residue const &
SimpleNode::get_current_ref() const
{
	debug_assert( current_residue_ );
	return *current_residue_;
}

/// @brief return the alternate state -- must not be 0
conformation::Residue const &
SimpleNode::get_alternate_ref() const
{
	debug_assert( alternate_residue_ );
	return *alternate_residue_;
}

Vector const &
SimpleNode::bb_centroid() const
{
	return bb_centroid_;
}

Real
SimpleNode::bb_radius() const
{
	return bb_radius_;
}

Vector const &
SimpleNode::curr_sc_centroid() const
{
	return curr_sc_centroid_;
}

Real
SimpleNode::curr_sc_radius() const
{
	return curr_sc_radius_;
}

Vector const &
SimpleNode::alt_sc_centroid() const
{
	return alt_sc_centroid_;
}

Real
SimpleNode::alt_sc_radius() const
{
	return alt_sc_radius_;
}

void
SimpleNode::setup_for_scoring_for_residue( conformation::ResidueCOP res, basic::datacache::BasicDataCache & residue_data_cache )
{
	//std::cout << "setup for scoring for residue " << res->seqpos() << " " << res->name() << std::endl;
	// call setup for scoring for residue on the energy methods that require it:
	for ( auto const & method : get_simple_ig_owner()->setup_for_scoring_for_residue_energy_methods() ) {
		method->setup_for_scoring_for_residue( *res, get_simple_ig_owner()->pose(),
			get_simple_ig_owner()->scorefunction(), residue_data_cache );
	}
}

void
SimpleNode::update_current_one_body_energy()
{
	auto const * ig( static_cast< SimpleInteractionGraph const * >(Node::get_owner()));
	runtime_assert( ig );
	EnergyMap emap;
	ig->scorefunction().eval_ci_intrares_energy( *current_residue_, ig->pose(), emap );
	ig->scorefunction().eval_cd_intrares_energy( *current_residue_, ig->pose(), emap );
	ig->scorefunction().eval_ci_1b( *current_residue_, ig->pose(), emap );
	ig->scorefunction().eval_cd_1b( *current_residue_, ig->pose(), emap );
	current_one_body_energy_ = ig->scorefunction().weights().dot(emap);
}

void
SimpleNode::update_alternate_one_body_energy()
{
	debug_assert( dynamic_cast< SimpleInteractionGraph const * >(Node::get_owner()) );

	auto const * ig( static_cast< SimpleInteractionGraph const * >(Node::get_owner()));
	EnergyMap emap;
	ig->scorefunction().eval_ci_intrares_energy( *alternate_residue_, ig->pose(), emap );
	ig->scorefunction().eval_cd_intrares_energy( *alternate_residue_, ig->pose(), emap );
	ig->scorefunction().eval_ci_1b( *alternate_residue_, ig->pose(), emap );
	ig->scorefunction().eval_cd_1b( *alternate_residue_, ig->pose(), emap );
	alternate_one_body_energy_ = ig->scorefunction().weights().dot( emap );
}

/*Vector
SimpleNode::calc_sc_centroid( conformation::Residue const & res ) const
{
Vector centroid( 0.0 );
Size count( 0 );
for ( Size ii = res.type().first_sidechain_atom(); ii <= res.type().nheavyatoms(); ++ii ) {
count += 1;
centroid += res.xyz( ii );
}
if ( count == 0 ) {
return bb_centroid_;
} else {
centroid /= count;
return centroid;
}
}

Real
SimpleNode::calc_sc_radius( conformation::Residue const & res, Vector const & centroid )
{
Real max_d2 = 0;
for ( Size ii = res.type().first_sidechain_atom(); ii <= res.type().nheavyatoms(); ++ii ) {
Real d2 = res.xyz(ii).distance_squared( centroid );
if ( d2 > max_d2 ) max_d2 = d2;
}
return std::sqrt( max_d2 );
}*/


/*SimpleEdge::SimpleEdge( Graph* owner ):
Edge( owner, 0, 0 ),
current_energy_(0),
proposed_energy_(0)
{
for ( Size ii = 0; ii < 3; ++ii ) {
for ( Size jj = 0; jj < 3; ++jj ) {
bb_bbE_calced_[ ii ][ jj ] = false;
bb_bbE_[ ii ][ jj ] = -1234; // sentinel; this should never be returned.
}
}
}*/


SimpleEdge::SimpleEdge( Graph* owner, Size res1, Size res2 ):
	Edge( owner, res1, res2 ),
	// long_range_energies_exist_( false ),
	current_energy_(0),
	proposed_energy_(0),
	calc_short_range_(true),
	calc_long_range_(true)
{
	// compute_energy();
	for ( Size ii = 0; ii < 3; ++ii ) {
		for ( Size jj = 0; jj < 3; ++jj ) {
			bb_bbE_calced_[ ii ][ jj ] = false;
			bb_bbE_[ ii ][ jj ] = -1234; // sentinel; this should never be returned.
		}
	}
}

SimpleEdge::~SimpleEdge() = default;

Real
SimpleEdge::get_current_energy() const
{
	return current_energy_;
}

Real
SimpleEdge::get_proposed_energy() const
{
	return proposed_energy_;
}

void
SimpleEdge::update_current_energy()
{
	compute_energy( true, true );
}

void
SimpleEdge::update_proposed_energy()
{
	auto * node1( static_cast< SimpleNode * > (this->get_node( (platform::Size) 0 )));
	auto * node2( static_cast< SimpleNode * > (this->get_node( (platform::Size) 1 )));
	bool use_current_node1 = true;
	bool use_current_node2 = true;
	if ( node1->moved() ) {
		use_current_node1 = false;
	}
	if ( node2->moved() ) {
		use_current_node2 = false;
	}
	compute_energy( use_current_node1, use_current_node2 );
}

void
SimpleEdge::commit_change()
{
	current_energy_ = proposed_energy_;
}

void
SimpleEdge::compute_energy( bool use_current_node1, bool use_current_node2 )
{
	//  TR.Debug << "num nodes " << (this->get_owner())->num_nodes() << std::endl;
	auto * node1( static_cast< SimpleNode * > (this->get_node( (platform::Size) 0 )));
	auto * node2( static_cast< SimpleNode * > (this->get_node( (platform::Size) 1 )));
	debug_assert( node1 && node2 );
	conformation::ResidueCOP res1;
	conformation::ResidueCOP res2;

	//Real   r1bb_rad( node1->bb_radius() + 0.5 ), r2bb_rad( node2->bb_radius() + 0.5 );
	Vector r1bb_centroid( node1->bb_centroid() ), r2bb_centroid( node2->bb_centroid() );
	Real   r1bb_rad( node1->bb_radius() ), r2bb_rad( node2->bb_radius() );
	Vector r1sc_centroid, r2sc_centroid;
	Real   r1sc_rad, r2sc_rad;

	//res1 = use_current_node1 ? node1->get_current() : node1->get_alternate();
	//res2 = use_current_node2 ? node2->get_current() : node2->get_alternate();

	if ( !use_current_node1 ) {
		res1 = node1->get_alternate();
		r1sc_centroid = node1->alt_sc_centroid();
		r1sc_rad = node1->alt_sc_radius();
		//node1_used= "alternate";
	} else {
		res1 = node1->get_current();
		r1sc_centroid = node1->curr_sc_centroid();
		r1sc_rad = node1->curr_sc_radius();
		//node1_used = "current";
	}

	if ( !use_current_node2 ) {
		//node2_used = "alternate";
		res2 = node2->get_alternate();
		r2sc_centroid = node2->alt_sc_centroid();
		r2sc_rad = node2->alt_sc_radius();
	} else {
		//node2_used = "current";
		res2 = node2->get_current();
		r2sc_centroid = node2->curr_sc_centroid();
		r2sc_rad = node2->curr_sc_radius();
	}

	// Add 1.0A fudge to the score-function cutoffs. There's a bug somewhere in either our thinking
	// about the meaning of the atomic-interaction-cutoff, or in the implementation of certain terms
	// that causes deltaE discrepancies in the absence of these
	r1sc_rad += 0.5;
	r2sc_rad += 0.5;
	r1bb_rad += 0.5;
	r2bb_rad += 0.5;

	// TR.Debug << res1->seqpos() << " using " << node1_used << " " << res2->seqpos() << " " << node2_used << std::endl;

	debug_assert( dynamic_cast< SimpleInteractionGraph const * >(this->get_owner()) );
	auto const * ig( static_cast< SimpleInteractionGraph const * >(this->get_owner()));
	pose::Pose const & pose = ig->pose();
	scoring::ScoreFunction const & sfxn = ig->scorefunction();
	EnergyMap emap; // APL Note: class TwoBodyEnergyMap has been deprecated / removed

	if ( calc_short_range_ ) {
		// Just call residue pair energy directly?
		//ig->scorefunction().eval_ci_2b( *res1, *res2, pose, emap );
		//ig->scorefunction().eval_cd_2b( *res1, *res2, pose, emap );

		// Or break up residue pair energies into sc/sc, sc/bb, and bb/bb interactions.
		scoring::eval_scsc_sr2b_energies( *res1, *res2, r1sc_centroid, r2sc_centroid, r1sc_rad, r2sc_rad, pose, sfxn, emap );
		scoring::eval_bbsc_sr2b_energies( *res1, *res2, r1bb_centroid, r2sc_centroid, r1bb_rad, r2sc_rad, pose, sfxn, emap );
		scoring::eval_bbsc_sr2b_energies( *res2, *res1, r2bb_centroid, r1sc_centroid, r2bb_rad, r1sc_rad, pose, sfxn, emap );

		// Instead of calculating bb/bb energies repeatedly, only evaluate bb/bb energies once,
		// allowing for three kinds of backbones: proline, glycine, and everything else. This
		// is handled by the get_bbE function below.
		//scoring::eval_bbbb_sr2b_energies( *res1, *res2, r1bb_centroid, r2bb_centroid, r1bb_rad, r2bb_rad, pose, sfxn, emap );

	}
	if ( calc_long_range_ ) {
		for ( auto it = sfxn.long_range_energies_begin(); it != sfxn.long_range_energies_end(); ++it ) {
			if ( (*it)->defines_residue_pair_energy( pose, get_first_node_ind(), get_second_node_ind()) ) {
				(*it)->residue_pair_energy( *res1, *res2, pose, sfxn, emap);
			}
		}
	}

	//if ( res1->seqpos() == 49 && res2->seqpos() == 72 ) {
	// std::cout << " SimpleIG::compute_energy " << res1->seqpos() << " " << res2->seqpos() << " ";
	// emap.show_weighted( std::cout, ig->scorefunction().weights() );
	// std::cout << std::endl;
	//}

	Real energy =
		ig->scorefunction().weights().dot( emap, ig->scorefunction().ci_2b_types() ) +
		ig->scorefunction().weights().dot( emap, ig->scorefunction().cd_2b_types() ) +
		ig->scorefunction().weights().dot( emap, ig->scorefunction().ci_lr_2b_types() ) +
		ig->scorefunction().weights().dot( emap, ig->scorefunction().cd_lr_2b_types() ) +
		get_bb_E( *res1, *res2 );

	if ( use_current_node1 && use_current_node2 ) {
		current_energy_ = energy;
	} else {
		proposed_energy_ = energy;
	}

}

Real
SimpleEdge::get_bb_E( conformation::Residue const & r1, conformation::Residue const & r2 )
{
	//return 0.0;

	Size const i1 = get_bb_index( r1 );
	Size const i2 = get_bb_index( r2 );

	if ( ! bb_bbE_calced( i1, i2 ) ) {
		auto const * ig( static_cast< SimpleInteractionGraph const * >(this->get_owner()));
		pose::Pose const & pose = ig->pose();
		EnergyMap emap;
		ig->scorefunction().eval_ci_2b_bb_bb( r1, r2, pose, emap );
		ig->scorefunction().eval_cd_2b_bb_bb( r1, r2, pose, emap );
		set_bb_bbE( i1, i2 , ig->scorefunction().weights().dot( emap, ig->scorefunction().ci_2b_types()) +
			ig->scorefunction().weights().dot( emap, ig->scorefunction().cd_2b_types() ) );
		set_bb_bbE_calced( i1, i2 );
		//std::cout << get_first_node_ind() << " " << get_second_node_ind() << " get bb bb E " << bb_bbE( i1, i2 ) << std::endl;
	}
	return bb_bbE( i1, i2 );
}

Size
SimpleEdge::get_bb_index( conformation::Residue const & r ) const
{
	using namespace chemical;
	switch ( r.aa() ) {
	case aa_pro : return 1;
	case aa_gly : return 2;
	default : return 0;
	}
}

bool
SimpleEdge::bb_bb_boundaries( Size ind1, Size ind2 ) const
{
	return ind1 <= 2 && ind2 <= 2; // 0 <= ind1 &&  0 <= ind2 always hold
}

bool
SimpleEdge::bb_bbE_calced( Size ind1, Size ind2 ) const
{
	debug_assert( bb_bb_boundaries( ind1, ind2 ));
	return bb_bbE_calced_[ ind1 ][ ind2 ];
}

void
SimpleEdge::set_bb_bbE_calced( Size ind1, Size ind2 )
{
	debug_assert( bb_bb_boundaries( ind1, ind2 ));
	bb_bbE_calced_[ ind1 ][ ind2 ] = true;
}

void
SimpleEdge::set_bb_bbE( Size ind1, Size ind2, Real val )
{
	debug_assert( bb_bb_boundaries( ind1, ind2 ));
	bb_bbE_[ ind1 ][ ind2 ] = val;
}


Real
SimpleEdge::bb_bbE( Size ind1, Size ind2 ) const
{
	debug_assert( bb_bb_boundaries( ind1, ind2 ));
	return bb_bbE_[ ind1 ][ ind2 ];
}

SimpleInteractionGraph::SimpleInteractionGraph() :
	Graph(),
	sfxn_(),
	pose_()
	// accumulated_ediff_(0)
{}

SimpleInteractionGraph::~SimpleInteractionGraph() = default;

void SimpleInteractionGraph::set_scorefunction( ScoreFunction const & sfxn )
{
	sfxn_ = sfxn.clone();
	if ( pose_ ) {
		for ( auto const & method : sfxn_->all_methods() ) {
			if ( method->requires_a_setup_for_scoring_for_residue_opportunity_during_regular_scoring( *pose_ ) )  {
				sfs_energy_methods_.push_back( method );
			}
		}
	}
}

pose::PoseCOP
SimpleInteractionGraph::pose_cop() const {
	return pose_;
}


void
SimpleInteractionGraph::initialize( pose::Pose const & pose){
	task::PackerTaskOP task( task::TaskFactory::create_packer_task(pose) );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		task->nonconst_residue_task( ii ).restrict_to_repacking();
	}
	initialize(pose, task);
}

void
SimpleInteractionGraph::initialize(
	pose::Pose const & pose,
	task::PackerTaskCOP task
)
{
	TR.Debug << "calling initialize on pose " << std::endl;
	set_pose_no_initialize( pose );

	//Prepare Edges For Short Range Interactions
	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, *sfxn_, task );
	copy_connectivity(*packer_neighbor_graph);
	for ( EdgeListIter iter = edge_list_begin(), iter_end = edge_list_end();
			iter != iter_end; ++iter ) {
		auto * existing_edge = static_cast< SimpleEdge * >( *iter );
		existing_edge->calc_short_range(true);
		existing_edge->calc_long_range(false);
	}


	//Prepare Edges For Long Range Interactions
	for ( auto it = sfxn_->long_range_energies_begin(); it != sfxn_->long_range_energies_end(); ++it ) {
		scoring::LREnergyContainerCOP lrec = pose.energies().long_range_container( (*it)->long_range_type() );

		if ( !lrec || lrec->empty() ) continue;
		for ( Size ii = 1; ii <= pose_->size(); ii++ ) {
			for ( scoring::ResidueNeighborConstIteratorOP rni = lrec->const_upper_neighbor_iterator_begin( ii ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii ); (*rni) != (*rniend); ++(*rni) ) {

				Size node1 = rni->lower_neighbor_id();
				Size node2 = rni->upper_neighbor_id();

				if ( Graph::get_edge_exists( node1, node2 ) ) {
					auto * existing_edge = static_cast< SimpleEdge * >( Graph::find_edge( node1, node2 ) );

					if ( existing_edge ) {
						existing_edge->calc_long_range(true);
					}
				} else {
					auto * new_edge = static_cast< SimpleEdge * > ( add_edge( node1, node2 ));
					new_edge->calc_short_range(false);
					new_edge->calc_long_range(true);
				}

			}//rni
		}//ii
	}//it

	for ( EdgeListIter iter = edge_list_begin(), iter_end = edge_list_end(); iter != iter_end; ++iter ) {
		auto * existing_edge = static_cast< SimpleEdge * >( *iter );
		existing_edge->compute_energy( true, true );
	}

	/*
	for ( Size ii = 2; ii <= pose_->size(); ii++ ) {
	for ( Size jj = 1; jj < ii; jj++ ) {

	bool calc_short_range_interactions = false;
	bool calc_long_range_interactions = false;

	//Short Range
	if ( sfxn_->are_they_neighbors( *pose_, jj, ii ) ) {
	calc_short_range_interactions = true;
	}

	//Long Range
	for( auto it = sfxn_->long_range_energies_begin(); it != sfxn_->long_range_energies_end(); ++it){
	if( (*it)->defines_residue_pair_energy(*pose_, ii, jj) ){
	calc_long_range_interactions = true;
	break;
	}
	}

	//Update Edge Settings and Compute Energy
	if ( calc_short_range_interactions || calc_long_range_interactions){
	if ( Graph::get_edge_exists( jj, ii ) ){
	SimpleEdge * existing_edge = static_cast< SimpleEdge * >( Graph::find_edge( jj, ii ) );

	if ( existing_edge ) {
	existing_edge->calc_short_range(calc_short_range_interactions);
	existing_edge->calc_long_range(calc_long_range_interactions);
	existing_edge->compute_energy( true, true );
	}
	}
	else{
	SimpleEdge * new_edge = static_cast< SimpleEdge * > ( add_edge( jj, ii ));
	new_edge->calc_short_range(calc_short_range_interactions);
	new_edge->calc_long_range(calc_long_range_interactions);
	new_edge->compute_energy( true, true );
	}
	}
	else {//Negate (remove?) any edges that do not represent any interactions
	if ( Graph::get_edge_exists( jj, ii ) ){
	SimpleEdge * existing_edge = static_cast< SimpleEdge * >( Graph::find_edge( jj, ii ) );

	if ( existing_edge ) {
	existing_edge->calc_short_range(false);
	existing_edge->calc_long_range(false);
	existing_edge->compute_energy( true, true );
	}
	}
	}

	}//jj
	}//ii*/
}

void
SimpleInteractionGraph::set_pose_no_initialize( pose::Pose const & pose )
{
	pose_ = pose::PoseCOP( pose::PoseOP( new pose::Pose( pose ) ) );
	if ( sfxn_ ) {
		for ( auto const & method : sfxn_->all_methods() ) {
			if ( method->requires_a_setup_for_scoring_for_residue_opportunity_during_regular_scoring( *pose_ ) )  {
				sfs_energy_methods_.push_back( method );
			}
		}
	}

	if ( num_nodes() != pose.size() ) {
		set_num_nodes( pose.size() );
	}
	for ( Size ii = 1; ii <= pose_->size(); ++ii ) {
		auto * newnode = static_cast< SimpleNode * >(get_node( ii ));
		runtime_assert( newnode );
		conformation::ResidueOP ii_res( new conformation::Residue( pose_->residue( ii ) ));
		newnode->set_current( ii_res, *ii_res->nonconst_data_ptr() );
	}
}


void SimpleInteractionGraph::setup_after_edge_addition()
{
	for ( Size ii = 1; ii <= pose_->total_residue(); ++ii ) {
		for ( auto iter = get_node( ii )->upper_edge_list_begin(), iter_end = get_node( ii )->upper_edge_list_end(); iter != iter_end; ++iter ) {
			(static_cast< SimpleEdge * > (*iter))->compute_energy( true, true );
		}
	}
}

//void SimpleInteractionGraph::setup_after_edge_addition() {}


void
SimpleInteractionGraph::commit_change( Size node_id ){
	static_cast< SimpleNode * >(get_node( node_id ))->commit_change();
	//SimpleNode * node = static_cast< SimpleNode * >(get_node( node_id ));
	//node->set_current( node->get_alternate() );
}

void
SimpleInteractionGraph::reject_change( Size node_id , conformation::ResidueCOP res, basic::datacache::BasicDataCache & residue_data_cache ){
	static_cast< SimpleNode * >( get_node( node_id ) )->reject_change( res, residue_data_cache );
}

/// @details Note, this function returns (currE - altE) which represents
/// the negative of the change in energy for the substition
Real
SimpleInteractionGraph::consider_substitution( Size node_id, conformation::ResidueCOP new_state, basic::datacache::BasicDataCache & residue_data_cache ){
	Real current_energy = 0.0;

	auto* thisnode = static_cast< SimpleNode *>(get_node( node_id ));

	for ( EdgeListIterator edge_itr = thisnode->edge_list_begin();
			edge_itr != thisnode->edge_list_end(); ++edge_itr ) {
		//update energies
		auto* current_edge(static_cast< SimpleEdge *>(*edge_itr));
		current_energy += current_edge->get_current_energy();
	}
	current_energy += thisnode->current_one_body_energy();

	thisnode->set_alternate( new_state, residue_data_cache );
	//iterate over all edges
	Real alternate_energy = 0.0;
	for ( EdgeListIterator edge_itr = thisnode->edge_list_begin();
			edge_itr != thisnode->edge_list_end(); ++edge_itr ) {
		//update energies
		auto* current_edge(static_cast< SimpleEdge *> (*edge_itr));
		alternate_energy += current_edge->get_proposed_energy();
	}
	alternate_energy += thisnode->proposed_one_body_energy();
	return (current_energy - alternate_energy); // this is -deltaE for the substitution
}


Real
SimpleInteractionGraph::total_energy() {
	//iterate over node ids
	Real total_energy = 0.0;
	for ( Size node_itr = 1; node_itr <= pose_->size(); node_itr++ ) {
		auto* thisnode = static_cast< SimpleNode *> (get_node( node_itr ));
		Real residue_energy = thisnode->one_body_energy();
		for ( EdgeListIter edge_iter = thisnode->upper_edge_list_begin();
				edge_iter != thisnode->upper_edge_list_end(); ++edge_iter ) {
			residue_energy += (static_cast< SimpleEdge *>( *edge_iter ))->get_current_energy();
		}
		//TR.Debug << node_itr << ": total residue energy is " << residue_energy << std::endl;
		total_energy += residue_energy;
	}
	return total_energy;

}

Node* SimpleInteractionGraph::create_new_node( platform::Size node_index ){
	return new SimpleNode( this, node_index );
}

Edge* SimpleInteractionGraph::create_new_edge( platform::Size index1, platform::Size index2 ){
	return new SimpleEdge( this, index1, index2 );
}

void SimpleInteractionGraph::delete_edge( Edge * edge )
{
	delete edge;
}

/**
Edge* SimpleInteractionGraph::create_new_edge( Edge* example_edge ){
static_cast < SimpleEdge * >
return new SimpleEdge( example_edge );
}
**/

/**
Real
SimpleInteractionGraph::total_alternate_energy() {
if( moved_ ){
//compute energies
moved_ = false;
}
return alternate_one_body_energy_ + alt_scsc_E + alt_bbsc_E + alt_scbb_E;
}
**/

/**
bool
SimpleInteractionGraph::node_altered( Real const node_id ){
return node(node_id).is_altered();
}
**/ //needed?

} //namespace interaction_graph
} //namespace pack
} //namespace core
