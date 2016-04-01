// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/TenANeighborGraph.cc
/// @brief  Neighbor graph to represent for each residue the number of other residues within 10 Angstroms
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/scoring/TenANeighborGraph.hh>

// Boost Headers
#include <core/graph/unordered_object_pool.hpp>

#include <utility/vector1.hh>
#include <boost/pool/pool.hpp>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {

///////////////////////////////
//   class TenANeighborNode ///
///////////////////////////////

TenANeighborNode::~TenANeighborNode() {}

TenANeighborNode::TenANeighborNode( graph::Graph * owner, Size node_id )
:
	parent( owner, node_id ),
	neighbor_mass_( 1.0 ), // Default -- 1 neighbor mass unit
	sum_of_neighbors_masses_( 0.0 ),
	since_last_sonm_update_( 0 )
{}

void
TenANeighborNode::copy_from( Node const * source )
{
	debug_assert( dynamic_cast< TenANeighborNode const * > (source) );
	TenANeighborNode const * tAsource( static_cast< TenANeighborNode const * > (source) );
	neighbor_mass_           = tAsource->neighbor_mass_;
	sum_of_neighbors_masses_ = tAsource->sum_of_neighbors_masses_;
	since_last_sonm_update_  = tAsource->since_last_sonm_update_;
}

/// @details triggers update of the neighbor mass sums for this nodes' neighbors
void
TenANeighborNode::neighbor_mass( Real mass ) {
	neighbor_mass_ = mass;
	for ( EdgeListConstIter eiter = const_edge_list_begin(),
			eiter_end = const_edge_list_end(); eiter != eiter_end; ++eiter ) {
		TenANeighborNode const * tenaneighb =
			static_cast< TenANeighborNode const * > ((*eiter)->get_other_node( get_node_index() ));
		tenaneighb->update_neighbor_mass_sum();
	}
}


Size
TenANeighborNode::count_static_memory() const {
	return sizeof( TenANeighborNode );
}


Size
TenANeighborNode::count_dynamic_memory() const {
	return parent::count_dynamic_memory();
}

void
TenANeighborNode::update_neighbor_mass_sum() const
{
	using namespace graph;
	sum_of_neighbors_masses_ = 0.0;
	for ( EdgeListConstIter eiter = const_edge_list_begin(),
			eiter_end = const_edge_list_end(); eiter != eiter_end; ++eiter ) {
		TenANeighborNode const * tenaneighb =
			static_cast< TenANeighborNode const * > ((*eiter)->get_other_node( get_node_index() ));
		sum_of_neighbors_masses_ += tenaneighb->neighbor_mass_;
	}
	since_last_sonm_update_ = 0;
}


///////////////////////////////
//   class TenANeighborEdge ///
///////////////////////////////

TenANeighborEdge::~TenANeighborEdge()
{
	get_TenANode( 0 )->subtract_neighbors_mass( get_TenANode( 1 )->neighbor_mass() );
	get_TenANode( 1 )->subtract_neighbors_mass( get_TenANode( 0 )->neighbor_mass() );
}

TenANeighborEdge::TenANeighborEdge(
	graph::Graph* owner,
	Size first_node_ind,
	Size second_node_ind )
:
	parent( owner, first_node_ind, second_node_ind )
{
	get_TenANode( 0 )->add_neighbors_mass( get_TenANode( 1 )->neighbor_mass() );
	get_TenANode( 1 )->add_neighbors_mass( get_TenANode( 0 )->neighbor_mass() );
}

void TenANeighborEdge::copy_from( Edge const * /*source*/ )
{
	//debug_assert( dynamic_cast< TenNeighborEdge const * > ( source ) );
}

Size TenANeighborEdge::count_static_memory() const
{
	return sizeof( TenANeighborEdge );
}

Size TenANeighborEdge::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}

///////////////////////////////
//  class TenANeighborGraph ///
///////////////////////////////

Distance        const TenANeighborGraph::tenA_( 10.0 );
DistanceSquared const TenANeighborGraph::tenA_squared_( tenA_ * tenA_);

TenANeighborGraph::~TenANeighborGraph() { delete_everything(); delete tenA_edge_pool_; tenA_edge_pool_ = 0; }

TenANeighborGraph::TenANeighborGraph()
:
	parent(),
	tenA_edge_pool_( new boost::unordered_object_pool< TenANeighborEdge > ( 256 ) )
{}

TenANeighborGraph::TenANeighborGraph( Size num_nodes)
:
	parent(),
	tenA_edge_pool_( new boost::unordered_object_pool< TenANeighborEdge > ( 256 ) )
{
	set_num_nodes( num_nodes );
}

TenANeighborGraph::TenANeighborGraph( TenANeighborGraph const & source )
:
	parent(),
	tenA_edge_pool_( new boost::unordered_object_pool< TenANeighborEdge > ( 256 ) )
{
	parent::operator = ( source );
}


TenANeighborGraph &
TenANeighborGraph::operator = ( TenANeighborGraph const & source )
{
	return static_cast< TenANeighborGraph & >  (parent::operator = ( source ));
}

Distance
TenANeighborGraph::neighbor_cutoff() const
{
	return tenA_;
}

void
TenANeighborGraph::conditionally_add_edge(
	Size lower_node_id,
	Size upper_node_id,
	DistanceSquared dsq
) {
	if ( dsq < tenA_squared_ ) add_edge( lower_node_id, upper_node_id );
}

ContextGraphOP
TenANeighborGraph::clone() const
{
	return ContextGraphOP( new TenANeighborGraph( *this ) );
}

void
TenANeighborGraph::update_from_pose(
	pose::Pose const & /*pose*/
)
{}

void TenANeighborGraph::delete_edge( graph::Edge * edge )
{
	//delete edge;
	debug_assert( dynamic_cast< TenANeighborEdge* > (edge) );
	tenA_edge_pool_->destroy( static_cast< TenANeighborEdge* > (edge) );

}

Size
TenANeighborGraph::count_static_memory() const {
	return sizeof ( TenANeighborGraph );
}

Size
TenANeighborGraph::count_dynamic_memory() const {
	return parent::count_dynamic_memory();
}

graph::Node*
TenANeighborGraph::create_new_node( Size node_index ) {
	return new TenANeighborNode( this, node_index );
}

graph::Edge*
TenANeighborGraph::create_new_edge( Size index1, Size index2)
{
	return tenA_edge_pool_->construct( this, index1, index2 );
}

graph::Edge*
TenANeighborGraph::create_new_edge( graph::Edge const * example_edge )
{
	return tenA_edge_pool_->construct(
		this,
		example_edge->get_first_node_ind(),
		example_edge->get_second_node_ind()
	);
}


} // scoring
} // core


#ifdef    SERIALIZATION

template < class Archive >
void
core::scoring::TenANeighborNode::save_to_archive( Archive & arc ) const {
	arc( neighbor_mass_, sum_of_neighbors_masses_, since_last_sonm_update_ );
}

template < class Archive >
void
core::scoring::TenANeighborNode::load_from_archive( Archive & arc ) {
	arc( neighbor_mass_, sum_of_neighbors_masses_, since_last_sonm_update_ );
}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::TenANeighborGraph::save( Archive & arc ) const
{
  arc( num_nodes() );
  for ( Size ii = 1; ii <= num_nodes(); ++ii ) {
    (static_cast< TenANeighborNode const * > ( get_node( ii ) ))->save_to_archive( arc );
  }
  arc( num_edges() );
  for ( EdgeListConstIter iter = const_edge_list_begin(), iter_end = const_edge_list_end(); iter != iter_end; ++iter ) {
    arc( (*iter)->get_first_node_ind(), (*iter)->get_second_node_ind() );
  }
	// Don't serialize the pool
	// EXEMPT tenA_edge_pool_

}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::TenANeighborGraph::load( Archive & arc )
{
  Size num_nodes(0); arc( num_nodes );
  set_num_nodes( num_nodes );

  for ( Size ii = 1; ii <= num_nodes; ++ii ) {
    (static_cast< TenANeighborNode * > ( get_node( ii ) ))->load_from_archive( arc );
  }

  Size num_edges(0); arc( num_edges );
  for ( Size ii = 1; ii <= num_edges; ++ii ) {
    Size node1(0), node2(0); arc( node1, node2 );
    /*Edge * new_edge =*/ add_edge( node1, node2 );
  }
	// Don't deserialize the pool, either
	// EXEMPT tenA_edge_pool_
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::TenANeighborGraph );
CEREAL_REGISTER_TYPE( core::scoring::TenANeighborGraph )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_TenANeighborGraph )
#endif // SERIALIZATION
