// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/TwelveANeighborGraph.cc
/// @brief  Neighbor graph to represent for each residue the number of other residues within 12 Angstroms
/// @author Mike Tyka mtyka@u.washington.edu

// Unit Headers
#include <core/scoring/TwelveANeighborGraph.hh>

#include <utility/vector1.hh>


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
//   class TwelveANeighborNode ///
///////////////////////////////

TwelveANeighborNode::~TwelveANeighborNode() {}

TwelveANeighborNode::TwelveANeighborNode( graph::Graph * owner, Size node_id )
:
	parent( owner, node_id )
{}

void
TwelveANeighborNode::copy_from( Node const * /*source*/ )
{
	//debug_assert( dynamic_cast< TwelveANeighborNode const * > (source) );
}

Size
TwelveANeighborNode::count_static_memory() const {
	return sizeof( TwelveANeighborNode );
}


Size
TwelveANeighborNode::count_dynamic_memory() const {
	return parent::count_dynamic_memory();
}

///////////////////////////////
//   class TwelveANeighborEdge ///
///////////////////////////////

TwelveANeighborEdge::~TwelveANeighborEdge() {}

TwelveANeighborEdge::TwelveANeighborEdge(
	graph::Graph* owner,
	Size first_node_ind,
	Size second_node_ind )
:
	parent( owner, first_node_ind, second_node_ind )
{}

void TwelveANeighborEdge::copy_from( Edge const * /*source*/ )
{
	//debug_assert( dynamic_cast< TwelveNeighborEdge const * > ( source ) );
}

Size TwelveANeighborEdge::count_static_memory() const
{
	return sizeof( TwelveANeighborEdge );
}

Size TwelveANeighborEdge::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}

///////////////////////////////
//  class TwelveANeighborGraph ///
///////////////////////////////

/// @details -- 12 A between nbr_atoms, + 6.12 A to the tip of arginine
Distance        const TwelveANeighborGraph::twelveA_( 18.21 );

DistanceSquared const TwelveANeighborGraph::twelveA_squared_( twelveA_ * twelveA_);

TwelveANeighborGraph::~TwelveANeighborGraph() { delete_everything(); }

TwelveANeighborGraph::TwelveANeighborGraph()
:
	parent()
{}

TwelveANeighborGraph::TwelveANeighborGraph( Size num_nodes )
:
	parent()
{
	set_num_nodes( num_nodes );
}

TwelveANeighborGraph::TwelveANeighborGraph( TwelveANeighborGraph const & source )
:
	parent( source )
{
	operator = ( source );
}


TwelveANeighborGraph &
TwelveANeighborGraph::operator = ( TwelveANeighborGraph const & source )
{
	return static_cast< TwelveANeighborGraph & >  (parent::operator = ( source ));
}

Distance
TwelveANeighborGraph::neighbor_cutoff() const
{
	return twelveA_;
}

void
TwelveANeighborGraph::conditionally_add_edge(
	Size lower_node_id,
	Size upper_node_id,
	DistanceSquared dsq
)
{
	if ( dsq < twelveA_squared_ ) add_edge( lower_node_id, upper_node_id );
}

ContextGraphOP
TwelveANeighborGraph::clone() const
{
	return ContextGraphOP( new TwelveANeighborGraph( *this ) );
}

void
TwelveANeighborGraph::update_from_pose(
	pose::Pose const & /*pose*/
)
{}

Size
TwelveANeighborGraph::count_static_memory() const {
	return sizeof ( TwelveANeighborGraph );
}

Size
TwelveANeighborGraph::count_dynamic_memory() const {
	return parent::count_dynamic_memory();
}

void TwelveANeighborGraph::delete_edge( graph::Edge * edge )
{
	delete edge;
}

graph::Node*
TwelveANeighborGraph::create_new_node( Size node_index ) {
	return new TwelveANeighborNode( this, node_index );
}

graph::Edge*
TwelveANeighborGraph::create_new_edge( Size index1, Size index2)
{
	return new TwelveANeighborEdge( this, index1, index2 );
}

graph::Edge*
TwelveANeighborGraph::create_new_edge( graph::Edge const * example_edge )
{
	return new TwelveANeighborEdge(
		this,
		example_edge->get_first_node_ind(),
		example_edge->get_second_node_ind()
	);
}


} // scoring
} // core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::TwelveANeighborGraph::save( Archive & arc ) const {
	arc( cereal::base_class< core::graph::Graph >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::TwelveANeighborGraph::load( Archive & arc ) {
	arc( cereal::base_class< core::graph::Graph >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::TwelveANeighborGraph );
CEREAL_REGISTER_TYPE( core::scoring::TwelveANeighborGraph )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_TwelveANeighborGraph )
#endif // SERIALIZATION
