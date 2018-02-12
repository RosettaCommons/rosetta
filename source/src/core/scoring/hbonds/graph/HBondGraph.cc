// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/graph/HBondGraph.cc
/// @brief HBondGraph, HBondNode, and HBondEdge classes
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <core/scoring/hbonds/graph/HBondGraph.hh>

#include <basic/Tracer.hh>
//#include <core/scoring/hbonds/HBondSet.hh>
#include <core/types.hh>

static basic::Tracer TR( "core.scoring.hbonds.graph.HBondGraph" );

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {

//dummy! please do not call these
HBondNode::HBondNode() :
	utility::graph::Node( nullptr, 0 ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_()
{ runtime_assert( false ); }

HBondNode::HBondNode( const HBondNode&  ) :
	utility::graph::Node( nullptr, 0 ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_()
{ runtime_assert( false ); }
///////////

//Constructor
HBondNode::HBondNode( utility::graph::Graph* owner, core::Size node_id ) :
	utility::graph::Node( owner, node_id ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_()
{}

HBondNode::HBondNode( utility::graph::Graph* owner, core::Size node_id, core::Size mres_id, core::Size rotamer_id ) :
	utility::graph::Node( owner, node_id ),
	mres_id_( mres_id ),
	rotamer_id_( rotamer_id ),
	ids_of_clashing_nodes_()
{}

//Destructor
HBondNode::~HBondNode() = default;

void HBondNode::copy_from( utility::graph::Node const * source ){
	utility::graph::Node::copy_from( source );

	auto const * src = dynamic_cast< HBondNode const * >( source );
	debug_assert( src );

	mres_id_ = src->mres_id_;
	rotamer_id_ = src->rotamer_id_;
	ids_of_clashing_nodes_ = src->ids_of_clashing_nodes_;
}

void HBondNode::print() const {
	TR << "HBondNode: rotamer_id:" << rotamer_id_ << " mres_id:" << mres_id_ << " num_edges: " << num_edges() << std::endl;
}

core::Size HBondNode::count_static_memory() const
{
	return sizeof( HBondNode );
}

core::Size HBondNode::count_dynamic_memory() const
{
	return utility::graph::Node::count_dynamic_memory() + ids_of_clashing_nodes_.size() * sizeof( unsigned int );
}


//Constructor
HBondEdge::HBondEdge( utility::graph::Graph* owner, core::Size first_node_ind, core::Size second_node_ind ):
	utility::graph::Edge( owner, first_node_ind, second_node_ind ),
	energy_( 0 )
{}


HBondEdge::HBondEdge( utility::graph::Graph* owner, core::Size first_node_ind, core::Size second_node_ind, core::Real energy ):
	utility::graph::Edge( owner, first_node_ind, second_node_ind ),
	energy_( energy )
{}

//Destructor
HBondEdge::~HBondEdge() = default;

void HBondEdge::copy_from( utility::graph::Edge const * source ){
	utility::graph::Edge::copy_from( source );

	auto const * src = dynamic_cast< HBondEdge const * >( source );
	debug_assert( src );
	energy_ = src->energy_;
}

core::Size HBondEdge::count_static_memory() const
{
	return sizeof( HBondEdge );
}

core::Size HBondEdge::count_dynamic_memory() const
{
	//Not sure how to count the HBond dynamic memory
	return utility::graph::Edge::count_dynamic_memory();

}


//Constructor
HBondGraph::HBondGraph() :
	AbstractHBondGraph(),
	hbond_edge_pool_( new boost::unordered_object_pool< HBondEdge > ( 1024 ) )
{}


HBondGraph::HBondGraph( core::Size num_nodes ) :
	AbstractHBondGraph(),
	hbond_edge_pool_( new boost::unordered_object_pool< HBondEdge > ( 1024 ) )
{
	set_num_nodes( num_nodes );
}

//Destructor

HBondGraph::~HBondGraph()
{
	delete_everything();
	delete hbond_edge_pool_;
	hbond_edge_pool_ = nullptr;
}

void
HBondGraph::set_num_nodes( platform::Size num_nodes ){
	all_nodes_.clear();
	all_nodes_.reserve( num_nodes );
	for ( core::Size ii = 1; ii <= num_nodes; ++ii ) {
		all_nodes_.emplace_back( this, ii );
	}
	utility::graph::Graph::set_num_nodes( num_nodes );
}

utility::graph::Node *
HBondGraph::create_new_node( platform::Size node_index ){
	return get_hbondnode( node_index );
	//return new HBondNode( this, node_index );
}


utility::graph::Edge *
HBondGraph::create_new_edge( core::Size index1, core::Size index2 ){
	return hbond_edge_pool_->construct( this, index1, index2 );
}


utility::graph::Edge *
HBondGraph::create_new_edge( utility::graph::Edge const * example_edge ){
	return hbond_edge_pool_->construct(
		this,
		example_edge->get_first_node_ind(),
		example_edge->get_second_node_ind()
	);
}


void
HBondGraph::delete_edge( utility::graph::Edge * edge )
{
	hbond_edge_pool_->destroy( static_cast< HBondEdge * >( edge ) );
}


core::Size
HBondGraph::count_static_memory() const
{
	return sizeof( HBondGraph );
}


core::Size
HBondGraph::count_dynamic_memory() const
{
	//so basically we are not really overriding this at the moment
	return utility::graph::Graph::count_dynamic_memory();
}

} //graph
} //hbonds
} //scoring
} //core
