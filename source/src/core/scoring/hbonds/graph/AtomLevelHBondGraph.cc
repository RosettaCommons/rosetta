// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/graph/AtomLevelHBondGraph.cc
/// @brief AtomLevelHBondGraph, AtomLevelHBondNode, and AtomLevelHBondEdge classes
/// @author Jack Maguire, jack@med.unc.edu

#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.hh>

#include <basic/Tracer.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/types.hh>

#include <utility/graph/unordered_object_pool.hpp>
#include <boost/pool/pool.hpp>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.hbonds.graph.AtomLevelHBondGraph" );

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {

//dummy! please do not call these
AtomLevelHBondNode::AtomLevelHBondNode() :
	HBondNode( ),
	polar_sc_atoms_not_satisfied_by_background_()
{ runtime_assert( false ); }

AtomLevelHBondNode::AtomLevelHBondNode( const AtomLevelHBondNode & ) :
	HBondNode( ),
	polar_sc_atoms_not_satisfied_by_background_()
{ runtime_assert( false ); }
///////////


//Constructor
AtomLevelHBondNode::AtomLevelHBondNode( utility::graph::Graph* owner, core::Size node_id ) :
	HBondNode( owner, node_id )
{}

AtomLevelHBondNode::AtomLevelHBondNode( utility::graph::Graph* owner, core::Size node_id, core::Size mres_id, core::Size rotamer_id ) :
	HBondNode( owner, node_id, mres_id, rotamer_id )
{}

//Destructor
AtomLevelHBondNode::~AtomLevelHBondNode()
{}

void AtomLevelHBondNode::copy_from( utility::graph::Node const * source ){
	HBondNode::copy_from( source );

	AtomLevelHBondNode const * src = dynamic_cast< AtomLevelHBondNode const * >( source );
	debug_assert( src );
	polar_sc_atoms_not_satisfied_by_background_ = src->polar_sc_atoms_not_satisfied_by_background_;
}

/*void AtomLevelHBondNode::print() const {
TR << "AtomLevelHBondNode: rotamer_id:" << rotamer_id_ << " mres_id:" << mres_id_ << " num_edges: " << num_edges() << std::endl;
}*/

core::Size AtomLevelHBondNode::count_static_memory() const
{
	return sizeof( AtomLevelHBondNode );
}

core::Size AtomLevelHBondNode::count_dynamic_memory() const
{
	return HBondNode::count_dynamic_memory() + polar_sc_atoms_not_satisfied_by_background_.size() * sizeof( AtomInfo );
}


//Constructor
AtomLevelHBondEdge::AtomLevelHBondEdge( utility::graph::Graph* owner, core::Size first_node_ind, core::Size second_node_ind ):
	HBondEdge( owner, first_node_ind, second_node_ind )
{}


AtomLevelHBondEdge::AtomLevelHBondEdge( utility::graph::Graph* owner, core::Size first_node_ind, core::Size second_node_ind, core::Real energy ):
	HBondEdge( owner, first_node_ind, second_node_ind, energy )
{}

//Destructor
AtomLevelHBondEdge::~AtomLevelHBondEdge()
{}

void AtomLevelHBondEdge::copy_from( utility::graph::Edge const * source ){
	HBondEdge::copy_from( source );

	AtomLevelHBondEdge const * src = dynamic_cast< AtomLevelHBondEdge const * >( source );
	debug_assert( src );
	hbonds_ = src->hbonds_;
}

core::Size AtomLevelHBondEdge::count_static_memory() const
{
	return sizeof( AtomLevelHBondEdge );
}

core::Size AtomLevelHBondEdge::count_dynamic_memory() const
{
	return HBondEdge::count_dynamic_memory() + hbonds_.size() * sizeof( HBondInfo );
}



//Constructor
AtomLevelHBondGraph::AtomLevelHBondGraph() :
	AbstractHBondGraph(),
	hbond_edge_pool_( new boost::unordered_object_pool< AtomLevelHBondEdge > ( 1024 ) )
{}


AtomLevelHBondGraph::AtomLevelHBondGraph( core::Size num_nodes ) :
	AbstractHBondGraph(),
	hbond_edge_pool_( new boost::unordered_object_pool< AtomLevelHBondEdge > ( 1024 ) )
{
	set_num_nodes( num_nodes );
}

//Destructor

AtomLevelHBondGraph::~AtomLevelHBondGraph()
{
	delete_everything();
	delete hbond_edge_pool_;
	hbond_edge_pool_ = nullptr;
}

void
AtomLevelHBondGraph::set_num_nodes( platform::Size num_nodes ){
	all_nodes_.clear();
	all_nodes_.reserve( num_nodes );
	for ( core::Size ii = 1; ii <= num_nodes; ++ii ) {
		all_nodes_.emplace_back( this, ii );
	}
	utility::graph::Graph::set_num_nodes( num_nodes );
}

utility::graph::Node *
AtomLevelHBondGraph::create_new_node( platform::Size node_index ){
	return get_hbondnode( node_index );
	//return new AtomLevelHBondNode( this, node_index );
}

utility::graph::Edge *
AtomLevelHBondGraph::create_new_edge( core::Size index1, core::Size index2 ){
	return hbond_edge_pool_->construct( this, index1, index2 );
}

utility::graph::Edge *
AtomLevelHBondGraph::create_new_edge( utility::graph::Edge const * example_edge ){
	return hbond_edge_pool_->construct(
		this,
		example_edge->get_first_node_ind(),
		example_edge->get_second_node_ind()
	);
}

void
AtomLevelHBondGraph::delete_edge( utility::graph::Edge * edge )
{
	hbond_edge_pool_->destroy( static_cast< AtomLevelHBondEdge * >( edge ) );
}

core::Size
AtomLevelHBondGraph::count_static_memory() const
{
	return sizeof( AtomLevelHBondGraph );
}

core::Size
AtomLevelHBondGraph::count_dynamic_memory() const
{
	//so basically we are not really overriding this at the moment
	return utility::graph::Graph::count_dynamic_memory();
}


} //graph
} //hbonds
} //scoring
} //core
