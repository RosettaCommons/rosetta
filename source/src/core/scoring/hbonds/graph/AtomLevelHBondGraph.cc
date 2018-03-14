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
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.hh>

#include <basic/Tracer.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/types.hh>

#include <utility/graph/unordered_object_pool.hpp>
#include <boost/pool/pool.hpp>

static basic::Tracer TR( "core.scoring.hbonds.graph.AtomLevelHBondGraph" );

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {

//dummy! please do not call these
AtomLevelHBondNode::AtomLevelHBondNode() :
	utility::graph::Node( nullptr, 0 ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_( 0 ),
	polar_sc_atoms_not_satisfied_by_background_()
{ runtime_assert( false ); }

AtomLevelHBondNode::AtomLevelHBondNode( const AtomLevelHBondNode & ) :
	utility::graph::Node( nullptr, 0 ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_( 0 ),
	polar_sc_atoms_not_satisfied_by_background_()
{ runtime_assert( false ); }
///////////


//Constructor
AtomLevelHBondNode::AtomLevelHBondNode( utility::graph::Graph* owner, Size node_id ) :
	utility::graph::Node( owner, node_id ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_(),
	polar_sc_atoms_not_satisfied_by_background_()
{}

AtomLevelHBondNode::AtomLevelHBondNode( utility::graph::Graph* owner, Size node_id, Size mres_id, Size rotamer_id ) :
	utility::graph::Node( owner, node_id ),
	mres_id_( mres_id ),
	rotamer_id_( rotamer_id ),
	ids_of_clashing_nodes_(),
	polar_sc_atoms_not_satisfied_by_background_()
{}

//Destructor
AtomLevelHBondNode::~AtomLevelHBondNode() = default;

void AtomLevelHBondNode::copy_from( utility::graph::Node const * source ){
	utility::graph::Node::copy_from( source );

	auto const * src = dynamic_cast< AtomLevelHBondNode const * >( source );
	debug_assert( src );

	mres_id_ = src->mres_id_;
	rotamer_id_ = src->rotamer_id_;
	ids_of_clashing_nodes_ = src->ids_of_clashing_nodes_;
	polar_sc_atoms_not_satisfied_by_background_ =
		src->polar_sc_atoms_not_satisfied_by_background_;
}

void AtomLevelHBondNode::print() const {
	TR << "AtomLevelHBondNode: rotamer_id:" << rotamer_id_ << " mres_id:" << mres_id_ << " num_edges: " << num_edges() << std::endl;
}

Size AtomLevelHBondNode::count_static_memory() const
{
	return sizeof( AtomLevelHBondNode );
}

Size AtomLevelHBondNode::count_dynamic_memory() const
{
	return utility::graph::Node::count_dynamic_memory()
		+ polar_sc_atoms_not_satisfied_by_background_.size() * sizeof( AtomInfo )
		+ ids_of_clashing_nodes_.size() * sizeof( unsigned int );
}


//Constructor
AtomLevelHBondEdge::AtomLevelHBondEdge( utility::graph::Graph* owner, Size first_node_ind, Size second_node_ind ):
	utility::graph::Edge( owner, first_node_ind, second_node_ind ),
	energy_( 0 ),
	hbonds_( 0 )
{}

AtomLevelHBondEdge::AtomLevelHBondEdge( utility::graph::Graph* owner, Size first_node_ind, Size second_node_ind, Real energy ):
	utility::graph::Edge( owner, first_node_ind, second_node_ind ),
	energy_( energy ),
	hbonds_( 0 )
{}

//Destructor
AtomLevelHBondEdge::~AtomLevelHBondEdge() = default;

void AtomLevelHBondEdge::copy_from( utility::graph::Edge const * source ){
	utility::graph::Edge::copy_from( source );

	auto const * src = dynamic_cast< AtomLevelHBondEdge const * >( source );
	debug_assert( src );

	energy_ = src->energy_;
	hbonds_ = src->hbonds_;
}

Size AtomLevelHBondEdge::count_static_memory() const
{
	return sizeof( AtomLevelHBondEdge );
}

Size AtomLevelHBondEdge::count_dynamic_memory() const
{
	return utility::graph::Edge::count_dynamic_memory() + hbonds_.size() * sizeof( HBondInfo );
}


//Constructor
AtomLevelHBondGraph::AtomLevelHBondGraph() :
	hbond_edge_pool_( new boost::unordered_object_pool< AtomLevelHBondEdge > ( 1024 ) ),
	hbond_node_pool_( new boost::unordered_object_pool< AtomLevelHBondNode > ( 1024 ) )
{}


AtomLevelHBondGraph::AtomLevelHBondGraph( Size num_nodes ) :
	hbond_edge_pool_( new boost::unordered_object_pool< AtomLevelHBondEdge > ( 1024 ) ),
	hbond_node_pool_( new boost::unordered_object_pool< AtomLevelHBondNode > ( num_nodes ) )
{
	set_num_nodes( num_nodes );
}

//Destructor

AtomLevelHBondGraph::~AtomLevelHBondGraph()
{
	delete_everything();
	delete hbond_edge_pool_;
	delete hbond_node_pool_;
	hbond_edge_pool_ = nullptr;
	hbond_node_pool_ = nullptr;
}

utility::graph::Node *
AtomLevelHBondGraph::create_new_node( platform::Size node_index ){
	return hbond_node_pool_->construct( this, node_index );
	//return new AtomLevelHBondNode( this, node_index );
}

utility::graph::Edge *
AtomLevelHBondGraph::create_new_edge( Size index1, Size index2 ){
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
AtomLevelHBondGraph::delete_node( utility::graph::Node * node )
{
	hbond_node_pool_->destroy( static_cast< AtomLevelHBondNode * >( node ) );
}

void
AtomLevelHBondGraph::delete_edge( utility::graph::Edge * edge )
{
	hbond_edge_pool_->destroy( static_cast< AtomLevelHBondEdge * >( edge ) );
}

Size
AtomLevelHBondGraph::count_static_memory() const
{
	return sizeof( AtomLevelHBondGraph );
}

Size
AtomLevelHBondGraph::count_dynamic_memory() const
{
	//so basically we are not really overriding this at the moment
	return utility::graph::Graph::count_dynamic_memory();
}

AtomLevelHBondEdge *
AtomLevelHBondGraph::register_hbond( Size rotamerA, Size rotamerB, Real score ) {
	AtomLevelHBondEdge * new_edge = static_cast< AtomLevelHBondEdge * >( add_edge( rotamerA, rotamerB ) );
	new_edge->set_energy( score );
	return new_edge;
}


} //graph
} //hbonds
} //scoring
} //core
