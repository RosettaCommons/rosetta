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
	utility::graph::LowMemNode( 0 ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_( 0 ),
	polar_sc_atoms_not_satisfied_by_background_()
{ runtime_assert( false ); }

AtomLevelHBondNode::AtomLevelHBondNode( const AtomLevelHBondNode & ) :
	utility::graph::LowMemNode( 0 ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_( 0 ),
	polar_sc_atoms_not_satisfied_by_background_()
{ runtime_assert( false ); }
///////////


//Constructor
AtomLevelHBondNode::AtomLevelHBondNode( Size node_id ) :
	utility::graph::LowMemNode( node_id ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_(),
	polar_sc_atoms_not_satisfied_by_background_()
{}

AtomLevelHBondNode::AtomLevelHBondNode( Size node_id, Size mres_id, Size rotamer_id ) :
	utility::graph::LowMemNode( node_id ),
	mres_id_( mres_id ),
	rotamer_id_( rotamer_id ),
	ids_of_clashing_nodes_(),
	polar_sc_atoms_not_satisfied_by_background_()
{}

//Destructor
AtomLevelHBondNode::~AtomLevelHBondNode() = default;

void AtomLevelHBondNode::print() const {
	TR << "AtomLevelHBondNode: rotamer_id:" << rotamer_id_ << " mres_id:" << mres_id_ << " num_edges: " << num_edges() << std::endl;
}

Size AtomLevelHBondNode::count_static_memory() const
{
	return sizeof( AtomLevelHBondNode );
}

Size AtomLevelHBondNode::count_dynamic_memory() const
{
	return utility::graph::LowMemNode::count_dynamic_memory()
		+ polar_sc_atoms_not_satisfied_by_background_.size() * sizeof( AtomInfo )
		+ ids_of_clashing_nodes_.size() * sizeof( unsigned int );
}


//Constructor
AtomLevelHBondEdge::AtomLevelHBondEdge( Size first_node_ind, Size second_node_ind ):
	utility::graph::LowMemEdge( first_node_ind, second_node_ind ),
	energy_( 0 ),
	hbonds_( 0 )
{}

AtomLevelHBondEdge::AtomLevelHBondEdge( Size first_node_ind, Size second_node_ind, Real energy ):
	utility::graph::LowMemEdge( first_node_ind, second_node_ind ),
	energy_( energy ),
	hbonds_( 0 )
{}

//Destructor
AtomLevelHBondEdge::~AtomLevelHBondEdge() = default;

Size AtomLevelHBondEdge::count_static_memory() const
{
	return sizeof( AtomLevelHBondEdge );
}

Size AtomLevelHBondEdge::count_dynamic_memory() const
{
	return utility::graph::LowMemEdge::count_dynamic_memory() + hbonds_.size() * sizeof( HBondInfo );
}


//Constructor
AtomLevelHBondGraph::AtomLevelHBondGraph()
{}


AtomLevelHBondGraph::AtomLevelHBondGraph( Size num_nodes ):
	PARENT( num_nodes )
{
}

//Destructor

AtomLevelHBondGraph::~AtomLevelHBondGraph()
{}


Size
AtomLevelHBondGraph::count_static_memory() const
{
	return sizeof( AtomLevelHBondGraph );
}

Size
AtomLevelHBondGraph::count_dynamic_memory() const
{
	//so basically we are not really overriding this at the moment
	return PARENT::count_dynamic_memory();
}

AtomLevelHBondEdge *
AtomLevelHBondGraph::register_hbond( Size rotamerA, Size rotamerB, Real score ) {
	AtomLevelHBondEdge * new_edge = add_edge( rotamerA, rotamerB );
	new_edge->set_energy( score );
	return new_edge;
}


} //graph
} //hbonds
} //scoring
} //core
