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
/// @author Jack Maguire

#include <core/scoring/hbonds/graph/HBondGraph.hh>

#include <basic/Tracer.hh>
#include <core/types.hh>


static basic::Tracer TR( "core.scoring.hbonds.graph.HBondGraph" );

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {


HBondNode & HBondNode_from_LowMemNode(utility::graph::LowMemNode & node) {
	return (HBondNode &)node;
}
HBondEdge & HBondEdge_from_LowMemEdge(utility::graph::LowMemEdge & edge) {
	return (HBondEdge &)(edge);
}



//dummy! please do not call these
HBondNode::HBondNode() :
	utility::graph::LowMemNode( 0 ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_(),
	polar_sc_atoms_not_satisfied_by_background_()
{ runtime_assert( false ); }

HBondNode::HBondNode( const HBondNode & ) :
	utility::graph::LowMemNode( 0 ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_(),
	polar_sc_atoms_not_satisfied_by_background_()
{ runtime_assert( false ); }
///////////


//Constructor
HBondNode::HBondNode( Size node_id ) :
	utility::graph::LowMemNode( node_id ),
	mres_id_( 0 ),
	rotamer_id_( 0 ),
	ids_of_clashing_nodes_(),
	polar_sc_atoms_not_satisfied_by_background_()
{}

HBondNode::HBondNode( NodeIDSize node_id, MResIDSize mres_id, RotamerIDSize rotamer_id ) :
	utility::graph::LowMemNode( node_id ),
	mres_id_( mres_id ),
	rotamer_id_( rotamer_id ),
	ids_of_clashing_nodes_(),
	polar_sc_atoms_not_satisfied_by_background_()
{}

//Destructor
HBondNode::~HBondNode() = default;

void HBondNode::print() const {
	TR << "HBondNode: rotamer_id:" << rotamer_id_ << " mres_id:" << mres_id_ << " num_edges: " << num_edges() << std::endl;
}

Size HBondNode::count_static_memory() const
{
	return sizeof( HBondNode );
}

Size HBondNode::count_dynamic_memory() const
{
	return utility::graph::LowMemNode::count_dynamic_memory()
		+ polar_sc_atoms_not_satisfied_by_background_.size() * sizeof( AtomInfo )
		+ ids_of_clashing_nodes_.size() * sizeof( unsigned int );
}

void
HBondNode::merge_data(
	HBondNode const & other,
	utility::vector1< Size > const & other_node_to_my_node,
	bool merge_with_OR_logic
) {
	debug_assert( other.moltenres() == moltenres() );
	debug_assert( other.local_rotamer_id() == local_rotamer_id() );

	for ( unsigned int other_new_clash_node : other.ids_of_clashing_nodes_ ) {
		Size my_new_clash_node = other_node_to_my_node[ other_new_clash_node ];
		if ( ! clashes( NodeIDSize( my_new_clash_node ) ) ) {
			//if ( ! clashes( my_new_clash_node ) ) {
			register_clash( NodeIDSize( my_new_clash_node ) );
		}
	}

	if ( merge_with_OR_logic ) {
		for ( AtomInfo const & info : other.polar_sc_atoms_not_satisfied_by_background() ) {
			add_polar_atom_if_doesnt_exist( info );
		}
	} else {
		//Merge with AND logic

		auto const & other_atoms = other.polar_sc_atoms_not_satisfied_by_background();

		using iter_type = AtomInfoSet::const_iterator;
		for ( iter_type iter = polar_sc_atoms_not_satisfied_by_background_.begin();
				iter != polar_sc_atoms_not_satisfied_by_background_.end(); ) {
			//remove if absent from other set
			if ( other_atoms.find( *iter ) == other_atoms.end() ) {
				iter = polar_sc_atoms_not_satisfied_by_background_.erase( iter );
			} else {
				++iter;
			}
		}

	}


}



//Constructor
HBondEdge::HBondEdge( Size first_node_ind, Size second_node_ind ):
	utility::graph::LowMemEdge( first_node_ind, second_node_ind ),
	energy_( 0 ),
	hbonds_( 0 )
{}

HBondEdge::HBondEdge( Size first_node_ind, Size second_node_ind, Real energy ):
	utility::graph::LowMemEdge( first_node_ind, second_node_ind ),
	energy_( energy ),
	hbonds_( 0 )
{}

//Destructor
HBondEdge::~HBondEdge() = default;

Size HBondEdge::count_static_memory() const
{
	return sizeof( HBondEdge );
}

Size HBondEdge::count_dynamic_memory() const
{
	return utility::graph::LowMemEdge::count_dynamic_memory() + hbonds_.size() * sizeof( HBondInfo );
}

void
HBondEdge::merge_data(
	HBondEdge const & other,
	utility::vector1< Size > const & other_node_to_my_node
) {
	const Size my_other_first_node_ind = other_node_to_my_node[ other.get_first_node_ind() ];

	bool swap_edge_direction = get_first_node_ind() != my_other_first_node_ind;

#ifndef NDEBUG
	if ( swap_edge_direction ) {
		debug_assert( get_second_node_ind() == my_other_first_node_ind );
		debug_assert( get_first_node_ind() == other_node_to_my_node[ other.get_second_node_ind() ] );
	} else {
		debug_assert( get_first_node_ind() == my_other_first_node_ind );
		debug_assert( get_second_node_ind() == other_node_to_my_node[ other.get_second_node_ind() ] );
	}
#endif

	for ( HBondInfo const & hb_info : other.hbonds() ) {
		HBondInfo new_info = hb_info;
		if ( swap_edge_direction ) new_info.first_node_is_donor( ! new_info.first_node_is_donor() );

		if ( std::find( hbonds_.begin(), hbonds_.end(), new_info) == hbonds_.end() ) {
			hbonds_.push_back( new_info );
		}
	}

	set_energy( std::max<float>( energy(), other.energy() ) );

}


//Constructor
HBondGraph::HBondGraph()
{}

HBondGraph::HBondGraph( Size num_nodes ):
	PARENT( num_nodes )
{
}


HBondGraph::~HBondGraph()
{}


Size
HBondGraph::count_static_memory() const
{
	return sizeof( HBondGraph );
}

Size
HBondGraph::count_dynamic_memory() const
{
	//so basically we are not really overriding this at the moment
	return PARENT::count_dynamic_memory();
}

HBondEdge *
HBondGraph::register_hbond( Size rotamerA, Size rotamerB, Real score ) {
	HBondEdge * new_edge = add_edge( rotamerA, rotamerB );
	new_edge->set_energy( score );
	return new_edge;
}

void
HBondGraph::merge( HBondGraph const & other, bool merge_nodes_with_OR_logic ) {

	// First we need to get a 1-1 mapping for the nodes
	utility::vector1< Size > other_node_to_my_node( other.num_nodes() );

	Size my_cur_node = 1;
	for ( Size other_cur_node = 1; other_cur_node <= other.num_nodes(); other_cur_node++ ) {

		HBondNode const * other_node = other.get_node( other_cur_node );
		debug_assert( other_node );
		const Size other_moltenres = other_node->moltenres();
		const Size other_rotamer_id = other_node->local_rotamer_id();

		while ( my_cur_node <= num_nodes() ) {

			HBondNode * my_node = get_node( my_cur_node );
			debug_assert( my_node );
			const Size my_moltenres = my_node->moltenres();
			const Size my_rotamer_id = my_node->local_rotamer_id();

			if ( my_moltenres != other_moltenres || my_rotamer_id != other_rotamer_id ) {
				my_cur_node++;
				continue;
			}

			// Now we have matching nodes
			break;

		}
		// This means we tried to merge a node we don't have
		runtime_assert( my_cur_node <= num_nodes() );

		other_node_to_my_node[ other_cur_node ] = my_cur_node;
		my_cur_node ++;
	}


	// Now that they're mapped, we can start merging

	for ( Size other_first_node = 1; other_first_node <= other.num_nodes(); other_first_node++ ) {

		Size my_first_node = other_node_to_my_node[ other_first_node ];
		HBondNode const * other_node = other.get_node( other_first_node );
		HBondNode * my_node = get_node( my_first_node );

		my_node->merge_data( *other_node, other_node_to_my_node, merge_nodes_with_OR_logic );

		for ( auto iter = other_node->const_edge_list_begin( other ), end = other_node->const_edge_list_end( other );
				iter != end;
				++iter ) {
			utility::graph::LowMemEdge const * other_edge = *iter;
			debug_assert( other_edge );
			const Size other_second_node = other_edge->get_other_ind( other_first_node );

			// We only want to add the edge in one-way
			if ( other_second_node < other_first_node ) continue;

			const Size my_second_node = other_node_to_my_node[ other_second_node ];

			utility::graph::LowMemEdge * my_equiv_edge = my_node->find_edge( my_second_node, *this );

			if ( ! my_equiv_edge ) {
				// Adding edges is always scary with the LowMemGraph. But we aren't holding any of our own Edge*
				my_equiv_edge = add_edge( my_first_node, my_second_node );
			}

			HBondEdge * my_equiv_casted_edge = static_cast< HBondEdge * >( my_equiv_edge );
			HBondEdge const * other_casted_edge = static_cast< HBondEdge const * >( other_edge );
			my_equiv_casted_edge->merge_data( *other_casted_edge, other_node_to_my_node );
		}
	}

}


} //graph
} //hbonds
} //scoring
} //core
