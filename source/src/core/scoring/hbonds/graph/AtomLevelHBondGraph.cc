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

void
AtomLevelHBondNode::merge_data(
	AtomLevelHBondNode const & other,
	utility::vector1< Size > const & other_node_to_my_node,
	bool merge_with_OR_logic
) {
	debug_assert( other.moltenres() == moltenres() );
	debug_assert( other.local_rotamer_id() == local_rotamer_id() );

	for ( unsigned int other_new_clash_node : other.ids_of_clashing_nodes_ ) {
		Size my_new_clash_node = other_node_to_my_node[ other_new_clash_node ];
		if ( ! clashes( my_new_clash_node ) ) {
			register_clash( my_new_clash_node );
		}
	}

	if ( merge_with_OR_logic ) {
		for ( AtomInfo const & info : other.polar_sc_atoms_not_satisfied_by_background() ) {
			add_polar_atom_if_doesnt_exist( info );
		}
	} else {
		utility::vector1< unsigned short int > local_atom_ids_to_remove;
		utility::vector1< AtomInfo > const & other_atoms = other.polar_sc_atoms_not_satisfied_by_background();

		auto other_iter = other_atoms.begin();
		auto other_end = other_atoms.end();
		auto my_iter = polar_sc_atoms_not_satisfied_by_background_.begin();
		auto my_end = polar_sc_atoms_not_satisfied_by_background_.end();

		while ( other_iter != other_end && my_iter != my_end ) {

			if ( *other_iter < *my_iter ) {
				other_iter++;
				continue;
			}
			if ( *my_iter < *other_iter ) {
				// If we get to here, this atom isn't in the other list
				local_atom_ids_to_remove.push_back( my_iter->local_atom_id() );
				my_iter++;
				continue;
			}
			// If we are here, the iters are pointing to the same atom
			other_iter++;
			my_iter++;
		}

		for ( unsigned short int local_atom_id : local_atom_ids_to_remove ) {
			remove_atom_info_stable( local_atom_id );
		}
	}


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

void
AtomLevelHBondEdge::merge_data(
	AtomLevelHBondEdge const & other,
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
AtomLevelHBondGraph::AtomLevelHBondGraph()
{}

AtomLevelHBondGraph::AtomLevelHBondGraph( Size num_nodes ):
	PARENT( num_nodes )
{
}


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

void
AtomLevelHBondGraph::merge( AtomLevelHBondGraph const & other, bool merge_nodes_with_OR_logic ) {

	// First we need to get a 1-1 mapping for the nodes
	utility::vector1< Size > other_node_to_my_node( other.num_nodes() );

	Size my_cur_node = 1;
	for ( Size other_cur_node = 1; other_cur_node <= other.num_nodes(); other_cur_node++ ) {

		AtomLevelHBondNode const * other_node = other.get_node( other_cur_node );
		debug_assert( other_node );
		const Size other_moltenres = other_node->moltenres();
		const Size other_rotamer_id = other_node->local_rotamer_id();

		while ( my_cur_node <= num_nodes() ) {

			AtomLevelHBondNode * my_node = get_node( my_cur_node );
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
		AtomLevelHBondNode const * other_node = other.get_node( other_first_node );
		AtomLevelHBondNode * my_node = get_node( my_first_node );

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

			AtomLevelHBondEdge * my_equiv_casted_edge = static_cast< AtomLevelHBondEdge * >( my_equiv_edge );
			AtomLevelHBondEdge const * other_casted_edge = static_cast< AtomLevelHBondEdge const * >( other_edge );
			my_equiv_casted_edge->merge_data( *other_casted_edge, other_node_to_my_node );
		}
	}

}


} //graph
} //hbonds
} //scoring
} //core
