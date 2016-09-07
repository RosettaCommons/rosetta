// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/EnergyGraph.cc
/// @brief  Energy graph class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/scoring/EnergyGraph.hh>

// Utility headers
#include <utility/vector1.hh>

// Boost Headers
#include <boost/pool/pool.hpp>
#include <core/graph/unordered_object_pool.hpp>

// C++ headers
#include <iostream>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

using namespace graph;

///////// Energy Node Class /////////////

EnergyNode::EnergyNode( Graph * owner, Size index ) :
	parent( owner, index )//, moved_( false )
{}

EnergyNode::~EnergyNode() = default;


void EnergyNode::print() const
{
	std::cout << "EnergyNode::print() deferring to parent::print()" << std::endl;
	parent::print();
}

/// @brief copy mmember data from source node
///
/// invoked by copy ctor and operator= methods from Graph base class
void EnergyNode::copy_from( parent const * source )
{
	//EnergyNode const * en_source = utility::down_cast< EnergyNode const * > ( source ); //nothing to copy, still -- want to assert the dynamic cast
	EnergyNode const * en_source = static_cast< EnergyNode const * > ( source );
	moved_ = en_source->moved_;
}

Size EnergyNode::count_static_memory() const
{
	return sizeof ( EnergyNode );
}

Size EnergyNode::count_dynamic_memory() const
{
	Size tot = 0;
	tot += parent::count_dynamic_memory(); //recurse to parent
	return tot;
}

bool EnergyNode::moved() const { return moved_; }
void EnergyNode::moved( bool setting ) { moved_ = setting; }


///////// Energy Edge Class /////////////

/// @brief Energy edge ctor
///
/// when first added to the graph, sets its enegies_not_yet_computed_ to true
/// the responsibility of marking the energies as having been computed falls on
/// class ScoringFunction
EnergyEdge::EnergyEdge(
	EnergyGraph * owner,
	Size n1,
	Size n2
) :
	parent( owner, n1, n2 ),
	energies_not_yet_computed_( true ),
	dsqr_( 0.0 ),
	array_( get_energy_owner()->array_pool().new_array() )
{}

EnergyEdge::EnergyEdge(
	EnergyGraph * owner,
	EnergyEdge const & example_edge
)
:
	parent( owner, example_edge.get_first_node_ind(), example_edge.get_second_node_ind() ),
	energies_not_yet_computed_( example_edge.energies_not_yet_computed_ ),
	//energy_map_( example_edge.energy_map_ )//,
	dsqr_( example_edge.dsqr_ ),
	array_( get_energy_owner()->array_pool().new_array() )
{
	array_.copy_array_contents( example_edge.array_ );
}

/// @brief virtual dstor; The EnergyEdge must free the array pool element it
/// holds before it disappears.
EnergyEdge::~EnergyEdge()
{
	get_energy_owner()->deallocate_arraypoolelement( array_ );
}

/// @brief copies data from EnergyEdge const * source;
///
/// called from the copy ctor and operator= methods defined in the Graph base class
void EnergyEdge::copy_from( parent const * source )
{
	EnergyEdge const * ee = static_cast< EnergyEdge const * > ( source );
	// down_cast is *supposed* to assert the dynamic cast in debug builds; doesn't work for some reason
	//EnergyEdge const * ee = utility::down_cast< EnergyEdge const * > ( source );
	energies_not_yet_computed_ = ee->energies_not_yet_computed_;
	dsqr_ = ee->dsqr_;
	array_.copy_array_contents( ee->array_ );
	return;
}


/// @brief set energies_computed_ to true for an edge
void EnergyEdge::mark_energies_computed() { energies_not_yet_computed_ = false; }

/// @brief set energies_computed_ to false for an edge
///
/// this function would allow edges to be marked as dirty without deleting
/// it from the graph.  Currently, edges are not marked as dirty; they are simply
/// deleted (and if still needed) returned to the graph.  If edge addition and
/// deletion starts to show up as a hotspot, we will examine if simply marking
/// edges as dirty is faster.  For now, this function is not used.
void EnergyEdge::mark_energies_uncomputed() { energies_not_yet_computed_ = true; }

/// @brief virtual call to determine the static size of an Edge object
/// dynamic memory use is counted through the recursive count_dynamic_memory()
/// calling path
Size EnergyEdge::count_static_memory() const
{
	return sizeof ( EnergyEdge );
}

/// @brief virtual call to determine the amount of dynamic memory allocated by
/// an edge; this function must recurse to the parent class to determine how
/// much memory the parent class is responsible for.  Do not account for
/// the size of the ArrayPool array here; instead, that is accounted for in
/// the EnergyGraph::count_dynamic_memory method.
Size EnergyEdge::count_dynamic_memory() const
{
	Size tot = 0;
	tot += parent::count_dynamic_memory(); //recurse to parent
	return tot;
}


///////// Energy Graph Class /////////////

EnergyEdge * EnergyGraph::find_energy_edge( Size n1, Size n2)
{
	Edge * edge( find_edge( n1, n2 ) );
	if ( edge ) {
		return static_cast< EnergyEdge * > ( edge );
		//return utility::down_cast< EnergyEdge * > ( edge );
	} else {
		return nullptr;
	}
}

EnergyEdge const * EnergyGraph::find_energy_edge( Size n1, Size n2) const
{
	Edge const * edge( find_edge( n1, n2 ) );
	if ( edge ) {
		return static_cast< EnergyEdge const * > ( edge );
		//return utility::down_cast< EnergyEdge const * > ( edge );
	} else {
		return nullptr;
	}
}


EnergyGraph::EnergyGraph()
:
	parent(),
	energy_edge_pool_( new boost::unordered_object_pool< EnergyEdge > ( 256 ) ),
	energy_array_pool_( 256 ),
	score_type_2_active_( n_shortranged_2b_score_types, -1 )
{}

/// @details This does not call the base class parent( Size ) constructor since
/// that produces calls to the polymorphic function create_new_node() and polymorphism
/// does not work during constructor intialization.
EnergyGraph::EnergyGraph( Size num_nodes )
:
	parent(),
	energy_edge_pool_( new boost::unordered_object_pool< EnergyEdge > ( 256 ) ),
	energy_array_pool_( 256 ),
	score_type_2_active_( n_shortranged_2b_score_types, -1 )
{
	set_num_nodes( num_nodes );
}

/// @details Notice that this does not call the parent( src ) copy constructor.
/// This is because the copy constructor relies on polymorphic functions which
/// are unavailable during the Graph constructor.  Instead, this function waits
/// until parent construction is complete, and relies on the assigmnent operator.
EnergyGraph::EnergyGraph( EnergyGraph const & src )
:
	parent( ),
	energy_edge_pool_( new boost::unordered_object_pool< EnergyEdge > ( 256 ) ),
	energy_array_pool_( 256 ),
	score_type_2_active_( n_shortranged_2b_score_types, -1 )
{
	active_score_types( src.active_2b_score_types_ );
	parent::operator = ( src );
}

EnergyGraph::~EnergyGraph() {
	delete_everything();
	delete energy_edge_pool_; energy_edge_pool_ = nullptr;
}

void
EnergyGraph::add_energy_edge(
	Size index1,
	Size index2,
	DistanceSquared dsq
)
{
	Edge* newedge = add_edge( index1, index2 );
	(static_cast< EnergyEdge * > (newedge))->square_distance( dsq );
}

//void
//EnergyGraph::add_energy_edge(
// Size index1,
// Size index2,
// EnergyMap const & emap
//)
//{
// Edge * newedge = add_edge( index1, index2 );
// (static_cast< EnergyEdge * > ( newedge ))->store_active_energies( emap );
//}


/// @brief assignment operator -- performs a deep copy
EnergyGraph &
EnergyGraph::operator = ( EnergyGraph const & rhs )
{
	active_score_types( rhs.active_2b_score_types_ );
	parent::operator = ( rhs );
	return *this;
}

/// @details The array pool only needs to be resized if the number
/// of active score types has changed; it doesn't have to be
/// resized if the actual active types have changed but the number
/// of active score types has remained constant.  However, if either
/// the number or the identity of the active score types have changed,
/// then all old edges in the EnergyGraph should be dropped.
/// The input array "active" must be sorted.
bool
EnergyGraph::active_score_types( ScoreTypes const & active )
{
	for ( Size ii = 1; ii < active.size(); ++ii ) {
		debug_assert( active[ ii ] < active[ ii + 1 ] );
	}
	bool clear_edges = false;
	if ( active.size() != active_2b_score_types_.size() ) clear_edges = true;

	bool const resize_array_pool = clear_edges;

	//ScoreTypes sactive = active;
	//std::sort( sactive.begin(), sactive.end() );
	if ( ! clear_edges ) {
		for ( Size ii = 1; ii <= active.size(); ++ii ) {
			if ( active_2b_score_types_[ ii ] != active[ ii ] ) {
				clear_edges = true;
				break;
			}
		}
	}

	if ( clear_edges ) {
		drop_all_edges();
		/// Update the mapping between active score types and their positions in the ArrayPoolElement arrays.
		active_2b_score_types_ = active;
		std::fill( score_type_2_active_.begin(), score_type_2_active_.end(), -1 );
		int count_active( 0 );
		for ( Size ii = 1; ii <= active_2b_score_types_.size(); ++ii ) {
			score_type_2_active_[ active_2b_score_types_[ ii ] ] = count_active;
			++count_active;
		}
	}

	if ( resize_array_pool ) {
		energy_array_pool_.set_array_size( active.size() );
	}

	return ! clear_edges;
}

void EnergyGraph::delete_edge( graph::Edge * edge )
{
	debug_assert( dynamic_cast< EnergyEdge* > (edge) );
	energy_edge_pool_->destroy( static_cast< EnergyEdge* > (edge) );
}

void EnergyGraph::deallocate_arraypoolelement(
	graph::ArrayPoolElement< Real > & element
)
{
	energy_array_pool_.deallocate_array( element );
}

Size EnergyGraph::count_static_memory() const
{
	return sizeof ( EnergyGraph );
}

Size EnergyGraph::count_dynamic_memory() const
{
	Size tot = energy_array_pool_.footprint();
	tot += parent::count_dynamic_memory(); //recurse to parent
	return tot;
}

/// @brief Factory method for node creation
Node*
EnergyGraph::create_new_node( Size index )
{
	return new EnergyNode( this, index );
}

/// @brief Factory method for edge creation
Edge*
EnergyGraph::create_new_edge( Size index1, Size index2 )
{
	return energy_edge_pool_->construct( this, index1, index2 );
}

/// @brief Factory copy-constructor method for edge creation
Edge*
EnergyGraph::create_new_edge( Edge const * example_edge )
{
	return energy_edge_pool_->construct(
		this,
		static_cast< EnergyEdge const & > (*example_edge)
	);
}

#ifdef    SERIALIZATION
template < class Archive >
void EnergyNode::save_to_archive( Archive & arc ) const
{
	arc( moved_ );
}

template < class Archive >
void EnergyNode::load_from_archive( Archive & arc )
{
	arc( moved_ );
}

template < class Archive >
void EnergyEdge::save_to_archive( Archive & archive ) const
{
  archive( energies_not_yet_computed_, dsqr_ );
  std::vector< Real > energies( array_.size() );
  for ( Size ii = 0; ii < array_.size(); ++ii ) {
    energies[ ii ] = array_[ ii ];
  }
  archive( energies );
}

template < class Archive >
void EnergyEdge::load_from_archive( Archive & archive )
{
  archive( energies_not_yet_computed_, dsqr_ );
  std::vector< Real > energies( array_.size() );
  archive( energies );
  for ( Size ii = 0; ii < array_.size(); ++ii ) {
    array_[ ii ] = energies[ ii ];
  }
}

template < class Archive >
void EnergyGraph::save( Archive & arc ) const
{
  arc( active_2b_score_types_ );

  arc( num_nodes() );
  for ( Size ii = 1; ii <= num_nodes(); ++ii ) {
    get_energy_node( ii )->save_to_archive( arc );
  }
  arc( num_edges() );
  for ( EdgeListConstIter iter = const_edge_list_begin(), iter_end = const_edge_list_end(); iter != iter_end; ++iter ) {
    arc( (*iter)->get_first_node_ind(), (*iter)->get_second_node_ind() );
		EnergyEdge const * eeptr = static_cast< EnergyEdge const * > (*iter);
    eeptr->save_to_archive( arc );
  }
	// The pools are not serialized;
	// EXEMPT energy_edge_pool_ energy_array_pool_
	// score_type_2_active_ is essentially derived from active_2b_score_types_
	// EXEMPT score_type_2_active_

}

template < class Archive >
void EnergyGraph::load( Archive & arc )
{
  ScoreTypes sts; arc( sts );
  active_score_types( sts );
	// The above line stores the active score types, and that in turn initializes score_type_2_active_
	// EXEMPT active_2b_score_types_ score_type_2_active_

  Size num_nodes(0); arc( num_nodes );
  set_num_nodes( num_nodes );

  for ( Size ii = 1; ii <= num_nodes; ++ii ) {
    get_energy_node( ii )->load_from_archive( arc );
  }

  Size num_edges(0); arc( num_edges );
  for ( Size ii = 1; ii <= num_edges; ++ii ) {
    Size node1(0), node2(0); arc( node1, node2 );
    Edge * new_edge = add_edge( node1, node2 );
		EnergyEdge * eeptr = static_cast< EnergyEdge * > ( new_edge );
    eeptr->load_from_archive( arc );
  }
	// The pools are not serialized;
	// EXEMPT energy_edge_pool_ energy_array_pool_
}

SAVE_AND_LOAD_SERIALIZABLE( EnergyGraph );
#endif // SERIALIZATION

} //namespace scoring
} //namespace core

#ifdef    SERIALIZATION
CEREAL_REGISTER_TYPE( core::scoring::EnergyGraph )
CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_EnergyGraph )
#endif // SERIALIZATION
