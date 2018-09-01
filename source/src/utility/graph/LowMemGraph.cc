// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/graph/LowMemGraph.cc
/// @brief A lower memory version of utility::graph::Graph with three key limitations
///        1. Due to std::vector::resize(), all of the Edge* can go invalid any time you call add_edge().
///        2. This doesn't have constant time deletion (a definite design goal for utility::graph::Graph)
///        3. Deleting an element doesn't actually delete it from memory
/// @author Brian Coventry (bcov@uw.edu)


#include <utility/graph/LowMemGraph.hh>

namespace utility {
namespace graph {

using platform::Size;




//////////////////////////// LowMemEdgeListIter //////////////////////////////

LowMemEdgeListIter::LowMemEdgeListIter()
: graph_(nullptr), owner_( nullptr ), cur_offset_( 0 ) {}

// This one is for iterating all the edges
LowMemEdgeListIter::LowMemEdgeListIter( LowMemGraphBase * graph, platform::Size cur_offset )
: graph_( graph ), owner_( nullptr ), cur_offset_( cur_offset )
{
	fast_forward();
}

// This one is for iterating the edges of a node
LowMemEdgeListIter::LowMemEdgeListIter( LowMemGraphBase * graph, LowMemNode * owner, platform::Size cur_offset )
: graph_( graph ), owner_( owner ), cur_offset_( cur_offset ) {}


LowMemEdgeListIter const &
LowMemEdgeListIter::operator ++ () {
	debug_assert( valid() );
	cur_offset_++;
	fast_forward();
	return *this;
}

LowMemEdgeListIter const &
LowMemEdgeListIter::operator -- () {
	cur_offset_--;
	fast_backward();
	debug_assert( valid() );
	return *this;
}


LowMemEdge *
LowMemEdgeListIter::operator * () const {
	debug_assert( valid() );
	if ( owner_ ) {
		return owner_->internal_get_edge( cur_offset_, *graph_ );
	} else {
		return graph_->internal_get_edge( cur_offset_ );
	}
}

LowMemEdge &
LowMemEdgeListIter::operator -> () const {
	debug_assert( valid() );
	if ( owner_ ) {
		return *(owner_->internal_get_edge( cur_offset_, *graph_ ));
	} else {
		return *(graph_->internal_get_edge( cur_offset_ ));
	}
}
bool
LowMemEdgeListIter::valid() const {
	if ( ! graph_ ) return false;
	if ( owner_ ) {
		return cur_offset_ < owner_->num_edges();
	} else {
		return cur_offset_ > 0 && cur_offset_ <= graph_->internal_edge_list_size();
	}
}

/// @brief The graph edge list can have deleted edges. This ensures we pass over them
void
LowMemEdgeListIter::fast_forward() {
	if ( owner_ ) return;       // Nodes can't have delete edges
	while ( cur_offset_ <= graph_->internal_edge_list_size() && (**this)->internal_deleted() ) cur_offset_++;
}

/// @brief This will allow cur_offset_ to go to 0. Which will throw an exception
///        later which is what we want if someone goes past begin.
void
LowMemEdgeListIter::fast_backward() {
	if ( owner_ ) return;       // Nodes can't have delete edges
	while ( cur_offset_ > 0 && (**this)->internal_deleted() ) cur_offset_--;
}

/////////////////////////// LowMemEdgeListConstIter //////////////////////////////


LowMemEdgeListConstIter::LowMemEdgeListConstIter()
: graph_(nullptr), owner_( nullptr ), cur_offset_( 0 ) {}

// This one is for iterating all the edges
LowMemEdgeListConstIter::LowMemEdgeListConstIter( LowMemGraphBase const * graph, platform::Size cur_offset )
: graph_( graph ), owner_( nullptr ), cur_offset_( cur_offset )
{
	fast_forward();
}

// This one is for iterating the edges of a node
LowMemEdgeListConstIter::LowMemEdgeListConstIter( LowMemGraphBase const * graph, LowMemNode const * owner, platform::Size cur_offset )
: graph_( graph ), owner_( owner ), cur_offset_( cur_offset ) {}



LowMemEdgeListConstIter const &
LowMemEdgeListConstIter::operator ++ () {
	debug_assert( valid() );
	cur_offset_++;
	fast_forward();
	return *this;
}
LowMemEdgeListConstIter const &
LowMemEdgeListConstIter::operator -- () {
	cur_offset_--;
	fast_backward();
	debug_assert( valid() );
	return *this;
}


LowMemEdge const *
LowMemEdgeListConstIter::operator * () const {
	debug_assert( valid() );
	if ( owner_ ) {
		return owner_->internal_get_edge( cur_offset_, *graph_ );
	} else {
		return graph_->internal_get_edge( cur_offset_ );
	}
}
LowMemEdge const &
LowMemEdgeListConstIter::operator -> () const {
	debug_assert( valid() );
	if ( owner_ ) {
		return *(owner_->internal_get_edge( cur_offset_, *graph_ ));
	} else {
		return *(graph_->internal_get_edge( cur_offset_ ));
	}
}
bool
LowMemEdgeListConstIter::valid() const {
	if ( ! graph_ ) return false;
	if ( owner_ ) {
		return cur_offset_ < owner_->num_edges();
	} else {
		return cur_offset_ > 0 && cur_offset_ <= graph_->internal_edge_list_size();
	}
}

/// @brief The graph edge list can have deleted edges. This ensures we pass over them
void
LowMemEdgeListConstIter::fast_forward() {
	if ( owner_ ) return;       // Nodes can't have delete edges
	while ( cur_offset_ <= graph_->internal_edge_list_size() && (**this)->internal_deleted() ) cur_offset_++;
}

/// @brief This will allow cur_offset_ to go to 0. Which will throw an exception
///        later which is what we want if someone goes past begin.
void
LowMemEdgeListConstIter::fast_backward() {
	if ( owner_ ) return;       // Nodes can't have delete edges
	while ( cur_offset_ > 0 && (**this)->internal_deleted() ) cur_offset_--;
}

///////////////////////// LowMemNode /////////////////////////////////////


void
LowMemNode::internal_add_edge( Size edge_offset ) {
	edge_vec_.push_back(edge_offset);
}


void
LowMemNode::drop_all_edges( LowMemGraphBase & graph ) {
	graph.drop_all_edges_for_node( get_node_index() );
}


LowMemEdge const *
LowMemNode::find_edge( uint32_t other_node_ind, LowMemGraphBase const & graph ) const {
	for ( Size offset : edge_vec_ ) {
		LowMemEdge const * edge = graph.internal_get_edge( offset );
		if ( edge->get_first_node_ind() == other_node_ind || edge->get_second_node_ind() == other_node_ind ) {
			return edge;
		}
	}
	return nullptr;
}
LowMemEdge *
LowMemNode::find_edge( uint32_t other_node_ind, LowMemGraphBase & graph ) {
	for ( Size offset : edge_vec_ ) {
		LowMemEdge * edge = graph.internal_get_edge( offset );
		if ( edge->get_first_node_ind() == other_node_ind || edge->get_second_node_ind() == other_node_ind ) {
			return edge;
		}
	}
	return nullptr;
}

// Returns local offset
Size
LowMemNode::internal_find_edge( uint32_t other_node_ind, LowMemGraphBase const & graph ) const {
	for ( Size local_offset = 0; local_offset < edge_vec_.size(); local_offset++ ) {   // starting at 0 on purpose
		Size offset = edge_vec_[ local_offset ];
		LowMemEdge const * edge = graph.internal_get_edge( offset );
		if ( edge->get_first_node_ind() == other_node_ind || edge->get_second_node_ind() == other_node_ind ) {
			return local_offset;
		}
	}
	runtime_assert( false );
}


LowMemEdge *
LowMemNode::internal_get_edge( Size local_offset, LowMemGraphBase & graph ) {
	return graph.internal_get_edge( edge_vec_[ local_offset ] );
}
LowMemEdge const *
LowMemNode::internal_get_edge( Size local_offset, LowMemGraphBase const & graph ) const {
	return graph.internal_get_edge( edge_vec_[ local_offset ] );
}


/// @brief memory accounting scheme
Size
LowMemNode::count_static_memory() const {
	return sizeof( LowMemNode );
}
/// @brief memory accounting scheme
Size
LowMemNode::count_dynamic_memory() const {
	return sizeof( Size ) * edge_vec_.capacity();
}

Size
LowMemNode::internal_drop_edge( LowMemEdge const * edge, LowMemGraphBase const & graph ) {
	uint32_t other_node_ind = edge->get_other_ind( get_node_index() );
	Size local_offset = internal_find_edge( other_node_ind, graph );
	Size global_offset = edge_vec_[ local_offset ];
	edge_vec_.erase( edge_vec_.begin() + local_offset );
	return global_offset;
}

///////////////////////////////// LowMemEdge ///////////////////////////////////////


uint32_t
LowMemEdge::get_other_ind( uint32_t node_ind ) const {
	if ( node_ind == first_node_ind_ ) return second_node_ind_;
	debug_assert( node_ind == second_node_ind_ );
	return first_node_ind_;
}

bool
LowMemEdge::same_edge( uint32_t node1, uint32_t node2 ) const {
	if ( node1 > node2 ) {
		uint32_t temp = node2;
		node2 = node1;
		node1 = temp;
	}
	return (node1 == first_node_ind_ && node2 == second_node_ind_);

}


LowMemNode const *
LowMemEdge::get_other_node( uint32_t node_index, LowMemGraphBase const & graph ) const {
	debug_assert( node_index == first_node_ind_ || node_index == second_node_ind_ );
	return first_node_ind_ == node_index ? graph.get_node( second_node_ind_ ) : graph.get_node( first_node_ind_ );
}

LowMemNode *
LowMemEdge::get_other_node( uint32_t node_index, LowMemGraphBase & graph ) {
	debug_assert( node_index == first_node_ind_ || node_index == second_node_ind_ );
	return first_node_ind_ == node_index ? graph.get_node( second_node_ind_ ) : graph.get_node( first_node_ind_ );
}

/// @brief memory accounting scheme
Size
LowMemEdge::count_static_memory() const {
	return sizeof( LowMemEdge );
}
/// @brief memory accounting scheme
Size
LowMemEdge::count_dynamic_memory() const {
	return 0;
}

bool
LowMemEdge::internal_deleted() const {
	debug_assert( ( first_node_ind_ == 0 || second_node_ind_ == 0 )
		? first_node_ind_ + second_node_ind_ == 0
		: true );
	return first_node_ind_ == 0;
}
void
LowMemEdge::internal_delete_self() {
	first_node_ind_ = 0;
	second_node_ind_ = 0;
}



template class LowMemGraph< LowMemNode, LowMemEdge >;  // Force instantiation of the default




}
}
