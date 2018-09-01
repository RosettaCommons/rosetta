// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/graph/LowMemGraph.hh
/// @brief A lower memory version of utility::graph::Graph with three key limitations
///        1. Due to std::vector::resize(), all of the LMEdge* can go invalid any time you call add_edge().
///        2. This doesn't have constant time deletion (a definite design goal for utility::graph::Graph)
///        3. Deleting an element doesn't actually delete it from memory
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_utility_graph_LowMemGraph_HH
#define INCLUDED_utility_graph_LowMemGraph_HH

#include <utility/graph/LowMemGraph.fwd.hh>


#include <utility/pointer/memory.hh>
#include <platform/types.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/backtrace.hh>
#include <utility/vector1.hh>

#include <vector>

namespace utility {
namespace graph {


/// @brief Non-const iterator for edges
class LowMemEdgeListIter {

public:

	/// @brief Default constructor. Don't use this
	LowMemEdgeListIter();

	/// @brief This one is for iterating all the edges
	LowMemEdgeListIter( LowMemGraphBase * graph, platform::Size cur_offset );

	// @brief This one is for iterating the edges of a node
	LowMemEdgeListIter( LowMemGraphBase * graph, LowMemNode * owner, platform::Size cur_offset );

	/// @brief increment operator.  Point this iterator at the next element in the list.
	LowMemEdgeListIter const & operator ++ ();

	/// @brief decrement operator.  Point this iterator at the previous element in the list.
	LowMemEdgeListIter const & operator -- ();

	/// @brief equality operator.  Do these elements point at the same list element?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	bool operator == ( LowMemEdgeListIter const & rhs ) const {
		debug_assert( owner_ == rhs.owner_ );
		debug_assert( graph_ == rhs.graph_ );
		return cur_offset_ == rhs.cur_offset_;
	}

	/// @brief inequality operator.  Do these elements point to different elements from the same list?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	bool operator != ( LowMemEdgeListIter const & rhs ) const {
		return ! ( operator == ( rhs ) );
	}

	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a non-const edge pointer, which defines this class
	/// as a non-const iterator.  There are no methods, though, to access the list element itself.
	LowMemEdge * operator * () const;

	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a non-const edge pointer, which defines this class
	/// as a non-const iterator.  There are no methods, though, to access the list element itself.
	LowMemEdge & operator -> () const;

	/// @brief check that this iterator is valid.  Will not guarantee that the iterator points at an element
	/// that has not been destroyed.
	bool valid() const;

protected:
	void fast_forward();
	void fast_backward();

private:
	LowMemGraphBase * graph_;
	LowMemNode * owner_;
	platform::Size cur_offset_;
};

/// @brief Const iterator for edges
class LowMemEdgeListConstIter {

public:

	/// @brief Default constructor. Don't use this
	LowMemEdgeListConstIter();

	/// @brief This one is for iterating all the edges
	LowMemEdgeListConstIter( LowMemGraphBase const * graph, platform::Size cur_offset );

	// @brief This one is for iterating the edges of a node
	LowMemEdgeListConstIter( LowMemGraphBase const * graph, LowMemNode const * owner, platform::Size cur_offset );

	/// @brief increment operator.  Point this iterator at the next element in the list.
	LowMemEdgeListConstIter const & operator ++ ();

	/// @brief decrement operator.  Point this iterator at the previous element in the list.
	LowMemEdgeListConstIter const & operator -- ();

	/// @brief equality operator.  Do these elements point at the same list element?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	bool operator == ( LowMemEdgeListConstIter const & rhs ) const {
		debug_assert( owner_ == rhs.owner_ );
		debug_assert( graph_ == rhs.graph_ );
		return cur_offset_ == rhs.cur_offset_;
	}

	/// @brief inequality operator.  Do these elements point to different elements from the same list?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	/// Elements need not be valid to be compared.
	bool operator != ( LowMemEdgeListConstIter const & rhs ) const {
		return ! ( operator == ( rhs ) );
	}

	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a const edge pointer, which defines this class
	/// as a const iterator.  There are no methods, of course, to access the list element itself.
	LowMemEdge const * operator * () const;

	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a const edge pointer, which defines this class
	/// as a const iterator.  There are no methods, of course, to access the list element itself.
	LowMemEdge const & operator -> () const;

	/// @brief check that this iterator is valid.  Will not guarantee that the iterator points at an element
	/// that has not been destroyed.
	bool valid() const;

protected:
	void fast_forward();
	void fast_backward();

private:
	LowMemGraphBase const * graph_;
	LowMemNode const * owner_;
	platform::Size cur_offset_;
};


/// @brief An Node class for LowMemGraph. Will often be overriden
/// @detail Be careful with this class! It doesn't use actual virtual
///         functions. Never do this:
///             LowMemNodeOP op = MyDerivedNodeOP()
///         Instead, if you wish to use an owning pointer, you must
///         do this:
///             MyDerivedNodeOP op = MyDerivedNodeOP()
///
///         Having any virtual functions increases a class's footprint by 8 bytes
///
/// @remarks If writing a derived class, "override" these four functions:
///            ~LowMemNode()
///            print()
///            count_static_memory()
///            count_dynamic_memory()
///
///          See core/scoring/hbonds/graph/AtomLevelHBondGraph.hh for example
///
class LowMemNode {

public:
	LowMemNode( uint32_t node_index )
	: node_index_( node_index ) {}

	/// @brief This class should have a virtual destructor. But adding virtual functions
	///        adds an extra 8 bytes per class. Be careful!!!
	/* virtual */ ~LowMemNode() {}

	/// @brief deletes all edges incident upon this node
	/// @details Although edges are deleted, the memory for the edges is not freed.
	void drop_all_edges( LowMemGraphBase & graph );

	/// @brief a "slow" (linear) search for an edge.
	LowMemEdge const * find_edge( uint32_t other_node_ind, LowMemGraphBase const & graph ) const;

	/// @brief a "slow" (linear) search for an edge.
	LowMemEdge * find_edge( uint32_t other_node_ind, LowMemGraphBase & graph );

	/// @brief send summaray data about this node to the screen
	/// @details This function isn't technically virtual, but if you override
	///           it, your code *WILL* be called by print_vertices()
	/* virtual */ void print() const {}


	/// @brief returns a non-const iterator to the beginning of its edge list
	LowMemEdgeListIter edge_list_begin( LowMemGraphBase & graph ) { return LowMemEdgeListIter( &graph, this, 0 ); }
	/// @brief returns a const iterator to the beginning of its edge list
	LowMemEdgeListConstIter const_edge_list_begin( LowMemGraphBase const & graph ) const { return LowMemEdgeListConstIter( &graph, this, 0 ); }

	/// @brief returns a non-const iterator to the end of its edge list
	LowMemEdgeListIter edge_list_end( LowMemGraphBase & graph ) { return LowMemEdgeListIter( &graph, this, num_edges() ); }
	/// @brief returns a const iterator to the end of its edge list
	LowMemEdgeListConstIter const_edge_list_end( LowMemGraphBase const & graph ) const { return LowMemEdgeListConstIter( &graph, this, num_edges() ); }

	/// @brief the index for this node
	uint32_t get_node_index() const { return node_index_; }

	/// @brief the number of edges incident on this node, which may include a loop edge
	platform::Size num_edges() const {
		return edge_vec_.size();
	}


	/// @brief how much memory is statically allocated by this node
	/// @details This function isn't technically virtual, but if you override
	///           it, your code *WILL* be called by getTotalMemoryUsage()
	/* virtual */ platform::Size count_static_memory() const;

	/// @brief how much memory is dynamically allocated by this node -- must be recursively invoked
	///        by a derived class.
	/// @details This function isn't technically virtual, but if you override
	///           it, your code *WILL* be called by getTotalMemoryUsage()
	/* virtual */ platform::Size count_dynamic_memory() const;


	/// @brief Don't call this!!!. Adds an edge from the edge list to this node
	/// @detail Only made public so we don't have to friend the templated LowMemGraph
	void internal_add_edge( platform::Size edge_offset );

	/// @brief Don't call this!!!. It only locally deletes the edge from this node.
	/// @detail Only made public so we don't have to friend the templated LowMemGraph.
	///         Returns offset in the global edge list
	platform::Size internal_drop_edge( LowMemEdge const * edge, LowMemGraphBase const & graph );

	/// @brief Don't call this!!!. Only to be used by LowMemGraph
	/// @detail Only made public so we don't have to friend the templated LowMemGraph
	std::vector<platform::Size> const & internal_get_edge_list( ) const { return edge_vec_; }

protected:

	LowMemEdge * internal_get_edge( platform::Size local_offset, LowMemGraphBase & graph );
	LowMemEdge const * internal_get_edge( platform::Size local_offset, LowMemGraphBase const & graph ) const;

	// Returns local_offset
	platform::Size internal_find_edge( uint32_t other_node_ind, LowMemGraphBase const & graph ) const;


public:

	friend class LowMemEdgeListIter;
	friend class LowMemEdgeListConstIter;

private:
	// Total size of member variables
	// 24 + 4 + (4) = 32

	std::vector<platform::Size> edge_vec_;   // sizeof(std::vector) == 24 ; sizeof(utility::vector1) == 32
	uint32_t node_index_;
	// Wasting 4 bytes here because classes have to be 8 byte aligned
};


/// @brief An Edge class for LowMemGraph. Will often be overriden.
/// @detail Be careful with this class! It doesn't use actual virtual
///         functions. Never do this:
///             LowMemEdgeOP op = MyDerivedEdgeOP()
///         Instead, if you wish to use an owning pointer, you must
///         do this:
///             MyDerivedEdgeOP op = MyDerivedEdgeOP()
///
///         Having any virtual functions increases a class's footprint by 8 bytes
///
/// @remarks If writing a derived class, "override" these three functions:
///            ~LowMemNode()
///            count_static_memory()
///            count_dynamic_memory()
///
///          See core/scoring/hbonds/graph/AtomLevelHBondGraph.hh for example
///
class LowMemEdge {

public:

	/// @brief Main edge constructor.  This should only be invoked by create_new_edge, which
	/// itself is only called by add_edge.  The ONLY way an edge should be added to a graph
	/// is through add_edge.  NOTE: edges should be only be deleted by a call to the Graph's
	/// delete_edge method
	LowMemEdge( uint32_t first_node_ind, uint32_t second_node_ind )
	: first_node_ind_( first_node_ind ), second_node_ind_( second_node_ind ) {}

	/// @brief This class should have a virtual destructor. But adding virtual functions
	///        adds an extra 8 bytes per class. Be careful!!!
	~LowMemEdge() {}

	/// @brief returns the index of the one node given the index of the other.
	///        node_index must be one of the two nodes that this edge is incident upon.
	uint32_t get_other_ind( uint32_t node_ind ) const;

	/// @brief returns a const pointer to one node given the index of the other.
	///        node_index must be one of the two nodes that this edge is incident upon.
	LowMemNode const * get_other_node( uint32_t node_index, LowMemGraphBase const & graph ) const;

	/// @brief returns a non-const pointer to one node given the index of the other.
	///        node_index must be one of the two nodes that this edge is incident upon.
	LowMemNode * get_other_node( uint32_t node_index, LowMemGraphBase & graph );

	/// @brief returns the index of the lower node
	uint32_t get_first_node_ind() const { return first_node_ind_; }

	/// @brief returns the index of the upper node
	uint32_t get_second_node_ind() const { return second_node_ind_; }

	/// @brief Is this the same edge as another edge (node1,node2)?  Note:
	///        this graph does not work for multi-graphs.  Edges must be unique.
	bool same_edge( uint32_t node1, uint32_t node2 ) const;


	/// @brief how much memory is statically allocated by this edge
	/// @details This function isn't technically virtual, but if you override
	///          it, your code *WILL* be called by getTotalMemoryUsage()
	/* virtual */ platform::Size count_static_memory() const;

	/// @brief how much memory is dynamically allocated by this edge -- must be recursively invoked
	///        by a derived class.
	/// @details This function isn't technically virtual, but if you override
	///          it, your code *WILL* be called by getTotalMemoryUsage()
	/* virtual */ platform::Size count_dynamic_memory() const;

	/// @brief Don't call this!!! Only to be called by LowMemGraph
	void internal_delete_self();

	friend class LowMemEdgeListIter;
	friend class LowMemEdgeListConstIter;

protected:
	bool internal_deleted() const;

private:
	// Total size of member variables
	// 4 + 4 = 8

	uint32_t first_node_ind_;
	uint32_t second_node_ind_;

};

/// @brief Pure virtual baseclass that was required to avoid templating Edges and Nodes
class LowMemGraphBase : public utility::pointer::ReferenceCount {

public:
	LowMemGraphBase() {}

	virtual ~LowMemGraphBase() {}

	virtual LowMemNode const * get_node( uint32_t index ) const = 0;

	virtual LowMemNode * get_node( uint32_t index ) = 0;

	virtual void drop_all_edges_for_node( uint32_t index ) = 0;

	virtual platform::Size internal_edge_list_size() const = 0;

	friend class LowMemEdgeListIter;
	friend class LowMemEdgeListConstIter;
	friend class LowMemNode;

protected:


	virtual LowMemEdge const * internal_get_edge( platform::Size offset ) const = 0;

	virtual LowMemEdge * internal_get_edge( platform::Size offset ) = 0;

};




// This typedef is in the forward declaration
//typedef LowMemGraph< LowMemNode, LowMemEdge > DefaultLowMemGraph;

/// @brief A graph with low memory use and constant time edge insertion. Extensible.
/// @detail For general use, use utility::graph::DefaultLowMemGraph.
//          Limited to 2^32 nodes and 2^64 edges.
///         Edges are 8 bytes and nodes are 32 bytes + 8 bytes for each connected edge.
///
///         Due to an implementation detail, adding an edge invalidates all
///         edge pointers. Be sure not to hold any Edge* when adding new edges.
///
///         Deleting edges does not decrease memory use. Memory may only be recaptured
///         by calling set_num_edges() or delete_everything()
///
///
///     If one wishes to inherit LowMemGraph, the correct way to inherit is:
///         class MyClass : public LowMemGraph<MyNode,MyEdge> {
///             typedef LowMemGraph<MyNode,MyEdge> PARENT;
///     Then only count_static_memory and count_dynamic_memory need to be overwritten.
///
///     See core/scoring/hbonds/graph/AtomLevelHBondGraph.hh for example on how to inherit.
///
/// @remarks This class saves an extra 8 bytes per edge by not storing a list of pointers
///           to a pool, but instead using templates. Templates cause all sorts
///           of problems and therefore their use here was limited. Unless absolutely
///           necessary, use the base classes instead of the template types. Under almost
///           all circumstances, it is possible to fully extend this class without introducing
///           any more templates. Look at LowMemGraphBase for instance. This class allows
///           LowMemEdge and LowMemNode to interact with LowMemGraph without themselves
///           needing to be templated.
///
///           Additionally, the entire implementation of LowMemGraph is left in the .hh file.
///           Templates are weird and putting it in the .cc file is likely to cause weird
///           linking errors.
///
template< class _LMNode, class _LMEdge >
class LowMemGraph : public LowMemGraphBase {

public:
	typedef LowMemGraph<_LMNode, _LMEdge> THIS;
	typedef _LMNode LMNode;
	typedef _LMEdge LMEdge;

	//https://stackoverflow.com/questions/44158904/static-assert-for-public-inheritance
	static_assert( std::is_convertible< _LMNode*, LowMemNode* >::value, "LMNode class must inherit LowMemNode as public" );
	static_assert( std::is_convertible< _LMEdge*, LowMemEdge* >::value, "LMEdge class must inherit LowMemEdge as public" );

	/// @brief ctor
	LowMemGraph() {}

	/// @brief num nodes ctor
	/// @details This does not suffer from the same limitations as utility::graph::Graph
	LowMemGraph( platform::Size num_nodes ) {
		set_num_nodes( num_nodes );
	}

	/// @brief the number of nodes in the graph
	platform::Size num_nodes() const {
		return nodes_.size();
	}

	/// @brief returns a pointer to the edge connecting nodes node1 and node2, if that edge exists
	/// in the graph, o.w. returns 0. Does not focus the graph on this edge
	LMEdge * find_edge( platform::Size node1, platform::Size node2 ) {
		if ( node1 > nodes_.size() ) return nullptr;
		return static_cast< LMEdge * > ( nodes_[ node1 ].find_edge( node2, *this ) );
	}

	/// @brief returns a const pointer to the edge connecting nodes node1 and node2, if that edge exists
	/// in the graph, o.w. returns 0. Does not focus the graph on this edge
	LMEdge const * find_edge( platform::Size node1, platform::Size node2 ) const {
		if ( node1 > nodes_.size() ) return nullptr;
		return static_cast< LMEdge const * > ( nodes_[ node1 ].find_edge( node2, *this ) );
	}

	/// @brief deallocate all nodes and edges from the graph
	void delete_everything() {
		nodes_.clear();
		edges_.clear();
		deleted_edges_offsets_.clear();
	}

	/// @brief Get a const pointer to a node
	LMNode const * get_node( uint32_t index ) const override {
		debug_assert( index > 0 &&  index <= nodes_.size() );
		return &nodes_[ index ];
	}

	/// @brief Get a non-const pointer to a node
	LMNode * get_node( uint32_t index ) override {
		debug_assert( index > 0 &&  index <= nodes_.size() );
		return &nodes_[ index ];
	}
	/// @brief calls print() on each of the nodes in the graph
	/// @details Will call "overriden" print in derived nodes
	void print_vertices() const {
		for ( platform::Size ii = 1; ii <= nodes_.size(); ii++ ) {
			// This calls the "virtual" functions becuase we aren't working with the base class
			nodes_[ii].print();
		}
	}

	/// @brief How many edges are there in the graph?
	platform::Size num_edges() const { return edges_.size() - deleted_edges_offsets_.size(); }

protected:


	uint32_t internal_create_new_node( uint32_t node_ind ) {
		nodes_.emplace_back( node_ind );
		return nodes_.size();
	}
	template< typename... Args >
	platform::Size internal_create_new_edge( uint32_t index1, uint32_t index2, Args&&... args ) {
		if ( deleted_edges_offsets_.size() == 0 ) {
			edges_.emplace_back( index1, index2, std::forward<Args>( args )... );
			return edges_.size();
		} else {
			Size offset = deleted_edges_offsets_.back();
			deleted_edges_offsets_.pop_back();
			LMEdge edge( index1, index2, std::forward<Args>( args )... );
			edges_[ offset ] = edge;
			return offset;
		}
	}

	LMEdge const * internal_get_edge( platform::Size offset ) const override {
		debug_assert( offset > 0 &&  offset <= edges_.size() );
		return &edges_[ offset ];
	}

	LMEdge * internal_get_edge( platform::Size offset ) override {
		debug_assert( offset > 0 &&  offset <= edges_.size() );
		return &edges_[ offset ];
	}


public:
	/// @brief How many slots are there internally. Used for unit testing.
	platform::Size internal_edge_list_size() const override {
		return edges_.size();
	}

	/// @brief set the number of nodes in the graph -- deletes any existing edges and nodes in the graph
	void set_num_nodes( uint32_t num_nodes ) {
		delete_everything();
		nodes_.reserve( num_nodes );
		for ( uint32_t ii = 1; ii <= num_nodes; ii++ ) internal_create_new_node( ii );
	}

	/// @brief add an edge between two vertices. Returns a pointer to the edge after its been added,
	///        allowing the calling function to immediately set data for this edge.
	/// @detail WARNING!!! After adding an edge, other edge pointers are not guarenteed to be valid!!
	template< typename... Args >
	LMEdge * add_edge( uint32_t node1, uint32_t node2, Args&&... args ) {
		debug_assert( ! get_edge_exists( node1, node2 ) );
		platform::Size temp = node1 < node2 ? node1 : node2;
		node2 = node1 < node2 ? node2 : node1;
		node1 = temp;
		platform::Size edge_offset = internal_create_new_edge( node1, node2, std::forward<Args>( args )... );
		nodes_[ node1 ].internal_add_edge( edge_offset );
		nodes_[ node2 ].internal_add_edge( edge_offset );
		return internal_get_edge( edge_offset );
	}

	/// @brief add an edge between two vertices. Returns a pointer to the edge after its been added,
	///        allowing the calling function to immediately set data for this edge.
	/// @detail WARNING!!! After adding an edge, other edge pointers are not guarenteed to be valid!!
	LMEdge * add_edge( LowMemEdge const * example_edge ) {
		return add_edge( example_edge->get_first_node_ind(), example_edge->get_second_node_ind() );
	}

	/// @brief is an edge already present in the graph? O(V) worst case.  O(1) iff all vertices have O(1) edges
	bool get_edge_exists( uint32_t node1, uint32_t node2) const {
		LowMemEdge const * edge = find_edge( node1, node2 );
		return ( edge != nullptr );
	}


	/// @brief returns a non-const iterator to the beginning of the (unordered) edge list for the graph.
	/// this edge list contains all the edges in the graph, not simply those for a particular vertex
	LowMemEdgeListIter edge_list_begin() { return LowMemEdgeListIter( this, 1 ); }

	/// @brief returns a const iterator to the beginning of the (unordered) edge list for the graph.
	/// this edge list contains all the edges in the graph, not simply those for a particular vertex
	LowMemEdgeListConstIter const_edge_list_begin() const { return LowMemEdgeListConstIter( this, 1 ); }

	/// @brief returns a non-const iterator to the end of the (unordered) edge list for the graph.
	/// this edge list contains all the edges in the graph, not simply those for a particular vertex
	LowMemEdgeListIter edge_list_end() { return LowMemEdgeListIter( this, internal_edge_list_size()+1 ); }

	/// @brief returns a const iterator to the end of the (unordered) edge list for the graph.
	/// this edge list contains all the edges in the graph, not simply those for a particular vertex
	LowMemEdgeListConstIter const_edge_list_end() const { return LowMemEdgeListConstIter( this, internal_edge_list_size()+1 ); }


	/// @brief remove an edge from the graph
	void delete_edge( LowMemEdge * edge ) {
		debug_assert( get_edge_exists( edge->get_first_node_ind(), edge->get_second_node_ind() ) );
		Size offset = nodes_[ edge->get_first_node_ind() ].internal_drop_edge( edge, *this );
		nodes_[ edge->get_second_node_ind() ].internal_drop_edge( edge, *this );
		edge->internal_delete_self();
		deleted_edges_offsets_.push_back( offset );
	}

	/// @brief delete all the edges for a single vertex in the graph
	/// @details This function will not free memory associated with the deleted edges
	void drop_all_edges_for_node( uint32_t index ) override {
		debug_assert( index > 0 &&  index <= nodes_.size() );
		std::vector<platform::Size> edge_list = nodes_[index].internal_get_edge_list();
		for ( platform::Size offset : edge_list ) delete_edge( &edges_[offset] );
	}

	/// @brief delete all the edges present in the graph
	/// @details This function will also free memory associated with the delete edges
	void drop_all_edges() {
		for ( platform::Size inode = 1; inode <= nodes_.size(); inode++ ) {
			drop_all_edges_for_node( inode );
		}
		edges_.clear();
	}

	/// @brief memory accounting scheme
	virtual platform::Size count_static_memory() const {
		return sizeof( THIS );
	}

	/// @brief memory accounting scheme
	virtual platform::Size count_dynamic_memory() const {
		return sizeof( LMNode ) * nodes_.capacity() + sizeof( LMEdge ) * edges_.capacity()
			+ sizeof( uint32_t ) * deleted_edges_offsets_.capacity();
	}

	/// @brief returns a count of all the memory used by every vertex and edge in a graph
	///        by invoking the "polymorphic" count_static_memory and count_dynamic_memory of each
	///        (possibly derived) node and edge object as well as for the (possibly derived) graph
	///        class.
	/// @detail As long as the derived LMEdge and LMNode classes implement count_dynamic_memory() and
	///         count_static_memory(), they will be called even though they aren't technically
	///         overriden
	platform::Size getTotalMemoryUsage() {
		platform::Size total_memory = 0;
		for ( platform::Size ii = 1; ii <= num_nodes(); ii++ ) {
			// These call the "virtual" functions becuase we aren't working with the base class
			total_memory += nodes_[ ii ].count_dynamic_memory();
			total_memory += nodes_[ ii ].count_static_memory();
		}
		for ( LowMemEdgeListConstIter iter = const_edge_list_begin(); iter != const_edge_list_end(); ++iter ) {
			// Performing a static cast to the templated type so we can call the "virtual" functions.
			//  They aren't actually virtual, so we can't call them through the base class!!!
			LMEdge const * edge = static_cast< LMEdge const * >( *iter );
			total_memory += edge->count_dynamic_memory();
			total_memory += edge->count_static_memory();
		}

		total_memory += count_dynamic_memory();
		total_memory += count_static_memory();

		return total_memory;
	}


private:
	// DO NOT CHANGE THESE TO std::vector!!!
	// Edges are deleted by setting their nodes = 0
	// 0 must be reserved!!!
	utility::vector1< LMNode > nodes_;
	utility::vector1< LMEdge > edges_;

	// Last in first out stack of edges for reuse
	utility::vector1< uint32_t > deleted_edges_offsets_;

};



}
}


#endif
