// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/graph/Digraph.hh
/// @brief  generic directed graph class header
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_graph_Digraph_hh
#define INCLUDED_utility_graph_Digraph_hh

// Unit Headers
#include <utility/graph/Digraph.fwd.hh>

// Package Headers
#include <utility/graph/unordered_object_pool.fwd.hpp>
#include <platform/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// STL Headers
#include <iosfwd>
#include <list>

#ifdef PYROSETTA
#include <utility/graph/unordered_object_pool.hpp>
#endif

#ifdef    SERIALIZATION
// Cereal Headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace utility {
namespace graph {

/// @brief An extensible directed graph class

/**
@li DirectedNodes are identified by positive integers -- each node is assigned a unique integer starting at 1.

@li DirectedNodes know about both their incoming and outgoing edges.  Iterators can range over all
edges, just the incoming edges or just the outgoing edges.

@li DirectedNodes store their incident edge information in edge lists; they provide const and non-const
iterators which return const and non-const pointers to Edge objects.

@li The Digraph object provides iterators for all of the edges in the graph; this list is unordered.

@li The Digraph base class is instantiatable (ie not a pure virtual class).  It was also built with extension
in mind.  See protocols/jd3/JobDigraph for an example on how to create an extension of
the DirectedNode, Edge and Digraph classes

@li Derived graph classes must define two factory methods: create_node and create_edge
which are called by the base class when adding graph elements.

@li Edges offer constant-time deletion -- edges remove themselves from the edge-lists
that their nodes maintain, and from the edge-list that the Digraph maintains without
having to iterate over those lists.

@li Digraphs own the vertices and edges they create.  Digraphs cannot share or give away edges.

@li The virtual functions Derived graph classes must override:
DirectedNode: copy_from, print, count_static_memory, count_dynamic_memory
Edge: copy_from, count_static_mmory, count_dynamic_memory
Digraph: create_node, create_edge, delete_edge, count_static_memory, count_dynamic_memory

**/


/// @brief Custom written edge list element class.  Little more than a struct.
/// Defined so that edge list memory management could rely on boost::pool like object
class DirectedEdgeListElement
{
public:
	DirectedEdgeListElement() : edge_( 0 ), previous_( 0 ), next_( 0 ) {}
	DirectedEdgeListElement( DirectedEdge * edge, DirectedEdgeListElement * previous, DirectedEdgeListElement * next )
	: edge_( edge ), previous_( previous ), next_( next ) {}


	~DirectedEdgeListElement() {}

	DirectedEdge * edge() {debug_assert( edge_ ); return edge_; }
	void edge( DirectedEdge * setting ) { edge_ = setting; }
	DirectedEdge const * const_edge() const {debug_assert( edge_ ); return edge_; }

	DirectedEdgeListElement * next() { return next_; }
	void next( DirectedEdgeListElement * setting ) { next_ = setting; }
	DirectedEdgeListElement const * const_next() const { return next_; }

	DirectedEdgeListElement * previous() { return previous_; }
	void previous( DirectedEdgeListElement * setting ) { previous_ = setting; }
	DirectedEdgeListElement const * const_previous() const { return previous_; }

	friend class DirectedEdgeList;

private:
	DirectedEdge * edge_;
	DirectedEdgeListElement * previous_;
	DirectedEdgeListElement * next_;

};

/// @brief Custom DirectedEdge list (non-const) iterator class, which can return non-const
/// DirectedEdge pointers.  This iterator cannot be used to change the structure of its list
/// without access to that list directly.  Customized since STL's const-iterator cannot
/// be prevented from giving non-const access to its data.  The former workaround to
/// this problem was to define two sets of edge lists on each vertex: a list< DirectedEdge * >
/// and a list< DirectedEdge const * >.
class DirectedEdgeListIterator
{
public:

	/// @brief default constructor, owner and element set to null
	DirectedEdgeListIterator()
	: owner_( 0 ), element_( 0 ) {}

	/// @brief owner constructor
	DirectedEdgeListIterator( DirectedEdgeList const * owner )
	: owner_( owner ), element_( 0 ) {}

	/// @brief owner and element constructor: points at a position in a list
	DirectedEdgeListIterator( DirectedEdgeList const * owner, DirectedEdgeListElement * element )
	: owner_( owner ), element_( element ) {}

	/// @brief copy constructor
	DirectedEdgeListIterator( DirectedEdgeListIterator const & src )
	: owner_( src.owner_ ), element_( src.element_ ) {}

	/// @brief non-virtual destructor, does nothing
	~DirectedEdgeListIterator() {}

	/// @brief assignmnet operator
	DirectedEdgeListIterator & operator = ( DirectedEdgeListIterator const & rhs )
	{
		owner_ = rhs.owner_;
		element_ = rhs.element_;
		return *this;
	}

	/// @brief increment operator.  Point this iterator at the next element in the list.
	inline
	DirectedEdgeListIterator const & operator ++ ()
	{debug_assert( valid() ); element_ = element_->next(); return *this; }

	/// @brief decrement operator.  Point this iterator at the previous element in the list.
	inline
	DirectedEdgeListIterator const & operator -- ()
	{debug_assert( valid() ); element_ = element_->previous(); return *this; }

	/// @brief equality operator.  Do these elements point at the same list element?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	inline
	bool operator == ( DirectedEdgeListIterator const & rhs ) const {
		debug_assert( owner_ == rhs.owner_ );
		return element_ == rhs.element_;
	}

	/// @brief inequality operator.  Do these elements point to different elements from the same list?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	inline
	bool operator != ( DirectedEdgeListIterator const & rhs ) const {
		return ! ( operator == ( rhs ) );
	}

	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a non-const edge pointer, which defines this class
	/// as a non-const iterator.  There are no methods, though, to access the list element itself.
	DirectedEdge * operator * () const {debug_assert( valid() ); return element_->edge(); }

	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a non-const edge pointer, which defines this class
	/// as a non-const iterator.  There are no methods, though, to access the list element itself.
	DirectedEdge & operator -> () const {debug_assert( valid() ); return * (element_->edge()) ; }

	/// @brief check that this iterator is valid.  Will not guarantee that the iterator points at an element
	/// that has not been destroyed.
	bool valid() const;

	friend class DirectedEdgeList;
	friend class DirectedEdgeListConstIterator;

private:
	DirectedEdgeList const * owner_;
	DirectedEdgeListElement * element_;

};

/// @brief Custom DirectedEdge list const iterator class, which returns only const
/// DirectedEdge pointers.  This iterator cannot be used to change the structure of its list
/// without access to that list directly.  Customized since STL's const-iterator cannot
/// be prevented from giving non-const access to its data.
class DirectedEdgeListConstIterator
{
public:

	/// @brief default constructor, owner and element set to null
	DirectedEdgeListConstIterator()
	: owner_( 0 ), element_( 0 ) {}

	/// @brief owner constructor
	DirectedEdgeListConstIterator( DirectedEdgeList const * owner )
	: owner_( owner ), element_( 0 ) {}

	/// @brief owner and element constructor: points at a position in a list
	DirectedEdgeListConstIterator( DirectedEdgeList const * owner, DirectedEdgeListElement const * element )
	: owner_( owner ), element_( element ) {}

	/// @brief copy constructor
	DirectedEdgeListConstIterator( DirectedEdgeListConstIterator const & src )
	: owner_( src.owner_ ), element_( src.element_ ) {}

	/// @brief const-cast constructor
	DirectedEdgeListConstIterator( DirectedEdgeListIterator const & src )
	: owner_( src.owner_ ), element_( src.element_ ) {}


	/// @brief non-virtual destructor, does nothing
	~DirectedEdgeListConstIterator() {}

	/// @brief assignmnet operator
	DirectedEdgeListConstIterator & operator = ( DirectedEdgeListConstIterator const & rhs )
	{
		owner_ = rhs.owner_;
		element_ = rhs.element_;
		return *this;
	}

	/// @brief increment operator.  Point this iterator at the next element in the list.
	inline
	DirectedEdgeListConstIterator const & operator ++ ()
	{debug_assert( valid() ); element_ = element_->const_next(); return *this; }

	/// @brief decrement operator.  Point this iterator at the previous element in the list.
	inline
	DirectedEdgeListConstIterator const & operator -- ()
	{debug_assert( valid() ); element_ = element_->const_previous(); return *this; }

	/// @brief equality operator.  Do these elements point at the same list element?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	inline
	bool operator == ( DirectedEdgeListConstIterator const & rhs ) const {
		debug_assert( owner_ == rhs.owner_ );
		return element_ == rhs.element_;
	}

	/// @brief const-cast equality operator.  Do these elements point at the same list element?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	inline
	bool operator == ( DirectedEdgeListIterator const & rhs ) const {
		debug_assert( owner_ == rhs.owner_ );
		return element_ == rhs.element_;
	}


	/// @brief inequality operator.  Do these elements point to different elements from the same list?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	/// Elements need not be valid to be compared.
	inline
	bool operator != ( DirectedEdgeListConstIterator const & rhs ) const {
		return ! ( operator == ( rhs ) );
	}

	/// @brief const-cast inequality operator.  Do these elements point to different elements from the same list?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	/// Elements need not be valid to be compared.
	inline
	bool operator != ( DirectedEdgeListIterator const & rhs ) const {
		return ! ( operator == ( rhs ) );
	}

	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a const edge pointer, which defines this class
	/// as a const iterator.  There are no methods, of course, to access the list element itself.
	DirectedEdge const * operator * () const {debug_assert( valid() ); return element_->const_edge(); }


	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a const edge pointer, which defines this class
	/// as a const iterator.  There are no methods, of course, to access the list element itself.
	DirectedEdge const & operator -> () const {debug_assert( valid() ); return * (element_->const_edge()); }

	/// @brief Is this a valid iterator?
	inline
	bool valid() const;

	friend class DirectedEdgeList;

private:
	DirectedEdgeList const * owner_;
	DirectedEdgeListElement const * element_;

};

/// @brief Custom edge list class.  Returns const-iterators which only return DirectedEdge const *'s
/// and non-const-iterators which can return either const or non-const DirectedEdge*'s.  Manages its
/// own memory using an unordered-object-pool for fast insertion and deletion of DirectedEdgeListElements.
/// Implemented as a doubly linked list, though there's no practical way to start at the
/// end of a list and work backward since decrementing the end iterator is not a valid operation.
class DirectedEdgeList
{
public:
	DirectedEdgeList( boost::unordered_object_pool< DirectedEdgeListElement > & edge_list_element_pool );
	~DirectedEdgeList();

	/// @brief create a new edge-list element and insert it at the front of the list
	void push_back( DirectedEdge * edgeptr );

	/// @brief create a new edge-list element and insert it at the end of the list
	void push_front( DirectedEdge * edgeptr );

	/// @brief insert a new edge-list element in the list infront of the input iterator and
	/// return an iterator to the position of the new element
	DirectedEdgeListIterator insert( DirectedEdgeListIterator const & element_to_insert_before, DirectedEdge * edgeptr );

	/// @brief returns a non-const iterator to the front of the list
	DirectedEdgeListIterator begin() { return DirectedEdgeListIterator( this, end_->next_ ); }

	/// @brief returns a const iterator to the front of the list
	DirectedEdgeListConstIterator begin() const { return DirectedEdgeListConstIterator( this, end_->next_ ); }

	/// @brief returns a const iterator to the front of the list
	DirectedEdgeListConstIterator const_begin() const { return DirectedEdgeListConstIterator( this, end_->next_ ); }

	/// @brief returns a non-const iterator to the last element in the list
	DirectedEdgeListIterator last() { return DirectedEdgeListIterator( this, end_->previous_ ); }

	/// @brief returns a const iterator to the last element in the list
	DirectedEdgeListConstIterator last() const { return DirectedEdgeListConstIterator( this, end_->previous_ ); }

	/// @brief returns a const iterator to the last element in the list
	DirectedEdgeListConstIterator const_last() const { return DirectedEdgeListConstIterator( this, end_->previous_ ); }

	/// @brief returns a non-const iterator to the end of the list
	DirectedEdgeListIterator end() { return DirectedEdgeListIterator( this, end_ ); }

	/// @brief returns a const iterator to the end of the list
	DirectedEdgeListConstIterator end() const { return DirectedEdgeListConstIterator( this, end_ ); }

	/// @brief returns a const iterator to the end of the list
	DirectedEdgeListConstIterator const_end() const { return DirectedEdgeListConstIterator( this, end_ ); }

	/// @brief removes an element from the list pointed to by the input iterator
	void erase( DirectedEdgeListIterator to_erase );

	/// @brief method invoked by the DirectedEdgeListIterator class: is an iterator the special
	/// end iterator for a list?
	bool is_end_element( DirectedEdgeListElement const * element ) const { return element == end_; }

	/// O(N)
	platform::Size size() const;

private:

	/// @brief Uncopyable -- private and unimplemented copy ctor
	DirectedEdgeList( DirectedEdgeList const & );
	/// @brief Uncopyable -- private and unimplemented assignment operator
	DirectedEdgeList const & operator = ( DirectedEdgeList const & );

	/// @brief this edge-list-element-pool reference is handed to the list
	/// to draw from.  This pool must outlive the edge-list itself.
	/// To guarantee this, for the case of class Digraph, the graph deletes its
	/// nodes (and their edge lists) before it deletes itself.
	boost::unordered_object_pool< DirectedEdgeListElement > & edge_list_element_pool_;

	/// @brief The special "end" position in the list.
	DirectedEdgeListElement * end_;
};


class DirectedNode
{
public:

	typedef DirectedEdgeListIterator      DirectedEdgeListIter;
	typedef DirectedEdgeListConstIterator DirectedEdgeListConstIter;

public:
	virtual ~DirectedNode();
	DirectedNode( Digraph*, platform::Size node_id );

	/// @brief invoked during graph assignment operators to copy any
	/// node data from one graph to another graph.  The source node must
	/// be the same type as this node.
	virtual void copy_from( DirectedNode const * source );

	/// @brief adds an edge pointer to node's edge list as an incoming edge
	/// i.e. this node is the head of the arrow; returns an iterator to the new
	/// list element.  Virtual so derived classes can observe when incoming
	/// edges are added.
	virtual
	void add_incoming_edge( DirectedEdge* edge_ptr, DirectedEdgeListIter & );

	/// @brief adds an edge pointer to node's edge list as an outgoing edge
	/// i.e. this node is the tail of the arrow; returns an iterator to the new
	/// list element. Virtual so derived classes can observe when incoming
	/// edges are added.
	virtual
	void add_outgoing_edge( DirectedEdge* edge_ptr, DirectedEdgeListIter & );

	/// @brief removes an edge iterator from the node's edge list.  Only called by DirectedEdge class.
	void drop_edge( DirectedEdgeListIter edge_iterator );

	/// @brief deletes all edges incident upon this node
	void drop_all_edges();

	/// @brief a "slow" (linear) search for an edge.
	DirectedEdge const * find_edge_to( platform::Size tail_node ) const;

	/// @brief a "slow" (linear) search for an edge.
	DirectedEdge * find_edge_to( platform::Size tail_node );

	/// @brief a "slow" (linear) search for an edge.
	DirectedEdge const * find_edge_from( platform::Size head_node ) const;

	/// @brief a "slow" (linear) search for an edge.
	DirectedEdge * find_edge_from( platform::Size head_node );

	/// @brief send summaray data about this node to the screen
	virtual void print() const;

	/// @brief returns a non-const iterator to the beginning of its edge list
	DirectedEdgeListIter edge_list_begin() { return incident_edge_list_.begin(); }
	/// @brief returns a const iterator to the beginning of its edge list
	DirectedEdgeListConstIter const_edge_list_begin() const { return incident_edge_list_.const_begin(); }

	/// @brief returns a non-const iterator to the end of its edge list
	DirectedEdgeListIter edge_list_end() { return incident_edge_list_.end(); }
	/// @brief returns a const iterator to the end of its edge list
	DirectedEdgeListConstIter const_edge_list_end() const { return incident_edge_list_.const_end(); }

	/// @brief returns a non-const iterator to the beginning of its incoming-edge list
	DirectedEdgeListIter incoming_edge_list_begin() { return incident_edge_list_.begin(); }
	/// @brief returns a const iterator to the beginning of its incoming-edge list
	DirectedEdgeListConstIter const_incoming_edge_list_begin() const { return incident_edge_list_.const_begin(); }

	/// @brief returns a non-const iterator to the end of its lower-edge list
	DirectedEdgeListIter incoming_edge_list_end() { return first_outgoing_edge_; }
	/// @brief returns a const iterator to the end of its lower-edge list
	DirectedEdgeListConstIter const_incoming_edge_list_end() const { return first_outgoing_edge_; }

	/// @brief returns a non-const iterator to the beginning of its outgoing-edge list
	DirectedEdgeListIter outgoing_edge_list_begin() { return first_outgoing_edge_; }
	/// @brief returns a const iterator to the beginning of its outgoing-edge list
	DirectedEdgeListConstIter const_outgoing_edge_list_begin() const { return first_outgoing_edge_; }

	/// @brief returns a non-const iterator to the end of its outgoing-edge list
	DirectedEdgeListIter outgoing_edge_list_end() { return incident_edge_list_.end(); }
	/// @brief returns a const iterator to the end of its outgoing-edge list
	DirectedEdgeListConstIter const_outgoing_edge_list_end() const { return incident_edge_list_.const_end(); }

	/// @brief the index for this node
	inline
	platform::Size get_node_index() const
	{
		return node_index_;
	}

	/// @brief the number of edges incident on this node, which may include a loop edge
	inline
	platform::Size num_edges() const
	{
		return num_incident_edges_;
	}

	/// @brief the number of neighbors counting "self" as a neighbor.
	inline
	platform::Size num_neighbors_counting_self() const
	{
		return num_incident_edges_ + 1;
	}

	/// @brief the number of incoming edges
	inline
	platform::Size indegree() const
	{
		return indegree_;
	}

	/// @brief the number of outgoing edges
	inline
	platform::Size outdegree() const
	{
		return outdegree_;
	}

	/// @brief memory accounting scheme
	virtual platform::Size count_static_memory() const;
	/// @brief memory accounting scheme
	virtual platform::Size count_dynamic_memory() const;

protected:

	/// @brief derived class access to the owner
	inline
	Digraph* get_owner() const
	{
		return owner_;
	}

private:

	platform::Size node_index_;
	platform::Size num_incident_edges_;
	platform::Size indegree_;
	platform::Size outdegree_;

	// edge list
	DirectedEdgeList incident_edge_list_;

	DirectedEdgeListIter first_outgoing_edge_;
	Digraph* owner_;

	//no default constructor, uncopyable
	DirectedNode();
	DirectedNode( DirectedNode const & );
	DirectedNode & operator = ( DirectedNode & );
};

class DirectedEdge
{
public:
	typedef DirectedEdgeListIterator      DirectedEdgeListIter;

public:
	virtual ~DirectedEdge();

	/// @brief Main edge constructor.  This should only be invoked by create_new_edge, which
	/// itself is only called by add_edge.  The ONLY way an edge should be added to a graph
	/// is through add_edge.  NOTE: edges should be only be deleted by a call to the Digraph's
	/// delete_edge method, and this method absolutely must be implemented by derived Digraph
	/// classes.
	DirectedEdge(Digraph* owner, platform::Size tail_node_ind, platform::Size head_node_ind);

	/// @brief copy-from for use in Digraph::operator= and copy ctors
	virtual void copy_from( DirectedEdge const * source );

	/// @brief returns the index of the one node given the index of the other.
	/// node_index must be one of the two nodes that this edge is incident upon.
	platform::Size get_other_ind(platform::Size node_index) const;

	/// @brief returns a const pointer to one node given the index of the other.
	/// node_index must be one of the two nodes that this edge is incident upon.
	DirectedNode const * get_other_node(platform::Size node_index) const;
	/// @brief returns a non-const pointer to one node given the index of the other.
	/// node_index must be one of the two nodes that this edge is incident upon.
	DirectedNode * get_other_node(platform::Size node_index);

	/// @brief returns the index of the tail node
	inline
	platform::Size get_tail_node_ind() const { return node_indices_[0]; }
	/// @brief returns the index of the upper node
	inline
	platform::Size get_head_node_ind() const { return node_indices_[1]; }

	/// @brief returns a pointer to the tail node
	inline
	DirectedNode const *
	get_tail_node() const { return get_node(0); }

	/// @brief returns a pointer to the tail node
	inline
	DirectedNode *
	get_tail_node() { return get_node(0); }

	/// @brief returns a pointer to the head node
	inline
	DirectedNode const *
	get_head_node() const { return get_node(1); }

	/// @brief returns a pointer to the head node
	inline
	DirectedNode *
	get_head_node() { return get_node(1); }

	inline
	bool is_tail_node( platform::Size node_index ) {
		debug_assert( node_indices_[ 0 ] == node_index || node_indices_[ 1 ] == node_index );
		return node_indices_[ 0 ] == node_index;
	}

	inline
	bool is_head_node( platform::Size node_index ) {
		debug_assert( node_indices_[ 0 ] == node_index || node_indices_[ 1 ] == node_index );
		return node_indices_[ 1 ] == node_index;
	}

	/// @brief called only by class Digraph, this function gives the DirectedEdge the data it needs
	/// to later delete itself from its owner's edge list in constant time.
	void set_pos_in_owners_list( DirectedEdgeListIter edge_iterator );

	/// @brief Is this the same edge as another edge (tail_node,head_node)?  Note:
	/// this graph does not work for multi-graphs.  DirectedEdges must be unique.
	bool same_edge(platform::Size tail_node, platform::Size head_node ) const;

	/// @brief how much memory is statically allocated by this edge
	virtual platform::Size count_static_memory() const;
	/// @brief how much memory is dynamically allocated by this edge -- must be recursively invoked
	/// by a derived class.
	virtual platform::Size count_dynamic_memory() const;

protected:

	//Read access to private data granted to derived classes

	/// @brief get the node index for one of the two nodes this edge is incident upon
	/// uses c-style index-from-0.  0 is the index of the tail node, 1 is the index
	/// of the head node.
	inline
	platform::Size get_node_index( platform::Size index ) const
	{
		debug_assert( index == 0 || index == 1 );
		return node_indices_[ index ];
	}

	/// @brief get a const * to one node that this edge is incident upon
	/// uses c-style index-from-0 for these two nodes. 0 is the index of the tail node,
	/// 1 is the index of the head node.
	inline
	DirectedNode const *
	get_node( platform::Size index ) const
	{
		debug_assert( index == 0 || index == 1 );
		return nodes_[ index ];
	}

	/// @brief get a non-const * to one node that this edge is incident upon
	/// uses c-style index-from-0 for these two nodes
	inline
	DirectedNode *
	get_node( platform::Size index )
	{
		debug_assert( index == 0 || index == 1 );
		return nodes_[ index ];
	}


	/// @brief get a const * to the owning graph
	inline
	Digraph const *
	get_owner() const
	{
		return owner_;
	}

	/// @brief get a non-const * to the owning graph
	inline
	Digraph *
	get_owner()
	{
		return owner_;
	}

private:
	platform::Size node_indices_[2];
	DirectedNode* nodes_[2];

	DirectedEdgeListIter pos_in_nodes_edge_list_[2];

	DirectedEdgeListIter pos_in_owners_edge_list_;
	Digraph* owner_;

	//no default constructor, uncopyable
	DirectedEdge();
	DirectedEdge( DirectedEdge const & );
	DirectedEdge & operator = ( DirectedEdge & );

};

/// @brief A Digraph with constant time edge insertion and deletion.  Extensible.
class Digraph : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< Digraph >
{
public:
	/// self pointers
	inline DigraphCOP get_self_ptr() const { return shared_from_this(); }
	inline DigraphOP get_self_ptr() { return shared_from_this(); }

public:
	typedef utility::vector1< DirectedNode* > DirectedNodeVector;

	typedef DirectedNode::DirectedEdgeListIter DirectedEdgeListIter;
	typedef DirectedNode::DirectedEdgeListConstIter DirectedEdgeListConstIter;

	typedef utility::pointer::ReferenceCount parent;

public:

	/// @brief virtual destructor.  Derived classes must ensure they've destroyed all their
	/// nodes and edges through a call to "destroy_everything" before this function is arrived at
	virtual ~Digraph();

	/// @brief ctor
	Digraph();

	/// @brief num nodes ctor
	Digraph(platform::Size num_nodes);

	/// @brief copy ctor.  Must not be called by derived class copy ctors.
	Digraph( Digraph const & source );

	/// @brief assignment operator.  source and this must have the same type.
	Digraph & operator = ( Digraph const & source );

	/// @brief copy the edge connectivity from a source graph with a potentially
	/// unknown type.
	void copy_connectivity( Digraph const & source );

	/// @brief the number of nodes in the graph
	inline
	platform::Size num_nodes() const
	{
		return num_nodes_;
	}

	/// @brief set the number of nodes in the graph -- deletes any existing edges in the graph
	void set_num_nodes( platform::Size num_nodes );

	/// @brief Add a new node to the graph; preserve all of the existing edges in the graph
	void add_node();

	/// @brief add a directed edge between two vertices.  Invokes "create_edge" from the derived class.
	/// Returns a pointer to the edge after its been added, allowing the calling function
	/// to immediately set data for this edge.
	DirectedEdge * add_edge(platform::Size tail_node, platform::Size head_node );
	/// @brief add an edge to this graph copying the data from an edge in another graph.
	/// Returns a pointer to the edge after its been added, allowing the calling function
	/// to immediately set data for this edge.
	DirectedEdge * add_edge( DirectedEdge const * example_edge );

	/// @brief is an edge already present in the graph? O(V) worst case.  O(1) iff all vertices have O(1) edges
	bool get_edge_exists(platform::Size tail_node, platform::Size head_node ) const;
	/// @brief delete all the edges present in the graph
	void drop_all_edges();
	/// @brief delete all the edges for a single vertex in the graph
	void drop_all_edges_for_node( platform::Size node );

	/// @brief send summary information to the screen for all vertices in the graph
	void print_vertices() const;

	/// @brief send an edge list to the stream os.
	void output_connectivity(std::ostream & os) const;
	/// @brief describe this graph in dimacs form to the stream os.
	void output_dimacs(std::ostream & os) const;


public:
	inline
	DirectedNode const * get_node( platform::Size index ) const
	{
		debug_assert( index > 0 && index <= num_nodes_ );
		return nodes_[ index ];
	}

	inline
	DirectedNode* get_node( platform::Size index )
	{
		debug_assert( index > 0 && index <= num_nodes_ );
		return nodes_[ index ];
	}


	inline
	platform::Size num_edges() const
	{
		return num_edges_;
	}

	/// @brief returns a const iterator to the beginning of the (unordered) edge list for the graph.
	/// this edge list contains all the edges in the graph, not simply those for a particular vertex
	inline
	DirectedEdgeListConstIter const_edge_list_begin() const
	{
		return edge_list_.const_begin();
	}

	/// @brief returns a non-const iterator to the beginning of the (unordered) edge list for the graph.
	/// this edge list contains all the edges in the graph, not simply those for a particular vertex
	inline
	DirectedEdgeListIter edge_list_begin()
	{
		return edge_list_.begin();
	}

	/// @brief returns a const iterator to the end of the (unordered) edge list for the graph.
	/// this edge list contains all the edges in the graph, not simply those for a particular vertex
	inline
	DirectedEdgeListConstIter const_edge_list_end() const
	{
		return edge_list_.const_end();
	}

	/// @brief returns a non-const iterator to the end of the (unordered) edge list for the graph.
	/// this edge list contains all the edges in the graph, not simply those for a particular vertex
	inline
	DirectedEdgeListIter edge_list_end()
	{
		return edge_list_.end();
	}

	/// @brief returns a pointer to the directed edge connecting nodes tail_node and head_node, if that edge exists
	/// in the graph, o.w. returns 0.  Focuses the graph on this edge for fast subsequent retrieval.
	DirectedEdge * find_edge( platform::Size tail_node, platform::Size head_node );

	/// @brief returns a const pointer to the directed edge connecting nodes tail_node and head_node, if that edge exists
	/// in the graph, o.w. returns 0.  Focuses the graph on this edge for fast subsequent retrieval.
	DirectedEdge const * find_edge( platform::Size tail_node, platform::Size head_node ) const;

	/// @brief returns a pointer to the focused edge
	DirectedEdge * focused_edge() { return focused_edge_;}
	/// @brief returns a const-pointer to the focused edge
	DirectedEdge const * focused_edge() const { return focused_edge_;}

	/// @brief remove an edge from the graph. (NEW AS OF 12/9/07) Never call C++'s
	/// "delete" function on an edge pointer directly.  Derived classes must implement this function.
	/// If they wish to use unordered_object_pools to manage their memory
	virtual void delete_edge( DirectedEdge * edge );

	/// @brief returns a count of all the memory used by every vertex and edge in a graph
	/// by invoking the polymorphic count_static_memory and count_dynamic_memory of each
	/// (possibly derived) node and edge object as well as for the (possibly derived) graph
	/// class.
	platform::Size getTotalMemoryUsage() const;

	friend class DirectedNode;
	friend class DirectedEdge;

#ifdef    SERIALIZATION
	template < class Archive >
	void save( Archive & archive ) const;

	template < class Archive >
	void load( Archive & archive );
#endif // SERIALIZATION

protected:

	virtual platform::Size count_static_memory() const;
	virtual platform::Size count_dynamic_memory() const;

	/// @brief remove an edge from the entire-graph edge list. Called only by class DirectedEdge
	/// during its destructor
	void drop_edge( DirectedEdgeListIter edge_iter );

	/// @brief deallocate all nodes and edges from the graph
	void delete_everything();

	/// @brief factory method for node creation, defined by derived graph
	/// classes, called by the base class
	virtual DirectedNode* create_new_node( platform::Size node_index );

	/// @brief factory method for edge creation, defined by derived graph
	/// classes, called by the base class
	virtual DirectedEdge* create_new_edge( platform::Size index1, platform::Size index2 );

	/// @brief factory method for edge copy-construction.  Derived class
	/// should downcast the example_edge pointer and may read that edge's data.
	virtual DirectedEdge* create_new_edge( DirectedEdge const * example_edge );

	/// @brief Used by class DirectedNode only, this is the pool from which edge lists are to
	/// allocate their edge lists elements from.
	boost::unordered_object_pool< DirectedEdgeListElement > &
	edge_list_element_pool() {
		return * edge_list_element_pool_;
	}

private:
	platform::Size num_nodes_;
	DirectedNodeVector nodes_;
	platform::Size num_edges_;

	/// @brief the pool from which edge lists are to allocate their edge list elements
	boost::unordered_object_pool< DirectedEdgeListElement > * edge_list_element_pool_;
	DirectedEdgeList edge_list_;

	/// @brief the pool from which class Digraph allocates DirectedEdge objects.
	/// Not used by derived classes
	boost::unordered_object_pool< DirectedEdge > * edge_pool_;

	/// @brief Quick-access to a frequently needed edge -- the most recently sought edge
	/// in a call to find_edge() or the most recently added edge
	mutable DirectedEdge* focused_edge_;

};

inline
bool DirectedEdgeListIterator::valid() const
{ return ( owner_ != 0 && element_ != 0 && ! owner_->is_end_element( element_ ) );}

inline
bool DirectedEdgeListConstIterator::valid() const
{ return ( owner_ != 0 && element_ != 0 && ! owner_->is_end_element( element_ ) ); }

/// @brief Construct a topological sort for the input directed graph, if it is a DAG,
/// and return whether or not the input graph is actually a DAG.
std::pair< std::list< platform::Size >, bool >
topological_sort( Digraph const & g );

/// @brief Return whether or not the input directed graph is a DAG -- this invokes
/// the topological_sort function and simply returns the second element in the
/// pair that it returns.
bool
digraph_is_a_DAG( Digraph const & g );


} //end namespace graph
} //end namespace utility

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( utility_graph_Digraph )
#endif // SERIALIZATION

#endif
