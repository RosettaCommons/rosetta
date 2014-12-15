// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/graph/Graph.hh
/// @brief  generic graph class header
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_graph_Graph_hh
#define INCLUDED_core_graph_Graph_hh

// Rosetta Headers
#include <core/graph/Graph.fwd.hh>

#include <platform/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

// STL Headers
#include <iosfwd>
// AUTO-REMOVED #include <list>
#include <cassert>

// Boost Headers
#include <core/graph/unordered_object_pool.fwd.hpp>

#ifdef PYROSETTA
#include <core/graph/unordered_object_pool.hpp>
#endif

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1_bool.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>


namespace core {
namespace graph {

	///@brief An extensible graph class

	/**
		@li Nodes are identified by positive integers -- each node is assigned a unique integer starting at 1.

		@li Nodes know about both their upper and their lower edges.  Iterators can range over all
			edges, just the upper edges or just the lower edges.

		@li Nodes store their incident edge information in edge lists; they provide const and non-const
			iterators which return const and non-const pointers to Edge objects.

		@li The Graph object provides iterators for all of the edges in the graph; this list is unordered.

		@li This graph is instantiatable (ie not a pure virtual class).  It was also built with extension
			in mind.  See core/scoring/EnergyGraph for an example on how to create an extension of
			the Node, Edge and Graph classes

		@li Derived graph classes must define two factory methods: create_node and create_edge
			which are called by the base class when adding graph elements.

		@li Edges offer constant-time deletion -- edges remove themselves from the edge-lists
			that their nodes maintain, and from the edge-list that the Graph maintains without
			having to iterate over those lists.

		@li Memory is carefully managed to make edge addition and deletion exceedingly fast.
			The management is through a class "unodered_object_pool" which is very much
			like boost::object_pool.  The base class ueses these unordered object pools for
			base class edges (when used) and for the edge lists elements.  Derived classes
			may take advantage of this by defining their own unordered_object_pool objects.

		@li Graphs own the vertices and edges they create.  Graphs cannot share or give away edges.

		@li To delete an edge pointed to by Edge* e in a graph g, call the virtual function
			"g->delete_edge( &e );"

		@li Derived classes must invoke the "destroy_everything" method in their virtual destructors so that
			they are empty by the time the base class destructor is called.  The destroy_everything method
			invokes the virtual delete_edge() method -- remember, virtual methods do not work in base class
			constructors or destructors.  If a derived class does not call destroy_everything in its destructor,
			~Graph() will call Graph::destroy_everything which will call (surprisingly) Graph::delete_edge
			and not Derived::delete_edge.

		@li Derived classes should not invoke the base class copy constructors for the same reason.
			The base class copy constructor will call the base class create_node and create_edge methods
			and not the desired derived class versions.

		@li The virtual functions Derived graph classes must override:
				Node: copy_from, print, count_static_memory, count_dynamic_memory
				Edge: copy_from, count_static_mmory, count_dynamic_memory
				Graph: create_node, create_edge, delete_edge, count_static_memory, count_dynamic_memory

	**/


/// @brief Custom written edge list element class.  Little more than a struct.
/// Defined so that edge list memory management could rely on boost::pool like object
class EdgeListElement
{
public:
	EdgeListElement() : edge_( 0 ), previous_( 0 ), next_( 0 ) {}
	EdgeListElement( Edge * edge, EdgeListElement * previous, EdgeListElement * next )
		: edge_( edge ), previous_( previous ), next_( next ) {}


	~EdgeListElement() {}

	Edge * edge() { assert( edge_ ); return edge_; }
	void edge( Edge * setting ) { edge_ = setting; }
	Edge const * const_edge() const { assert( edge_ ); return edge_; }

	EdgeListElement * next() { return next_; }
	void next( EdgeListElement * setting ) { next_ = setting; }
	EdgeListElement const * const_next() const { return next_; }

	EdgeListElement * previous() { return previous_; }
	void previous( EdgeListElement * setting ) { previous_ = setting; }
	EdgeListElement const * const_previous() const { return previous_; }

	friend class EdgeList;

private:
	Edge * edge_;
	EdgeListElement * previous_;
	EdgeListElement * next_;

};

/// @brief Custom Edge list (non-const) iterator class, which can return non-const
/// Edge pointers.  This iterator cannot be used to change the structure of its list
/// without access to that list directly.  Customized since STL's const-iterator cannot
/// be prevented from giving non-const access to its data.  The former workaround to
/// this problem was to define two sets of edge lists on each vertex: a list< Edge * >
/// and a list< Edge const * >.
class EdgeListIterator
{
public:

	/// @brief default constructor, owner and element set to null
	EdgeListIterator()
		: owner_( 0 ), element_( 0 ) {}

	/// @brief owner constructor
	EdgeListIterator( EdgeList const * owner )
		: owner_( owner ), element_( 0 ) {}

	/// @brief owner and element constructor: points at a position in a list
	EdgeListIterator( EdgeList const * owner, EdgeListElement * element )
		: owner_( owner ), element_( element ) {}

	/// @brief copy constructor
	EdgeListIterator( EdgeListIterator const & src )
		: owner_( src.owner_ ), element_( src.element_ ) {}

	/// @brief non-virtual destructor, does nothing
	~EdgeListIterator() {}

	/// @brief assignmnet operator
	EdgeListIterator const & operator = ( EdgeListIterator const & rhs )
	{
		owner_ = rhs.owner_;
		element_ = rhs.element_;
		return *this;
	}

	/// @brief increment operator.  Point this iterator at the next element in the list.
	inline
	EdgeListIterator const & operator ++ ()
	{ assert( valid() ); element_ = element_->next(); return *this; }

	/// @brief decrement operator.  Point this iterator at the previous element in the list.
	inline
	EdgeListIterator const & operator -- ()
	{ assert( valid() ); element_ = element_->previous(); return *this; }

	/// @brief equality operator.  Do these elements point at the same list element?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	inline
	bool operator == ( EdgeListIterator const & rhs ) const {
		assert( owner_ == rhs.owner_ );
		return element_ == rhs.element_;
	}

	/// @brief inequality operator.  Do these elements point to different elements from the same list?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	inline
	bool operator != ( EdgeListIterator const & rhs ) const {
		return ! ( operator == ( rhs ) );
	}

	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a non-const edge pointer, which defines this class
	/// as a non-const iterator.  There are no methods, though, to access the list element itself.
	Edge * operator * () const { assert( valid() ); return element_->edge(); }

	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a non-const edge pointer, which defines this class
	/// as a non-const iterator.  There are no methods, though, to access the list element itself.
	Edge & operator -> () const { assert( valid() ); return * (element_->edge()) ; }

	/// @brief check that this iterator is valid.  Will not guarantee that the iterator points at an element
	/// that has not been destroyed.
	bool valid() const;

	friend class EdgeList;
	friend class EdgeListConstIterator;

private:
	EdgeList const * owner_;
	EdgeListElement * element_;

};

/// @brief Custom Edge list const iterator class, which returns only const
/// Edge pointers.  This iterator cannot be used to change the structure of its list
/// without access to that list directly.  Customized since STL's const-iterator cannot
/// be prevented from giving non-const access to its data.  The former workaround to
/// this problem was to define two sets of edge lists on each vertex: a list< Edge * >
/// and a list< Edge const * >.
class EdgeListConstIterator
{
public:

	/// @brief default constructor, owner and element set to null
	EdgeListConstIterator()
		: owner_( 0 ), element_( 0 ) {}

	/// @brief owner constructor
	EdgeListConstIterator( EdgeList const * owner )
		: owner_( owner ), element_( 0 ) {}

	/// @brief owner and element constructor: points at a position in a list
	EdgeListConstIterator( EdgeList const * owner, EdgeListElement const * element )
		: owner_( owner ), element_( element ) {}

	/// @brief copy constructor
	EdgeListConstIterator( EdgeListConstIterator const & src )
		: owner_( src.owner_ ), element_( src.element_ ) {}

	/// @brief const-cast constructor
	EdgeListConstIterator( EdgeListIterator const & src )
		: owner_( src.owner_ ), element_( src.element_ ) {}


	/// @brief non-virtual destructor, does nothing
	~EdgeListConstIterator() {}

	/// @brief assignmnet operator
	EdgeListConstIterator const & operator = ( EdgeListConstIterator const & rhs )
	{
		owner_ = rhs.owner_;
		element_ = rhs.element_;
		return *this;
	}

	/// @brief increment operator.  Point this iterator at the next element in the list.
	inline
	EdgeListConstIterator const & operator ++ ()
	{ assert( valid() ); element_ = element_->const_next(); return *this; }

	/// @brief decrement operator.  Point this iterator at the previous element in the list.
	inline
	EdgeListConstIterator const & operator -- ()
	{ assert( valid() ); element_ = element_->const_previous(); return *this; }

	/// @brief equality operator.  Do these elements point at the same list element?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	inline
	bool operator == ( EdgeListConstIterator const & rhs ) const {
		assert( owner_ == rhs.owner_ );
		return element_ == rhs.element_;
	}

	/// @brief const-cast equality operator.  Do these elements point at the same list element?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	inline
	bool operator == ( EdgeListIterator const & rhs ) const {
		assert( owner_ == rhs.owner_ );
		return element_ == rhs.element_;
	}


	/// @brief inequality operator.  Do these elements point to different elements from the same list?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	/// Elements need not be valid to be compared.
	inline
	bool operator != ( EdgeListConstIterator const & rhs ) const {
		return ! ( operator == ( rhs ) );
	}

	/// @brief const-cast inequality operator.  Do these elements point to different elements from the same list?  Asserts
	/// that they belong to the same list.  Never compare elements from separate lists.
	/// Elements need not be valid to be compared.
	inline
	bool operator != ( EdgeListIterator const & rhs ) const {
		return ! ( operator == ( rhs ) );
	}

	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a const edge pointer, which defines this class
	/// as a const iterator.  There are no methods, of course, to access the list element itself.
	Edge const * operator * () const { assert( valid() ); return element_->const_edge(); }


	/// @brief dereference operator: return the edge pointer pointed to by the list element this iterator
	/// is hovering over.  This method returns a const edge pointer, which defines this class
	/// as a const iterator.  There are no methods, of course, to access the list element itself.
	Edge const & operator -> () const { assert( valid() ); return * (element_->const_edge()); }

	/// @brief Is this a valid iterator?
	inline
	bool valid() const;

	friend class EdgeList;

private:
	EdgeList const * owner_;
	EdgeListElement const * element_;

};

/// @brief Custom edge list class.  Returns const-iterators which only return Edge const *'s
/// and non-const-iterators which can return either const or non-const Edge*'s.  Manages its
/// own memory using an unordered-object-pool for fast insertion and deletion of EdgeListElements.
/// Implemented as a doubly linked list, though there's no practical way to start at the
/// end of a list and work backward since decrementing the end iterator is not a valid operation.
class EdgeList
{
public:
	EdgeList( boost::unordered_object_pool< EdgeListElement > & edge_list_element_pool );
	~EdgeList();

	/// @brief create a new edge-list element and insert it at the front of the list
	void push_back( Edge * edgeptr );

	/// @brief create a new edge-list element and insert it at the end of the list
	void push_front( Edge * edgeptr );

	/// @brief insert a new edge-list element in the list infront of the input iterator and
	/// return an iterator to the position of the new element
	EdgeListIterator insert( EdgeListIterator const & element_to_insert_before, Edge * edgeptr );

	/// @brief returns a non-const iterator to the front of the list
	EdgeListIterator begin() { return EdgeListIterator( this, end_->next_ ); }

	/// @brief returns a const iterator to the front of the list
	EdgeListConstIterator begin() const { return EdgeListConstIterator( this, end_->next_ ); }

	/// @brief returns a const iterator to the front of the list
	EdgeListConstIterator const_begin() const { return EdgeListConstIterator( this, end_->next_ ); }

	/// @brief returns a non-const iterator to the last element in the list
	EdgeListIterator last() { return EdgeListIterator( this, end_->previous_ ); }

	/// @brief returns a const iterator to the last element in the list
	EdgeListConstIterator last() const { return EdgeListConstIterator( this, end_->previous_ ); }

	/// @brief returns a const iterator to the last element in the list
	EdgeListConstIterator const_last() const { return EdgeListConstIterator( this, end_->previous_ ); }

	/// @brief returns a non-const iterator to the end of the list
	EdgeListIterator end() { return EdgeListIterator( this, end_ ); }

	/// @brief returns a const iterator to the end of the list
	EdgeListConstIterator end() const { return EdgeListConstIterator( this, end_ ); }

	/// @brief returns a const iterator to the end of the list
	EdgeListConstIterator const_end() const { return EdgeListConstIterator( this, end_ ); }

	/// @brief removes an element from the list pointed to by the input iterator
	void erase( EdgeListIterator to_erase );

	/// @brief method invoked by the EdgeListIterator class: is an iterator the special
	/// end iterator for a list?
	bool is_end_element( EdgeListElement const * element ) const { return element == end_; }

	/// O(N)
	platform::Size size() const;

private:

	///@brief Uncopyable -- private and unimplemented copy ctor
	EdgeList( EdgeList const & );
	///@brief Uncopyable -- private and unimplemented assignment operator
	EdgeList const & operator = ( EdgeList const & );

	/// @brief this edge-list-element-pool reference is handed to the list
	/// to draw from.  This pool must outlive the edge-list itself.
	/// To guarantee this, for the case of class Graph, the graph deletes its
	/// nodes (and their edge lists) before it deletes itself.
	boost::unordered_object_pool< EdgeListElement > & edge_list_element_pool_;

	/// @brief The special "end" position in the list.
	EdgeListElement * end_;
};


class Node
{
public:

	typedef EdgeListIterator      EdgeListIter;
	typedef EdgeListConstIterator EdgeListConstIter;

public:
	virtual ~Node();
	Node( Graph*, platform::Size node_id );

	/// @brief invoked during graph assignment operators to copy any
	/// node data from one graph to another graph.  The source node must
	/// be the same type as this node.
	virtual void copy_from( Node const * source );

	/// @brief adds edge pointer to edge list; returns an iterator to the new
	/// list element
	void add_edge( Edge* edge_ptr, EdgeListIter & );

	/// @brief removes an edge iterator from the node's edge list.  Only called by Edge class.
	void drop_edge( EdgeListIter edge_iterator );

	/// @brief deletes all edges incident upon this node
	void drop_all_edges();

	/// @details manually change the number of neighbors for a Node. Used
	/// for symmetry scoring
	void set_num_neighbors_counting_self_static( platform::Size neighbor );

	/// @brief a "slow" (linear) search for an edge.
	Edge const * find_edge(platform::Size other_node_index) const;
	Edge * find_edge(platform::Size other_node_index);

	/// @brief send summaray data about this node to the screen
	virtual void print() const;

	/// @brief returns a non-const iterator to the beginning of its edge list
	EdgeListIter edge_list_begin() { return incident_edge_list_.begin(); }
	/// @brief returns a const iterator to the beginning of its edge list
	EdgeListConstIter const_edge_list_begin() const { return incident_edge_list_.const_begin(); }

	/// @brief returns a non-const iterator to the end of its edge list
	EdgeListIter edge_list_end() { return incident_edge_list_.end(); }
	/// @brief returns a const iterator to the end of its edge list
	EdgeListConstIter const_edge_list_end() const { return incident_edge_list_.const_end(); }

	/// @brief returns a non-const iterator to the beginning of its lower-edge list
	EdgeListIter lower_edge_list_begin() { return incident_edge_list_.begin(); }
	/// @brief returns a const iterator to the beginning of its lower-edge list
	EdgeListConstIter const_lower_edge_list_begin() const { return incident_edge_list_.const_begin(); }

	/// @brief returns a non-const iterator to the end of its lower-edge list
	EdgeListIter lower_edge_list_end() { return first_upper_edge_; }
	/// @brief returns a const iterator to the end of its lower-edge list
	EdgeListConstIter const_lower_edge_list_end() const { return first_upper_edge_; }

	/// @brief returns a non-const iterator to the beginning of its upper-edge list
	EdgeListIter upper_edge_list_begin() { return first_upper_edge_; }
	/// @brief returns a const iterator to the beginning of its upper-edge list
	EdgeListConstIter const_upper_edge_list_begin() const { return first_upper_edge_; }

	/// @brief returns a non-const iterator to the end of its upper-edge list
	EdgeListIter upper_edge_list_end() { return incident_edge_list_.end(); }
	/// @brief returns a const iterator to the end of its upper-edge list
	EdgeListConstIter const_upper_edge_list_end() const { return incident_edge_list_.const_end(); }

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
		if ( loop_incident_ ) { return num_incident_edges_; }
		else { return num_incident_edges_ + 1; }
	}

	/// @brief the number of neighbors counting "self" as neighbor. Defaults to
  /// num_neighbors_counting_self() but can be set to other values as well.
  /// Useful in calculation of symmetrical structures.
  inline
  platform::Size num_neighbors_counting_self_static() const
  {
    return num_neighbors_counting_self_static_;
  }

	/// @brief the number of lower neighbors
	inline
	platform::Size get_num_edges_to_smaller_indexed_nodes() const
	{
		return num_edges_to_smaller_indexed_nodes_;
	}

	/// @brief the number of upper neighbors -- which "self" neighborness is counted if a loop edge
	/// is present
	inline
	platform::Size get_num_edges_to_larger_indexed_nodes() const
	{
		return num_edges_to_larger_indexed_nodes_;
	}

	/// @brief memory accounting scheme
	virtual platform::Size count_static_memory() const;
	/// @brief memory accounting scheme
	virtual platform::Size count_dynamic_memory() const;

	/// NOTE TO SELF: remove loop support
	bool loop_incident() const { return loop_incident_; }
protected:

	/// @brief derived class access to the owner
	inline
	Graph* get_owner() const
	{
		return owner_;
	}

private:

	platform::Size node_index_;
	platform::Size num_incident_edges_;
	platform::Size num_neighbors_counting_self_static_;
	bool loop_incident_; /// NOTE TO SELF: remove loop support
	platform::Size num_edges_to_smaller_indexed_nodes_;
	platform::Size num_edges_to_larger_indexed_nodes_;

	// edge list
	EdgeList incident_edge_list_;

	EdgeListIter first_upper_edge_;
	Graph* owner_;

	//no default constructor, uncopyable
	Node();
	Node( Node const & );
	Node & operator = ( Node & );
};

class Edge
{
public:
	typedef EdgeListIterator      EdgeListIter;

public:
	virtual ~Edge();

	/// @brief Main edge constructor.  This should only be invoked by create_new_edge, which
	/// itself is only called by add_edge.  The ONLY way an edge should be added to a graph
	/// is through add_edge.  NOTE: edges should be only be deleted by a call to the Graph's
	/// delete_edge method, and this method absolutely must be implemented by derived Graph
	/// classes.
	Edge(Graph* owner, platform::Size first_node_ind, platform::Size second_node_ind);

	/// @brief copy-from for use in Graph::operator= and copy ctors
	virtual void copy_from( Edge const * source );

	/// @brief returns the index of the one node given the index of the other.
	/// node_index must be one of the two nodes that this edge is incident upon.
	platform::Size get_other_ind(platform::Size node_index) const;

	/// @brief returns a const pointer to one node given the index of the other.
	/// node_index must be one of the two nodes that this edge is incident upon.
	Node const * get_other_node(platform::Size node_index) const;
	/// @brief returns a non-const pointer to one node given the index of the other.
	/// node_index must be one of the two nodes that this edge is incident upon.
	Node * get_other_node(platform::Size node_index);

	/// @brief returns the index of the lower node
	inline
	platform::Size get_first_node_ind() const { return node_indices_[0]; }
	/// @brief returns the index of the upper node
	inline
	platform::Size get_second_node_ind() const { return node_indices_[1]; }

	/// @brief called only by class Graph, this function gives the Edge the data it needs
	/// to later delete itself from its owner's edge list in constant time.
	void set_pos_in_owners_list( EdgeListIter edge_iterator );

	/// @brief Is this the same edge as another edge (node1,node2)?  Note:
	/// this graph does not work for multi-graphs.  Edges must be unique.
	bool same_edge(platform::Size node1, platform::Size node2) const;

	/// @brief Is this edge a loop? In Pseudographs, loop edges are incident twice on a single vertex.
	bool is_loop() const { return node_indices_[ 0 ] == node_indices_[ 1 ];}

	/// @brief how much memory is statically allocated by this edge
	virtual platform::Size count_static_memory() const;
	/// @brief how much memory is dynamically allocated by this edge -- must be recursively invoked
	/// by a derived class.
	virtual platform::Size count_dynamic_memory() const;

protected:

	//Read access to private data granted to derived classes

	/// @brief get the node index for one of the two nodes this edge is incident upon
	/// uses c-style index-from-0.
	inline
	platform::Size get_node_index( platform::Size index ) const
	{
		assert( index == 0 || index == 1 );
		return node_indices_[ index ];
	}

	/// @brief get a const * to one node that this edge is incident upon
	/// uses c-style index-from-0 for these two nodes
	inline
	Node const *
	get_node( platform::Size index ) const
	{
		assert( index == 0 || index == 1 );
		return nodes_[ index ];
	}

	/// @brief get a non-const * to one node that this edge is incident upon
	/// uses c-style index-from-0 for these two nodes
	inline
	Node *
	get_node( platform::Size index )
	{
		assert( index == 0 || index == 1 );
		return nodes_[ index ];
	}


	/// @brief get a const * to the owning graph
	inline
	Graph const *
	get_owner() const
	{
		return owner_;
	}

	/// @brief get a non-const * to the owning graph
	inline
	Graph *
	get_owner()
	{
		return owner_;
	}

private:
	platform::Size node_indices_[2];
	Node* nodes_[2];

	EdgeListIter pos_in_nodes_edge_list_[2];

	EdgeListIter pos_in_owners_edge_list_;
	Graph* owner_;

	//no default constructor, uncopyable
	Edge();
	Edge( Edge const & );
	Edge & operator = ( Edge & );

};

/// @brief A Graph with constant time edge insertion and deletion.  Extensible.
class Graph : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< Graph >
{
public:
	/// self pointers
	inline GraphCOP get_self_ptr() const { return shared_from_this(); }
	inline GraphOP get_self_ptr() { return shared_from_this(); }

public:
	typedef utility::vector1< Node* > NodeVector;

	typedef Node::EdgeListIter EdgeListIter;
	typedef Node::EdgeListConstIter EdgeListConstIter;

	typedef utility::pointer::ReferenceCount parent;

public:

	/// @brief virtual destructor.  Derived classes must ensure they've destroyed all their
	/// nodes and edges through a call to "destroy_everything" before this function is arrived at
	virtual ~Graph();

	/// @brief ctor
	Graph();

	/// @brief num nodes ctor
	Graph(platform::Size num_nodes);

	/// @brief copy ctor.  Must not be called by derived class copy ctors.
	Graph( Graph const & source );

	/// @brief assignment operator.  source and this must have the same type.
	Graph & operator = ( Graph const & source );

	/// @brief copy the edge connectivity from a source graph with a potentially
	/// unknown type.
	void copy_connectivity( Graph const & source );

	/// @brief the number of nodes in the graph
	inline
	platform::Size num_nodes() const
	{
		return num_nodes_;
	}

	/// @brief set the number of nodes in the graph -- deletes any existing edges in the graph
	void set_num_nodes( platform::Size num_nodes );

	/// @brief add an edge between two vertices.  Invokes "create_edge" from the derived class.
	/// Returns a pointer to the edge after its been added, allowing the calling function
	/// to immediately set data for this edge.
	Edge * add_edge(platform::Size node1, platform::Size node2);
	/// @brief add an edge to this graph copying the data from an edge in another graph.
	/// Returns a pointer to the edge after its been added, allowing the calling function
	/// to immediately set data for this edge.
	Edge * add_edge( Edge const * example_edge );

	/// @brief is an edge already present in the graph? O(V) worst case.  O(1) iff all vertices have O(1) edges
	bool get_edge_exists(platform::Size node1, platform::Size node2) const;
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

	/// @brief O(V^3).  Computes all pairs shortest paths using Warshall's algorithm
	/// and writes all the path distances to the two-dimensional table.
	ObjexxFCL::FArray2D_int all_pairs_shortest_paths() const;

public:
	inline
	Node const * get_node( platform::Size index ) const
	{
		assert( index > 0 && index <= num_nodes_ );
		return nodes_[ index ];
	}

	inline
	Node* get_node( platform::Size index )
	{
		assert( index > 0 && index <= num_nodes_ );
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
	EdgeListConstIter const_edge_list_begin() const
	{
		return edge_list_.const_begin();
	}

	/// @brief returns a non-const iterator to the beginning of the (unordered) edge list for the graph.
	/// this edge list contains all the edges in the graph, not simply those for a particular vertex
	inline
	EdgeListIter edge_list_begin()
	{
		return edge_list_.begin();
	}

	/// @brief returns a const iterator to the end of the (unordered) edge list for the graph.
	/// this edge list contains all the edges in the graph, not simply those for a particular vertex
	inline
	EdgeListConstIter const_edge_list_end() const
	{
		return edge_list_.const_end();
	}

	/// @brief returns a non-const iterator to the end of the (unordered) edge list for the graph.
	/// this edge list contains all the edges in the graph, not simply those for a particular vertex
	inline
	EdgeListIter edge_list_end()
	{
		return edge_list_.end();
	}

	/// @brief returns a pointer to the edge connecting nodes node1 and node2, if that edge exists
	/// in the graph, o.w. returns 0.  Focuses the graph on this edge for fast subsequent retrieval.
	Edge * find_edge( platform::Size node1, platform::Size node2 );
	/// @brief returns a const pointer to the edge connecting nodes node1 and node2, if that edge exists
	/// in the graph, o.w. returns 0.  Focuses the graph on this edge for fast subsequent retrieval.
	Edge const * find_edge(platform::Size node1, platform::Size node2) const;

	/// @brief returns a pointer to the focused edge
	Edge * focused_edge() { return focused_edge_;}
	/// @brief returns a const-pointer to the focused edge
	Edge const * focused_edge() const { return focused_edge_;}

	/// @brief remove an edge from the graph. (NEW AS OF 12/9/07) Never call C++'s
	/// "delete" function on an edge pointer directly.  Derived classes must implement this function.
	/// If they wish to use unordered_object_pools to manage their memory
	virtual void delete_edge( Edge * edge );

	/// @brief returns a count of all the memory used by every vertex and edge in a graph
	/// by invoking the polymorphic count_static_memory and count_dynamic_memory of each
	/// (possibly derived) node and edge object as well as for the (possibly derived) graph
	/// class.
	platform::Size getTotalMemoryUsage() const;

	friend class Node;
	friend class Edge;

protected:

	virtual platform::Size count_static_memory() const;
	virtual platform::Size count_dynamic_memory() const;

	/// @brief remove an edge from the entire-graph edge list. Called only by class Edge
	/// during its destructor
	void drop_edge( EdgeListIter edge_iter );

	/// @brief deallocate all nodes and edges from the graph
	void delete_everything();

	/// @brief factory method for node creation, defined by derived graph
	/// classes, called by the base class
	virtual Node* create_new_node( platform::Size node_index );

	/// @brief factory method for edge creation, defined by derived graph
	/// classes, called by the base class
	virtual Edge* create_new_edge( platform::Size index1, platform::Size index2 );

	/// @brief factory method for edge copy-construction.  Derived class
	/// should downcast the example_edge pointer and may read that edge's data.
	virtual Edge* create_new_edge( Edge const * example_edge );

	/// @brief Used by class Node only, this is the pool from which edge lists are to
	/// allocate their edge lists elements from.
	boost::unordered_object_pool< EdgeListElement > &
	edge_list_element_pool() {
		return * edge_list_element_pool_;
	}

private:
	platform::Size num_nodes_;
	NodeVector nodes_;
	platform::Size num_edges_;

	/// @brief the pool from which edge lists are to allocate their edge list elements
	boost::unordered_object_pool< EdgeListElement > * edge_list_element_pool_;
	EdgeList edge_list_;

	/// @brief the pool from which class Graph allocates Edge objects.
	/// Not used by derived classes
	boost::unordered_object_pool< Edge > * edge_pool_;

	/// @brief Quick-access to a frequently needed edge -- the most recently sought edge
	/// in a call to find_edge() or the most recently added edge
	mutable Edge* focused_edge_;

};

inline
bool EdgeListIterator::valid() const
{ return ( owner_ != 0 && element_ != 0 && ! owner_->is_end_element( element_ ) );}

inline
bool EdgeListConstIterator::valid() const
{ return ( owner_ != 0 && element_ != 0 && ! owner_->is_end_element( element_ ) ); }


} //end namespace graph
} //end namespace core

#endif
