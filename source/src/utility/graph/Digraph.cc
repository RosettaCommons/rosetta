// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/graph/Digraph.cc
/// @brief  directed graph base classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <utility/graph/Digraph.hh>

// Package headers
#include <utility/graph/unordered_object_pool.hpp>

// Utility headers
#include <utility/assert.hh>
#include <utility/graph/unordered_object_pool.hpp>

//STL Headers
#include <iostream>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

// Boost Headers
#include <boost/pool/pool.hpp>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

using namespace ObjexxFCL;

namespace utility {
namespace graph {

DirectedEdgeList::DirectedEdgeList(
	boost::unordered_object_pool< DirectedEdgeListElement > & edge_list_element_pool
)
:
	edge_list_element_pool_( edge_list_element_pool ),
	end_( new DirectedEdgeListElement )
{
	end_->next_ = end_;
	end_->previous_ = end_;
}

DirectedEdgeList::~DirectedEdgeList()
{
	for ( DirectedEdgeListIterator
			eiter = begin(), eiter_end = end();
			eiter != eiter_end; /* no increment */ ) {
		DirectedEdgeListIterator eiter_next( eiter );
		++eiter_next;

		// remove the list element
		//delete eiter.element_;
		edge_list_element_pool_.destroy( eiter.element_ );
		eiter = eiter_next;
	}
	delete end_;
}

void
DirectedEdgeList::push_back( DirectedEdge * edgeptr )
{
	debug_assert( edgeptr ); // Do not accept null-pointing edges

	//DirectedEdgeListElement * new_node = new DirectedEdgeListElement( edgeptr, end_->previous_, end_ );
	DirectedEdgeListElement * new_node =  edge_list_element_pool_.construct( edgeptr, end_->previous_, end_ );

	end_->previous_->next_ = new_node; // may be end_->next_ if 0-element list
	end_->previous_ = new_node;
}

void
DirectedEdgeList::push_front( DirectedEdge * edgeptr )
{
	debug_assert( edgeptr ); // Do not accept null-pointing edges

	//DirectedEdgeListElement * new_node = new DirectedEdgeListElement( edgeptr, end_, end_->next_ );
	DirectedEdgeListElement * new_node =  edge_list_element_pool_.construct( edgeptr, end_, end_->next_ );

	end_->next_->previous_ = new_node; // may be end_->next_ if 0-element list
	end_->next_ = new_node;

}


DirectedEdgeListIterator
DirectedEdgeList::insert(
	DirectedEdgeListIterator const & element_to_insert_before,
	DirectedEdge * edgeptr
)
{
	debug_assert( edgeptr );
	debug_assert( element_to_insert_before.owner_ == this );

	DirectedEdgeListElement * next_node = element_to_insert_before.element_;
	DirectedEdgeListElement * prev_node = element_to_insert_before.element_->previous_;

	//DirectedEdgeListElement * new_node = new DirectedEdgeListElement( edgeptr, prev_node, next_node );
	DirectedEdgeListElement * new_node =  edge_list_element_pool_.construct( edgeptr, prev_node, next_node );

	next_node->previous_ = new_node;
	prev_node->next_     = new_node;
	return DirectedEdgeListIterator( this, new_node );
}

void DirectedEdgeList::erase( DirectedEdgeListIterator to_erase )
{
	debug_assert( to_erase.owner_ == this );
	DirectedEdgeListElement * next_node = to_erase.element_->next_;
	DirectedEdgeListElement * prev_node = to_erase.element_->previous_;

	debug_assert( next_node->previous_ == to_erase.element_ );
	debug_assert( prev_node->next_ == to_erase.element_ );

	// deallocate the list element
	//delete to_erase.element_;
	edge_list_element_pool_.destroy( to_erase.element_ );


	next_node->previous_ = prev_node;
	prev_node->next_     = next_node;

}

/// O(N)
platform::Size
DirectedEdgeList::size() const
{
	platform::Size nelements( 0 );
	for ( DirectedEdgeListConstIterator
			eiter = begin(), eiter_end = end();
			eiter != eiter_end; ++eiter ) {
		++nelements;
	}
	return nelements;
}


//----------------------------------------------------------------------------//
//------------------------------  Digraph DirectedNode Class ---------------------------//
//----------------------------------------------------------------------------//


// what is this? bool output_interaction_graph_memory_usage = false;


/// @brief
/// virtual destructor
DirectedNode::~DirectedNode()
{}

/// @brief
/// Main constructor, no default constructor nor copy constructor
DirectedNode::DirectedNode(Digraph * owner, platform::Size node_id ) :
	node_index_(node_id),
	num_incident_edges_(0),
	indegree_(0),
	outdegree_(0),
	incident_edge_list_( owner->edge_list_element_pool() ),
	first_outgoing_edge_( incident_edge_list_.end() ),
	owner_(owner)
{}

/// @brief copy-from for use in Digraph::operator= and copy ctors;
/// derived classes must define their own version of this function
//  to copy any data stored on nodes
void DirectedNode::copy_from( DirectedNode const * ) {}


/// @details Add an "incoming" edge (i.e. this is the tail node)
///
/// @param edge_ptr - [in] - the new edge
///
void
DirectedNode::add_incoming_edge( DirectedEdge* edge_ptr, DirectedEdgeListIter & eiter )
{
	++num_incident_edges_;
	++indegree_;
	eiter = incident_edge_list_.insert( incident_edge_list_.begin(), edge_ptr);
}

/// @details Add an "outgoing" edge (i.e. this is the head node)
///
/// @param edge_ptr - [in] - the new edge
///
void
DirectedNode::add_outgoing_edge( DirectedEdge* edge_ptr, DirectedEdgeListIter & eiter )
{
	++num_incident_edges_;
	++outdegree_;
	eiter = incident_edge_list_.insert( incident_edge_list_.end(), edge_ptr );
	if ( outdegree_ == 1 ) {
		first_outgoing_edge_ = eiter;
	}
}


/// @details edges efficiently delete themselves from the edge lists of the nodes they
/// are incident upon by keeping a pair of iterators.  DirectedEdges request nodes
/// delete them by handing the iterator back to the node.
///
/// @param
/// edge - [in] - the iterator for this node's edge list that points at the
///               edge which is trying to delete itself
///
void DirectedNode::drop_edge( DirectedEdgeListIter eiter )
{
	if ( first_outgoing_edge_ == eiter ) { ++first_outgoing_edge_; }

	if ( (*eiter)->is_tail_node( node_index_ ) ) {
		--outdegree_;
	} else {
		--indegree_;
	}
	incident_edge_list_.erase( eiter );
	--num_incident_edges_;
}


/// @details As edges delete themselves, they invalidate any iterators
/// that point to their (former) positions in the node and graph edge lists.
/// Therefore, before calling delete on an edge object, one must grab the next
/// iterator in a list.  Below, nextiter copies iter and is incremented before
/// iter's edge is deleted.  Note also that "++iter" does not appear in the
/// for loop.
void DirectedNode::drop_all_edges()
{
	for ( DirectedEdgeListIter iter = incident_edge_list_.begin();
			iter != incident_edge_list_.end(); /*no increment statement*/ ) {

		DirectedEdgeListIter nextiter = iter;
		++nextiter;
		owner_->delete_edge( *iter ); iter = nextiter;
	}
}


/// @details Constant time if each vertex has a constant number of edges.  DirectedEdges are
/// identified by the index of the node to which the edge connects this node.
/// Returns NULL when there is no such connecting edge.
///
/// @param
/// other_node - [in] - the index of the node that the desired
///   edge connects this node to
///
DirectedEdge const * DirectedNode::find_edge_to(platform::Size other_node) const
{
	//call non-const version of this function, which in fact does
	//not change any data, but simply returns a DirectedEdge * instead of
	//the desired DirectedEdge const *
	return const_cast< DirectedNode * > (this)->find_edge_to( other_node );
}

/// @brief non-const edge finding method; changes no data, but returns a non-const pointer
DirectedEdge * DirectedNode::find_edge_to(platform::Size other_node)
{
	DirectedEdgeListIter start, end;
	start = first_outgoing_edge_;
	end = incident_edge_list_.end();

	//iterate over range of edges
	for ( DirectedEdgeListIter iter = start; iter != end; ++iter ) {
		if ( (*iter)->same_edge( node_index_, other_node) ) {
			return (*iter);
		}
	}
	return 0;
}

/// @details Constant time if each vertex has a constant number of edges.  DirectedEdges are
/// identified by the index of the node to which the edge connects this node.
/// Returns NULL when there is no such connecting edge.
///
/// @param
/// other_node - [in] - the index of the node that the desired
///   edge connects this node to
///
DirectedEdge const * DirectedNode::find_edge_from( platform::Size other_node ) const
{
	//call non-const version of this function, which in fact does
	//not change any data, but simply returns a DirectedEdge * instead of
	//the desired DirectedEdge const *
	return const_cast< DirectedNode * > (this)->find_edge_from( other_node );
}

/// @brief non-const edge finding method; changes no data, but returns a non-const pointer
DirectedEdge * DirectedNode::find_edge_from( platform::Size other_node )
{
	DirectedEdgeListIter start, end;
	start = incident_edge_list_.begin();
	end = first_outgoing_edge_;

	//iterate over range of edges
	for ( DirectedEdgeListIter iter = start; iter != end; ++iter ) {
		if ( (*iter)->same_edge( other_node, node_index_ ) ) {
			return (*iter);
		}
	}
	return 0;
}


/// @brief virtual function to print node to standard out
void DirectedNode::print() const
{
	std::cout << "DirectedNode " << get_node_index() << " attached to edges: " << std::endl;
	for ( DirectedEdgeListConstIter
			iter = incident_edge_list_.const_begin(),
			iter_end = incident_edge_list_.const_end();
			iter != iter_end; ++iter ) {
		std::cout << "   DirectedEdge( " << (*iter)->get_tail_node_ind() << " --> ";
		std::cout << (*iter)->get_head_node_ind() << " )" << std::endl;
	}
}


/// @details called on most-derived class.  The most-derived class should NOT recursively call this method
/// on its parent class.  The sizeof function will handle the whole DirectedNode (or DerivedDirectedNode).
platform::Size DirectedNode::count_static_memory() const
{
	return sizeof( DirectedNode );
}


/// @details recursively descend through heirarchy accounting for heap memory usage.  Each derived
/// class in the heirarchy should recursively add the amount of dynamic memory its parent
/// allocates by calling parent::count_dynamic_memory
platform::Size DirectedNode::count_dynamic_memory() const
{
	platform::Size tot = 0;
	tot += sizeof( DirectedEdgeListElement ) * ( num_incident_edges_ + 1 ); // edge list
	return tot;
}

//----------------------------------------------------------------------------//
//------------------------------ Digraph DirectedEdge Class ----------------------------//
//----------------------------------------------------------------------------//

/// @brief destructor
///
/// @details removes all record of this edge from edge-lists of
/// the 1) nodes this edge is incident upon and 2) the owning
/// interaction graph
DirectedEdge::~DirectedEdge()
{
	nodes_[0]->drop_edge( pos_in_nodes_edge_list_[0] );
	nodes_[1]->drop_edge( pos_in_nodes_edge_list_[1] );
	owner_->drop_edge( pos_in_owners_edge_list_ );
}

/// @brief main constructor for edge, no default nor copy constructors
///
/// @details edge adds itself to the edge list of the two nodes its set to be
/// incident upon, and stores the list-iterators that the nodes return.
///
/// @param owner - [in] - owning InteractionDigraph
/// @param tail_node_ind - [in] - the index of the tail node
/// @param head_node_ind - [in] - the index of the head node
///
DirectedEdge::DirectedEdge
(
	Digraph* owner,
	platform::Size tail_node_ind,
	platform::Size head_node_ind
)
: owner_(owner)
{
	debug_assert( tail_node_ind != head_node_ind );
	node_indices_[0]    = tail_node_ind;
	node_indices_[1]    = head_node_ind;
	nodes_[0]           = owner->nodes_[ node_indices_[0] ];
	nodes_[1]           = owner->nodes_[ node_indices_[1] ];

	nodes_[0]->add_outgoing_edge( this, pos_in_nodes_edge_list_[0] );
	nodes_[1]->add_incoming_edge( this, pos_in_nodes_edge_list_[1] );

	return;
}


/// @details derived classes should recursively call the copy_from method to ensure all parent class
/// data is copied.  It just so happens that this method does nothing, but that could change
/// and the derived class should include a call to this function for that reason.
void DirectedEdge::copy_from( DirectedEdge const * ) {}

/// @brief returns the index of the other node that the edge is incident upon
platform::Size DirectedEdge::get_other_ind(platform::Size node_ind) const
{
	debug_assert( node_ind == node_indices_[0] || node_ind == node_indices_[1]);
	return node_indices_[0] == node_ind ? node_indices_[1] : node_indices_[0];
}

/// @brief returns a pointer to the other node that the edge is incident upon
DirectedNode * DirectedEdge::get_other_node(platform::Size node_ind)
{
	debug_assert( node_ind == node_indices_[0] || node_ind == node_indices_[1]);
	return node_indices_[0] == node_ind ? nodes_[1] : nodes_[0];
}

/// @brief return a const pointer to the other node that the edge is incident upon
DirectedNode const * DirectedEdge::get_other_node(platform::Size node_ind) const
{
	return const_cast< DirectedEdge * >(this)->get_other_node( node_ind );
}

/// @brief sets the iterator for this edge's position in its owner's edge list
void DirectedEdge::set_pos_in_owners_list( DirectedEdgeListIter iter )
{
	debug_assert( this == *iter );
	pos_in_owners_edge_list_ = iter;
	return;
}


/// @brief returns true if this edge connects nodes of index tail_node and head_node
/// the order of tail_node and head_node is not important
///
/// @param tail_node - [in] - index of the tail node
/// @param head_node - [in] - index of the head node
bool DirectedEdge::same_edge(platform::Size tail_node, platform::Size head_node) const
{
	return (tail_node == node_indices_[0] && head_node == node_indices_[1]);
}

/// @brief memory accouting scheme
///
/// @details This is called non-recursively on the most-derived class
platform::Size DirectedEdge::count_static_memory() const
{
	return sizeof( DirectedEdge );
}

/// @brief memory accounting scheme
///
/// @details This method should be called recursively by derived classes -- that is, each class should
/// recurse to its parent.
platform::Size DirectedEdge::count_dynamic_memory() const
{
	platform::Size tot = 0;
	//no dynamic memory
	return tot;
}

//----------------------------------------------------------------------------//
//---------------------------------  Digraph Class -----------------------------//
//----------------------------------------------------------------------------//

/// @brief destructor
Digraph::~Digraph()
{
	delete_everything();
	delete edge_list_element_pool_; edge_list_element_pool_ = 0;
	delete edge_pool_; edge_pool_ = 0;
}

/// @brief default constructor; creates an empty graph (no nodes, no edges)
Digraph::Digraph() :
	parent(),
	num_nodes_( 0 ),
	nodes_(),
	num_edges_( 0 ),
	edge_list_element_pool_( new boost::unordered_object_pool< DirectedEdgeListElement > ( 256 ) ),
	edge_list_( *edge_list_element_pool_ ),
	edge_pool_( new boost::unordered_object_pool< DirectedEdge > ( 256 ) ),
	focused_edge_( 0 )
{}

/// @details Do not call this constructor from a derived class in the initialization list,
/// since this constructor calls the polymorphic function create_new_node, and polymorphism
/// does not work during constructors or destructors.
///
/// @param num_ig_nodes - [in] - number of nodes that this graph will contain
Digraph::Digraph( platform::Size num_nodes ) :
	parent(),
	num_nodes_(num_nodes),
	nodes_(num_nodes, (DirectedNode*) 0),
	num_edges_( 0 ),
	edge_list_element_pool_( new boost::unordered_object_pool< DirectedEdgeListElement > ( 256 ) ),
	edge_list_( *edge_list_element_pool_ ),
	edge_pool_( new boost::unordered_object_pool< DirectedEdge > ( 256 ) ),
	focused_edge_( 0 )
{
	for ( platform::Size ii = 1; ii <= num_nodes; ++ii ) {
		nodes_[ ii ] = create_new_node( ii );
	}
}


/// @brief copy constructor relies on factory methods and virtual "copy_from" methods
///
/// @details This copy constructor should NOT be used by derived classes.  At the time
/// this is called, the identity of this has not yet been resolved -- this constructor will
/// produce DirectedNode and DirectedEdge objects when it calls "create_new_node" and not
/// DerivedDirectedNode or DerivedDirectedEdge objects.  Derived class copy constructors should call
/// the base class assignment operator once the initial construction has been completed.
Digraph::Digraph( Digraph const & source ) :
	parent(),
	utility::pointer::enable_shared_from_this< Digraph >(),
	num_nodes_( source.num_nodes_ ),
	nodes_( num_nodes_, (DirectedNode *) 0 ),
	num_edges_( 0 ),
	edge_list_element_pool_( new boost::unordered_object_pool< DirectedEdgeListElement > ( 256 ) ),
	edge_list_( *edge_list_element_pool_ ),
	edge_pool_( new boost::unordered_object_pool< DirectedEdge > ( 256 ) ),
	focused_edge_( 0 )
{
	for ( platform::Size ii = 1; ii <= num_nodes_; ++ii ) {
		nodes_[ ii ] = create_new_node( ii );
		nodes_[ ii ]->copy_from( source.nodes_[ii] );
	}

	for ( DirectedEdgeListConstIter
			iter = source.const_edge_list_begin(),
			iter_end = source.const_edge_list_end();
			iter != iter_end; ++iter ) {
		platform::Size const n1((*iter)->get_tail_node_ind()), n2((*iter)->get_head_node_ind());
		add_edge( n1, n2 );
		DirectedEdge * the_edge = find_edge( n1, n2 ); // O(1) op since focused_edge_ points to it already.
		the_edge->copy_from( *iter );
	}

}

/// @brief operator = ().  Relies on factory methods and virtual "copy_from" methods
///
/// @details operator= must only be performed on graphs of the same type e.g.
/// an EnergyDigraph may be copied from another EnergyDigraph, but should
/// not be copied from a Digraph.
Digraph &
Digraph::operator = ( Digraph const & source )
{
	if ( this == &source ) return *this;

	if ( num_nodes_ != source.num_nodes_ ) { set_num_nodes( source.num_nodes_ );}
	drop_all_edges();

	for ( platform::Size ii = 1; ii <= num_nodes_; ++ii ) {
		nodes_[ ii ]->copy_from( source.nodes_[ii] );
	}

	for ( DirectedEdgeListConstIter
			iter = source.const_edge_list_begin(),
			iter_end = source.const_edge_list_end();
			iter != iter_end; ++iter ) {
		add_edge( *iter ); // no longer calling copy_from method!
	}
	return *this;
}

/// @brief copy the connectivity of the source graph, but do not copy the data stored in the
/// nodes and edges of the source graph.  Useful for when copying a graph of a different type
/// (e.g. from an EnergyDigraph into a Digraph)
void Digraph::copy_connectivity( Digraph const & source )
{
	if ( num_nodes_ != source.num_nodes_ ) { set_num_nodes( source.num_nodes_ );}
	drop_all_edges();

	for ( DirectedEdgeListConstIter
			iter = source.const_edge_list_begin(),
			iter_end = source.const_edge_list_end();
			iter != iter_end; ++iter ) {
		platform::Size const n1((*iter)->get_tail_node_ind()), n2((*iter)->get_head_node_ind());
		add_edge( n1, n2 );
	}
}

/// @brief creates a new edge between nodes tail_index and head_index.  DirectedNodes do
/// not have to be listed in order.  For speed, does NOT check to see if
/// edge already exists -- except in debug mode.
///
/// @details uses factory method create_new_edge and adds the created edge to the graph's edge list.
/// Not threadsafe.  Only a single thread should add edges to the graph at a time.
///
/// @param tail_index - [in] - index of one of the two nodes the edge is to connect
/// @param head_index - [in] - index of the second of the two nodes the edge is to connect
DirectedEdge *
Digraph::add_edge(platform::Size tail_index, platform::Size head_index)
{
	debug_assert( ! get_edge_exists( tail_index, head_index ) );

	DirectedEdge* new_edge = create_new_edge(tail_index, head_index);
	edge_list_.push_back( new_edge );
	++num_edges_;
	new_edge->set_pos_in_owners_list( edge_list_.last() );
	focused_edge_ = new_edge;
	return new_edge;
}

/// @brief for use in Digraph operator=
///
/// @details Uses the edge copy constructor so that data stored
/// on edges of one graph may be placed rapidly into the new edge
/// of this graph.
DirectedEdge *
Digraph::add_edge( DirectedEdge const * example_edge )
{
	debug_assert( ! get_edge_exists(
		example_edge->get_tail_node_ind(),
		example_edge->get_head_node_ind() ) );

	DirectedEdge* new_edge = create_new_edge( example_edge );
	edge_list_.push_front( new_edge );
	++num_edges_;
	new_edge->set_pos_in_owners_list( edge_list_.begin() );
	focused_edge_ = new_edge;
	return new_edge;
}


/// @brief alternative to integer constructor; first create an empty graph and
/// later tell the graph how many nodes it has.  If the graph is not already
/// empty, it will delete everything its holding.
void Digraph::set_num_nodes( platform::Size num_nodes )
{
	delete_everything();
	num_nodes_ = num_nodes;
	nodes_.resize( num_nodes_ );
	for ( platform::Size ii = 1; ii <= num_nodes_; ++ii ) nodes_[ ii ] = create_new_node( ii );
}

void Digraph::add_node() {
	++num_nodes_;
	nodes_.resize( num_nodes_ );
	nodes_[ num_nodes_ ] = create_new_node( num_nodes_ );
}

/// @brief returns true if an edge between tail_node and head_node exists
///
/// @param tail_node - [in] - index of the one of the nodes
/// @param head_node - [in] - index of the other node
bool Digraph::get_edge_exists(platform::Size tail_node, platform::Size head_node) const
{
	DirectedEdge const * edge = find_edge( tail_node, head_node );
	return (edge != NULL);
}

/// @brief deletes all edges adjacent to the node specified
///
/// @param node - [in] - index of the node
void Digraph::drop_all_edges_for_node( platform::Size node )
{
	DirectedNode* nodeptr = get_node( node );
	nodeptr->drop_all_edges();
}

/// @brief deletes all edges in the graph
void Digraph::drop_all_edges()
{
	for ( DirectedEdgeListIter iter = edge_list_.begin(), iter_end = edge_list_.end();
			iter != iter_end; /*no increment*/ ) {

		DirectedEdgeListIter iter_next = iter;
		++iter_next;
		delete_edge(*iter);
		iter = iter_next;
	}
}


/// @brief calls print() on each of the nodes in the graph
void Digraph::print_vertices() const
{
	for ( platform::Size ii = 1; ii <= num_nodes_; ii++ ) {
		nodes_[ii]->print();
	}
	return;
}

/// @brief writes out a list of all the edges in the graph
///
/// @param os - [in] - the output stream to write to
void Digraph::output_connectivity(std::ostream & os) const
{
	platform::Size counter = 1;
	for ( DirectedEdgeListConstIter iter = edge_list_.begin(); iter != edge_list_.end(); ++iter ) {
		os << "edge " << counter << " between " << (*iter)->get_tail_node_ind()
			<< " " << (*iter)->get_head_node_ind() << std::endl;
		counter++;
	}
	return;
}

/// @brief writes out a connectivity description of the graph in the famous
/// dimacs format. (where the first column "DIMACS:" should be sed'ed out)
///
/// @param os - [in] - the output stream to write to
void Digraph::output_dimacs(std::ostream & os) const
{
	platform::Size num_edges = edge_list_.size();
	os << "DIMACS: " << "p edges " << num_nodes_ << " " ;
	os << num_edges << std::endl;
	for ( DirectedEdgeListConstIter iter = edge_list_.begin(); iter != edge_list_.end(); ++iter ) {
		os << "DIMACS: " << "e " << (*iter)->get_tail_node_ind();
		os << " " << (*iter)->get_head_node_ind() << std::endl;
	}

	return;
}


/// @brief removes edge from edge list at iterator iter
///
/// @details each edge keeps track of its position in its owner's graph's edge list
/// so it can efficiently delete itself should it need to.
///
/// @param iter - [in] - the iterator in the non-const edge list pointing at the edge that's deleting itself
/// @param citer - [in] - the iterator in the const edge list pointing at the edge that's deleting itself
void Digraph::drop_edge( DirectedEdgeListIter iter )
{
	if ( *iter == focused_edge_ ) focused_edge_ = NULL; //invalidate focused_edge_

	--num_edges_;
	edge_list_.erase(iter);

	return;
}

/// @brief deletes each edge in the graph and then deletes each node
///
/// @details its important to note that nodes must outlive their incident edges
void Digraph::delete_everything()
{
	for ( DirectedEdgeListIter iter = edge_list_.begin();
			iter != edge_list_.end(); /*no increment*/ ) {
		DirectedEdgeListIter next_iter = iter;
		++next_iter;
		delete_edge( *iter );
		iter = next_iter;
	}
	for ( platform::Size ii = 1; ii <= num_nodes_; ii++ ) { delete nodes_[ii]; nodes_[ii] = 0; }
	num_nodes_ = 0;
	nodes_.resize( 0 );
	focused_edge_ = 0;
}

/// @brief
/// returns the edge connecting tail_node and head_node (const version)
///
/// @details
/// graph keeps a pointer to the last edge that was accessed to that search is
/// fairly efficient.
///
/// @param
/// tail_node - [in] - index of the tail node
/// @param
/// head_node - [in] - index of the head node
DirectedEdge const * Digraph::find_edge(platform::Size tail_node, platform::Size head_node) const
{
	if ( focused_edge_ == NULL || !( focused_edge_->same_edge(tail_node, head_node)) ) {
		if ( nodes_[ tail_node ]->outdegree() < nodes_[ head_node ]->indegree() ) {
			focused_edge_ = nodes_[tail_node]->find_edge_to(head_node);
		} else {
			focused_edge_ = nodes_[head_node]->find_edge_from(tail_node);
		}
	}
	return focused_edge_;
}

/// @brief
/// returns the edge connecting tail_node and head_node
///
/// @details graph keeps a pointer to the last edge that was accessed to that search is
/// fairly efficient.
///
/// @param
/// tail_node - [in] - index of the first node
/// @param
/// head_node - [in] - index of the second node
DirectedEdge * Digraph::find_edge(platform::Size tail_node, platform::Size head_node)
{
	if ( focused_edge_ == NULL || !( focused_edge_->same_edge(tail_node, head_node)) ) {
		if ( nodes_[ tail_node ]->outdegree() < nodes_[ head_node ]->indegree() ) {
			focused_edge_ = nodes_[tail_node]->find_edge_to(head_node);
		} else {
			focused_edge_ = nodes_[head_node]->find_edge_from(tail_node);
		}
	}
	return focused_edge_;
}


void Digraph::delete_edge( DirectedEdge * edge )
{
	edge_pool_->destroy( edge );
}

platform::Size
Digraph::getTotalMemoryUsage() const
{
	platform::Size total_memory = 0;
	for ( platform::Size ii = 1; ii <= num_nodes(); ++ii ) {
		total_memory += nodes_[ ii ]->count_dynamic_memory();
		total_memory += nodes_[ ii ]->count_static_memory();
	}
	for ( DirectedEdgeListConstIter iter = edge_list_.const_begin();
			iter != edge_list_.const_end(); ++iter ) {
		total_memory += (*iter)->count_dynamic_memory();
		total_memory += (*iter)->count_static_memory();
	}

	total_memory += count_dynamic_memory();
	total_memory += count_static_memory();

	return total_memory;
}

platform::Size Digraph::count_static_memory() const
{
	return sizeof( Digraph );
}

platform::Size Digraph::count_dynamic_memory() const
{
	platform::Size tot = 0;
	tot += sizeof( DirectedNode* ) * num_nodes_;
	tot += sizeof( DirectedEdgeListElement ) * ( num_edges_ + 1 ); // edge list
	return tot;
}


/// @brief factory method for node creation
///   Should be overriden in derived classes
DirectedNode* Digraph::create_new_node( platform::Size index )
{
	return new DirectedNode( this, index );
}

/// @brief factory method for edge creation
///   Should be overriden in derived classes
DirectedEdge* Digraph::create_new_edge( platform::Size tail_index, platform::Size head_index )
{
	return edge_pool_->construct( this, tail_index, head_index );
}

DirectedEdge* Digraph::create_new_edge( DirectedEdge const * example_edge )
{
	return edge_pool_->construct(
		this,
		example_edge->get_tail_node_ind(),
		example_edge->get_head_node_ind()
	);
}

#ifdef    SERIALIZATION
template < class Archive >
void Digraph::save( Archive & archive ) const
{
  archive( num_nodes_ );

	// DirectedNodes and edges will be freshly created when this graph is deserialized
	// EXEMPT nodes_ edge_pool_ edge_list_element_pool_ edge_list_ focused_edge_

  // for ( Size ii = 1; ii <= num_nodes_; ++ii ) {
  //   nodes_[ ii ]->save( archive );
  // }
  archive( num_edges_ );
  for ( DirectedEdgeListConstIter iter = const_edge_list_begin(), iter_end = const_edge_list_end(); iter != iter_end; ++iter ) {
    archive( (*iter)->get_tail_node_ind(), (*iter)->get_head_node_ind() );
    // (*iter)->save( archive );
  }
}

template < class Archive >
void Digraph::load( Archive & archive )
{
  Size num_nodes(0); archive( num_nodes );
  set_num_nodes( num_nodes );
	// EXEMPT num_nodes_
	// The nodes will be freshly instantiated on this end
	// so they will not be deserialized.  Same for edges.
	// EXEMPT nodes_ num_edges_ edge_pool_ edge_list_element_pool_ edge_list_ focused_edge_

  // for ( Size ii = 1; ii <= num_nodes; ++ii ) {
  //   nodes_[ ii ]->load( archive );
  // }

  Size num_edges(0); archive( num_edges );
  for ( Size ii = 1; ii <= num_edges; ++ii ) {
    Size tail_node(0), head_node(0); archive( tail_node, head_node );
    /* DirectedEdge * new_edge = */ add_edge( tail_node, head_node );
    // new_edge->load( archive );
  }
}

SAVE_AND_LOAD_SERIALIZABLE( Digraph );
#endif // SERIALIZATION


// Variables for the topological sort algorithm
platform::Size const NOT_VISITED = 0;
platform::Size const TEMPORARY_VISITED = 1;
platform::Size const VISITED = 2;

void
visit(
	Digraph const & g,
	platform::Size node_index,
	utility::vector1< platform::Size > & visited_status,
	std::list< platform::Size > & toposort_order,
	bool & is_DAG
)
{
	if ( visited_status[ node_index ] == TEMPORARY_VISITED ) {
		is_DAG = false;
		return;
	}
	if ( visited_status[ node_index ] == NOT_VISITED ) {
		visited_status[ node_index ] = TEMPORARY_VISITED;
		DirectedNode const & node = *g.get_node( node_index );
		for ( DirectedNode::DirectedEdgeListConstIter iter = node.const_outgoing_edge_list_begin();
				iter != node.const_outgoing_edge_list_end(); ++iter ) {
			visit( g, (*iter)->get_head_node_ind(), visited_status, toposort_order, is_DAG );
			if ( ! is_DAG ) { return; }
		}
		visited_status[ node_index ] = VISITED;
		toposort_order.push_front( node_index );
	}
}

platform::Size
find_unmarked_node(
	utility::vector1< platform::Size > & visited_status,
	platform::Size last_descend_from
)
{
	for ( platform::Size ii = last_descend_from+1; ii <= visited_status.size(); ++ii ) {
		if ( visited_status[ ii ] == NOT_VISITED ) {
			return ii;
		}
	}
	return 0;
}

std::pair< std::list< platform::Size >, bool >
topological_sort( Digraph const & g )
{
	using platform::Size;

	utility::vector1< Size > visited_status( g.num_nodes(), NOT_VISITED );
	std::list< Size > toposort_order;

	bool is_dag = true;
	Size last_descend_from( 0 );
	while ( true ) {

		// find a starting node for this iteration
		Size descend_from = find_unmarked_node( visited_status, last_descend_from );
		if ( descend_from == 0 ) break;
		last_descend_from = descend_from;

		visit( g, descend_from, visited_status, toposort_order, is_dag );
		if ( ! is_dag ) {
			return std::make_pair( std::list< Size >(), false );
		}
	}
	return std::make_pair( toposort_order, true );
}


bool
digraph_is_a_DAG( Digraph const & g )
{
	std::pair< std::list< platform::Size >, bool > result = topological_sort( g );
	return result.second;
}


} //end namespace graph
} //end namespace utility

#ifdef    SERIALIZATION
CEREAL_REGISTER_TYPE( utility::graph::Digraph )
CEREAL_REGISTER_DYNAMIC_INIT( utility_graph_Digraph )
#endif // SERIALIZATION
