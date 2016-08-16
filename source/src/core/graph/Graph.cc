// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/graph/Graph.cc
/// @brief  graph base classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/graph/Graph.hh>

// Package headers
#include <core/graph/unordered_object_pool.hpp>

//STL Headers
#include <iostream>
#include <utility/assert.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

// Boost Headers
#include <core/graph/unordered_object_pool.hpp>

#include <boost/pool/pool.hpp>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


// #include <basic/Tracer.hh>

//static THREAD_LOCAL basic::Tracer TR( "core.graph.Graph" );
using namespace ObjexxFCL;

namespace core {
namespace graph {

EdgeList::EdgeList(
	boost::unordered_object_pool< EdgeListElement > & edge_list_element_pool
)
:
	edge_list_element_pool_( edge_list_element_pool ),
	end_( new EdgeListElement )
{
	end_->next_ = end_;
	end_->previous_ = end_;
}

EdgeList::~EdgeList()
{
	for ( EdgeListIterator
			eiter = begin(), eiter_end = end();
			eiter != eiter_end; /* no increment */ ) {
		EdgeListIterator eiter_next( eiter );
		++eiter_next;

		// remove the list element
		//delete eiter.element_;
		edge_list_element_pool_.destroy( eiter.element_ );
		eiter = eiter_next;
	}
	delete end_;
}

void
EdgeList::push_back( Edge * edgeptr )
{
	debug_assert( edgeptr ); // Do not accept null-pointing edges

	//EdgeListElement * new_node = new EdgeListElement( edgeptr, end_->previous_, end_ );
	EdgeListElement * new_node =  edge_list_element_pool_.construct( edgeptr, end_->previous_, end_ );

	end_->previous_->next_ = new_node; // may be end_->next_ if 0-element list
	end_->previous_ = new_node;
}

void
EdgeList::push_front( Edge * edgeptr )
{
	debug_assert( edgeptr ); // Do not accept null-pointing edges

	//EdgeListElement * new_node = new EdgeListElement( edgeptr, end_, end_->next_ );
	EdgeListElement * new_node =  edge_list_element_pool_.construct( edgeptr, end_, end_->next_ );

	end_->next_->previous_ = new_node; // may be end_->next_ if 0-element list
	end_->next_ = new_node;

}


EdgeListIterator
EdgeList::insert(
	EdgeListIterator const & element_to_insert_before,
	Edge * edgeptr
)
{
	debug_assert( edgeptr );
	debug_assert( element_to_insert_before.owner_ == this );

	EdgeListElement * next_node = element_to_insert_before.element_;
	EdgeListElement * prev_node = element_to_insert_before.element_->previous_;

	//EdgeListElement * new_node = new EdgeListElement( edgeptr, prev_node, next_node );
	EdgeListElement * new_node =  edge_list_element_pool_.construct( edgeptr, prev_node, next_node );

	next_node->previous_ = new_node;
	prev_node->next_     = new_node;
	return EdgeListIterator( this, new_node );
}

void EdgeList::erase( EdgeListIterator to_erase )
{
	debug_assert( to_erase.owner_ == this );
	EdgeListElement * next_node = to_erase.element_->next_;
	EdgeListElement * prev_node = to_erase.element_->previous_;

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
EdgeList::size() const
{
	platform::Size nelements( 0 );
	for ( EdgeListConstIterator
			eiter = begin(), eiter_end = end();
			eiter != eiter_end; ++eiter ) {
		++nelements;
	}
	return nelements;
}


//----------------------------------------------------------------------------//
//------------------------------  Graph Node Class ---------------------------//
//----------------------------------------------------------------------------//


bool output_interaction_graph_memory_usage = false;


/// @brief
/// virtual destructor
Node::~Node()
{}

/// @brief
/// Main constructor, no default constructor nor copy constructor
Node::Node(Graph * owner, platform::Size node_id ) :
	node_index_(node_id), num_incident_edges_(0),
	num_neighbors_counting_self_static_(1),
	loop_incident_( false ),
	num_edges_to_smaller_indexed_nodes_(0),
	num_edges_to_larger_indexed_nodes_(0),
	incident_edge_list_( owner->edge_list_element_pool() ),
	first_upper_edge_( incident_edge_list_.end() ),
	owner_(owner)
{}

/// @brief copy-from for use in Graph::operator= and copy ctors;
/// derived classes must define their own version of this function
//  to copy any data stored on nodes
void Node::copy_from( Node const * ) {}


/// @details If the other node this node is attached to by edge_ptr has a higher index
/// then the edge is added to the end of its edge list; if the node has a
/// smaller index, the edge pointer is added to the front of the edge list.
/// The presence of a new edge means the edge vector is not up to date.
///
/// @param edge_ptr - [in] - the new edge
///
void
Node::add_edge( Edge* edge_ptr, EdgeListIter & eiter )
{
	++num_incident_edges_;
	platform::Size other_node_index = edge_ptr->get_other_ind( node_index_);
	if ( other_node_index <  node_index_ ) {
		++num_edges_to_smaller_indexed_nodes_;
		eiter = incident_edge_list_.insert( incident_edge_list_.begin(), edge_ptr);
	} else {
		if ( !( edge_ptr->is_loop() && loop_incident_) ) {
			eiter = incident_edge_list_.insert(incident_edge_list_.end(), edge_ptr);
			if ( edge_ptr->is_loop() ) loop_incident_ = true;
		} else {
			//loop already attached; return 0 eiter/ceiter as a dummy
			// fixing odd iterator behavior with g++ -v 4.1.1
			debug_assert( num_edges_to_larger_indexed_nodes_ != 0 );
			eiter = incident_edge_list_.end(); //eiter = 0;
		}

		++num_edges_to_larger_indexed_nodes_;
		if ( num_edges_to_larger_indexed_nodes_ == 1 ) {
			first_upper_edge_ = eiter;
		}
	}
	num_neighbors_counting_self_static_ = num_neighbors_counting_self();
}


/// @details edges efficiently delete themselves from the edge lists of the nodes they
/// are incident upon by keeping a pair of iterators.  Edges request nodes
/// delete them by handing the iterator back to the node.
///
/// @param
/// edge - [in] - the iterator for this node's edge list that points at the
///               edge which is trying to delete itself
///
void Node::drop_edge( EdgeListIter eiter )
{
	if ( first_upper_edge_ == eiter ) { ++first_upper_edge_; }

	if ( ! (*eiter)->is_loop() ) {

		platform::Size other_node_index = (*eiter)->get_other_ind( node_index_ );
		if ( node_index_ < other_node_index ) {
			--num_edges_to_larger_indexed_nodes_;
		} else {
			--num_edges_to_smaller_indexed_nodes_;
		}
		incident_edge_list_.erase( eiter );
		--num_incident_edges_;
	} else {
		if ( loop_incident_ ) {
			loop_incident_ = false;

			incident_edge_list_.erase(eiter);

			--num_edges_to_larger_indexed_nodes_;
			--num_incident_edges_;
		}
	}
	num_neighbors_counting_self_static_ = num_neighbors_counting_self();
}


/// @details As edges delete themselves, they invalidate any iterators
/// that point to their (former) positions in the node and graph edge lists.
/// Therefore, before calling delete on an edge object, one must grab the next
/// iterator in a list.  Below, nextiter copies iter and is incremented before
/// iter's edge is deleted.  Note also that "++iter" does not appear in the
/// for loop.
void Node::drop_all_edges()
{
	for ( EdgeListIter iter = incident_edge_list_.begin();
			iter != incident_edge_list_.end(); /*no increment statement*/ ) {

		EdgeListIter nextiter = iter;
		++nextiter;
		owner_->delete_edge( *iter ); iter = nextiter;
	}
}

/// @details manually change the number of neighbors for a Node. Used
/// for symmetry scoring
void Node::set_num_neighbors_counting_self_static( platform::Size neighbor )
{
	num_neighbors_counting_self_static_ = neighbor;
}

/// @details Constant time if each vertex has a constant number of edges.  Edges are
/// identified by the index of the node to which the edge connects this node.
/// Returns NULL when there is no such connecting edge.
///
/// @param
/// other_node - [in] - the index of the node that the desired
///   edge connects this node to
///
Edge const * Node::find_edge(platform::Size other_node) const
{
	//call non-const version of this function, which in fact does
	//not change any data, but simply returns a Edge * instead of
	//the desired Edge const *
	return const_cast< Node * > (this)->find_edge( other_node );
}

/// @brief non-const edge finding method; changes no data, but returns a non-const pointer
Edge * Node::find_edge(platform::Size other_node)
{
	EdgeListIter start, end;
	if ( other_node > get_node_index() ) {
		start = first_upper_edge_;
		end = incident_edge_list_.end();
	} else {
		start = incident_edge_list_.begin();
		end = first_upper_edge_;
	}

	//iterate over range of edges
	for ( EdgeListIter iter = start; iter != end; ++iter ) {
		if ( (*iter)->same_edge( node_index_, other_node) ) {
			return (*iter);
		}
	}
	return 0;
}


/// @brief virtual function to print node to standard out
void Node::print() const
{
	// TR.Debug << "Node " << get_node_index() << " attached to edges: " << std::endl;
	for ( EdgeListConstIter
			iter = incident_edge_list_.const_begin(),
			iter_end = incident_edge_list_.const_end();
			iter != iter_end; ++iter ) {
		//  TR.Debug << "   Edge( " << (*iter)->get_first_node_ind() << ", ";
		//  TR.Debug << (*iter)->get_second_node_ind() << ")" << std::endl;
	}
}


/// @details called on most-derived class.  The most-derived class should NOT recursively call this method
/// on its parent class.  The sizeof function will handle the whole Node (or DerivedNode).
platform::Size Node::count_static_memory() const
{
	return sizeof( Node );
}


/// @details recursively descend through heirarchy accounting for heap memory usage.  Each derived
/// class in the heirarchy should recursively add the amount of dynamic memory its parent
/// allocates by calling parent::count_dynamic_memory
platform::Size Node::count_dynamic_memory() const
{
	platform::Size tot = 0;
	tot += sizeof( EdgeListElement ) * ( num_incident_edges_ + 1 ); // edge list
	return tot;
}

//----------------------------------------------------------------------------//
//------------------------------ Graph Edge Class ----------------------------//
//----------------------------------------------------------------------------//

/// @brief destructor
///
/// @details removes all record of this edge from edge-lists of
/// the 1) nodes this edge is incident upon and 2) the owning
/// interaction graph
Edge::~Edge()
{
	nodes_[0]->drop_edge( pos_in_nodes_edge_list_[0] );
	if ( !is_loop() ) {
		nodes_[1]->drop_edge( pos_in_nodes_edge_list_[1] );
	}
	owner_->drop_edge( pos_in_owners_edge_list_ );
}

/// @brief main constructor for edge, no default nor copy constructors
///
/// @details edge adds itself to the edge list of the two nodes its set to be
/// incident upon, and stores the list-iterators that the nodes return.
///
/// @param owner - [in] - owning InteractionGraph
/// @param first_node_ind - [in] - the index of the first node
/// @param second_node_ind - [in] - the index of the second node
///
Edge::Edge
(
	Graph* owner,
	platform::Size first_node_ind,
	platform::Size second_node_ind
)
: owner_(owner)
{
	debug_assert( first_node_ind <= second_node_ind );
	node_indices_[0]    = first_node_ind;
	node_indices_[1]    = second_node_ind;
	nodes_[0]           = owner->nodes_[ node_indices_[0] ];
	nodes_[1]           = owner->nodes_[ node_indices_[1] ];

	nodes_[0]->add_edge( this, pos_in_nodes_edge_list_[0] );
	nodes_[1]->add_edge( this, pos_in_nodes_edge_list_[1] );

	return;
}


/// @details derived classes should recursively call the copy_from method to ensure all parent class
/// data is copied.  It just so happens that this method does nothing, but that could change
/// and the derived class should include a call to this function for that reason.
void Edge::copy_from( Edge const * ) {}

/// @brief returns the index of the other node that the edge is incident upon
platform::Size Edge::get_other_ind(platform::Size node_ind) const
{
	debug_assert( node_ind == node_indices_[0] || node_ind == node_indices_[1]);
	return node_indices_[0] == node_ind ? node_indices_[1] : node_indices_[0];
}

/// @brief returns a pointer to the other node that the edge is incident upon
Node * Edge::get_other_node(platform::Size node_ind)
{
	debug_assert( node_ind == node_indices_[0] || node_ind == node_indices_[1]);
	return node_indices_[0] == node_ind ? nodes_[1] : nodes_[0];
}

/// @brief return a const pointer to the other node that the edge is incident upon
Node const * Edge::get_other_node(platform::Size node_ind) const
{
	return const_cast< Edge * >(this)->get_other_node( node_ind );
}

/// @brief sets the iterator for this edge's position in its owner's edge list
void Edge::set_pos_in_owners_list( EdgeListIter iter )
{
	debug_assert( this == *iter );
	pos_in_owners_edge_list_ = iter;
	return;
}


/// @brief returns true if this edge connects nodes of index node1 and node2
/// the order of node1 and node2 is not important
///
/// @param node1 - [in] - index of one of the two nodes
/// @param node2 - [in] - index of the other of the two nodes
bool Edge::same_edge(platform::Size node1, platform::Size node2) const
{
	if ( node1 > node2 ) {
		platform::Size temp = node2;
		node2 = node1;
		node1 = temp;
	}
	return (node1 == node_indices_[0] && node2 == node_indices_[1]);
}

/// @brief memory accouting scheme
///
/// @details This is called non-recursively on the most-derived class
platform::Size Edge::count_static_memory() const
{
	return sizeof( Edge );
}

/// @brief memory accounting scheme
///
/// @details This method should be called recursively by derived classes -- that is, each class should
/// recurse to its parent.
platform::Size Edge::count_dynamic_memory() const
{
	platform::Size tot = 0;
	//no dynamic memory
	return tot;
}

//----------------------------------------------------------------------------//
//---------------------------------  Graph Class -----------------------------//
//----------------------------------------------------------------------------//

/// @brief destructor
Graph::~Graph()
{
	delete_everything();
	delete edge_list_element_pool_; edge_list_element_pool_ = 0;
	delete edge_pool_; edge_pool_ = 0;
}

/// @brief default constructor; creates an empty graph (no nodes, no edges)
Graph::Graph() :
	parent(),
	num_nodes_( 0 ),
	nodes_(),
	num_edges_( 0 ),
	edge_list_element_pool_( new boost::unordered_object_pool< EdgeListElement > ( 256 ) ),
	edge_list_( *edge_list_element_pool_ ),
	edge_pool_( new boost::unordered_object_pool< Edge > ( 256 ) ),
	focused_edge_( 0 )
{}

/// @details Do not call this constructor from a derived class in the initialization list,
/// since this constructor calls the polymorphic function create_new_node, and polymorphism
/// does not work during constructors or destructors.
///
/// @param num_ig_nodes - [in] - number of nodes that this graph will contain
Graph::Graph( platform::Size num_nodes ) :
	parent(),
	num_nodes_(num_nodes),
	nodes_(num_nodes, (Node*) 0),
	num_edges_( 0 ),
	edge_list_element_pool_( new boost::unordered_object_pool< EdgeListElement > ( 256 ) ),
	edge_list_( *edge_list_element_pool_ ),
	edge_pool_( new boost::unordered_object_pool< Edge > ( 256 ) ),
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
/// produce Node and Edge objects when it calls "create_new_node" and not
/// DerivedNode or DerivedEdge objects.  Derived class copy constructors should call
/// the base class assignment operator once the initial construction has been completed.
Graph::Graph( Graph const & source ) :
	parent(),
	utility::pointer::enable_shared_from_this< Graph >(),
	num_nodes_( source.num_nodes_ ),
	nodes_( num_nodes_, (Node *) 0 ),
	num_edges_( 0 ),
	edge_list_element_pool_( new boost::unordered_object_pool< EdgeListElement > ( 256 ) ),
	edge_list_( *edge_list_element_pool_ ),
	edge_pool_( new boost::unordered_object_pool< Edge > ( 256 ) ),
	focused_edge_( 0 )
{
	for ( platform::Size ii = 1; ii <= num_nodes_; ++ii ) {
		nodes_[ ii ] = create_new_node( ii );
		nodes_[ ii ]->copy_from( source.nodes_[ii] );
	}

	for ( EdgeListConstIter
			iter = source.const_edge_list_begin(),
			iter_end = source.const_edge_list_end();
			iter != iter_end; ++iter ) {
		platform::Size const n1((*iter)->get_first_node_ind()), n2((*iter)->get_second_node_ind());
		add_edge( n1, n2 );
		Edge * the_edge = find_edge( n1, n2 ); // O(1) op since focused_edge_ points to it already.
		the_edge->copy_from( *iter );
	}

}

/// @brief operator = ().  Relies on factory methods and virtual "copy_from" methods
///
/// @details operator= must only be performed on graphs of the same type e.g.
/// an EnergyGraph may be copied from another EnergyGraph, but should
/// not be copied from a Graph.
Graph &
Graph::operator = ( Graph const & source )
{
	if ( this == &source ) return *this;

	if ( num_nodes_ != source.num_nodes_ ) { set_num_nodes( source.num_nodes_ );}
	drop_all_edges();

	for ( platform::Size ii = 1; ii <= num_nodes_; ++ii ) {
		nodes_[ ii ]->copy_from( source.nodes_[ii] );
	}

	for ( EdgeListConstIter
			iter = source.const_edge_list_begin(),
			iter_end = source.const_edge_list_end();
			iter != iter_end; ++iter ) {
		add_edge( *iter ); // no longer calling copy_from method!
	}
	return *this;
}

/// @brief copy the connectivity of the source graph, but do not copy the data stored in the
/// nodes and edges of the source graph.  Useful for when copying a graph of a different type
/// (e.g. from an EnergyGraph into a Graph)
void Graph::copy_connectivity( Graph const & source )
{
	if ( num_nodes_ != source.num_nodes_ ) { set_num_nodes( source.num_nodes_ );}
	drop_all_edges();

	for ( EdgeListConstIter
			iter = source.const_edge_list_begin(),
			iter_end = source.const_edge_list_end();
			iter != iter_end; ++iter ) {
		platform::Size const n1((*iter)->get_first_node_ind()), n2((*iter)->get_second_node_ind());
		add_edge( n1, n2 );
	}
}

/// @brief creates a new edge between nodes index1 and index2.  Nodes do
/// not have to be listed in order.  For speed, does NOT check to see if
/// edge already exists -- except in debug mode.
///
/// @details uses factory method create_new_edge and adds the created edge to the graph's edge list.
/// Not threadsafe.  Only a single thread should add edges to the graph at a time.
///
/// @param index1 - [in] - index of one of the two nodes the edge is to connect
/// @param index2 - [in] - index of the second of the two nodes the edge is to connect
Edge *
Graph::add_edge(platform::Size index1, platform::Size index2)
{
	debug_assert( ! get_edge_exists( index1, index2 ) );

	//swap so that index1 < index2
	platform::Size temp = index1 < index2 ? index1 : index2;
	index2 = index1 < index2 ? index2 : index1;
	index1 = temp;

	Edge* new_edge = create_new_edge(index1, index2);
	edge_list_.push_back( new_edge );
	++num_edges_;
	new_edge->set_pos_in_owners_list( edge_list_.last() );
	focused_edge_ = new_edge;
	return new_edge;
}

/// @brief for use in Graph operator=
///
/// @details Uses the edge copy constructor so that data stored
/// on edges of one graph may be placed rapidly into the new edge
/// of this graph.
Edge *
Graph::add_edge( Edge const * example_edge )
{
	debug_assert( ! get_edge_exists(
		example_edge->get_first_node_ind(),
		example_edge->get_second_node_ind() ) );

	Edge* new_edge = create_new_edge( example_edge );
	edge_list_.push_front( new_edge );
	++num_edges_;
	new_edge->set_pos_in_owners_list( edge_list_.begin() );
	focused_edge_ = new_edge;
	return new_edge;
}


/// @brief alternative to integer constructor; first create an empty graph and
/// later tell the graph how many nodes it has.  If the graph is not already
/// empty, it will delete everything its holding.
void Graph::set_num_nodes( platform::Size num_nodes )
{
	delete_everything();
	num_nodes_ = num_nodes;
	nodes_.resize( num_nodes_ );
	for ( platform::Size ii = 1; ii <= num_nodes_; ++ii ) nodes_[ ii ] = create_new_node( ii );
}

/// @brief returns true if an edge between node1 and node2 exists
///
/// @param node1 - [in] - index of the one of the nodes
/// @param node2 - [in] - index of the other node
bool Graph::get_edge_exists(platform::Size node1, platform::Size node2) const
{
	Edge const * edge = find_edge( node1, node2 );
	return (edge != NULL);
}

/// @brief deletes all edges adjacent to the node specified
///
/// @param node - [in] - index of the node
void Graph::drop_all_edges_for_node( platform::Size node )
{
	Node* nodeptr = get_node( node );
	nodeptr->drop_all_edges();
}

/// @brief deletes all edges in the graph
void Graph::drop_all_edges()
{
	for ( EdgeListIter iter = edge_list_.begin(), iter_end = edge_list_.end();
			iter != iter_end; /*no increment*/ ) {

		EdgeListIter iter_next = iter;
		++iter_next;
		delete_edge(*iter);
		iter = iter_next;
	}
}


/// @brief calls print() on each of the nodes in the graph
void Graph::print_vertices() const
{
	for ( platform::Size ii = 1; ii <= num_nodes_; ii++ ) {
		nodes_[ii]->print();
	}
	return;
}

/// @brief writes out a list of all the edges in the graph
///
/// @param os - [in] - the output stream to write to
void Graph::output_connectivity(std::ostream & os) const
{
	platform::Size counter = 1;
	for ( EdgeListConstIter iter = edge_list_.begin(); iter != edge_list_.end(); ++iter ) {
		os << "edge " << counter << " between " << (*iter)->get_first_node_ind()
			<< " " << (*iter)->get_second_node_ind() << std::endl;
		counter++;
	}
	return;
}

/// @brief writes out a connectivity description of the graph in the famous
/// dimacs format. (where the first column "DIMACS:" should be sed'ed out)
///
/// @param os - [in] - the output stream to write to
void Graph::output_dimacs(std::ostream & os) const
{
	platform::Size num_edges = edge_list_.size();
	os << "DIMACS: " << "p edges " << num_nodes_ << " " ;
	os << num_edges << std::endl;
	for ( EdgeListConstIter iter = edge_list_.begin(); iter != edge_list_.end(); ++iter ) {
		os << "DIMACS: " << "e " << (*iter)->get_first_node_ind();
		os << " " << (*iter)->get_second_node_ind() << std::endl;
	}

	return;
}


/// @brief compute all node distances and return FArray2D with that information
///
/// O(n^3)
FArray2D_int
Graph::all_pairs_shortest_paths() const
{
	platform::Size const inf( 12345678 ); //assumption: fewer than 12 million nodes in the graph.
	debug_assert( num_nodes_ < inf );

	FArray2D_int distance_table( num_nodes_, num_nodes_, inf);
	for ( platform::Size ii = 1; ii <= num_nodes_; ++ii ) distance_table( ii, ii ) = 0; //nodes are 0 distance from themselves.

	for ( EdgeListConstIter iter = edge_list_.begin();
			iter != edge_list_.end(); ++iter ) {
		platform::Size n1 = (*iter)->get_first_node_ind();
		platform::Size n2 = (*iter)->get_second_node_ind();

		if ( ! (*iter)->is_loop() ) {
			distance_table( n2, n1 ) = 1;
			distance_table( n1, n2 ) = 1;
		}

	}

	// Warshall algorithm
	// symmetry makes this marginally inefficient, but easy to read
	// if this shows up in a hotspot, it can be made more efficient
	for ( platform::Size ii = 1; ii <= num_nodes_; ++ii ) {
		for ( platform::Size jj = 1; jj <= num_nodes_; ++jj ) {
			for ( platform::Size kk = 1; kk <= num_nodes_; ++kk ) {
				int const jj_2_kk = distance_table( jj, kk );
				int const jj_2_ii = distance_table( jj, ii );
				int const ii_2_kk = distance_table( ii, kk );

				int const jj_2_ii_2_kk = jj_2_ii + ii_2_kk;

				if ( jj_2_kk > jj_2_ii_2_kk ) {
					distance_table( jj, kk ) =  jj_2_ii_2_kk;
					distance_table( kk, jj ) =  jj_2_ii_2_kk;
				}
			}
		}
	}

	return distance_table;
}


/// @brief removes edge from edge list at iterator iter
///
/// @details each edge keeps track of its position in its owner's graph's edge list
/// so it can efficiently delete itself should it need to.
///
/// @param iter - [in] - the iterator in the non-const edge list pointing at the edge that's deleting itself
/// @param citer - [in] - the iterator in the const edge list pointing at the edge that's deleting itself
void Graph::drop_edge( EdgeListIter iter )
{
	if ( *iter == focused_edge_ ) focused_edge_ = NULL; //invalidate focused_edge_

	--num_edges_;
	edge_list_.erase(iter);

	return;
}

/// @brief deletes each edge in the graph and then deletes each node
///
/// @details its important to note that nodes must outlive their incident edges
void Graph::delete_everything()
{
	for ( EdgeListIter iter = edge_list_.begin();
			iter != edge_list_.end(); /*no increment*/ ) {
		EdgeListIter next_iter = iter;
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
/// returns the edge connecting node1 and node2 (const version)
///
/// @details
/// graph keeps a pointer to the last edge that was accessed to that search is
/// fairly efficient.
///
/// @param
/// node1 - [in] - index of the first node
/// @param
/// node2 - [in] - index of the second node
Edge const * Graph::find_edge(platform::Size node1, platform::Size node2) const
{
	if ( focused_edge_ == NULL || !( focused_edge_->same_edge(node1, node2)) ) {
		focused_edge_ = nodes_[node1]->find_edge(node2);
	}
	return focused_edge_;
}

/// @brief
/// returns the edge connecting node1 and node2
///
/// @details graph keeps a pointer to the last edge that was accessed to that search is
/// fairly efficient.
///
/// @param
/// node1 - [in] - index of the first node
/// @param
/// node2 - [in] - index of the second node
Edge * Graph::find_edge(platform::Size node1, platform::Size node2)
{
	if ( focused_edge_ == NULL || !( focused_edge_->same_edge(node1, node2)) ) {
		focused_edge_ = nodes_[node1]->find_edge(node2);
	}
	return focused_edge_;
}


void Graph::delete_edge( Edge * edge )
{
	edge_pool_->destroy( edge );
}

platform::Size
Graph::getTotalMemoryUsage() const
{
	platform::Size total_memory = 0;
	for ( platform::Size ii = 1; ii <= num_nodes(); ++ii ) {
		total_memory += nodes_[ ii ]->count_dynamic_memory();
		total_memory += nodes_[ ii ]->count_static_memory();
	}
	for ( EdgeListConstIter iter = edge_list_.const_begin();
			iter != edge_list_.const_end(); ++iter ) {
		total_memory += (*iter)->count_dynamic_memory();
		total_memory += (*iter)->count_static_memory();
	}

	total_memory += count_dynamic_memory();
	total_memory += count_static_memory();

	return total_memory;
}

platform::Size Graph::count_static_memory() const
{
	return sizeof( Graph );
}

platform::Size Graph::count_dynamic_memory() const
{
	platform::Size tot = 0;
	tot += sizeof( Node* ) * num_nodes_;
	tot += sizeof( EdgeListElement ) * ( num_edges_ + 1 ); // edge list
	return tot;
}


/// @brief factory method for node creation
///   Should be overriden in derived classes
Node* Graph::create_new_node( platform::Size index )
{
	return new Node( this, index );
}

/// @brief factory method for edge creation
///   Should be overriden in derived classes
Edge* Graph::create_new_edge( platform::Size index1, platform::Size index2 )
{
	return edge_pool_->construct( this, index1, index2 );
}

Edge* Graph::create_new_edge( Edge const * example_edge )
{
	return edge_pool_->construct(
		this,
		example_edge->get_first_node_ind(),
		example_edge->get_second_node_ind()
	);
}

#ifdef    SERIALIZATION
template < class Archive >
void Graph::save( Archive & archive ) const
{
  archive( num_nodes_ );

	// Nodes and edges will be freshly created when this graph is deserialized
	// EXEMPT nodes_ edge_pool_ edge_list_element_pool_ edge_list_ focused_edge_

  // for ( Size ii = 1; ii <= num_nodes_; ++ii ) {
  //   nodes_[ ii ]->save( archive );
  // }
  archive( num_edges_ );
  for ( EdgeListConstIter iter = const_edge_list_begin(), iter_end = const_edge_list_end(); iter != iter_end; ++iter ) {
    archive( (*iter)->get_first_node_ind(), (*iter)->get_second_node_ind() );
    // (*iter)->save( archive );
  }
}

template < class Archive >
void Graph::load( Archive & archive )
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
    Size node1(0), node2(0); archive( node1, node2 );
    /* Edge * new_edge = */ add_edge( node1, node2 );
    // new_edge->load( archive );
  }
}

SAVE_AND_LOAD_SERIALIZABLE( Graph );
#endif // SERIALIZATION

} //end namespace graph
} //end namespace core

#ifdef    SERIALIZATION
CEREAL_REGISTER_TYPE( core::graph::Graph )
CEREAL_REGISTER_DYNAMIC_INIT( core_graph_Graph )
#endif // SERIALIZATION
