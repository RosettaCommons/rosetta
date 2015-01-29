// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/InteractionGraphBase.cc
/// @brief  Interaction graph base class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/pack/interaction_graph/InteractionGraphBase.hh>

/// Utility headers
#include <utility/exit.hh>

//STL Headers
#include <list>
// AUTO-REMOVED #include <vector>
#include <algorithm>
#include <iostream>
#include <utility/assert.hh>

//ObjexxFCL Headers
#include <ObjexxFCL/FArray1A.hh>

using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace interaction_graph {

//----------------------------------------------------------------------------//
//---------------------- Interaction Graph Node Base Class -------------------//
//----------------------------------------------------------------------------//


//bool output_interaction_graph_memory_usage = false;

////////////////////////////////////////////////////////////////////////////////
/// @begin NodeBase::~NodeBase
///
/// @brief
/// virtual destructor
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
NodeBase::~NodeBase()
{}

////////////////////////////////////////////////////////////////////////////////
/// @begin NodeBase::NodeBase(InteractionGraphBase *, int, int)
///
/// @brief
/// Main constructor, no default constructor nor copy constructor
///
/// @detailed
///
/// @param
/// owner - [in] - the owning interaction graph
/// node_id - [in] - the index for this node amongst its owners set
/// num_states - [in] - the number of states for this node
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
NodeBase::NodeBase(InteractionGraphBase * owner, int node_id, int num_states) :
	node_index_(node_id), num_states_(num_states), num_incident_edges_(0),
	num_edges_to_smaller_indexed_nodes_(0),
	num_edges_to_larger_indexed_nodes_(0),
	edge_vector_up_to_date_(false), owner_(owner)
{}


////////////////////////////////////////////////////////////////////////////////
/// @begin NodeBase::~get_num_states
///
/// @brief
/// returns the number of states for this node
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
int  NodeBase::get_num_states() const {return num_states_;}

////////////////////////////////////////////////////////////////////////////////
/// @begin NodeBase::add_edge
///
/// @brief adds edge pointer to edge list; returns an iterator to the new
/// list element
///
/// @detailed
/// If the other node this node is attached to by edge_ptr has a higher index
/// then the edge is added to the end of its edge list; if the node has a
/// smaller index, the edge pointer is added to the front of the edge list.
/// The presence of a new edge means the edge vector is not up to date.
///
/// @param
/// edge_ptr - [in] - the new edge
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
std::list< EdgeBase* >::iterator NodeBase::add_edge(EdgeBase* edge_ptr)
{
	++num_incident_edges_;
	edge_vector_up_to_date_ = false;
	int other_node_index = edge_ptr->get_other_ind( node_index_);
	if (other_node_index <  node_index_) {
		++num_edges_to_smaller_indexed_nodes_;
		return incident_edge_list_.insert( incident_edge_list_.begin(), edge_ptr);
	} else {
		++num_edges_to_larger_indexed_nodes_;
		return incident_edge_list_.insert(incident_edge_list_.end(), edge_ptr);
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @begin NodeBase::drop_edge
///
/// @brief removes an edge iterator from the node's edge list
///
/// @detailed
/// edges efficiently delete themselves from the edge lists of the nodes they
/// are incident upon by keeping a pair of iterators.  Edges request nodes
/// delete them by handing the iterator back to the node.
///
/// @param
/// edge - [in] - the iterator for this node's edge list that points at the
///               edge which is trying to delete itself
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void NodeBase::drop_edge(std::list< EdgeBase* >::iterator edge)
{
	int other_node_index = (*edge)->get_other_ind( node_index_ );
	if (node_index_ < other_node_index) {
		--num_edges_to_larger_indexed_nodes_;
	} else {
		--num_edges_to_smaller_indexed_nodes_;
	}
	incident_edge_list_.erase(edge);
	--num_incident_edges_;
	edge_vector_up_to_date_ = false;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin NodeBase::drop_all_edges
///
/// @brief
/// deletes all edges incident upon this node
///
/// @detailed
///
/// @param
//
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void NodeBase::drop_all_edges()
{
	for (std::list< EdgeBase* >::iterator iter = incident_edge_list_.begin();
			iter != incident_edge_list_.end(); ) {
		std::list< EdgeBase* >::iterator nextiter = iter;
		++nextiter;
		delete *iter; iter = nextiter;
	}
	edge_vector_up_to_date_ = false;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin NodeBase::find_edge
///
/// @brief a slow (linear) search for an edge.  The edge is identified by the
/// index of the node to which the edge connects this node. Returns NULL when
/// there is no such connecting edge.
///
/// @detailed
///
/// @param
/// other_node - [in] - the index of the node that the desired
///   edge connects this node to
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
EdgeBase* NodeBase::find_edge(int other_node) const
{
	for (std::list< EdgeBase* >::const_iterator iter = incident_edge_list_.begin();
			iter != incident_edge_list_.end(); iter++) {
		if ( (*iter)->same_edge( node_index_, other_node) ) return (*iter);
	}
	return NULL;
}


////////////////////////////////////////////////////////////////////////////////
/// @begin NodeBase::depth_first_connected_component_counting
///
/// @brief
/// performs a depth first traversal of the graph.  Each node informs
/// the graph that the traversal resulted in arriving at the node.
///
/// @detailed
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void NodeBase::depth_first_connected_component_counting()
{
	if ( owner_->vertex_already_reached( node_index_) ) return;

	//std::cerr << "Arrived at node: " << node_index_ << std::endl;
	owner_->note_vertex_reached( node_index_ );

	for (std::list< EdgeBase* >::const_iterator iter = incident_edge_list_.begin();
			iter != incident_edge_list_.end(); ++iter) {
		((*iter)->get_other_node( node_index_ ))->
			depth_first_connected_component_counting();
	}

}


std::list< EdgeBase * >::const_iterator
NodeBase::edge_list_begin()
{
	return incident_edge_list_.begin();
}

std::list< EdgeBase * >::const_iterator
NodeBase::edge_list_end()
{
	return incident_edge_list_.end();
}


////////////////////////////////////////////////////////////////////////////////
/// @begin NodeBase::update_edge_vector
///
/// @brief converts edge-list to edge-vector representation
///
/// @detailed
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void NodeBase::update_edge_vector()
{
	incident_edge_vector_.resize( num_incident_edges_ + 1);
	adjacent_node_ind_.resize( num_incident_edges_ + 1);
	adjacent_node_.resize(num_incident_edges_ + 1);

	std::vector< EdgeBase* >::iterator position1 = incident_edge_vector_.begin();
	++position1;

	std::copy( incident_edge_list_.begin(), incident_edge_list_.end(), position1 );

	for (int ii = 1; ii <= num_incident_edges_; ++ii) {
		incident_edge_vector_[ii]->set_pos_in_node_edgevector( node_index_, ii );
		adjacent_node_ind_[ ii ] = incident_edge_vector_[ii]->get_other_ind( node_index_);

	debug_assert( (adjacent_node_ind_[ii] < node_index_ &&
			ii <= num_edges_to_smaller_indexed_nodes_)
			||
			( adjacent_node_ind_[ii] > node_index_ &&
			ii > num_edges_to_smaller_indexed_nodes_ ) );

		adjacent_node_[ ii ] = incident_edge_vector_[ii]->get_other_node( node_index_ );
	}

	edge_vector_up_to_date_ = true;
	return;
}

/// @brief memory accounting scheme
unsigned int
NodeBase::count_dynamic_memory() const
{
	unsigned int total_memory = 0;
	total_memory += 4 * incident_edge_list_.size() * sizeof( EdgeBase* );
	total_memory += incident_edge_vector_.size() * sizeof( EdgeBase* );
	total_memory += adjacent_node_ind_.size() * sizeof( int );
	total_memory += adjacent_node_.size() * sizeof( NodeBase* );
	return total_memory;
}


////////////////////////////////////////////////////////////////////////////////
/// @begin NodeBase::NodeBase( NodeBase const & rhs)
///
/// @brief copy constructor, do not use
///
/// @detailed
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
//NodeBase::NodeBase( NodeBase const & ) {}

//----------------------------------------------------------------------------//
//--------------------- Interaction Graph Edge Base Class --------------------//
//----------------------------------------------------------------------------//

////////////////////////////////////////////////////////////////////////////////
/// @begin EdgeBase::~EdgeBase
///
/// @brief destructor
///
/// @detailed removes all record of this edge from edge-lists of
/// the 1) nodes this edge is incident upon and 2) the owning
/// interaction graph
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
EdgeBase::~EdgeBase()
{
	//std::cerr << "~EdgeBase(): " << this << " node 1 " << node_indices_[0]
	//<< " node 2 " << node_indices_[1] << " owner " <<
	//*pos_in_owners_edge_list_ << " ";
//debug_assert( this == *pos_in_owners_edge_list_);
	nodes_[0]->drop_edge(pos_in_nodes_edge_list_[0]);
	nodes_[1]->drop_edge(pos_in_nodes_edge_list_[1]);
	owner_->drop_edge(pos_in_owners_edge_list_);
}

////////////////////////////////////////////////////////////////////////////////
/// @begin EdgeBase::EdgeBase(InteractionGraphBase, int, int)
///
/// @brief main constructor for edge, no default nor copy constructors
///
/// @detailed edge adds itself to the edge list of the two nodes its set to be
/// incident upon, and stores the list-iterators that the nodes return.
///
/// @param
/// owner - [in] - owning InteractionGraph
/// @param
/// first_node_ind - [in] - the index of the first node
/// @param
/// second_node_ind - [in] - the index of the second node
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
EdgeBase::EdgeBase
(
	InteractionGraphBase* owner,
	int first_node_ind,
	int second_node_ind
)
	: owner_(owner),
	edge_weight_( 1.0 )
{  //pre condition: first_node_ind < second_node_ind.
	node_indices_[0]    = first_node_ind;
	node_indices_[1]    = second_node_ind;
	nodes_[0]           = owner->ig_nodes_[ node_indices_[0] ];
	nodes_[1]           = owner->ig_nodes_[ node_indices_[1] ];
	num_node_states_[0] = nodes_[0]->get_num_states();
	num_node_states_[1] = nodes_[1]->get_num_states();

	pos_in_nodes_edge_list_[0] = nodes_[0]->add_edge(this);
	pos_in_nodes_edge_list_[1] = nodes_[1]->add_edge(this);

	return;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin EdgeBase::get_other_ind
///
/// @brief returns the index of the other node that the edge is incident upon
///
/// @detailed
///
/// @param
/// node_ind - [in] - the node index of the node whose index is already known
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
int EdgeBase::get_other_ind(int node_ind) const
{      debug_assert( node_ind == node_indices_[0] || node_ind == node_indices_[1]);
	return node_indices_[0] == node_ind ? node_indices_[1] : node_indices_[0];
}

////////////////////////////////////////////////////////////////////////////////
/// @begin EdgeBase::get_other_node
///
/// @brief returns a pointer to the other node that the edge is incident upon
///
/// @detailed
///
/// @param
/// node_ind - [in] - the node index of the node whose index is already known
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
NodeBase* EdgeBase::get_other_node(int node_ind) const
{ debug_assert( node_ind == node_indices_[0] || node_ind == node_indices_[1]);
	return node_indices_[0] == node_ind ? nodes_[1] : nodes_[0];
}

////////////////////////////////////////////////////////////////////////////////
/// @begin EdgeBase::get_first_node_ind
///
/// @brief returns the index of the smaller-indexed node
///
/// @detailed
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
int EdgeBase::get_first_node_ind() const
{       return node_indices_[0]; }

////////////////////////////////////////////////////////////////////////////////
/// @begin EdgeBase::get_second_node_ind()
///
/// @brief returns the index of the larger-indexed node
///
/// @detailed
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
int EdgeBase::get_second_node_ind() const
{       return node_indices_[1]; }

////////////////////////////////////////////////////////////////////////////////
/// @begin EdgeBase::set_pos_in_owners_list
///
/// @brief edge keeps iterator to its position in it's owner's edge list
///
/// @detailed
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void EdgeBase::set_pos_in_owners_list( std::list< EdgeBase* >::iterator iter )
{
debug_assert( this == *iter);
	pos_in_owners_edge_list_ = iter;
	return;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin EdgeBase::set_pos_in_node_edgevector
///
/// @brief edge keeps index it has in node_ind's edge vector
///
/// @detailed
///
/// @param
/// node_ind - [in] - the index of the node calling this method
/// @param
/// vect_position - [in] - the position for this edge in the node's edge vector
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void EdgeBase::set_pos_in_node_edgevector(int node_ind, int vect_position)
{
debug_assert( node_ind == node_indices_[0] || node_ind == node_indices_[1]);
	int node_pos = (node_ind == node_indices_[0] ? 0 : 1 );
	pos_in_nodes_edge_vector_[node_pos] = vect_position;
	return;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin EdgeBase::same_edge
///
/// @brief returns true if this edge connects nodes of index node1 and node2
/// the order of node1 and node2 is not important
///
/// @detailed
///
/// @param
/// node1 - [in] - index of one of the two nodes
/// @param
/// node2 - [in] - index of the other of the two nodes
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
bool EdgeBase::same_edge(int node1, int node2) const
{
	//swap
	if (node1 > node2)
	{
			int temp = node2;
			node2 = node1;
			node1 = temp;
	}
	return (node1 == node_indices_[0] && node2 == node_indices_[1]);
}


unsigned int
EdgeBase::count_dynamic_memory() const
{
	return 0;
}

void
EdgeBase::edge_weight( Real new_weight )
{
	edge_weight_ = new_weight;
}


//----------------------------------------------------------------------------//
//------------------------- Interaction Graph Base Class ---------------------//
//----------------------------------------------------------------------------//

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::~InteractionGraphBase
///
/// @brief destructor
///
/// @detailed deletes each edge in the graph and deletes each node
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
InteractionGraphBase::~InteractionGraphBase()
{

	for (std::list< EdgeBase* >::iterator iter = ig_edge_list_.begin();
			iter != ig_edge_list_.end(); )
	{
			std::list< EdgeBase* >::iterator next_iter = iter;
			next_iter++;
			delete (*iter);
			iter = next_iter;
	}
	for (int ii = 1; ii <= num_ig_nodes_; ii++)
			delete ig_nodes_[ii];
}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::InteractionGraphBase
///
/// @brief main constructor
///
/// @detailed no default or copy constructors provided.
///
/// @param
/// num_ig_nodes - [in] - number of nodes that this graph will contain
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
InteractionGraphBase::InteractionGraphBase(int num_ig_nodes) :
	num_ig_nodes_(num_ig_nodes),
	ig_nodes_(num_ig_nodes + 1, (NodeBase*) NULL),
	node_state_offsets_( num_ig_nodes + 1, 0 ),
	num_total_states_(0),
	focused_edge_(NULL),
	num_energy_sum_groups_( -1 )
{}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::set_num_states_for_node
///
/// @brief sets the number of states for a node of a particular index
/// NEW REQUIREMENT: Nodes must have their num-states set in ascending order by
/// node index; that is, node 1 must go first, node 2 next, and so on.
///
/// @detailed once the graph knows how many states a node has, it instantiates
/// a new node using the NodeBase(int) constructor through the graph's
/// factory method create_new_node()
///
/// @param
/// node_index - [in] - the index of the node
/// @param
/// num_states - [in] - the number of states for that node
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void InteractionGraphBase::set_num_states_for_node
(
	int node_index,
	int num_states
)
{
debug_assert (ig_nodes_[node_index] == NULL);
	ig_nodes_[node_index] = create_new_node( node_index, num_states);
	num_total_states_ += num_states;
	if ( node_index != num_ig_nodes_ )
	{
			node_state_offsets_[ node_index + 1 ] =
					node_state_offsets_[ node_index ] + num_states;
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::get_num_states_for_node
///
/// @brief returns the number of states for a particular node
///
/// @detailed
///
/// @param
/// node_index - [in] - the index of the node in question
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
int  InteractionGraphBase::get_num_states_for_node(int node_index) const
{
debug_assert( ig_nodes_[node_index] );
	return ig_nodes_[node_index]->get_num_states();
}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::add_edge
///
/// @brief creates a new edge between nodes index1 and index2.  Nodes do
/// not have to be listed in order
///
/// @detailed uses factory method create_new_edge and adds
/// the created edge to the graph's edge list.
///
/// @param
/// index1 - [in] - index of one of the two nodes the edge is to connect
/// @param
/// index2 - [in] - index of the second of the two nodes the edge is to connect
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void InteractionGraphBase::add_edge(int index1, int index2)
{
	//swap so that index1 < index2
	int temp = index1 < index2 ? index1 : index2;
	index2 = index1 < index2 ? index2 : index1;
	index1 = temp;

debug_assert( index1 != index2 );

	EdgeBase* new_edge = create_new_edge(index1, index2);
	ig_edge_list_.push_front( new_edge );
	new_edge->set_pos_in_owners_list( ig_edge_list_.begin() );
	focused_edge_ = new_edge;
	return;
}


////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::get_edge_exists
///
/// @brief returns true if an edge between node1 and node2 exists
///
/// @detailed
///
/// @param
/// node1 - [in] - index of the one of the nodes
/// node2 - [in] - index of the other node
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
bool InteractionGraphBase::get_edge_exists(int node1, int node2)
{
	EdgeBase* edge = find_edge( node1, node2 );
	return (edge != NULL);
}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::drop_all_edges_for_node
///
/// @brief
/// deletes all edges adjacent to the node specified
///
/// @detailed
///
/// @param
/// node - [in] - index of the node
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void InteractionGraphBase::drop_all_edges_for_node( int node )
{
	NodeBase* nodeptr = get_node( node );
	nodeptr->drop_all_edges();
}

/// @details Edges may decide to delete themselves during this subroutine; therefore
/// edges are prepared first.  Afterwards, the nodes must update their
/// edge vector representation.
void InteractionGraphBase::prepare_for_simulated_annealing()
{
	for (std::list< EdgeBase* >::iterator iter = get_edge_list_begin();
			iter != get_edge_list_end();
			/* note: no increment statement here */ ) {
		std::list< EdgeBase* >::iterator next_iter = iter;
		next_iter++;
		//edges sometimes delete themselves, invalidating iterators, so
		//get the next iterator before calling prepare_for_simulated_annealing
		(*iter)->prepare_for_simulated_annealing();
		iter = next_iter;
	}

	for (int ii = 1; ii <= get_num_nodes(); ++ii) {
		get_node(ii)->prepare_for_simulated_annealing();
	}
	return;
}


////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::print_vertices
///
/// @brief calls print() on each of the nodes in the graph
///
/// @detailed
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void InteractionGraphBase::print_vertices() const
{
	for (int ii = 1; ii <= num_ig_nodes_; ii++)
			ig_nodes_[ii]->print();
	return;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::output_connectivity
///
/// @brief writes out a list of all the edges in the graph
///
/// @detailed
///
/// @param
/// os - [in] - the output stream to write to
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void InteractionGraphBase::output_connectivity(std::ostream & os) const
{
	int counter = 1;
	for (std::list< EdgeBase* >::const_iterator iter = ig_edge_list_.begin();
			iter != ig_edge_list_.end(); iter++)
	{  os << "edge " << counter << " between " << (*iter)->get_first_node_ind()
					<< " " << (*iter)->get_second_node_ind() << std::endl;
			counter++;
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::output_dimacs
///
/// @brief writes out a connectivity description of the graph in the famous
/// dimacs format. (where the first column "DIMACS:" should be sed'ed out)
///
/// @detailed
///
/// @param
/// os - [in] - the output stream to write to
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void InteractionGraphBase::output_dimacs(std::ostream & os) const
{
	int num_edges = ig_edge_list_.size();
	os << "DIMACS: " << "p edges " << num_ig_nodes_ << " " ;
	os << num_edges << std::endl;
	for (std::list< EdgeBase* >::const_iterator iter = ig_edge_list_.begin();
			iter != ig_edge_list_.end(); iter++)
	{
			os << "DIMACS: " << "e " << (*iter)->get_first_node_ind();
			os << " " << (*iter)->get_second_node_ind() << std::endl;
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::output_dimacs
///
/// @brief
/// Returns true if any node in the graph is in state 0, the unassigned state.
///
/// @detailed
/// Useful for debugging.  If simulated annealing completes, and any vertex
/// remains in state 0, then the state assignment is not meaningful nor
/// is the energy for that assignment.  The cases in which state-0 problems
/// turn up are from passing an annealer a list of states to restrict itself to
/// wherein some vertex has no states listed.  Such cases are bugs; this
/// subroutine helps identify them.
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
bool
InteractionGraphBase::any_vertex_state_unassigned() const
{
	for (int ii = 1; ii <= get_num_nodes(); ++ii)
	{
		if ( get_node( ii )->state_unassigned() ) return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::add_to_one_body_energies
///
/// @brief
/// takes one FArray of energies -- one energy for each state for each node
///
/// @detailed
/// The input array should have \sum_{v\inV} |S_v| entries, where S_v is the
/// state-space for vertex v.  If the graph has two vertices, with 15 and 20
/// states respectively, then entries 1 to 15 correspond to vertex 1, and
/// entries 16 to 36 correspond to vertex 2.  This compact enumeration scheme
/// is the same as the enumeration scheme used by rotindex, rotcoord, etc
///
/// @param
/// one_body_energies - [in] - the array of one body energies
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void InteractionGraphBase::add_to_one_body_energies
(
	FArray1< core::PackerEnergy > & one_body_energies
)
{
	for (int ii = 1; ii <= num_ig_nodes_; ++ii )
	{
			int ii_num_states;
			if ( ii == num_ig_nodes_ )
			{
					ii_num_states = num_total_states_ - node_state_offsets_[ii];
			}
			else
			{
					ii_num_states = node_state_offsets_[ ii + 1 ] - node_state_offsets_[ ii ];
			}
			FArray1A< core::PackerEnergy > ii_one_body_energies( one_body_energies( node_state_offsets_[ii] + 1), ii_num_states );
			ig_nodes_[ ii ]->add_to_one_body_energies( ii_one_body_energies );
	}
}


////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::add_to_one_body_energies
///
/// @brief
/// decrements the one body energies by the values held in old_energy1b,
/// increments the one body energies by the values held in new_energy1b,
/// and copies new_energy1b into old_energy1b.
///
/// @detailed
///
/// @param
/// old_energy1b - [in/out] - the one body energies representing interactions
///   with portions of the background that are no longer valid
/// new_energy1b - [in] - the one body energies representing interactions
///   with the background after the background has changed
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// useful if you want to move something in the background like a ligand.
/// moving the background does not invalidate the two-body energies.
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void InteractionGraphBase::update_one_body_energies
(
	FArray1< core::PackerEnergy > & old_energy1b,
	FArray1< core::PackerEnergy > & new_energy1b
)
{
	old_energy1b *= -1;
	add_to_one_body_energies( old_energy1b );
	add_to_one_body_energies( new_energy1b );
	old_energy1b = new_energy1b;
}

void InteractionGraphBase::zero_one_body_energies_for_node( int node )
{
	ig_nodes_[ node ]->zero_one_body_energies();
}


/// @param node_ind - [in] - the node in question
/// @param one_body_energies - [in] - the energies to be added to the one-body energies
/// 	on that node.  One entry per state.
void InteractionGraphBase::add_to_nodes_one_body_energy
(
	int node_ind,
	utility::vector1< core::PackerEnergy > const & one_body_energies
)
{
	int num_states_for_node = get_node( node_ind )->get_num_states();
	for (int ii = 1; ii <= num_states_for_node; ++ii ) {
		get_node( node_ind )->add_to_one_body_energy( ii, one_body_energies[ ii ] );
	}
	return;
}

/// @param node_ind - [in] - the node in question
/// @param one_body_energies - [in] - the energies to be added to the one-body energies
/// 	on that node.  One entry per state.
void InteractionGraphBase::add_to_nodes_one_body_energy
(
	int node_ind,
	FArray1< core::PackerEnergy > const & one_body_energies
)
{
	int num_states_for_node = get_node( node_ind )->get_num_states();
	for (int ii = 1; ii <= num_states_for_node; ++ii ) {
		get_node( node_ind )->add_to_one_body_energy( ii, one_body_energies( ii ) );
	}
	return;
}

/// @param node_ind - [in] - the index of the node in question
/// @param state_id - [in] - the state in question
/// @param one_body_energy - [in] - the energy to be added
void InteractionGraphBase::add_to_nodes_one_body_energy
(
	int node_ind,
	int state_id,
	core::PackerEnergy const one_body_energy
)
{
	get_node( node_ind )->add_to_one_body_energy( state_id, one_body_energy );
	return;
}

/// @details invokes polymorphic method of EdgeBase class; no-op if edge is not part
/// of the graph.
void
InteractionGraphBase::set_edge_weight(
	int node1,
	int node2,
	Real edge_weight
)
{
	if ( edge_weight == 0.0 ) {
		utility_exit_with_message( "Error: set edge weight to 0 not a legal operation.  Consider an edge deletion. ");
	}

	EdgeBase * edge = find_edge( node1, node2 );
	if ( edge ) {
		edge->set_edge_weight( edge_weight );
	}
}

/// @details returns 0 if edge is not part of the graph
Real
InteractionGraphBase::get_edge_weight(
	int node1,
	int node2
) const
{
	EdgeBase const * edge = find_edge( node1, node2 );
	if ( edge ) {
		return edge->edge_weight();
	} else {
		return 0.0;
	}
}


// @brief Memory account algorithm.  Each node/edge/graph class in the heirarchy reports the dynamic
// memory it uses, and "recurses" on it's parent class to have the parent report its dynamic memory use.
// (NOTE: class writer must actually implement the recursive call, it will not happen automatically.)
// The static memory is reported by the "final" class or the most-derived class, in a straitforward sizeof
// call.  This way, all the dynamically allocated memory is reported and the statically allocated memory
// is reported without being double counted.
unsigned int
InteractionGraphBase::getTotalMemoryUsage() const
{
	//std::cout << "calling InteractionGraphBase::getTotalMemoryUsage() const" << std::endl;

	unsigned int total_memory = 0;
	for (int ii = 1; ii <= num_ig_nodes_; ++ii) {
		total_memory += ig_nodes_[ ii ]->count_dynamic_memory();
		total_memory += ig_nodes_[ ii ]->count_static_memory();
	}
	for (std::list< EdgeBase* >::const_iterator iter = ig_edge_list_.begin();
			iter != ig_edge_list_.end(); iter++) {
		total_memory += (*iter)->count_dynamic_memory();
		total_memory += (*iter)->count_static_memory();
	}

	total_memory += count_dynamic_memory();
	total_memory += count_static_memory();

	return total_memory;
}

/// @brief set the Graph's (single) edge list iterator to the beginning of the edge list
/// for a particular node
void InteractionGraphBase::reset_edge_list_iterator_for_node( int node_index ) const
{
	focused_edge_iterator_ = ig_nodes_[ node_index ]->edge_list_begin();
	focused_edge_iterator_end_ = ig_nodes_[ node_index ]->edge_list_end();
}

/// @brief increment the (single) edge list iterator to the next element
void InteractionGraphBase::increment_edge_list_iterator() const
{
debug_assert( focused_edge_iterator_ != focused_edge_iterator_end_ );
	++focused_edge_iterator_;
}

/// @brief test: have we arrived at the edge list end?
bool InteractionGraphBase::edge_list_iterator_at_end() const
{
	return focused_edge_iterator_ == focused_edge_iterator_end_;
}

/// @brief return a const reference to an edge pointed at by the list iterator
EdgeBase const & InteractionGraphBase::get_edge() const
{
	return **focused_edge_iterator_;
}

// @brief IGBase dynamic memory: memory cost is for a vector of node pointers and for a list
// of edge pointers -- stl list elements cost four pointers (I think), so the math is 4 times
// the number of edges in the graph.  Since this is not a time-sensitive function, I use the
// list::size() method which is O(N).
//
// There is unaccounted memory in the base class ATM: ANDREW come back to this later and account
// for everything.
unsigned int
InteractionGraphBase::count_dynamic_memory() const
{
	unsigned int tot = 0;
	tot += sizeof( NodeBase* ) * num_ig_nodes_;
	tot += 4 * sizeof( EdgeBase* ) * ig_edge_list_.size();
	return tot;
}


////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::set_number_of_energy_sum_vertex_groups
///
/// @brief
/// a user may define subsets of the vertex set for which they would like to
/// know the internal energy sum.  For instance in a graph with 6 vertices,
/// {a,b,c,d,e,f}
/// a user may be interested in the sum of the one- and two-body energies
/// for vertices {a,b,c}.  The graph will return sum of the one body energies
/// for vertices a b and c and also any two-body energies for the edges in the
/// subgraph induced by a,b, and c.  (In this case, edges {a,b}, {a,c} and {b,c}
/// if these edges are part of the graph.  The edge {a,d} will not be counted
/// if it is part of the graph.)
///
/// First you must declare how many groups you are interested in.  Do that
/// with this method.
/// Second you must declare which node is a member of each group.  Only
/// tell the graph which node is a member, do not tell the graph if a node is
/// not a member.
/// Third, when you want to know the energy sum for the group in the graph's
/// current state assignment, call get_energy_sum_for_vertex_group( group_id)
///
/// @detailed
///
/// @param
/// num_groups - [in] - the number of groups; set this at most once.
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void
InteractionGraphBase::set_number_of_energy_sum_vertex_groups( int num_groups )
{
debug_assert( num_energy_sum_groups_ == -1 && num_groups > 0 );
	num_energy_sum_groups_ = num_groups;
	energy_sum_group_membership_.dimension(
			num_ig_nodes_, num_energy_sum_groups_ );
	energy_sum_group_membership_ = false;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::
///   count_connected_components_and_initialize_vertex_groups
///
/// @brief
/// makes a depth first traversal of the graph, counting the number of
/// connected components, and initializes the vertex group memberships
/// to reflect the connected components.  Returns the number of connected
/// components in the graph.
///
/// @detailed
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
int
InteractionGraphBase::count_connected_components_and_initialize_vertex_groups()
{
	component_membership_.dimension( num_ig_nodes_) = 0;
	num_energy_sum_groups_ = 0;
	for (int ii = 1; ii <= num_ig_nodes_; ++ii)
	{
			if ( vertex_already_reached( ii ) ) continue;
			++num_energy_sum_groups_;
			//std::cerr << "Starting depth first search at node: " << ii << " of group # " << num_energy_sum_groups_ << std::endl;

			ig_nodes_[ ii ]->depth_first_connected_component_counting();
	}

	energy_sum_group_membership_.dimension( num_ig_nodes_, num_energy_sum_groups_ );
	energy_sum_group_membership_ = false;
	for ( int ii = 1; ii <= num_ig_nodes_; ++ii )
	{
			energy_sum_group_membership_( ii, component_membership_(ii)) = true;
	}
	component_membership_.dimension( 0 );
	return num_energy_sum_groups_;
}


////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::note_vertex_reached
///
/// @brief
/// marks a vertex as belonging to the connected component currently being
/// traversed in the depth first traversal.
///
/// @detailed
///
/// @param
/// node_index - [in] - the index of the node invoking this method.
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// This method should be used by the NodeBase class only.
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void
InteractionGraphBase::note_vertex_reached( int node_index )
{
debug_assert( component_membership_( node_index ) == 0 );
	component_membership_( node_index ) = num_energy_sum_groups_;
	//std::cerr << "Marked node " << node_index << " in group " << num_energy_sum_groups_ << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::
///   vertex_already_reached
///
/// @brief
/// used by class NodeBase during the depth-first traversal to determine the
/// number of connected components in the graph. returns true if the dft has
/// already reached the node.
///
/// @detailed
///
/// @param
/// node_index - [in] - the index of the node calling the method.
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// This method should be used by the NodeBase class only.
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
bool
InteractionGraphBase::vertex_already_reached( int node_index )
{
	return component_membership_( node_index ) != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::set_vertex_member_of_group
///
/// @brief
/// a user may define subsets of the vertex set for which they would like to
/// know the internal energy sum.  For instance in a graph with 6 vertices,
/// {a,b,c,d,e,f}
/// a user may be interested in the sum of the one- and two-body energies
/// for vertices {a,b,c}.  The graph will return sum of the one body energies
/// for vertices a b and c and also any two-body energies for the edges in the
/// subgraph induced by a,b, and c.  (In this case, edges {a,b}, {a,c} and {b,c}
/// if these edges are part of the graph.  The edge {a,d} will not be counted
/// if it is part of the graph.)
///
/// @detailed
/// tell the graph which vertices you want to be part of which groups
///
/// @param
/// vertex - [in] - the index of the vertex you wish to include in the group
/// group - [in] - the group index you wish to add the vertex to.
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void
InteractionGraphBase::set_vertex_member_of_group( int vertex, int group )
{
	energy_sum_group_membership_( vertex, group ) = true;
}

void
InteractionGraphBase::print_vertex_groups()
{
	for (int ii = 1; ii <= get_num_nodes(); ++ii )
	{
			std::cerr << "Node " << ii << ": ";
			for (int jj = 1; jj <= num_energy_sum_groups_; ++jj)
			{
					std::cerr << energy_sum_group_membership_( ii, jj ) << " ";
			}
			std::cerr << std::endl;
	}
}
/*
unsigned int
InteractionGraphBase::getTotalMemoryUsage() const
{
	unsigned int total_memory = 0;
	for (int ii = 1; ii <= get_num_nodes(); ++ii)
	{
			total_memory += ig_nodes_[ ii ]->getMemoryUsageInBytes();
	}
	for (std::list< EdgeBase* >::const_iterator iter = ig_edge_list_.begin();
			iter != ig_edge_list_.end(); ++iter )
	{
			total_memory += (*iter)->getMemoryUsageInBytes();
	}

	total_memory += getMemoryUsageInBytes();
	return total_memory;
}

unsigned int
InteractionGraphBase::getMemoryUsageInBytes() const
{
	unsigned int total_memory = 0;
	total_memory += ig_nodes_.size() * sizeof ( NodeBase * );
	total_memory += 4 * ig_nodes_.size() * sizeof ( EdgeBase * );
	total_memory += node_state_offsets_.size() * sizeof ( int );
	return total_memory;
}
*/

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::drop_edge
///
/// @brief removes edge from edge list at iterator iter
///
/// @detailed
/// each edge keeps track of its position in its owner's graph's edge list
/// so it can efficiently delete itself should it need to.
///
/// @param
/// iter - [in] - the iterator pointing at the edge that's deleting itself
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void InteractionGraphBase::drop_edge(std::list< EdgeBase* >::iterator iter)
{
	//std::cerr << "ig.drop_edge: " << *iter << " ";
	ig_edge_list_.erase(iter);
	//std::cerr << *ig_edge_list_.begin() << " " << *(++ig_edge_list_.begin()) << std::endl;

	//invalidate focused_edge_
	focused_edge_ = NULL;
	return;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin InteractionGraphBase::find_edge
///
/// @brief
/// returns the edge connecting node1 and node2
///
/// @detailed
/// graph keeps a pointer to the last edge that was accessed to that search is
/// fairly efficient.
///
/// @param
/// node1 - [in] - index of the first node
/// @param
/// node2 - [in] - index of the second node
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
EdgeBase const * InteractionGraphBase::find_edge(int node1, int node2) const
{
	if (focused_edge_ == NULL || !( focused_edge_->same_edge(node1, node2)) ) {
		focused_edge_ = ig_nodes_[node1]->find_edge(node2);
	}
	return focused_edge_;
}

EdgeBase * InteractionGraphBase::find_edge(int node1, int node2)
{
	if ( focused_edge_ == NULL || !( focused_edge_->same_edge(node1, node2)) ) {
		focused_edge_ = ig_nodes_[node1]->find_edge(node2);
	}
	return focused_edge_;
}

} //end namespace interaction_graph
} //end namespace pack
} //end namespace core

