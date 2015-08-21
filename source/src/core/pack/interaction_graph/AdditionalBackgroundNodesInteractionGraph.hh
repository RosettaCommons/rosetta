// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/AdditionalBackgroundNodesInteractionGraph.hh
/// @brief  Interaction graph which can handle some residues not changing during packing.
/// @author Andrew Leaver-Fay
/// @author Ron Jacak (ron.jacak@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_AdditionalBackgroundNodesInteractionGraph_HH
#define INCLUDED_core_pack_interaction_graph_AdditionalBackgroundNodesInteractionGraph_HH


//STL Headers
#include <list>

#include <iostream>

#include <utility/vector1_bool.hh>


namespace core {
namespace pack {
namespace interaction_graph {


template < typename V, typename E, typename G > class FirstClassNode;
template < typename V, typename E, typename G > class FirstClassEdge;
template < typename V, typename E, typename G > class BackgroundNode;
template < typename V, typename E, typename G > class BackgroundToFirstClassEdge;
template < typename V, typename E, typename G > class AdditionalBackgroundNodesInteractionGraph;


//----------------------------------------------------------------------------//
//------------------------- First Class Node Class ---------------------------//
//----------------------------------------------------------------------------//

///
/// @brief
///
/// @details
/// A background node is a node which is not changing type or rotamer throughout a simulation.
///
///
/// @remarks
/// No public default constructor makes this class uncopyable.
///
template < typename V, typename E, typename G >
class FirstClassNode : public V {

public:
	typedef typename std::list< BackgroundToFirstClassEdge< V, E, G >* >             BackgroundEdgeList;
	typedef typename std::list< BackgroundToFirstClassEdge< V, E, G >* >::iterator   BackgroundEdgeListIter;

	typedef typename std::vector< BackgroundNode< V, E, G >* >                       BackgroundNodeVector;
	typedef typename std::vector< BackgroundToFirstClassEdge< V, E, G >* >           BackgroundEdgeVector;
	typedef typename std::vector< BackgroundToFirstClassEdge< V, E, G >* >::iterator BackgroundEdgeVectorIter;

public:
	virtual ~FirstClassNode();
	FirstClassNode( G * owner, int node_id, int num_states);
	BackgroundEdgeListIter add_background_edge(BackgroundToFirstClassEdge< V, E, G >* edge );
	void drop_background_edge( BackgroundEdgeListIter edge );
	virtual unsigned int count_dynamic_memory() const;

protected:
	inline
	int get_num_edges_to_background_nodes() const {
		return num_edges_to_bg_nodes_;
	}

	inline
	BackgroundToFirstClassEdge< V, E, G >* get_edge_to_bg_node( int index ) const {
		debug_assert( bg_edge_vector_up_to_date_ );
		return bg_edge_vector_[ index ];
	}

	inline
	int get_index_of_adjacent_background_node( int index ) const {
		debug_assert( bg_edge_vector_up_to_date_ );
		return adjacent_bg_node_indices_[ index ];
	}

	inline
	BackgroundNode< V, E, G >* get_adjacent_background_node( int index ) const {
		debug_assert( bg_edge_vector_up_to_date_ );
		return bg_node_vector_[ index ];
	}

	inline
	bool get_bg_edge_vector_up_to_date_() const {
		return bg_edge_vector_up_to_date_;
	}

	void update_bg_edge_vector();

private:
	int num_edges_to_bg_nodes_;
	BackgroundEdgeList bg_edge_list_;
	BackgroundEdgeVector bg_edge_vector_;
	std::vector< int > adjacent_bg_node_indices_;
	BackgroundNodeVector bg_node_vector_;
	bool bg_edge_vector_up_to_date_;

	//no default constructor, uncopyable
	FirstClassNode();
	FirstClassNode( FirstClassNode< V, E, G > const & );
	FirstClassNode< V, E, G > & operator = ( FirstClassNode< V, E, G > const & );
};


//----------------------------------------------------------------------------//
//------------------------- First Class Edge Class ---------------------------//
//----------------------------------------------------------------------------//

///
/// @brief
///
/// @details
///
/// @remarks
/// Defines only a constructor and count_dynamic_memory().
///
template < typename V, typename E, typename G >
class FirstClassEdge : public E {

public:
	virtual ~FirstClassEdge();
	FirstClassEdge( G * owner, int first_node_ind, int second_node_ind );
	virtual unsigned int count_dynamic_memory() const { return E::count_dynamic_memory(); }

private:
	//no default constructor, uncopyable
	FirstClassEdge();
	FirstClassEdge( FirstClassEdge< V, E, G > const & );
	FirstClassEdge< V, E, G > & operator = ( FirstClassEdge< V, E, G > const & );

};


//----------------------------------------------------------------------------//
//------------------- Second Class (Background) Node Class -------------------//
//----------------------------------------------------------------------------//

///
///
/// @brief
/// A node which is not changing type or rotamer throughout a simulation.
///
/// @details
/// In the case of solvent-accessible surface area (SASA) scoring, a background node would be one whose SASA score can
/// change due to neighboring residues being changed.  It itself will not change residue type or rotamer during a
/// simulation, but the scores on this residue can change.
///
/// @remarks
/// Non-instantiable due to pure virtual methods, prepare_for_simulated_annealing(), print(), and count_static_memory().
///
///
template < typename V, typename E, typename G >
class BackgroundNode {

public:
	typedef typename std::list< BackgroundToFirstClassEdge< V, E, G >* >::iterator   BackgroundEdgeListIter;
	typedef typename std::list< BackgroundToFirstClassEdge< V, E, G >* >             BackgroundEdgeList;
	typedef typename std::vector< BackgroundToFirstClassEdge< V, E, G >* >           BackgroundEdgeVector;
	typedef typename std::vector< BackgroundToFirstClassEdge< V, E, G >* >::iterator BackgroundEdgeVectorIter;
	typedef typename std::vector< FirstClassNode< V, E, G >* >                       FirstClassNodeVector;

public:
	virtual ~BackgroundNode();
	BackgroundNode( AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner, int node_index );

	BackgroundEdgeListIter add_edge( BackgroundToFirstClassEdge< V, E, G >* edge_ptr);
	void drop_edge( BackgroundEdgeListIter edge);
	BackgroundToFirstClassEdge< V, E, G >* find_edge( int fc_node_index );

	virtual void prepare_for_simulated_annealing() = 0;
	virtual void print() const = 0;
	virtual unsigned int count_static_memory() const = 0;
	virtual unsigned int count_dynamic_memory() const;

protected:
	void update_edge_vector();

	inline
	int get_node_index() const { return node_index_; }

	inline
	int get_num_incident_edges() const { return num_incident_edges_; }

	inline
	BackgroundToFirstClassEdge< V, E, G >* get_incident_edge( int index ) const {
		debug_assert( edge_vector_up_to_date_ );
		return edge_vector_[ index ];
	}

	inline
	int get_index_of_adjacent_first_class_node( int index ) const {
		debug_assert( edge_vector_up_to_date_ );
		return adjacent_first_class_node_indices_[ index ];
	}

	inline
	FirstClassNode< V, E, G>* get_adjacent_first_class_node( int index ) const {
		debug_assert( edge_vector_up_to_date_ );
		return adjacent_nodes_[ index ];
	}

	inline
	bool get_edge_vector_up_to_date() const { return edge_vector_up_to_date_; }

	inline
	AdditionalBackgroundNodesInteractionGraph< V, E, G >* get_owner() const { return owner_; }


private:
	int node_index_;
	int num_incident_edges_;
	BackgroundEdgeList edge_list_;
	BackgroundEdgeVector edge_vector_;
	std::vector< int > adjacent_first_class_node_indices_;
	FirstClassNodeVector adjacent_nodes_;
	bool edge_vector_up_to_date_;
	AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner_;

	//no default constructor, uncopyable
	BackgroundNode();
	BackgroundNode( BackgroundNode< V, E, G > const & );
	BackgroundNode< V, E, G > & operator = ( BackgroundNode< V, E, G > const & );
};


//----------------------------------------------------------------------------//
//----------------- Second Class To First Class Edge Class -------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// An edge between a background node and a first class node.
///
/// @details
///
/// @remarks
/// Only derived classes can get non-const access to the FirstClassNode and BackgroundNode members.  Const access is
/// provided as a public method.  Non-instantiable due to pure virtual methods.
/// Not sure what the 'pos_in_fc_nodes_edge_list_' and other Iterator returning methods are used for.
///
///
template < typename V, typename E, typename G >
class BackgroundToFirstClassEdge {

public:
	typedef typename std::list< BackgroundToFirstClassEdge< V, E, G >* >::iterator BackgroundEdgeListIter;

public:
	virtual ~BackgroundToFirstClassEdge();
	BackgroundToFirstClassEdge( AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner, int fc_node_index, int bg_node_index);

	inline
	int get_first_class_node_index() const { return first_class_node_index_; }

	inline
	int get_background_node_index() const { return background_node_index_; }

	int get_other_ind( FirstClassNode< V, E, G >* ) const;
	int get_other_ind( BackgroundNode< V, E, G >* ) const;

	FirstClassNode< V, E, G >* get_other_node( BackgroundNode< V, E, G >*  ) const;
	BackgroundNode< V, E, G >* get_other_node( FirstClassNode< V, E, G >*  ) const;

	void set_pos_in_owners_list( BackgroundEdgeListIter );
	void set_pos_in_node_edgevector( FirstClassNode< V, E, G >* caller, int pos );
	void set_pos_in_node_edgevector( BackgroundNode< V, E, G >* caller, int pos);

	bool same_edge( int fc_node_index, int bg_node_index ) const;

	virtual void prepare_for_simulated_annealing() = 0;

	virtual unsigned int count_static_memory() const = 0;
	virtual unsigned int count_dynamic_memory() const;

protected:
	inline
	FirstClassNode< V, E, G >* get_first_class_node() const { return first_class_node_; }

	inline
	BackgroundNode< V, E, G >* get_background_node() const { return background_node_; }

	inline
	AdditionalBackgroundNodesInteractionGraph< V, E, G >* get_owner() const { return owner_; }


private:
	int first_class_node_index_;
	int background_node_index_;

	FirstClassNode< V, E, G >* first_class_node_;
	BackgroundNode< V, E, G >* background_node_;

	int pos_in_fc_edge_vector_;
	int pos_in_bg_edge_vector_;
	BackgroundEdgeListIter pos_in_fc_nodes_edge_list_;
	BackgroundEdgeListIter pos_in_bg_nodes_edge_list_;
	BackgroundEdgeListIter pos_in_owners_edge_list_;

	AdditionalBackgroundNodesInteractionGraph< V, E, G > * owner_;

	//no default constructor, uncopyable
	BackgroundToFirstClassEdge();
	BackgroundToFirstClassEdge(BackgroundToFirstClassEdge< V, E, G > const &);
	BackgroundToFirstClassEdge< V, E, G > & operator = (BackgroundToFirstClassEdge< V, E, G > const & );

};


//----------------------------------------------------------------------------//
//----------- Additional Background Nodes Interaction Graph Class ------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// Function declarations for the AdditionalBackgroundNodesInteractionGraph.
///
template < typename V, typename E, typename G >
class AdditionalBackgroundNodesInteractionGraph : public G {

public:
	typedef typename std::vector< BackgroundNode< V, E, G > * >                           BackgroundNodeVector;
	typedef typename std::list< BackgroundToFirstClassEdge< V, E, G > * >                 BackgroundEdgeList;
	typedef typename std::list< BackgroundToFirstClassEdge< V, E, G > * >::iterator       BackgroundEdgeListIter;
	typedef typename std::list< BackgroundToFirstClassEdge< V, E, G > * >::const_iterator BackgroundEdgeListConstIter;

public:
	AdditionalBackgroundNodesInteractionGraph( int num_nodes );
	virtual ~AdditionalBackgroundNodesInteractionGraph();

	friend class BackgroundToFirstClassEdge< V, E, G >;

	virtual unsigned int getTotalMemoryUsage() const;
	virtual unsigned int count_dynamic_memory() const;

	inline
	int get_num_background_nodes() const { return num_bg_nodes_; }

protected:
	void drop_background_edge( BackgroundEdgeListIter iter);

	// Factory methods for background nodes and edges
	virtual BackgroundNode< V, E, G >* create_background_node( int bg_node_index ) = 0;
	virtual BackgroundToFirstClassEdge< V, E, G >* create_background_edge( int fc_node_index, int bg_node_index ) = 0;

	// call at most once
	void set_num_background_nodes( int num_bg_nodes );

	void add_background_edge( int first_class_node_index, int bg_node_index );

	BackgroundToFirstClassEdge< V, E, G >* find_background_edge( int first_class_node_index, int bg_node_index ) const;

	inline
	int get_num_bg_edges() const { return bg_to_fc_edge_list_.size(); }

	inline
	BackgroundNode< V, E, G >* get_background_node( int index ) const {
		debug_assert( index > 0 && index <= num_bg_nodes_ );
		return bg_nodes_[ index ];
	}

	void prepare_for_simulated_annealing();

	inline
	BackgroundEdgeListConstIter get_bg_edge_list_begin() const {
		return bg_to_fc_edge_list_.begin();
	}

	inline
	BackgroundEdgeListConstIter get_bg_edge_list_end() const {
		return bg_to_fc_edge_list_.end();
	}

private:
	int num_bg_nodes_;
	BackgroundNodeVector bg_nodes_;
	BackgroundEdgeList bg_to_fc_edge_list_;

	mutable BackgroundToFirstClassEdge< V, E, G >* focused_bg_edge_;

	//no default constructor, uncopyable
	AdditionalBackgroundNodesInteractionGraph();
	AdditionalBackgroundNodesInteractionGraph( AdditionalBackgroundNodesInteractionGraph< V, E, G > const & );
	AdditionalBackgroundNodesInteractionGraph< V, E, G > & operator = (AdditionalBackgroundNodesInteractionGraph< V, E, G > const & );
};


//----------------------------------------------------------------------------//
//------------------------- First Class Node Class ---------------------------//
//----------------------------------------------------------------------------//


/// @brief
/// FirstClassNode constructor
///
/// @param
/// owner - [in] - the owning interaction graph
/// node_id - [in] - the index for this node amongst its owners set
/// num_states - [in] - the number of states for this node
///
template < typename V, typename E, typename G >
FirstClassNode< V, E, G > ::FirstClassNode( G* owner, int node_id, int num_states ) :
V( owner, node_id, num_states ),
num_edges_to_bg_nodes_(0),
bg_edge_vector_up_to_date_(false)
{}


/// @brief FirstClassNode destructor
///
template < typename V, typename E, typename G >
FirstClassNode< V, E, G > ::~FirstClassNode()
{}


/// @brief
/// Adds a BackgroundToFirstClassEdge to the node's list of such edges and returns an iterator to that list position.
///
/// @remarks
/// Edge addition invalidates the bg_edge_vector_up_to_date_ flag
///
/// @param
/// edge_ptr - [in] - a pointer to the BackgroundToFirstClassEdge
///
template < typename V, typename E, typename G >
typename std::list< BackgroundToFirstClassEdge< V, E, G >* >::iterator
FirstClassNode< V, E, G >::add_background_edge( BackgroundToFirstClassEdge< V, E, G >* edge_ptr ) {
	++num_edges_to_bg_nodes_;
	bg_edge_vector_up_to_date_ = false;
	return bg_edge_list_.insert( bg_edge_list_.end(), edge_ptr );
}


/// @brief
/// Removes an edge from the node's BackgroundToFirstClassEdge list
///
/// @remarks
/// Invalidates bg_edge_vector_up_to_date_ flag
///
/// @param
/// edge - [in] - the iterator to the edge; in order to guarantee constant
///   time edge deletion, each edge keeps iterators to its position in the
///   edge lists of the nodes it's incident upon.  It hands these iterators
///   back to the nodes when it wants to delete itself.
///
template < typename V, typename E, typename G >
void FirstClassNode< V, E, G >::drop_background_edge( BackgroundEdgeListIter edge ) {
	--num_edges_to_bg_nodes_;
	bg_edge_list_.erase( edge );
	bg_edge_vector_up_to_date_ = false;
}


/// @brief
/// Returns an int representing the amount of memory in bytes used by this node
///
/// @remarks
/// ronj lists consume more memory because they maintain pointers to the previous and next nodes
/// ronj current calculation uses 4 * # elements in list * sizeof(pointer), but I think this is wrong
/// ronj isn't it only 3 * # elements in list * sizeof(pointer)?
///
/// @param
/// total_memory - [out] - the amount of memory used by this node
///
template < typename V, typename E, typename G >
unsigned int FirstClassNode< V, E, G >::count_dynamic_memory() const {

	//ronj calculate how much dynamic memory the base class is using
	unsigned int total_memory = V::count_dynamic_memory();

	total_memory += 4 * bg_edge_list_.size() * sizeof( BackgroundToFirstClassEdge< V, E, G >* );
	total_memory += bg_edge_vector_.size() * sizeof( BackgroundToFirstClassEdge< V, E, G >* );
	total_memory += adjacent_bg_node_indices_.size() * sizeof( int );
	total_memory += bg_node_vector_.size() * sizeof( BackgroundNode< V, E, G >* );

	return total_memory;
}


/// @brief
/// Syncs the background edge vector with the background edge list
///
/// @details
/// Updates the vector of pointers to background edges to correspond to the list of background edges. Informs its
/// incident edges of their position in its vector.
///
/// @remarks
/// We are resizing these vectors to num_edges_to_bg_nodes PLUS ONE. This makes the vector 1-based.
///
template < typename V, typename E, typename G >
void FirstClassNode< V, E, G >::update_bg_edge_vector() {

	bg_edge_vector_.resize( num_edges_to_bg_nodes_ + 1 );
	adjacent_bg_node_indices_.resize( num_edges_to_bg_nodes_ + 1 );
	bg_node_vector_.resize( num_edges_to_bg_nodes_ + 1 );

	BackgroundEdgeVectorIter position1 = bg_edge_vector_.begin();
	++position1;

	std::copy( bg_edge_list_.begin(), bg_edge_list_.end(), position1 );

	for ( int ii = 1; ii <= num_edges_to_bg_nodes_; ++ii ) {

		//ronj set_pos_in_node_edgevector updates the edge vector in this FirstClassNode
		//ronj this FirstClassNode is telling the BackgroundToFirstClass edge connected to it what position in the
		//ronj edge vector (stored in this class) that edge is in. Hey, edge, you're fourth in my vector of incident edges
		bg_edge_vector_[ii]->set_pos_in_node_edgevector( this, ii );

		//ronj get_other_ind is passed 'this' FirstClassNode, so it will return the index of the BackgroundNode that has
		//ronj an edge with this FirstClassNode
		adjacent_bg_node_indices_[ii] = bg_edge_vector_[ii]->get_other_ind( this );
		bg_node_vector_[ii] = bg_edge_vector_[ii]->get_other_node( this );
	}

	bg_edge_vector_up_to_date_ = true;
	return;
}


//----------------------------------------------------------------------------//
//------------------------- First Class Edge Class ---------------------------//
//----------------------------------------------------------------------------//


/// @brief
/// FirstClassEdge constructor
///
/// @details
/// This class is almost empty; a FirstClassEdge connects two first class vertices, and is ignorant of the
/// existence of both second class nodes (synonymous with BackgroundNode) and second-to-first-class edges
/// (synonymous with BackgroundToFirstClassEdge).  I (apl) cannot think of any data that a FirstClassEdge should
/// hold that is not held in a PrecomputedPairEnergiesEdge.  Nor can I think of any methods that the edge should
/// implement.  However, I include the class in the inheritance hierarchy just in case I come up with something
/// later.
///
/// @param
/// owner - [in] - the owning interaction graph
/// first_node_index - [in] - the index of the lower-indexed first class node
/// second_node_index - [in] - the index of the higher-indexed first class node
///
template < typename V, typename E, typename G >
FirstClassEdge< V, E, G >::FirstClassEdge( G* owner, int first_node_ind, int second_node_ind) :
E ( owner, first_node_ind, second_node_ind ) {}


/// @brief
/// FirstClassEdge destructor
///
template < typename V, typename E, typename G >
FirstClassEdge< V, E, G >::~FirstClassEdge() {}


//----------------------------------------------------------------------------//
//------------------- Second Class (Background) Node Class -------------------//
//----------------------------------------------------------------------------//


/// @brief
/// BackgroundNode constructor - no default or copy constructors; no operator =
///
/// @details
/// I use "background node" synonymously with "second class node".
/// This kind of node contributes to the energy of the system.  It does not
/// have any assignable states -- the energy it contributes is a function of
/// the states of the first class vertices that surround it.
/// A simulated annealer is unaware of the presence of this kind of node;
/// the owning interaction graph obscures it from view of the annealer.
///
/// @param
/// owner - [in] - the owning graph
/// node_index - [in] - the index of the background node
///
template < typename V, typename E, typename G >
BackgroundNode< V, E, G >::BackgroundNode( AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner, int node_index ) :
node_index_( node_index ),
num_incident_edges_( 0 ),
edge_vector_up_to_date_( false ),
owner_( owner ) {
}


/// @brief
/// BackgroundNode destructor
///
template < typename V, typename E, typename G >
BackgroundNode< V, E, G >::~BackgroundNode() {}


/// @brief
/// Adds a BackgroundToFirstClassEdge to the edge list for this node.  Returns an iterator to the new list element.
///
/// @remarks
/// Invalidates edge_vector_up_to_date_ flag
///
/// @param
/// edge_ptr - [in] - pointer to the edge being added
///
template < typename V, typename E, typename G >
typename std::list< BackgroundToFirstClassEdge< V, E, G >* >::iterator
BackgroundNode< V, E, G > ::add_edge( BackgroundToFirstClassEdge< V, E, G > * edge_ptr ) {
	++num_incident_edges_;
	edge_vector_up_to_date_ = false;
	return edge_list_.insert( edge_list_.end(),  edge_ptr );
}


/// @brief
/// removes an edge from the nodes edge list in constant time
///
/// @remarks
/// Invalidates edge_vector_up_to_date_ flag
///
/// @param
/// edge - [in] - the iterator to the edge's position in node's edge list
///
template < typename V, typename E, typename G >
void BackgroundNode< V, E, G >::drop_edge( BackgroundEdgeListIter edge ) {
	--num_incident_edges_;
	edge_vector_up_to_date_ = false;
	edge_list_.erase( edge );
}


/// @brief
/// Syncs the edge vector with the edge list
///
/// @details
/// Updates the vector of pointers to background edges to correspond to the list of background edges. Informs its
/// incident edges of their position in its vector.
///
/// @remarks
/// sets edge_vector_up_to_date_ flag to true
///
template < typename V, typename E, typename G >
void BackgroundNode< V, E, G >::update_edge_vector() {

	edge_vector_.resize( num_incident_edges_ + 1);
	adjacent_first_class_node_indices_.resize( num_incident_edges_ + 1);
	adjacent_nodes_.resize(num_incident_edges_ + 1);

	BackgroundEdgeVectorIter position1 = edge_vector_.begin();
	++position1;

	std::copy( edge_list_.begin(), edge_list_.end(), position1 );

	for ( int ii = 1; ii <= num_incident_edges_; ++ii ) {
		//ronj set_pos_in_node_edgevector updates the edge vector in this BackgroundNode.
		//ronj this BackgroundNode is telling the BackgroundToFirstClass edge connected to it what position in the
		//ronj edge vector (stored in this class) that edge is in. Hey, edge, you're fourth in my vector of incident edges.
		edge_vector_[ii]->set_pos_in_node_edgevector( this, ii );

		adjacent_first_class_node_indices_[ ii ] = edge_vector_[ii]->get_other_ind( this );
		adjacent_nodes_[ ii ] = edge_vector_[ii]->get_other_node( this );
	}

	edge_vector_up_to_date_ = true;
}


/// @brief
/// Linear time edge lookup function
///
/// @param
/// fc_node_index - [in] - the index of the first class node that is on the other end of the sought edge
///
template < typename V, typename E, typename G >
BackgroundToFirstClassEdge< V, E, G >* BackgroundNode< V, E, G >::find_edge( int fc_node_index ) {

	for ( BackgroundEdgeListIter iter = edge_list_.begin(); iter != edge_list_.end(); ++iter ) {
		if ( (*iter)->same_edge( fc_node_index, node_index_ ) ) {
			return *iter;
		}
	}
	return 0;
}


/// @brief
/// Returns an int representing the amount of memory in bytes used by this node
///
/// @remarks
/// ronj lists consume more memory because they maintain pointers to the previous and next nodes
/// ronj current calculation uses 4 * # elements in list * sizeof(pointer), but I think this is wrong
/// ronj isn't it only 3 * # elements in list * sizeof(pointer)?
///
/// @param
/// total_memory - [out] - the amount of memory used by this node
///
template < typename V, typename E, typename G >
unsigned int BackgroundNode< V, E, G >::count_dynamic_memory() const {

	unsigned int total_memory = 0;
	total_memory += adjacent_first_class_node_indices_.size() * sizeof( int );
	total_memory += adjacent_nodes_.size() * sizeof( FirstClassNode< V, E, G >* );
	total_memory += 4 * edge_list_.size() * sizeof( BackgroundToFirstClassEdge< V, E, G >* );
	total_memory += edge_vector_.size() * sizeof( BackgroundToFirstClassEdge< V, E, G >* );

	return total_memory;
}


//----------------------------------------------------------------------------//
//----------------- Second Class To First Class Edge Class -------------------//
//----------------------------------------------------------------------------//


/// @brief
/// BackgroundToFirstClassEdge constructor - no default or copy constructors; no operator =.
///
/// @details
/// This class of edge connects first class and second class nodes.  This class of edge is very unlikely to hold the
/// same kind of data in a (concrete) derived class.  e.g. The SASABackgroundEdge holds pre-computed rotamer sphere
/// overlaps for each rotamer-on-a-SASANode and the single-(fixed)-rotamer-on-a-background residue.
/// The first class edge does not attempt to store precomputed sphere overlaps as they would be prohibitively expensive.
///
/// @param
/// owner - [in] - the owning graph
/// fc_node_index - [in] - the index of the first class node
/// bg_node_index - [in] - the index of the second class node
///
template < typename V, typename E, typename G >
BackgroundToFirstClassEdge< V, E, G >::BackgroundToFirstClassEdge
( AdditionalBackgroundNodesInteractionGraph< V, E, G > * owner, int fc_node_index, int bg_node_index ) :
first_class_node_index_( fc_node_index ),
background_node_index_( bg_node_index),
owner_( owner )
{
	first_class_node_ = ( FirstClassNode< V, E, G >* )(owner_->get_node( first_class_node_index_ ));
	background_node_ = owner_->get_background_node( background_node_index_);

	pos_in_fc_nodes_edge_list_ = first_class_node_->add_background_edge(this);
	pos_in_bg_nodes_edge_list_ = background_node_->add_edge(this);
}


/// @brief
/// virtual destructor.  The edge removes itself from the graph by informing the two vertices its incident upon
/// to drop it from their edge lists.
/// Constant time edge removal.
///
template < typename V, typename E, typename G >
BackgroundToFirstClassEdge< V, E, G >::~BackgroundToFirstClassEdge() {
	first_class_node_->drop_background_edge( pos_in_fc_nodes_edge_list_ );
	background_node_->drop_edge( pos_in_bg_nodes_edge_list_ );
	owner_->drop_background_edge( pos_in_owners_edge_list_ );
}


/// @brief
/// returns the index of the second class node
///
/// @details
/// A first class vertex may request a BackgroundToFirstClassEdge for the index of the background node the edge connects
/// it to by invoking edge->get_other_ind( this );  The this pointer simply tells the compiler which of the two
/// overloaded get_other_ind() methods to invoke.
///
/// @param
/// (unnamed parameter) - [in] - First Class nodes want information about second class nodes when they refer to
/// "the other node".  The compiler resolves which of the two overloaded methods to invoke.
///
template < typename V, typename E, typename G >
int BackgroundToFirstClassEdge< V, E, G >::get_other_ind( FirstClassNode< V, E, G >* ) const {
	return background_node_index_;
}


/// @brief
/// returns the index of the first class node
///
/// @details
/// A second class vertex may request a BackgroundToFirstClassEdge for the index of the first class node the edge connects
/// it to by invoking edge->get_other_ind( this );  The this pointer simply tells the compiler which of the two
/// overloaded get_other_ind() methods to invoke.
///
/// @param
/// (unnamed parameter) - [in] - second class nodes want information about first class nodes when they refer to
/// "the other node".  The compiler resolves which of the two overloaded methods to invoke.
///
template < typename V, typename E, typename G >
int BackgroundToFirstClassEdge< V, E, G >::get_other_ind(  BackgroundNode< V, E, G > * ) const {
	return first_class_node_index_;
}


/// @brief
/// returns a pointer to the second class node
///
/// @details
/// unnamed parameter is usually the this pointer where the method is invoked
/// inside a method of a first class node.  the pointer type must be resolved
/// at compile time.  Same ideas here as in get_other_ind()
///
template < typename V, typename E, typename G >
BackgroundNode< V, E, G >* BackgroundToFirstClassEdge< V, E, G >::get_other_node( FirstClassNode< V, E, G >* ) const {
	return background_node_;
}


/// @brief
/// returns a pointer to the first class node
///
/// @details
/// unnamed parameter is usually the this pointer where the method is invoked
/// inside a method of a second class node.  the pointer type must be resolved
/// at compile time.  Same ideas here as in get_other_ind()
///
template < typename V, typename E, typename G >
FirstClassNode< V, E, G >* BackgroundToFirstClassEdge< V, E, G > ::get_other_node( BackgroundNode< V, E, G >* ) const {
	return first_class_node_;
}


/// @brief
/// stores the iterator to this edge in the owning graph's list of background-to-first-class edges.
///
/// @details
/// required for constant time edge deletion
///
template < typename V, typename E, typename G >
void
BackgroundToFirstClassEdge< V, E, G >::set_pos_in_owners_list( BackgroundEdgeListIter iter ) {
	pos_in_owners_edge_list_ = iter;
}


/// @brief
/// stores the index of this edge in its first class node's edge vector
///
/// @details
/// again, first (unnamed) parameter is the 'this' pointer where the method has
/// been invoked inside a method of the first class node.
///
/// @param
/// (unnamed parameter) - [in] - pointer identifying the type of the node that invoked this overloaded method.
/// pos - [in] - the position of 'this' edge in the edge vector of the first class node
///
template < typename V, typename E, typename G >
void BackgroundToFirstClassEdge< V, E, G >::set_pos_in_node_edgevector( FirstClassNode< V, E, G >*, int pos ) {
	pos_in_fc_edge_vector_ = pos;
}


/// @brief
/// stores the index of this edge in its second class node's edge vector
///
/// @details
/// again, first (unnamed) parameter is the this pointer where the method has
/// been invoked inside a method of the second class node.
///
/// @param
/// (unnamed parameter) - [in] - pointer identifying the type of the node that invoked this overloaded method.
/// pos - [in] - the position of this node in the edge vector of the second class node
///
template < typename V, typename E, typename G >
void BackgroundToFirstClassEdge< V, E, G >::set_pos_in_node_edgevector( BackgroundNode< V, E, G >*, int pos ) {
	pos_in_bg_edge_vector_ = pos;
}


/// @brief
/// returns true if this node is incident upon the two input vertex indices
///
/// @param
/// fc_node_index - [in] - the index of the first class node
/// bg_node_index - [in] - the index of the second class node
///
/// @remarks
/// This graph setup will not support multigraphs: graphs that include multiple edges to the same vertices.
/// In my (apl) original dynamic programming algorithm I didn't want to support multigraphs either; but in the
/// adaptive dynamic programming implementation, multiple copies of "the same" edge ended up saving time -- though
/// maybe not memory.  So, this framework here will not readily support adaptive dynamic programming; but that's ok.
///
template < typename V, typename E, typename G >
bool BackgroundToFirstClassEdge< V, E, G >::same_edge( int fc_node_index, int bg_node_index ) const {
	return ( (first_class_node_index_ == fc_node_index) && (background_node_index_ == bg_node_index) );
}


/// @brief
/// Returns an int representing the amount of memory in bytes used by this node
///
/// @remarks
/// return 0 because no memory is dynamically allocated for these edges.
///
template < typename V, typename E, typename G >
unsigned int BackgroundToFirstClassEdge< V, E, G >::count_dynamic_memory() const {
	return 0;
}


//----------------------------------------------------------------------------//
//----------- Additional Background Nodes Interaction Graph Class ------------//
//----------------------------------------------------------------------------//


/// @brief
/// AdditionalBackgroundNodesInteractionGraph constructor; no default or copy constructors; no operator =
///
/// @param
/// num_nodes - [in] - the number of first class nodes
///
template < typename V, typename E, typename G >
AdditionalBackgroundNodesInteractionGraph< V, E, G >::AdditionalBackgroundNodesInteractionGraph( int num_nodes ) :
G( num_nodes ),
num_bg_nodes_( -1 ),
focused_bg_edge_( 0 )
{}


/// @brief
/// AdditionalBackgroundNodesInteractionGraph destructor
///
/// @details
/// deallocates each BackgroundToFirstClassEdge. Then deallocates each background
/// node.  Order is important.  The edges assume that its vertices still exist
/// at the time it is removed.  This destructor enforces that.
///
template < typename V, typename E, typename G >
AdditionalBackgroundNodesInteractionGraph< V, E, G >::~AdditionalBackgroundNodesInteractionGraph() {

	for ( BackgroundEdgeListIter iter = bg_to_fc_edge_list_.begin(); iter != bg_to_fc_edge_list_.end(); ) { //note: no increment statement here
		BackgroundEdgeListIter next_iter = iter;
		next_iter++;
		//edges delete themselves from this list, invalidating iterators, so
		//get the next iterator before calling delete
		delete (*iter);
		iter = next_iter;
	}

	for ( int ii = 1; ii <= num_bg_nodes_; ++ii ) {
		delete bg_nodes_[ ii ];
	}
}


/// @brief
/// Constant time edge removal.
///
/// @details
///
/// @param
/// iter - [in] - the iterator to the position in the graph's edge list for the edge that is removing itself.
///
template < typename V, typename E, typename G >
void AdditionalBackgroundNodesInteractionGraph< V, E, G > ::drop_background_edge( BackgroundEdgeListIter iter ) {
	bg_to_fc_edge_list_.erase( iter );
}


/// @brief
/// sets the number of background nodes in the graph.  Should be called no more than once.
/// Some problem instances do not require background nodes.
///
/// @details
/// Allocates the background nodes using a factory method.  Their indices start from 1. The create_background_node()
/// method is virtual.
///
template < typename V, typename E, typename G >
void AdditionalBackgroundNodesInteractionGraph< V, E, G >::set_num_background_nodes( int num_bg_nodes ) {

	debug_assert( num_bg_nodes_ == -1 && num_bg_nodes >= 0); //call this method at most once
	num_bg_nodes_ = num_bg_nodes;

	if ( num_bg_nodes_ == 0 ) {
		return;
	}

	bg_nodes_.resize( num_bg_nodes_ + 1 );
	for ( int ii = 1; ii <= num_bg_nodes_; ++ii ) {
		bg_nodes_[ii] = create_background_node(ii);
	}
}


/// @brief
/// adds a BackgroundToFirstClassEdge to the graph and performs the requisite bookkeepking.
///
/// @param
/// first_class_node_index - [in] - the index of the first class node
/// bg_node_index - [in] - the index of the second class node
///
template < typename V, typename E, typename G >
void AdditionalBackgroundNodesInteractionGraph< V, E, G >::add_background_edge( int first_class_node_index, int bg_node_index ) {

	BackgroundToFirstClassEdge< V, E, G >* new_edge = create_background_edge( first_class_node_index, bg_node_index );
	bg_to_fc_edge_list_.push_front( new_edge );
	new_edge->set_pos_in_owners_list( bg_to_fc_edge_list_.begin() );
	focused_bg_edge_ = new_edge;
	return;
}


/// @brief
/// returns a pointer to the background edge, identified by the indices of the first and second class nodes it connects.
/// Returns the null pointer if the edge does not exist.  Possibly linear time operation (in the number of
/// BackgroundToFirstClassEdges).
///
/// @details
/// The graph keeps a pointer to a "focused edge" so that repeated requests for the same edge take constant time.
///
/// @param
/// first_class_node_index - [in] - the index of the first class node
/// bg_node_index - [in] - the index of the second class node
///
template < typename V, typename E, typename G >
BackgroundToFirstClassEdge< V, E, G >*
AdditionalBackgroundNodesInteractionGraph< V, E, G >::find_background_edge( int first_class_node_index, int bg_node_index ) const {

	if ( (focused_bg_edge_ != 0) && !( focused_bg_edge_->same_edge(first_class_node_index, bg_node_index) ) ) {
		focused_bg_edge_ = bg_nodes_[ bg_node_index ]->find_edge( first_class_node_index );
	}
	return focused_bg_edge_;
}


/// @brief
/// invokes prepare_for_simulated_annealing on each of the BackgroundToFirstClassEdges and then invokes
/// prepare_for_simulated_annealing on each of the BackgroundNodes.
///
/// @details
/// A BackgroundToFirstClassEdges may decide to delete itself before sim annealing begins. Since the second class node
/// will likely update is edge vector in its (virtual) call to prepare_for_sim_annealing, any edges that should be
/// should be removed before the second class nodes update their edge vectors.
///
/// @remarks
/// does not invoke InteractionGraphBase::prepare_for_simulated_annealing()
///
template < typename V, typename E, typename G >
void AdditionalBackgroundNodesInteractionGraph< V, E, G >::prepare_for_simulated_annealing() {

	for ( BackgroundEdgeListIter iter = bg_to_fc_edge_list_.begin(); iter != bg_to_fc_edge_list_.end(); ) { //note: no increment statement here

		BackgroundEdgeListIter next_iter = iter;
		next_iter++;

		// edges sometimes delete themselves, invalidating iterators, so
		// get the next iterator before calling prepare_for_simulated_annealing

		// BackgroundToFirstClassEdge::prepare_for_simulated_annealing() is a pure virtual
		(*iter)->prepare_for_simulated_annealing();
		iter = next_iter;
	}

	for ( int ii = 1; ii <= num_bg_nodes_; ++ii ) {
		// BackgroundNode::prepare_for_simulated_annealing() is a pure virtual
		bg_nodes_[ ii ]->prepare_for_simulated_annealing();
	}

	return;
}


/// @brief
/// Returns an int representing the total amount of memory in bytes used by this graph
///
template < typename V, typename E, typename G >
unsigned int AdditionalBackgroundNodesInteractionGraph< V, E, G >::getTotalMemoryUsage() const {

	unsigned int total_memory = G::getTotalMemoryUsage();

	for ( int ii = 1; ii <= get_num_background_nodes(); ++ii ) {
		total_memory += bg_nodes_[ii]->count_dynamic_memory();
		total_memory += bg_nodes_[ii]->count_static_memory();
	}

	for ( BackgroundEdgeListConstIter iter = bg_to_fc_edge_list_.begin(); iter != bg_to_fc_edge_list_.end(); ++iter ) {
		total_memory += (*iter)->count_dynamic_memory();
		total_memory += (*iter)->count_static_memory();
	}

	return total_memory;
}


/// @brief
/// Returns an int representing the amount of dynamic memory in bytes used by this graph
///
template < typename V, typename E, typename G >
unsigned int AdditionalBackgroundNodesInteractionGraph< V, E, G >::count_dynamic_memory() const {

	unsigned int total_memory = G::count_dynamic_memory();

	total_memory += bg_nodes_.size() * sizeof( BackgroundNode< V, E, G >* );
	total_memory += 4 * bg_to_fc_edge_list_.size() * sizeof( BackgroundToFirstClassEdge< V, E, G >* );

	return total_memory;
}


}
}
} //end namespace pack

#endif

