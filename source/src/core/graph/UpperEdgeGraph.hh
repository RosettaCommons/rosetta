// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/graph/UpperEdgeGraph.hh
/// @brief  templated graph for fast edge additions
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_graph_UpperEdgeGraph_hh
#define INCLUDED_core_graph_UpperEdgeGraph_hh

// Unit Headers
#include <core/graph/UpperEdgeGraph.fwd.hh>

// Package headers
#include <platform/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.fwd.hh>


namespace core {
namespace graph {

// The point of this graph is to provide fast edge additions: new() is somewhat slow, whereas
// push_back with a vector that has reserved space for its edges ( through the method reserve())
// doesn't require new.

// Vertices contain a data member of class V so that their size is known at compile time and so that data can be stored on vertices.
// Edges contain a data member of class E so that their size is known at compile time and so that data can be stored on edges

// e.g. class EnergyEdgeData { public: EnergyMap emap; }; could be held inside a UEEdge object

// This class does not allow tight control over data integrity in the V and E classes. For instance, you can imagine
// wanting the vertex objects to update the sum of their two body energies if you ever changed the two body energy of
// one of their edges.  If you were to derive from class Graph (as scoring::EnergyGraph does) then you could have that graph
// enforce that bit of data integrity (EnergyGraph does not do that, currently, but it could).
// However, that is not an option with this graph class.

// The other major weakness is that each vertex is unaware of who its lower neighbors are; its not possible to
// iterate over all of the edges incident on a vertex.  This is also its major strength -- the cost of
// maintaining upper and lower edge lists is much greater than the cost of maintaining only upper edge lists.
// Each vertex is aware of how many lower edges it has.

// This graph does not (yet) support loops (edges incident on a single vertex)

// If the edge data is quite large, this graph will occupy more space than class Graph would, which stores edge data
// for only those edges that actually exist.  (With reserve(), vertices allocate space for edges that do not exist).

template < class V, class E >
class UEVertex : public utility::pointer::ReferenceCount
{
public:
	typedef UEVertex< V, E >                VertexClass;
 	typedef UEEdge< V, E >                  EdgeClass;
	typedef UpperEdgeGraph< V, E >          GraphClass;
	typedef typename utility::vector1< EdgeClass >   UpperEdgeVector;
	typedef typename utility::vector1< EdgeClass >::iterator       UpperEdgeListIter;
	typedef typename utility::vector1< EdgeClass >::const_iterator UpperEdgeListConstIter;

public:

	/// @brief standard constructor; though if vertices are to be allocated
	/// in a single vector1, then this ctor requires creating and destroying
	/// an anonymous vertex to add each vertex to the graph... wasteful
	UEVertex( GraphClass * owner, platform::Size index ) :
		index_( index ),
		upper_edges_( V::NUM_EDGES_TO_RESERVE ),
		num_upper_edges_( 0 ),
		num_lower_edges_( 0 ),
		lazily_deleted_edges_present_( false ),
		owner_( owner )
	{
		upper_edges_.resize( 0 );
	}

	/// @brief default constructor; useful for allocating all vertices at once
	/// in a single vector1 using the resize() method.  Note that both
	/// index_ and owner_ are set to 0.  This state must be corrected before the vertex
	/// is of any use.
	UEVertex() :
		index_( 0 ),
		upper_edges_( V::NUM_EDGES_TO_RESERVE ),
		num_upper_edges_( 0 ),
		num_lower_edges_( 0 ),
		lazily_deleted_edges_present_( false ),
		owner_( 0 )
	{
		upper_edges_.reserve( 0 );
	}

	/// @brief dstor, should be non-virtual
	~UEVertex()
	{
	}

private:

	/// @brief method for completing the construction of a vertex if the defaul ctor was used.
	///
	/// @detailed to be called by class UEGraph only
	void
	set_index_and_owner(
		platform::Size index,
		GraphClass * owner
	)
	{
		assert( index_ == 0 && owner_ == 0 ); // this can be called at most once
		index_ = index;
		owner_ = owner;
	}

public:

	/// @brief reserves space for the edges.  Due to the nature of vertex construction
	/// in the UEGraph, this method should be called separately from the constructor.
	///
	/// @details class V must define a static member "NUM_EDGES_TO_RESERVE()"; the
	/// UEVertex constructor allocates space for that many edges -- if fewer than
	/// that many edges are added to a sigle node, then new() need never be called.
	/// if more are added, the class does not crash -- it simply grows the vector out
	void
	reserve_edge_space() {
		upper_edges_.reserve( V::NUM_EDGES_TO_RESERVE );
	}

	EdgeClass * get_edge( platform::Size upper_vertex_id )
	{
		make_edge_vector_current();
		for ( platform::Size ii = 1; ii <= upper_edges_.size(); ++ii )
		{
			if ( upper_edges_[ ii ].upper_vertex() == upper_vertex_id )
			{
				return & upper_edges_[ ii ];
			}
		}
		return 0;
	}

	EdgeClass const * get_edge( platform::Size upper_vertex_id ) const
	{
		//pretend this is a non-const vertex
		return const_cast< VertexClass * > (this)->get_edge( upper_vertex_id );
	}

	bool edge_exists( platform::Size upper_vertex_id )
	{
		make_edge_vector_current();
		for ( platform::Size ii = 1; ii <= upper_edges_.size(); ++ii )
		{
			if ( upper_edges_[ ii ].upper_vertex() == upper_vertex_id )
			{
				return true;
			}
		}
		return false;
	}

	UpperEdgeListIter
	upper_edge_list_begin()
	{
		make_edge_vector_current();
		return upper_edges_.begin();
	}

	UpperEdgeListIter
	upper_edge_list_end()
	{
		make_edge_vector_current();
		return upper_edges_.end();
	}

	UpperEdgeListConstIter
	const_upper_edge_list_begin() const
	{
		make_edge_vector_current();
		return upper_edges_.begin();
	}

	UpperEdgeListConstIter
	const_upper_edge_list_end() const
	{
		make_edge_vector_current();
		return upper_edges_.end();
	}

	platform::Size num_upper_neighbors() const { return num_upper_edges_; }
	platform::Size num_lower_neighbors() const { return num_lower_edges_; }
	platform::Size num_neighbors() const { return num_upper_edges_ + num_lower_edges_; }
	platform::Size num_neighbors_counting_self() const { return num_upper_edges_ + num_lower_edges_ + 1; }

	void drop_all_edges() {
		num_upper_edges_ = 0;
		num_lower_edges_ = 0;
		lazily_deleted_edges_present_ = false;
		upper_edges_.resize( 0 ); // does not deallocate space
	}

	// accessors for data held by this vertex
	V & data() { return data_; }
	V const & data() const { return data_; }

	friend class UEEdge< V, E >;
	friend class UpperEdgeGraph< V, E >;

private:

	/// @brief add an edge
	///
	/// @detailed do not try to use this method or the other add_edge method to add edges to
	/// a vertex; they are for class UpperEdgeGraph only.  Invoke add_edge on the graph itself.
	void
	add_edge( platform::Size upper_vertex_index ) // called by UpperEdgeGraph
	{
		//if ( upper_edges_.capacity() == upper_edges_.size() ) {
		//	std::cout << this << " " << index_  << " About to resize edge vector on vertex " <<  index_ << " adding edge to " << upper_vertex_index << " upper_edges_.size() " << upper_edges_.size()  << std::endl;
		//}
		upper_edges_.push_back( EdgeClass( owner_, index_, upper_vertex_index ) );
		add_edge_common();
	}

	void
	add_edge( platform::Size upper_vertex_index, E const & edge_data ) // called by UpperEdgeGraph
	{
		//if ( upper_edges_.capacity() == upper_edges_.size() ) {
		//	std::cout << this << " " << index_ << " About to resize edge vector on vertex " <<  index_ << " adding edge to " << upper_vertex_index << " upper_edges_.size() " << upper_edges_.size()  << std::endl;
		//}
		upper_edges_.push_back( EdgeClass( owner_, index_, upper_vertex_index, edge_data ) );
		add_edge_common();
	}


	// called by UEEdge
	void note_upper_edge_deleted()
	{
		lazily_deleted_edges_present_ = true;
		--num_upper_edges_;
	}

	// called by UEEdge
	void note_lower_edge_deleted() { --num_lower_edges_; }

	// called by UEEdge
 	void note_lower_edge_added() { ++num_lower_edges_; }


	void make_edge_vector_current() const
	{
		if ( ! lazily_deleted_edges_present_ ) return;

		platform::Size shift( 0 );
		for ( platform::Size ii = 1; ii <= upper_edges_.size(); ++ii )
		{
			if ( upper_edges_[ii].deleted() )
			{
				++shift;
			}
			else
			{
				if ( shift > 0 )
				{
					upper_edges_[ ii - shift ] = upper_edges_[ ii ];
				}
			}

			assert ( ii != upper_edges_.size() || (ii - shift == num_upper_edges_) );
		}
		upper_edges_.resize( num_upper_edges_ );
		lazily_deleted_edges_present_ = false;
	}

	inline
	void
	add_edge_common()
	{
		++num_upper_edges_;
	}

private:
	// Data

	platform::Size index_;
	mutable UpperEdgeVector upper_edges_;

	platform::Size num_upper_edges_;
	platform::Size num_lower_edges_; // each vertex knows how many lower edges it has; it just doesn't know who they are

	mutable bool lazily_deleted_edges_present_;

	GraphClass * owner_;
	V data_;


};

template < class V, class E >
class UEEdge
{
public:
	typedef UpperEdgeGraph< V, E > GraphClass;
	typedef UEVertex< V, E >       VertexClass;

public:

	UEEdge() :
		owner_( 0 ),
		lower_vertex_( 0 ),
		upper_vertex_( 0 ),
		deleted_( true ),
		upper_vertex_index_( 0 )
	{}

	UEEdge( GraphClass * owner, int lower_node, int upper_node )
	:
		owner_( owner ),
		lower_vertex_( owner->get_vertex_ptr( lower_node )),
		upper_vertex_( owner->get_vertex_ptr( upper_node )),
		deleted_( false ),
		upper_vertex_index_( upper_node )
	{
		upper_vertex_->note_lower_edge_added();
	}

	~UEEdge()
	{
		//std::cout << "UEEdge dstor" << std::endl;
	}


	UEEdge( GraphClass * owner, int lower_node, int upper_node, E const & data )
	:
		owner_( owner ),
		lower_vertex_( owner->get_vertex_ptr( lower_node )),
		upper_vertex_( owner->get_vertex_ptr( upper_node )),
		deleted_( false ),
		upper_vertex_index_( upper_node ),
		data_( data )
	{
		upper_vertex_->note_lower_edge_added();
	}


	platform::Size upper_vertex() const { return upper_vertex_index_; }

	void delete_edge() {
		if ( deleted_ ) return;
		deleted_ = true;
		lower_vertex_->note_upper_edge_deleted();
		upper_vertex_->note_lower_edge_deleted();
		owner_->note_edge_deleted();
	}

	bool deleted() const { return deleted_; } // vertex needs read info, no one else should

	// accessors for data held by this edge
	E & data() { return data_; }
	E const & data() const { return data_; }

	friend class UpperEdgeGraph< V, E >;

private:

	GraphClass * owner_;
	VertexClass * lower_vertex_;
	VertexClass * upper_vertex_;

	bool deleted_;
	platform::Size upper_vertex_index_;

	E data_;

};

// Copy ctor and operator = not yet implemented
// do not use!
template < class V, class E >
class UpperEdgeGraph : public utility::pointer::ReferenceCount
{
public:
	typedef UEVertex< V, E >                         VertexClass;
	typedef UEEdge< V, E >                           EdgeClass;
	typedef utility::vector1< utility::pointer::owning_ptr< VertexClass > > VertexVector;
	typedef typename utility::vector1< UEEdge< V, E > >::iterator      UpperEdgeListIter;
	typedef typename utility::vector1< UEEdge< V, E > >::const_iterator UpperEdgeListConstIter;


public:

	UpperEdgeGraph() : num_vertices_( 0 ), num_edges_( 0 ), vertices_( 0 ) {}
	UpperEdgeGraph( platform::Size nverts ) : num_vertices_( nverts ), num_edges_( 0 ), vertices_( nverts, 0 )
	{
		create_vertices();
	}

	UpperEdgeGraph( UpperEdgeGraph< V, E > const & other ) :
		utility::pointer::ReferenceCount()
	{
		copy_from( other );
	}

	virtual ~UpperEdgeGraph() {};

	UpperEdgeGraph< V, E > const &
	operator = ( UpperEdgeGraph< V, E > const & other )
	{
		copy_from( other );
		return *this;
	}

	//clears all edge data
	void
	set_num_vertices( platform::Size num_vertices )
	{
		if ( num_vertices_ != num_vertices ) {
			num_vertices_ = num_vertices;
			create_vertices();
		} else {
			drop_all_edges();
		}
	}

	VertexClass &
	get_vertex( platform::Size index ) { return *vertices_[ index ]; }

	VertexClass const &
	get_vertex( platform::Size index ) const { return *vertices_[ index ]; }

	// add an edge
	void add_edge( platform::Size lower_vertex, platform::Size upper_vertex )
	{
		//assert( lower_vertex < upper_vertex );
		assert( ! edge_exists( lower_vertex, upper_vertex ));
		assert( ! edge_exists( upper_vertex, lower_vertex ));  //fpd
		vertices_[ lower_vertex ]->add_edge( upper_vertex );
		++num_edges_;
	}

	// add an edge and set its data
	void
	add_edge( platform::Size lower_vertex, platform::Size upper_vertex, E const & edge_data )
	{
		//assert( lower_vertex < upper_vertex );
		assert( ! edge_exists( lower_vertex, upper_vertex ));
		assert( ! edge_exists( upper_vertex, lower_vertex ));  //fpd
		vertices_[ lower_vertex ]->add_edge( upper_vertex, edge_data );
		++num_edges_;
	}

	// slow; O(V)
	bool
	edge_exists( platform::Size lower_vertex, platform::Size upper_vertex )
	{
		return vertices_[ lower_vertex ]->edge_exists( upper_vertex );
	}

	EdgeClass * get_edge( platform::Size lower_vertex, platform::Size upper_vertex )
	{
		return vertices_[ lower_vertex ]->get_edge( upper_vertex );
	}

	platform::Size num_edges() const { return num_edges_; }
	platform::Size num_vertices() const { return num_vertices_; }

	void drop_all_edges() {
		for ( platform::Size ii = 1; ii <= vertices_.size(); ++ii ) vertices_[ ii ]->drop_all_edges();
		num_edges_ = 0;
	}

	friend class UEEdge< V, E >;
	friend class UEVertex< V, E >;

private:

	void
	create_vertices()
	{
		if ( vertices_.size() != 0 ) vertices_.clear();
		if ( vertices_.size() != num_vertices_ ) vertices_.resize( num_vertices_ );
		for ( platform::Size ii = 1; ii <= num_vertices_; ++ii )
		{
			vertices_[ ii ] = utility::pointer::owning_ptr< UEVertex< V, E > >( new UEVertex< V, E >( this, ii ));
			//vertices_[ ii ]->set_index_and_owner( ii, this );
			//vertices_[ ii ]->reserve_edge_space();
		}
		num_edges_ = 0;
	}

	void
	copy_from( UpperEdgeGraph< V, E > const & other )
	{
		num_vertices_ = other.num_vertices_;
		create_vertices();
		for ( platform::Size ii = 1; ii <= other.num_vertices_; ++ii )
		{
			//std::cout << "vertices_.size() " << vertices_.size() << std::endl;
			//std::cout << "vertices_[ ii ].index_ " << vertices_[ ii ].index_ << std::endl;
			vertices_[ ii ]->data() = other.vertices_[ ii ]->data();
			//std::cout << "vertex data copied; vertices_.size() " << vertices_.size() << std::endl;

			//std::cout << "vertices_[ ii ].index_ " << vertices_[ ii ].index_ << std::endl;
			for ( UpperEdgeListConstIter
				iter = other.get_vertex( ii ).const_upper_edge_list_begin(),
				eiter = other.get_vertex( ii ).const_upper_edge_list_end();
				iter != eiter; ++iter )
			{
				//std::cout << "Adding edge.vertices_.size() " << vertices_.size() << std::endl;
				add_edge( ii, iter->upper_vertex(), iter->data() );
			}
		}
	}

	// called by edge class
	VertexClass * get_vertex_ptr( int index ) { return vertices_[ index ].get(); }

	// called by edge class
	void note_edge_deleted() { --num_edges_; }

private:
	// Data

	platform::Size num_vertices_;
	platform::Size num_edges_;
	VertexVector vertices_;

};



}
}

#endif
