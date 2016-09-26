// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobDigraph.cc
/// @brief  graph base classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/jd3/JobDigraph.hh>

// Package headers
#include <utility/graph/unordered_object_pool.hpp>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/assert.hh>
#include <utility/excn/Exceptions.hh>

//STL Headers
#include <iostream>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

// Boost Headers
#include <utility/graph/unordered_object_pool.hpp>
#include <boost/pool/pool.hpp>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


// #include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.jd3.JobDigraph" );
using namespace ObjexxFCL;

namespace protocols {
namespace jd3 {


//--------------------------------------------------------------------------------------------//
//------------------------------  JobDigraph JobDirectedNode Class ---------------------------//
//--------------------------------------------------------------------------------------------//


/// @brief
/// virtual destructor
JobDirectedNode::~JobDirectedNode()
{}

/// @brief
/// Main constructor, no default constructor nor copy constructor
JobDirectedNode::JobDirectedNode( utility::graph::Digraph * owner, platform::Size node_id ) :
	parent( owner, node_id ),
	all_jobs_completed_( false ),
	all_jobs_started_( false ),
	n_predecessors_w_outstanding_jobs_( 0 )
{}

bool JobDirectedNode::all_jobs_completed() const { return all_jobs_completed_; }
bool JobDirectedNode::all_jobs_started() const { return all_jobs_started_; }
core::Size JobDirectedNode::n_predecessors_w_outstanding_jobs() const
{
	return n_predecessors_w_outstanding_jobs_;
}

/// @details Declaring that a Job node has completed, the Node visits all of its
/// downstream neighbors (i.e. all neighbors connected by a directed edge
/// originating on this node), decrementing their n_predecessors_w_outstanding_jobs
/// counters.
void JobDirectedNode::all_jobs_completed( bool setting )
{
	all_jobs_completed_ = setting;
	if ( all_jobs_completed_ ) {
		for ( DirectedEdgeListIter iter = outgoing_edge_list_begin();
					iter != outgoing_edge_list_end(); ++iter ) {
			JobDirectedNode * neighb = static_cast< JobDirectedNode * >
				( (*iter)->get_head_node() );
			debug_assert( neighb->n_predecessors_w_outstanding_jobs_ > 0 );
			--neighb->n_predecessors_w_outstanding_jobs_;
		}
	}
}

void JobDirectedNode::all_jobs_started( bool setting ) { all_jobs_started_ = setting; }


/// @brief copy-from for use in JobDigraph::operator= and copy ctors;
/// derived classes must define their own version of this function
//  to copy any data stored on nodes
void JobDirectedNode::copy_from( DirectedNode const * source )
{
	debug_assert( dynamic_cast< JobDirectedNode const * > ( source ) );
	parent::copy_from( source );
	JobDirectedNode const * job_source = static_cast< JobDirectedNode const * > ( source );
	all_jobs_completed_ = job_source->all_jobs_completed_;
	all_jobs_started_   = job_source->all_jobs_started_;
}

void
JobDirectedNode::add_incoming_edge(
	utility::graph::DirectedEdge * edge_ptr,
	DirectedEdgeListIter & new_iterator
)
{
	// parent does all the work
	parent::add_incoming_edge( edge_ptr, new_iterator );

	// ask -- is the other node complete? If not, then increment
	// the n_predecessors_w_outstanding_jobs_ counter
	// go to the owner directly, sine the edge we are adding is not yet
	// complete (though, it must have already stored which is the tail node
	// and which is the head node)
	if ( ! get_job_owner()->get_job_node( edge_ptr->get_tail_node_ind() )->all_jobs_completed() ) {
		++n_predecessors_w_outstanding_jobs_;
	}
}

/// @details called on most-derived class.  The most-derived class should NOT recursively call this method
/// on its parent class.  The sizeof function will handle the whole JobDirectedNode (or DerivedJobDirectedNode).
platform::Size JobDirectedNode::count_static_memory() const
{
	return sizeof( JobDirectedNode );
}


/// @details recursively descend through heirarchy accounting for heap memory usage.  Each derived
/// class in the heirarchy should recursively add the amount of dynamic memory its parent
/// allocates by calling parent::count_dynamic_memory
platform::Size JobDirectedNode::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}

//--------------------------------------------------------------------------------------------//
//------------------------------ JobDigraph JobDirectedEdge Class ----------------------------//
//--------------------------------------------------------------------------------------------//

/// @brief destructor
JobDirectedEdge::~JobDirectedEdge() {}

/// @brief main constructor for edge, no default nor copy constructors
///
/// @param owner - [in] - owning InteractionDigraph
/// @param tail_node_ind - [in] - the index of the tail node
/// @param head_node_ind - [in] - the index of the head node
///
JobDirectedEdge::JobDirectedEdge
(
	utility::graph::Digraph* owner,
	platform::Size tail_node_ind,
	platform::Size head_node_ind
) :
	parent( owner, tail_node_ind, head_node_ind )
{}


/// @details derived classes should recursively call the copy_from method to ensure all parent class
/// data is copied.  It just so happens that this method does nothing, but that could change
/// and the derived class should include a call to this function for that reason.
void JobDirectedEdge::copy_from( utility::graph::DirectedEdge const * source )
{
	parent::copy_from( source );
}

/// @brief memory accouting scheme
///
/// @details This is called non-recursively on the most-derived class
platform::Size JobDirectedEdge::count_static_memory() const
{
	return sizeof( JobDirectedEdge );
}

/// @brief memory accounting scheme
///
/// @details This method should be called recursively by derived classes -- that is, each class should
/// recurse to its parent.
platform::Size JobDirectedEdge::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}

//----------------------------------------------------------------------------//
//---------------------------------  JobDigraph Class -----------------------------//
//----------------------------------------------------------------------------//

/// @brief destructor
JobDigraph::~JobDigraph()
{
	delete_everything();
	delete job_edge_pool_; job_edge_pool_ = 0;
}

/// @brief default constructor; creates an empty graph (no nodes, no edges)
JobDigraph::JobDigraph() :
	parent(),
	job_edge_pool_( new boost::unordered_object_pool< JobDirectedEdge > ( 256 ) )
{}

/// @details Do not call the base class's num-nodes constructor in the initialization list,
/// since that constructor calls the polymorphic function create_new_node, and polymorphism
/// does not work during constructors or destructors.
///
/// @param num_ig_nodes - [in] - number of nodes that this graph will contain
JobDigraph::JobDigraph( platform::Size num_nodes ) :
	parent(),
	job_edge_pool_( new boost::unordered_object_pool< JobDirectedEdge > ( 256 ) )
{
	set_num_nodes( num_nodes );
}


/// @brief copy constructor relies on factory methods and virtual "copy_from" methods
/// by calling the assignment operator
JobDigraph::JobDigraph( JobDigraph const & source ) :
	parent(),
	job_edge_pool_( new boost::unordered_object_pool< JobDirectedEdge > ( 256 ) )
{
	parent::operator = ( source );
}

/// @brief operator = ().  Relies on factory methods and virtual "copy_from" methods
///
/// @details operator= must only be performed on graphs of the same type e.g.
/// an EnergyJobDigraph may be copied from another EnergyJobDigraph, but should
/// not be copied from a JobDigraph.
JobDigraph &
JobDigraph::operator = ( JobDigraph const & source )
{
	parent::operator = ( source );
	return *this;
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
JobDirectedEdge const *
JobDigraph::find_job_edge(platform::Size tail_node, platform::Size head_node) const
{
	utility::graph::DirectedEdge const * edge = find_edge( tail_node, head_node );
	return static_cast< JobDirectedEdge const * > ( edge );
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
JobDirectedEdge * JobDigraph::find_job_edge(platform::Size tail_node, platform::Size head_node)
{
	utility::graph::DirectedEdge * edge = find_edge( tail_node, head_node );
	return static_cast< JobDirectedEdge * > ( edge );
}


void JobDigraph::delete_edge( utility::graph::DirectedEdge * edge )
{
	JobDirectedEdge * job_edge = static_cast< JobDirectedEdge * > ( edge );
	job_edge_pool_->destroy( job_edge );
}

platform::Size JobDigraph::count_static_memory() const
{
	return sizeof( JobDigraph );
}

platform::Size JobDigraph::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}


utility::graph::DirectedNode *
JobDigraph::create_new_node( platform::Size index )
{
	return new JobDirectedNode( this, index );
}

utility::graph::DirectedEdge *
JobDigraph::create_new_edge( platform::Size tail_index, platform::Size head_index )
{
	return job_edge_pool_->construct( this, tail_index, head_index );
}

utility::graph::DirectedEdge *
JobDigraph::create_new_edge( utility::graph::DirectedEdge const * example_edge )
{
	debug_assert( dynamic_cast< JobDirectedEdge const * > ( example_edge ) );
	return job_edge_pool_->construct(
		this,
		example_edge->get_tail_node_ind(),
		example_edge->get_head_node_ind()
	);
}

#ifdef    SERIALIZATION
template < class Archive >
void JobDirectedNode::save_to_archive( Archive & arc ) const
{
	// EXEMPT n_predecessors_w_outstanding_jobs_
	arc( all_jobs_completed_, all_jobs_started_ );
}

template < class Archive >
void JobDirectedNode::load_from_archive( Archive & arc )
{
	// EXEMPT n_predecessors_w_outstanding_jobs_
	arc( all_jobs_completed_, all_jobs_started_ );
}

template < class Archive >
void JobDigraph::save( Archive & archive ) const
{
  archive( num_nodes() );

	// JobDirectedNodes and edges will be freshly created when this graph is deserialized
	// EXEMPT job_edge_pool_

  for ( Size ii = 1; ii <= num_nodes(); ++ii ) {
		get_job_node( ii )->save_to_archive( archive );
  }
  archive( num_edges() );
  for ( DirectedEdgeListConstIter iter = const_edge_list_begin(), iter_end = const_edge_list_end(); iter != iter_end; ++iter ) {
    archive( (*iter)->get_tail_node_ind(), (*iter)->get_head_node_ind() );
  }
}

template < class Archive >
void JobDigraph::load( Archive & archive )
{
  Size num_nodes(0); archive( num_nodes );
  set_num_nodes( num_nodes );
	// EXEMPT job_edge_pool_

  for ( Size ii = 1; ii <= num_nodes; ++ii ) {
		get_job_node( ii )->load_from_archive( archive );
  }

  Size num_edges(0); archive( num_edges );
  for ( Size ii = 1; ii <= num_edges; ++ii ) {
    Size tail_node(0), head_node(0); archive( tail_node, head_node );
    add_edge( tail_node, head_node );
  }
}

SAVE_AND_LOAD_SERIALIZABLE( JobDigraph );
#endif // SERIALIZATION

JobDigraphUpdater::JobDigraphUpdater( JobDigraphOP job_digraph ) :
	job_digraph_( job_digraph ),
	orig_num_nodes_( job_digraph->num_nodes() )
{}

JobDigraphCOP
JobDigraphUpdater::job_dag() const
{
	return job_digraph_;
}

void JobDigraphUpdater::add_node()
{
	job_digraph_->add_node();
}

void JobDigraphUpdater::add_edge_to_new_node( core::Size tail_node, core::Size head_node )
{
	if ( head_node <= orig_num_nodes_ ) {
		throw utility::excn::EXCN_Msg_Exception( "Error: JobQueen has tried to add an edge that lands on a node that hasn't been freshly added to the JobDigraph" );
	}
	job_digraph_->add_edge( tail_node, head_node );
}


} //end namespace jd3
} //end namespace protocols

#ifdef    SERIALIZATION
CEREAL_REGISTER_TYPE( protocols::jd3::JobDigraph )
CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_JobDigraph )
#endif // SERIALIZATION
