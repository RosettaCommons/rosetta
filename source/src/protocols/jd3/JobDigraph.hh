// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobDigraph.hh
/// @brief  Directed graph for representing the relationship between batches of jobs.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_JobDigraph_HH
#define INCLUDED_protocols_jd3_JobDigraph_HH

// Unit Headers
#include <protocols/jd3/JobDigraph.fwd.hh>
#include <utility/graph/Digraph.hh>

// Package Headers
#include <utility/graph/unordered_object_pool.fwd.hpp>
#include <core/types.hh>

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

namespace protocols {
namespace jd3 {

/// @brief A node to use in JobDigraphs that holds information about
/// how much work has completed for the set of jobs that it represents.
class JobDirectedNode : public utility::graph::DirectedNode
{
public:
	typedef utility::graph::DirectedNode parent;
public:
	virtual ~JobDirectedNode();
	JobDirectedNode( utility::graph::Digraph*, platform::Size node_id );

	bool all_jobs_completed() const;
	bool all_jobs_started() const;
	core::Size n_predecessors_w_outstanding_jobs() const;

	void all_jobs_completed( bool setting );
	void all_jobs_started( bool setting );
	/// @brief decrement the n-incomplete-predecessors count by one
	//void note_predecessor_completed();

	/// @brief invoked during graph assignment operators to copy any
	/// node data from one graph to another graph.  The source node must
	/// be the same type as this node.
	virtual void copy_from( DirectedNode const * source );

	virtual
	void add_incoming_edge( utility::graph::DirectedEdge* edge_ptr, DirectedEdgeListIter & );

	/// @brief memory accounting scheme
	virtual platform::Size count_static_memory() const;
	/// @brief memory accounting scheme
	virtual platform::Size count_dynamic_memory() const;

#ifdef SERIALIZATION
	/// @brief Serialization function, but one that is called by the graph so that
	/// the JobDirectedNode object can save and then later restore its data.
	/// The archive is not expected to construct the JobDirectedNode iself, so instead
	/// of calling this function "save" as the archive would need if it were to
	/// perform the construction itself, it's named "save_to_archive".
	template < class Archive >
	void save_to_archive( Archive & archive ) const;

	/// @brief Deserialization function, but one that is called on the already-constructed
	/// JobDirectedNode object instead of being called by the archive directly.
	template < class Archive >
	void load_from_archive( Archive & archive );
#endif // SERIALIZATION

protected:

	/// @brief derived class access to the owner
	inline
	JobDigraph* get_job_owner() const;

private:

	bool all_jobs_completed_;
	bool all_jobs_started_;

	core::Size n_predecessors_w_outstanding_jobs_;

	//no default constructor, uncopyable
	JobDirectedNode();
	JobDirectedNode( JobDirectedNode const & );
	JobDirectedNode & operator = ( JobDirectedNode & );
};

/// @brief A run-of-the-mill directed edge to use in JobDigraphs
class JobDirectedEdge : public utility::graph::DirectedEdge
{
public:
	typedef utility::graph::DirectedEdge parent;
public:

	virtual ~JobDirectedEdge();

	/// @brief Main edge constructor.  This should only be invoked by create_new_edge, which
	/// itself is only called by add_edge.  The ONLY way an edge should be added to a graph
	/// is through add_edge.  NOTE: edges should be only be deleted by a call to the Digraph's
	/// delete_edge method, and this method absolutely must be implemented by derived Digraph
	/// classes.
	JobDirectedEdge( utility::graph::Digraph* owner, platform::Size tail_node_ind, platform::Size head_node_ind);

	/// @brief copy-from for use in Digraph::operator= and copy ctors.  The source node
	/// must be a JobDirectedEdge
	virtual void copy_from( utility::graph::DirectedEdge const * source );

	/// @brief how much memory is statically allocated by this edge
	virtual platform::Size count_static_memory() const;
	/// @brief how much memory is dynamically allocated by this edge -- must be recursively invoked
	/// by a derived class.
	virtual platform::Size count_dynamic_memory() const;

protected:

	//Read access to private data granted to derived classes

	/// @brief get a const * to one node that this edge is incident upon
	/// uses c-style index-from-0 for these two nodes. 0 is the index of the tail node,
	/// 1 is the index of the head node.
	inline
	JobDirectedNode const *
	get_job_node( platform::Size index ) const
	{
		return static_cast< JobDirectedNode const * > (get_node( index ));
	}

	/// @brief get a non-const * to one node that this edge is incident upon
	/// uses c-style index-from-0 for these two nodes
	inline
	JobDirectedNode *
	get_job_node( platform::Size index )
	{
		return static_cast< JobDirectedNode * > (get_node( index ));
	}


	/// @brief get a const * to the owning graph
	JobDigraph const *
	get_job_owner() const;

	/// @brief get a non-const * to the owning graph
	JobDigraph *
	get_job_owner();


private:

	//no default constructor, uncopyable
	JobDirectedEdge();
	JobDirectedEdge( JobDirectedEdge const & );
	JobDirectedEdge & operator = ( JobDirectedEdge & );

};

/// @brief A Digraph for storing the relationship between groups of jobs, so that
/// each node in this graph represents a group of jobs, and if the outputs from one
/// group are inputs to another group, then a directed edge from the first group
/// to the second group will be in the graph.  The JobDistributor will run all of the
/// jobs in the first group before any of the jobs in the second group will start.
class JobDigraph : public utility::graph::Digraph
{
public:
	typedef utility::graph::Digraph parent;

public:

	/// @brief virtual destructor.  Derived classes must ensure they've destroyed all their
	/// nodes and edges through a call to "destroy_everything" before this function is arrived at
	virtual ~JobDigraph();

	/// @brief ctor
	JobDigraph();

	/// @brief num nodes ctor
	JobDigraph(platform::Size num_nodes);

	/// @brief copy ctor.  Must not be called by derived class copy ctors.
	JobDigraph( JobDigraph const & source );

	/// @brief assignment operator.  source and this must have the same type.
	JobDigraph & operator = ( JobDigraph const & source );

public:
	inline
	JobDirectedNode const * get_job_node( platform::Size index ) const
	{
		return static_cast< JobDirectedNode const * > (get_node( index ));
	}

	inline
	JobDirectedNode* get_job_node( platform::Size index )
	{
		return static_cast< JobDirectedNode * > (get_node( index ));
	}

	/// @brief returns a pointer to the directed edge connecting nodes tail_node and head_node, if that edge exists
	/// in the graph, o.w. returns 0.  Focuses the graph on this edge for fast subsequent retrieval.
	JobDirectedEdge * find_job_edge( platform::Size tail_node, platform::Size head_node );

	/// @brief returns a const pointer to the directed edge connecting nodes tail_node and head_node, if that edge exists
	/// in the graph, o.w. returns 0.  Focuses the graph on this edge for fast subsequent retrieval.
	JobDirectedEdge const * find_job_edge( platform::Size tail_node, platform::Size head_node ) const;

	/// @brief remove an edge from the graph. (NEW AS OF 12/9/07) Never call C++'s
	/// "delete" function on an edge pointer directly.  Derived classes must implement this function.
	/// If they wish to use unordered_object_pools to manage their memory
	virtual void delete_edge( utility::graph::DirectedEdge * edge );

#ifdef    SERIALIZATION
	template < class Archive >
	void save( Archive & archive ) const;

	template < class Archive >
	void load( Archive & archive );
#endif // SERIALIZATION

protected:

	virtual platform::Size count_static_memory() const;
	virtual platform::Size count_dynamic_memory() const;

	/// @brief factory method for job-directed node creation
	virtual utility::graph::DirectedNode* create_new_node( platform::Size node_index );

	/// @brief factory method for job-directed edge creation
	virtual utility::graph::DirectedEdge* create_new_edge( platform::Size index1, platform::Size index2 );

	/// @brief factory method for edge copy-construction
	virtual utility::graph::DirectedEdge* create_new_edge( utility::graph::DirectedEdge const * example_edge );

private:

	/// @brief the pool from which class JobDigraph allocates JobDirectedEdge objects.
	boost::unordered_object_pool< JobDirectedEdge > * job_edge_pool_;

};

/// @brief This class defines the set of operations that a JobQueen can perform to update the JobDigraph that
/// she originally gave to the JobDistributor.  In particular, the JobQueen may only update the JobDigraph
/// by adding new nodes to the graph (which will receive incrementally larger indexes) and then adding
/// edges to the graph such that the head node for the edge must land on one of the newly added nodes in the
/// graph.  The JobDistributor will hand the JobQueen a JobDigraphUpdater through the JobQueen's
/// update_job_digraph method, and in this call, the JobQueen may add as many nodes as she wishes, and as many
/// edges that land on those new nodes, but when the method exits, those nodes petrify: the JobQueen may
/// not add any more edges that land on those new nodes.
class JobDigraphUpdater
{
public:
	JobDigraphUpdater( JobDigraphOP );

	/// @brief Read access to the JobDigraph that the updater holds
	JobDigraphCOP job_dag() const;

	/// @brief Add a new node to the JobDigraph that this class holds.
	void add_node();

	/// @brief Add an edge to the graph where the head_node for this edge has to be one of the nodes added
	/// to the graph since the creation of this updater
	void add_edge_to_new_node( core::Size tail_node, core::Size head_node );

	core::Size orig_num_nodes() const { return orig_num_nodes_; }

private:
	JobDigraphOP job_digraph_;
	core::Size orig_num_nodes_; // used to ensure that all added edges land on one of the newly added nodes

};

inline
JobDigraph* JobDirectedNode::get_job_owner() const
{
	return static_cast< JobDigraph * > ( get_owner() );
}

/// @brief get a const * to the owning graph
inline
JobDigraph const *
JobDirectedEdge::get_job_owner() const
{
	return static_cast< JobDigraph const * > ( get_owner() );
}

/// @brief get a non-const * to the owning graph
inline
JobDigraph *
JobDirectedEdge::get_job_owner()
{
	return static_cast< JobDigraph * > ( get_owner() );
}


} //end namespace jd3
} //end namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_JobDigraph )
#endif // SERIALIZATION

#endif
