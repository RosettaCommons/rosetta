// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobGenealogist.hh
/// @brief  class declaration for JobGenealogist
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jack Maguire, jackmaguire1444@gmail.com
/// @details See here for more information: https://www.rosettacommons.org/docs/latest/development_documentation/tutorials/jd3_derived_jq/classes/job_genealogist
/*!

This class was designed to help job queen developers
keep track of the sources for their jobs.

The underlying data in this class are represented
by a directed acyclic graph alternating between
job nodes and results nodes. Consider this example:

Job1 -|-- Result 1 -|-- Job 3 --- Result 1
......|             |
......|             |-- Job 4 -|- Result 1
......|                        |
......|-- Result 2             |- Result 2


Job2 -|-- Result 1 -|-- Job 5
....................|
....................|-- Job 6


Let's say Nobs 1 & 2 come from the first node in the job DAG
and the rest come from the second job DAG node. We can use this to:

1) Track the original input source (<Job> tag, most likely) for each job.

2) Find job results to discard. Maybe we originally wanted to keep result { Job1, Result2 }
but now we can delete all of the jobs from DAG node 1 that do not have offspring such as { Job1, Result2 }.
The garbage_collection() method helps us find those job results.


Jobs are represented using JGJobNodes and results with JGResultNodes.
This does not have to be a tree; you can have a JGJobNode with multiple parents (JGResultNodes).


To keep this up to date (only on node 0, of course), make sure to:

1) call JobGenealogist::register_new_job() inside JobQueen::determine_job_list()

2) call JobGenealogist::note_job_completed() inside JobQueen::note_job_completed()

3) call JobGenealogist::discard_job_result() for every job result id
returned in JobQueen::job_results_that_should_be_discarded()
*/
#ifndef INCLUDED_protocols_jd3_JobGenealogist_hh
#define INCLUDED_protocols_jd3_JobGenealogist_hh

// Project headers
#include <core/types.hh>
#include <protocols/jd3/JobGenealogist.fwd.hh>
#include <protocols/jd3/CompletedJobOutput.fwd.hh>

// Utility headers
#include <utility/graph/Digraph.hh>
#include <utility/vector1.hh>
#include <utility/graph/unordered_object_pool.hpp>
#include <boost/pool/pool.hpp>

// Numeric headers
#include <numeric/DiscreteIntervalEncodingTree.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>
#include <set>
#include <utility>

namespace protocols {
namespace jd3 {

class JGJobNode{

public:

	JGJobNode();

	JGJobNode(
		core::Size global_job,
		unsigned int job_dag_node,
		JGResultNode * parent_node,
		unsigned int input_source
	);

	~JGJobNode() = default;

	///@brief This is purely for the sake of sorting
	bool operator < ( JGJobNode const & rhs) const {
		return global_job_id_ < rhs.global_job_id_;
	}

	///@brief Like humans, the JGJobNode can have more than one parent.
	///@details please only set tell_parent_to_add_child to false if you plan on adding this job node to p's vector of children. The parent/child relationship is very complicated when only one of the two parties knows about the relationship.
	void add_parent( JGResultNode * p, bool tell_parent_to_add_child = true );

	///@details please only set tell_parent_to_add_child to false if you plan on removing this job node from p's vector of children. The parent/child relationship is very complicated when only one of the two parties knows that the relationship has ended.
	bool remove_parent( JGResultNode * p, bool tell_parent_to_remove_child = true );

public: //getters and setters

	core::Size global_job_id() const {
		return global_job_id_;
	}

	void global_job_id( core::Size global_job_id ) {
		global_job_id_ = global_job_id;
	}

	unsigned int node() const {
		return node_;
	}

	///@brief "Node" is the most over-used word in this file. This method calls node() but has a more descriptive name so that we know what type of node we are talking about.
	unsigned int job_dag_node() const {
		return node_;
	}

	void node( unsigned int node ) {
		node_ = node;
	}

	///@brief "Node" is the most over-used word in this file. This method calls node() but has a more descriptive name so that we know what type of node we are talking about.
	void job_dag_node( unsigned int node ) {
		node_ = node;
	}

	unsigned int input_source_id() const {
		return input_source_id_;
	}

	void input_source_id( unsigned int input_source_id ) {
		input_source_id_ = input_source_id;
	}

	utility::vector1< JGResultNode * > & parents() {
		return parents_;
	}

	///@brief Be careful here! The vector is const but the elements are not
	utility::vector1< JGResultNode * > const & parents() const {
		return parents_;
	}

	utility::vector1< JGResultNode * > & children() {
		return children_;
	}

	///@brief Be careful here! The vector is const but the elements are not
	utility::vector1< JGResultNode * > const & children() const {
		return children_;
	}


private:
	core::Size global_job_id_;
	unsigned int node_;
	unsigned int input_source_id_;

	utility::vector1< JGResultNode * > parents_;
	utility::vector1< JGResultNode * > children_;
};


class JGResultNode {

public:

	friend class JGJobNode;
	friend class JobGenealogist;

	JGResultNode();

	JGResultNode( unsigned int result, JGJobNode * par );

	~JGResultNode() = default;

	///@brief Some people like to sort. This is for sorting.
	bool operator < ( JGResultNode const & rhs) const {
		if ( parent_->node() != rhs.parent_->node() ) return parent_->node() < rhs.parent_->node();
		return result_id_ < rhs.result_id_;
	}

	///@brief Like humans, the JGResultNode can have more than one child.
	///@details please only set tell_child_to_add_parent to false if you plan on adding this result node to c's vector of parents. The parent/child relationship is very complicated when only one of the two parties knows about the relationship.
	void add_child( JGJobNode * c, bool tell_child_to_add_parent = true );

	///@details please only set tell_child_to_remove_parent to false if you plan on removing this result node from c's vector of parents. The parent/child relationship is very complicated when only one of the two parties knows that the relationship has ended.
	bool remove_child( JGJobNode * c, bool tell_child_to_remove_parent = true );
public://getters and setters

	unsigned int result_id() const {
		return result_id_;
	}

	void result_id( unsigned int result_id ) {
		result_id_ = result_id;
	}

	JGJobNode const * parent() const {
		return parent_;
	}

	void parent( JGJobNode * parent ) {
		parent_ = parent;
	}

	///@brief Be careful here! The vector is const but the elements are not
	utility::vector1< JGJobNode * > const & children() const {
		return children_;
	}

protected: //nonconst getters

	JGJobNode * parent() {
		return parent_;
	}

	utility::vector1< JGJobNode * > & children() {
		return children_;
	}

private:

	unsigned int result_id_;
	JGJobNode * parent_;
	utility::vector1< JGJobNode * > children_;
};

struct compare_job_nodes : public std::binary_function< JGJobNode const *, JGJobNode const *, bool >{
	bool operator()( JGJobNode const * const a, JGJobNode const * const b) const {
		return a->global_job_id() < b->global_job_id();
	}

	bool operator()( JGJobNode * const & a, const JGJobNode & b ) const {
		return a->global_job_id() < b.global_job_id();
	}
};


class JobGenealogist : public utility::pointer::ReferenceCount{

public:
	JobGenealogist(
		core::Size num_job_dag_nodes,
		core::Size num_input_sources
	);

	~JobGenealogist();

public:
	///@brief register a new job that does not depend on a previous job result but rather takes a pose directly from an input source
	JGJobNode * register_new_job(
		core::Size job_dag_node_id,
		core::Size global_job_id,
		core::Size input_source_id
	);

	///@brief register a new job that depends on a single parent job result
	JGJobNode * register_new_job(
		core::Size job_dag_node_id,
		core::Size global_job_id,
		core::Size job_dag_node_id_of_parent,
		core::Size global_job_id_of_parent,
		core::Size result_id_of_parent
	);

	///@brief register a new job that depends on a single parent job result
	JGJobNode * register_new_job(
		core::Size job_dag_node_id,
		core::Size global_job_id,
		core::Size job_dag_node_id_of_parent,
		jd3::JobResultID const & id
	) {
		return register_new_job( job_dag_node_id, global_job_id, job_dag_node_id_of_parent, id.first, id.second );
	}

	///@brief register a new job that depends on a single parent job result
	JGJobNode * register_new_job(
		core::Size job_dag_node_id,
		core::Size global_job_id,
		JGResultNode * parent
	);

	///@brief register a new job that depends on multiple parent job results
	JGJobNode * register_new_job(
		core::Size job_dag_node_id,
		core::Size global_job_id,
		utility::vector1< JGResultNode * > const & parents
	);

	///@brief Creates nresults new JGResultNodes for this job_node
	void note_job_completed( JGJobNode * job_node, core::Size nresults );

	///@brief wrapper for the other overload. This one is designed to more closely match the argument provided to JobQueen::note_job_completed
	void note_job_completed( core::Size dag_node_id, core::Size global_job_id, core::Size nresults ) {
		JGJobNode * job_node = get_job_node( dag_node_id, global_job_id );
		debug_assert( job_node );
		note_job_completed( job_node, nresults );
	}

	///@brief At the end of job_results_that_should_be_discarded(), call this funciton for every JobResultID in the list
	void discard_job_result( core::Size job_dag_node, core::Size global_job_id, core::Size result_id );

	///@brief At the end of job_results_that_should_be_discarded(), call this funciton for every JobResultID in the list
	void discard_job_result( core::Size job_dag_node, jd3::JobResultID const & id ) {
		discard_job_result( job_dag_node, id.first, id.second );
	}

	///@brief return every job result id in this job_dag_node that does not have any jobs spawned after it.
	///@param[out] container_for_discarded_result_ids
	///@details Consider the patter Job1 -> Result1 -> Job5 and Job5 has no children (results). Should Job1/Result1 be garbage collected? If so, set delete_downstream_job_if_it_has_no_results to true.
	void garbage_collection(
		core::Size job_dag_node,
		bool delete_downstream_job_if_it_has_no_results,
		std::list< jd3::JobResultID > & container_for_discarded_result_ids
	);

	///@brief This method populates the list with every JobResultID for this job dag node
	///@param[out] container_for_output
	void all_job_results_for_node(
		core::Size job_dag_node,
		std::list< jd3::JobResultID > & container_for_output
	) const;

	///@brief This returns the connectivity of the graph as a newick tree. If you graph is not a tree, this will have duplicate elements.
	std::string newick_tree() const;

	core::Size input_source_for_job( core::Size job_dag_node, core::Size global_job_id ) const;

	JGJobNode const * get_const_job_node( core::Size job_dag_node, core::Size global_job_id ) const;

	JGResultNode const * get_const_result_node( core::Size node, core::Size global_job_id, core::Size result_id ) const;

	///@brief This is more for debugging than anything. Print the global_job_ids for every job dag node to the screen
	void print_all_nodes();

protected:

	JGJobNode * get_job_node( core::Size job_dag_node, core::Size global_job_id );

	JGResultNode * get_result_node( core::Size node, core::Size global_job_id, core::Size result_id );

	///@brief Are you all done with this JGJobNode? If so, this method deletes it and removes all traces of it.
	///THIS DOES NOT DELETE ANY OF THE CHILDREN
	void delete_node( JGJobNode * job_node, unsigned int job_dag_node, bool delete_from_vec = true );

	///@brief Are you all done with this JGJobNode? If so, this method deletes it and removes all traces of it.
	///THIS DOES NOT DELETE ANY OF THE CHILDREN
	void delete_node( JGJobNode * job_node, bool delete_from_vec = true ){
		debug_assert( job_node->node() );
		delete_node( job_node, job_node->node(), delete_from_vec );
	}

	///@brief Are you all done with this JGResultNode? If so, this method deletes it and removes all traces of it.
	///THIS DOES NOT DELETE ANY OF THE CHILDREN
	void delete_node( JGResultNode * result_node ){
		result_node_pool_.destroy( result_node );
	}

	///@brief This is currently just a wrapper for JGJobNode::input_source_id(). I am leaving it as a method because it may become more complicated in the future.
	unsigned int input_source_for_node( JGJobNode const * job_node ) const {
		debug_assert( job_node );
		return job_node->input_source_id();
	}

private:
	///@brief recursive utility function for std::string newick_tree();
	void add_newick_tree_for_node( JGResultNode const *, std::stringstream & ) const;

private:
	core::Size num_input_sources_;

	utility::vector1< utility::vector1< JGJobNode * > > job_nodes_for_dag_node_;

	boost::unordered_object_pool< JGJobNode > job_node_pool_;
	boost::unordered_object_pool< JGResultNode > result_node_pool_;

	compare_job_nodes sorter;
};


inline core::Size JobGenealogist::input_source_for_job( core::Size job_dag_node, core::Size global_job_id ) const {
	JGJobNode const * const job_node = get_const_job_node( job_dag_node, global_job_id );
	debug_assert( job_node );
	return input_source_for_node( job_node );
}

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_JobGenealogist_HH
