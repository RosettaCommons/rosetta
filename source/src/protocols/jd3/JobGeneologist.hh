// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobGeneologist.hh
/// @brief  class declaration for JobGeneologist
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_JobGeneologist_hh
#define INCLUDED_protocols_jd3_JobGeneologist_hh

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/graph/Digraph.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/DiscreteIntervalEncodingTree.hh>

// C++ headers
#include <map>
#include <set>
#include <utility>

namespace protocols {
namespace jd3 {

/// @brief The JobGeneologist is meant to manage the annoying / complex / delicate
/// task of keeping track of what jobs belong to which node and what their parents
/// and whether a job has any undiscarded descendants. This task would otherwise fall
/// to the JobQueen, and there are already enough tasks for the JobQueen to manage.
///
/// @details The typical usage scenario is that the JobQueen tells the JobGeneologist
/// I have N nodes and for each node, here is an upper bound on the number of jobs
/// that will be run for this node. It's alright if the JobQueen doesn't have all of
/// that information to start; she can tell the JobGeneologist about new nodes in
/// the JobDAG as she gets to them. (The JobQueen does not need to tell the
/// JobGeneologist about the edges in the JobDAG -- only the nodes). The JG will
/// alot job-index ranges for each node based on the upper bounds. These are the
/// "target ranges."
///
/// For example, if you had a 5 stage protocol where you were going to take the
/// best 10 structures out of those that pass some quality filter that resulted
/// from each stage and use them as seeds for 1000 new jobs (100 replicates for
/// each one), then you would say that each of the 5 nodes has a projected 1000
/// jobs -- but you don't know if there are going to actually be 1000 jobs,
/// because it's possible that only 9 strucures will pass the filters.
///
/// When it is time to begin preparing for stage 2, e.g., the JQ will tell the
/// JG, jobs 150, 221, 388, ..., and 658 will be the parent to 100 replicates
/// using the "append_parents_and_n_replicate_jobs_for_node" function. As results
/// for round 2 trickle in, and the JobQueen filters out the worst ones, she
/// will tell the JG that those jobs have been discarded using the note_job_discarded()
/// function.
///
/// The geneologist needs to be informed of the actual job range for a node if you
/// are not using the "append_parents_..." method. If you have a set of jobs that
/// do not have parents, you can either call the "append_parents_..." function
/// where you pass in a vector of empty vectors, OR you can call
/// set_actual_njobs_for_node.
///
/// The JobQueen can even ask the JG to keep track of the jobs she has handed out
/// for each node using the jobs_remain_for_node and get_next_job_for_node functions
///
/// Perhaps the most important functionality of the JobGeneologist is keeping track
/// of ancestries. If you want to keep hold of intermediate structures so that you
/// can output the structures that led to the most interesting / lowest-energy
/// structures at the end of protocol, then JG will manage the tricky task of telling
/// you which intermidiate structures its safe to discard given that they do not
/// have any descendants that have not been discarded. (If you do go and discard those
/// jobs, make sure to tell the JobGeneologist about it.) The JG can also track which
/// jobs have and have not been discarded. It uses a Discrete Interval Encoding Tree (DIET)
/// to efficiently do so.
class JobGeneologist
{
public:
	typedef core::Size Size;
	typedef utility::vector1< Size > Sizes;

public:
	JobGeneologist();
	JobGeneologist( JobGeneologist const & );
	~JobGeneologist();
	JobGeneologist & operator = ( JobGeneologist const & rhs );

	void set_num_nodes( Size num_nodes );
	void add_node();
	void set_target_num_jobs_for_node( Size node_id, Size num_jobs );

	void append_parents_and_n_replicate_jobs_for_node(
		Size node_id,
		Sizes const & parent_job_index,
		Sizes const & n_replicate_jobs_for_each_parent_group
	);

	void append_parents_and_n_replicate_jobs_for_node(
		Size node_id,
		Sizes const & parent_job_index,
		Size n_replicate_jobs_for_all_parent_group
	);

	void append_parents_and_n_replicate_jobs_for_node(
		Size node_id,
		utility::vector1< Sizes > const & parent_job_indices,
		Sizes const & n_replicate_jobs_for_each_parent_group
	);

	void append_parents_and_n_replicate_jobs_for_node(
		Size node_id,
		utility::vector1< Sizes > const & parents_job_indices,
		Size n_replicate_jobs_for_all_parent_group
	);

	/// @brief Inform the geneologist of how many jobs for a given node there
	/// are in the event that the the node will not have
	/// append_ancestors_and_n_replicate_jobs_for_node called for it.
	void set_actual_njobs_for_node( Size node_id, Size num_jobs );

	void note_job_discarded( Size job_id );
	bool job_has_been_discarded( Size job_id ) const;

	bool jobs_remain_for_node( Size node_id ) const;
	Size get_next_job_for_node( Size node_id );

	Size get_num_nodes() const;
	Size get_node_target_range_begin( Size node_id ) const;
	Size get_node_target_range_end(   Size node_id ) const;
	Size get_node_actual_range_begin( Size node_id ) const;
	Size get_node_actual_range_end(   Size node_id ) const;

	bool job_has_any_parents( Size job_id ) const;
	bool job_from_node_has_any_parents( Size job_id, Size node_id ) const;
	Size get_parent_for_job( Size job_id ) const;
	Sizes const & get_parents_for_job( Size job_id ) const;
	Sizes const & get_parents_for_job_for_node( Size job_id, Size node_id ) const;
	Sizes get_all_ancestors_for_job( Size job_id ) const;

	std::list< Size >
	find_descendentless_jobs_backwards_from_node( Size job_id );

	Size node_for_jobid( Size job_id ) const;

	Size replicate_id_for_jobid( Size job_id ) const;

	std::pair< Size, Size > node_and_replicate_id_for_jobid( Size job_id ) const;

	Size replicate_id_for_jobid_from_node( Size job_id, Size node_id ) const;

private:

	Size pjg_ind_for_job_from_node( Size job_id, Size node_id ) const;

	bool job_is_from_node( Size job_id, Size node_id ) const;

	void
	add_upstream_nodes_to_bfs_queue(
		Size const target_node,
		std::list< core::Size > & bfs_node_list,
		utility::vector1< char > & node_in_bfs_list
	) const;

private:

	Size num_nodes_;
	Size committed_job_range_max_;
	Sizes target_njobs_for_node_;
	Sizes actual_njobs_for_node_;
	//utility::vector1< bool > nodes_finished_;
	utility::vector1< std::pair< Size, Size > > target_jobid_ranges_for_node_;
	utility::vector1< std::pair< Size, Size > > actual_jobid_ranges_for_node_;

	utility::vector1< utility::vector1< Sizes > > parent_job_groups_for_node_;
	utility::vector1< Sizes > n_replicate_jobs_for_pjg_for_node_;
	utility::vector1< Sizes > start_job_index_for_pjg_for_node_;
	utility::vector1< std::set< Size > > all_parental_jobs_for_node_;
	utility::vector1< std::map< Size, Sizes > > pjgs_for_parental_job_for_node_;

	utility::vector1< std::set< Size > > living_jobs_w_descendants_for_node_;

	Sizes last_delivered_job_for_node_;

	Sizes job_lookup_jobid_start_;
	Sizes job_lookup_node_index_;

	numeric::DiscreteIntervalEncodingTree< Size > discarded_jobs_;
	utility::graph::Digraph job_dag_;

};


} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_JobGeneologist_HH
