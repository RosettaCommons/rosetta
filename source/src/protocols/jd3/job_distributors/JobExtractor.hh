// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/job_distributors/JobExtractor.hh
/// @brief  Class to handle the process of requesting jobs from the JobQueen, keeping
///         track of which nodes have their work completed, and which nodes each jobs
///         belong to
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_job_distributors_JobExtractor_HH
#define INCLUDED_protocols_jd3_job_distributors_JobExtractor_HH

// Unit headers
#include <protocols/jd3/job_distributors/JobExtractor.fwd.hh>

// Package headers
#include <protocols/jd3/LarvalJob.fwd.hh>
#include <protocols/jd3/JobDigraph.fwd.hh>
#include <protocols/jd3/JobQueen.fwd.hh>

// Project headers
#include <core/types.hh>

// Numeric headers
#include <numeric/DiscreteIntervalEncodingTree.hh>

// Utility headers
#include <utility/heap.hh>

// C++ headers
#include <list>
#include <map>
#include <set>
#include <string>

namespace protocols {
namespace jd3 {
namespace job_distributors {


class JobExtractor : public utility::pointer::ReferenceCount
{
public:
	typedef std::list< core::Size > SizeList;
	typedef std::map< Size, LarvalJobOP > JobMap;
	typedef utility::pointer::shared_ptr< JobMap > JobMapOP;
	typedef std::set< core::Size > JobSet;
	typedef utility::pointer::shared_ptr< JobSet > JobSetOP;
	typedef std::map< core::Size, JobSetOP > OutstandingJobsForDigraphNodeMap;
	typedef std::map< core::Size, core::Size > DigraphNodeForJobMap;
	typedef std::map< core::Size, core::Size > WorkerNodeForJobMap;

public:

	JobExtractor();

	virtual ~JobExtractor();

	/// @brief dummy for master/slave version
	void
	set_job_queen( JobQueenOP queen );

	void set_maximum_jobs_to_hold_in_memory( core::Size max_njobs_at_once );

	JobDigraphOP create_initial_job_dag();

	bool job_queue_empty() const;
	LarvalJobOP pop_job_from_queue();

	void
	note_job_no_longer_running( Size job_id );

	/// @brief Did we just declare a node complete? Returns true if so, and
	/// sets the internal tracking variable to false.
	bool retrieve_and_reset_node_recently_completed();

	/// @brief Should the JobDistributor keep going based on there being jobs in the job queue,
	/// or outstanding jobs that have not completed, or nodes that have not yet been marked
	/// as completed, or the JobQueen providing new nodes in the JobDAG
	bool
	not_done();

	bool
	jobs_remain();

	LarvalJobOP running_job( core::Size job_index ) const;

	bool
	complete() const;

private:

	void
	query_job_queen_for_more_jobs_for_current_node();

	void
	mark_node_as_complete( Size digraph_node );

	void
	find_jobs_for_next_node();

	void
	queue_initial_digraph_nodes_and_jobs();

	//bool
	//jobs_ready_to_go();

private:

	JobQueenOP job_queen_;
	JobDigraphOP job_dag_;
	SizeList digraph_nodes_ready_to_be_run_;

	SizeList worker_nodes_waiting_for_jobs_;

	/// @brief The digraph node for which we are currently
	/// assigning new jobs -- it is possible for multiple
	/// digraph nodes to have their jobs running concurrently.
	Size current_digraph_node_;
	LarvalJobs jobs_for_current_digraph_node_;

	// The set of all currently running jobs
	JobMap running_jobs_;

	// the job-digraph node that each job was spawned as part of
	DigraphNodeForJobMap digraph_node_for_job_;

	// the list of jobs for each digraph node that have not yet completed
	OutstandingJobsForDigraphNodeMap jobs_running_for_digraph_nodes_;

	// Node 0 keeps track of which nodes have receieved which deallocation messages
	// and when all of the remote nodes have received a particular deallocation message,
	// then that deallocation message is deleted
	bool first_call_to_determine_job_list_;

	// Deallocation messages will be requested from the JQ when the last job for a node
	// completes -- I need to keep track of that in the JobExtractor
	bool node_recently_completed_;

	numeric::DiscreteIntervalEncodingTree< core::Size > job_indices_seen_;

	bool complete_;

	core::Size maximum_jobs_to_hold_in_memory_;
};

}//job_distributors
}//jd3
}//protocols

#endif // INCLUDED_protocols_jd3_job_distributors_JobExtractor_HH
