// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/job_distributors/JobExtractor.cc
/// @brief  A class for managing the complexity of requesting jobs from the JobQueen
///         and keeping track of what JobNodes in the JobDAG still have outstanding jobs
///         and for updating the JobDAG.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <protocols/jd3/job_distributors/JobExtractor.hh>

// Package headers
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/JobDigraph.hh>

// Basic headers
#include <basic/Tracer.hh>
//#include <basic/options/option.hh>
//#include <basic/options/keys/jd3.OptionKeys.gen.hh>
//#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility headers
#include <utility/assert.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
//#include <utility/mpi_util.hh>

static basic::Tracer TR( "protocols.jd3.job_distributors.JobExtractor" );

namespace protocols {
namespace jd3 {
namespace job_distributors {

using core::Size;

JobExtractor::JobExtractor() :
	current_digraph_node_( 0 ),
	first_call_to_determine_job_list_( true ),
	node_recently_completed_( false ),
	complete_( false ),
	maximum_jobs_to_hold_in_memory_( 1000 )
{}

JobExtractor::~JobExtractor() {}

void
JobExtractor::set_job_queen( JobQueenOP queen )
{
	job_queen_ = queen;
}

void JobExtractor::set_maximum_jobs_to_hold_in_memory( core::Size max_njobs_at_once )
{
	maximum_jobs_to_hold_in_memory_ = max_njobs_at_once;
}

JobDigraphOP
JobExtractor::create_initial_job_dag()
{
	job_dag_ = job_queen_->initial_job_dag();
	queue_initial_digraph_nodes_and_jobs();
	return job_dag_;
}

bool
JobExtractor::job_queue_empty() const
{
	return jobs_for_current_digraph_node_.empty();
}

LarvalJobOP
JobExtractor::pop_job_from_queue()
{
	// Invariant that the JobExtractor maintains: either the jobs_for_current_digraph_node_
	// queue has jobs in it, or the digraph_nodes_ready_to_be_run_ list is empty.

	// After a job is popped off the queue, the JobExtractor needs to check whether the
	// queue is empty -- if it is empty, then we should ask the JobQueen for more jobs
	// from the current_digraph_node_, and if she does not provide any, then to look
	// at the other job nodes that should be ready.

	debug_assert( ! jobs_for_current_digraph_node_.empty() );

	LarvalJobOP job = jobs_for_current_digraph_node_.front();
	jobs_for_current_digraph_node_.pop_front();
	digraph_node_for_job_[ job->job_index() ] = current_digraph_node_;
	if ( jobs_running_for_digraph_nodes_.count( current_digraph_node_ ) == 0 ) {
		jobs_running_for_digraph_nodes_[ current_digraph_node_ ] = JobSetOP( new JobSet );
	}
	jobs_running_for_digraph_nodes_[ current_digraph_node_ ]->insert( job->job_index() );
	running_jobs_[ job->job_index() ] = job;

	if ( jobs_for_current_digraph_node_.empty() ) {
		query_job_queen_for_more_jobs_for_current_node();
	}

	debug_assert( digraph_nodes_ready_to_be_run_.empty() || ! jobs_for_current_digraph_node_.empty() );

	return job;
}

void
JobExtractor::note_job_no_longer_running( Size job_id )
{
	// ok, now remove the completed/failed job from the maps keeping track of
	// outstanding jobs

	//TR << "Job no longer running: " << job_id << std::endl;

	running_jobs_.erase( job_id );
	Size digraph_node = digraph_node_for_job_[ job_id ];
	digraph_node_for_job_.erase( job_id );

	JobSetOP digraph_nodes_remaining_jobs =
		jobs_running_for_digraph_nodes_[ digraph_node ];
	debug_assert( digraph_nodes_remaining_jobs );
	digraph_nodes_remaining_jobs->erase( job_id );

	if ( digraph_nodes_remaining_jobs->empty() &&
			( current_digraph_node_ != digraph_node || jobs_for_current_digraph_node_.empty() ) ) {
		// previously, when we ran out of jobs in the jobs_for_current_digraph_node_
		// list to dispense for this node, we queried the JobQueen for more jobs, and
		// she told us there were none.
		// Now, we have completed the last job for this node, and the JobQueen has
		// seen the JobSummary for that last job.
		// At this point, it is appropriate to ask the JobQueen if there should
		// be any new nodes added to the JobDigraph, and to otherwise update
		// the internal set of job nodes which are ready to be run.
		mark_node_as_complete( digraph_node );
	}
}

bool
JobExtractor::retrieve_and_reset_node_recently_completed() {
	bool return_val = node_recently_completed_;
	node_recently_completed_ = false;
	return return_val;
}

bool
JobExtractor::not_done()
{
	// we aren't done if
	// 1. There are jobs in the jobs_for_current_digraph_node_ list
	// 2. The JobQueen adds new nodes to the JobDigraph.
	// 3. There are jobs that have been started but have not completed
	//    which would be found in the jobs_running_for_digraph_nodes_ map

	// Precondition:
	// There should never be a situation in which the jobs_for_current_digraph_node_
	// list is empty, but the digraph_nodes_ready_to_be_run_ list has elements.
	// find_jobs_for_next_node() should enforce that.

	debug_assert( digraph_nodes_ready_to_be_run_.empty() ||
		! jobs_for_current_digraph_node_.empty() );

	if ( ! jobs_for_current_digraph_node_.empty() ) return true;


	// ok -- in here, the JE asks the JQ to update the JobDigraph by adding
	// new nodes and edges to those new nodes.  If she does add new nodes, then
	// this function should look at those nodes to see if they are ready to be run.
	JobDigraphUpdater updater( job_dag_ );
	job_queen_->update_job_dag( updater );
	if ( job_dag_->num_nodes() != updater.orig_num_nodes() ) {
		bool found_digraph_node_ready_to_run = false;
		for ( Size ii = updater.orig_num_nodes() + 1;
				ii <= job_dag_->num_nodes(); ++ii ) {
			if ( job_dag_->get_job_node( ii )->n_predecessors_w_outstanding_jobs() == 0 ) {
				digraph_nodes_ready_to_be_run_.push_back( ii );
				found_digraph_node_ready_to_run = true;
			}
		}
		if ( found_digraph_node_ready_to_run ) {
			find_jobs_for_next_node();
		}
		return true;
	}

	// ok -- the JQ is perhaps not ready to add new nodes to the JobDigraph
	// and we should guarantee her the chance to do so, as long as there are
	// jobs outstanding.  Moreover, we must keep listening to worker nodes
	// that have not finished, so the listener loop must continue.  But if the
	// JQ has not added any new nodes, and we have no running jobs, then
	// that's it -- we're done.
	bool any_jobs_running = ! running_jobs_.empty();

	if ( ! any_jobs_running ) {
		// ok! we're completely done:
		// there were no still-running jobs, and yet the JobQueen did not
		// put any new nodes into the JobDigraph
		complete_ = true;
	}

	return any_jobs_running;
}

/// @brief Are there any jobs that have not yet been executed, but perhaps are not
/// ready to be submitted because the JobDirectedNode they belong to is downstream
/// of another Node whose jobs are still running?  The JobQueen would not have
/// told the JobDistributor about these jobs yet.  Perhaps the JobQueen has not even
/// told the JobDistributor about the JobDirectedNodes yet.  Basically, we must say
/// "yes" as long as there are jobs that have not yet completed unless we've emptied
/// the digraph_nodes_ready_to_be_run_ queue and then asked the JobQueen to update
/// the job DAG, and she has declined to add any new nodes.
///
/// @details This function relies on the not_done() function to have asked the
/// JobQueen to update the JobDigraph, and then to check the
/// jobs_running_for_digraph_nodes_ map to see if it's still empty.
bool
JobExtractor::jobs_remain()
{
	return ! complete_;
}

LarvalJobOP
JobExtractor::running_job( core::Size job_index ) const
{
	auto iter = running_jobs_.find( job_index );
	debug_assert( iter != running_jobs_.end() );
	return iter->second;
}

bool
JobExtractor::complete() const
{
	return complete_;
}



/// @details We have run out of jobs for the digraph node indicated; ask the
/// JobQueen for more jobs, and if she doesn't give us any, then consider
/// the node in the digraph exhausted.
void
JobExtractor::query_job_queen_for_more_jobs_for_current_node()
{
	debug_assert( current_digraph_node_ );

	first_call_to_determine_job_list_ = false;

	LarvalJobs jobs_for_current_node = job_queen_->determine_job_list(
		current_digraph_node_, maximum_jobs_to_hold_in_memory_ );

	// Make sure that we haven't encountered any previous jobs that have the same job index as this job.
	for ( LarvalJobs::const_iterator iter = jobs_for_current_node.begin();
			iter != jobs_for_current_node.end(); ++iter ) {
		if ( ! *iter ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "determine_job_list has returned a null-pointer" );
		}
		if ( job_indices_seen_.member( (*iter)->job_index() ) ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "determine_job_list has returned two jobs with the same job index: " +
				utility::to_string( (*iter)->job_index() ) + ". The job distributor requires that all jobs are given a unique job index." );
		}
		job_indices_seen_.insert( (*iter)->job_index() );
	}

	if ( jobs_for_current_node.empty() ) {
		// mark this node as complete
		job_dag_->get_job_node( current_digraph_node_ )->all_jobs_started();
		// recursive call here:
		// walk through the digraph_nodes_ready_to_be_run_ list to find a node
		// that the JobQueen has jobs for
		find_jobs_for_next_node();
	} else {
		jobs_for_current_digraph_node_.splice( jobs_for_current_digraph_node_.end(), jobs_for_current_node );
	}
}

/// @details Once the JobQueen has informed the %JobExtractor that no more
/// jobs remain for a particular node, then we are in a position where we
/// need to check the Job DAG to see if there were nodes waiting for this particular
/// node to complete.  So we iterate across all of the edges leaving the completed
/// node, and for each node downstream, we look at all of its upstream parents.
/// If each of the upstream parents has completed (all of its jobs have completed),
/// then the node is ready to be queued.
void
JobExtractor::mark_node_as_complete( Size digraph_node )
{
	node_recently_completed_ = true;

	// we no longer have to look at this node in the call to not_done(), so
	// delete it from the jobs_running_for_digraph_ map.
	// scope the running_jobs_iter
	{
		// If we never got any jobs for this node, then its index was
		// never inserted into the jobs_running_for_digraph_nodes_ map,
		// so look to see if it's there before calling erase.
		OutstandingJobsForDigraphNodeMap::iterator running_jobs_iter =
			jobs_running_for_digraph_nodes_.find( digraph_node );
		if ( running_jobs_iter != jobs_running_for_digraph_nodes_.end () ) {
			jobs_running_for_digraph_nodes_.erase( running_jobs_iter );
		}
	}

	JobDirectedNode * done_node = job_dag_->get_job_node( digraph_node );
	// trigger the update of the n_precessors_w_outstanding_jobs counter
	// for all nodes downstream of this node.
	done_node->all_jobs_completed( true );

	// bool job_node_queue_was_empty = digraph_nodes_ready_to_be_run_.empty();

	// look at all children of the node that just completed and ask them
	// if all of their parent's have completed, and if so, then mark those
	// nodes as ready to launch. The JobDirectedNode maintains the invariant
	// that its n_predecessors_w_outstanding_jobs() counter exactly matches
	// the number of incoming edges that connect the node to any upstream
	// node with its "all_jobs_complete_" status set as false; this gives
	// us an O( E ) expense of visiting nodes in the JobDigraph.
	bool found_any_nodes_ready_to_run = false;
	for ( auto done_child_iter = done_node->const_outgoing_edge_list_begin();
			done_child_iter != done_node->const_outgoing_edge_list_end();
			++done_child_iter ) {
		JobDirectedNode const * done_child = dynamic_cast< JobDirectedNode const * >
			((*done_child_iter)->get_head_node());
		if ( done_child->n_predecessors_w_outstanding_jobs() == 0 ) {
			digraph_nodes_ready_to_be_run_.push_back( done_child->get_node_index() );
			found_any_nodes_ready_to_run = true;
		}
	}

	if ( found_any_nodes_ready_to_run ) {
		find_jobs_for_next_node();
	}
}

/// @details We need to find the next set of jobs to run, and so we'll look at the
/// nodes in the digraph_nodes_ready_to_be_run_queue_. Pop one of the nodes off
/// and ask the job queen if there are any nodes for this job.  It is entirely
/// possible that the job queen will return an empty list of jobs for this node,
/// (which we'll detect by looking at the current_digraph_node_ index), in which
/// case, we need to mark the node as complete, which could in turn repopulate
/// the digraph_nodes_read_to_be_run_ queue.
///
/// The call to query_job_queen_for_more_jobs_for_current_node function itself
/// may call find_jobs_for_next_node: infinite recursion is avoided by the following
/// two facts: 1) if the digraph_nodes_ready_to_be_run_ queue is not empty, then
/// we decrease its size by one by popping an element off of it, and 2) each node
/// is only put into the digraph_nodes_ready_to_be_run_ queue a single time.
///
void
JobExtractor::find_jobs_for_next_node()
{
	if ( ! digraph_nodes_ready_to_be_run_.empty() ) {
		Size next_node_to_check = digraph_nodes_ready_to_be_run_.front();
		current_digraph_node_ = next_node_to_check;
		digraph_nodes_ready_to_be_run_.pop_front();

		// Possibly recursive call to end up back at find_jobs_for_next_node.
		// After this exits, we will have to check whether we actually received
		// any jobs from the JobQueen for this node.
		query_job_queen_for_more_jobs_for_current_node();

		if ( current_digraph_node_ != next_node_to_check ) {
			// ok -- then we know that the JobQueen did not request any
			// new jobs for this node!  HA!
			// go ahead and mark all of this node's jobs as having completed.
			mark_node_as_complete( next_node_to_check );
		}
	} else {
		// OK! We have exhausted the full set of digraph_nodes_waiting_to_be_run
		// without having found any jobs to run.  Set the current_digraph_node_
		// to 0 (which is not a legal digraph node index) so that recursive
		// calls to this function, as they recover control of flow, will know
		// that the nodes they queued for checking turned up no jobs.
		current_digraph_node_ = 0;
	}

	// post condition: if there are any nodes that say they are ready
	// to have their jobs run, then the jobs_for_current_digraph_node_
	// queue is not empty.
	debug_assert( digraph_nodes_ready_to_be_run_.empty() ||
		! jobs_for_current_digraph_node_.empty() );

}

void
JobExtractor::queue_initial_digraph_nodes_and_jobs()
{
	runtime_assert( digraph_is_a_DAG( *job_dag_ ) );

	// iterate across all nodes in the job_digraph, and note every
	// node with an indegree of 0
	for ( Size ii = 1; ii <= job_dag_->num_nodes(); ++ii ) {
		if ( job_dag_->get_node(ii)->indegree() == 0 ) {
			digraph_nodes_ready_to_be_run_.push_back( ii );
		}
	}
	debug_assert( ! digraph_nodes_ready_to_be_run_.empty() );
	find_jobs_for_next_node();
}

///// @brief With the invariant that the jobs_for_current_digraph_node_ queue should
///// never have anything in it if the jobs_for_current_digraph_node_ list is
///// empty, then we can simply check whether the jobs_for_current_digraph_node_ queue
///// is empty in order to decide if there are any jobs that are ready to be launched.
//bool
//JobExtractor::jobs_ready_to_go()
//{
// return ! jobs_for_current_digraph_node_.empty();
//}



} // namespace job_distributors
} // namespace jd3
} // namespace protocols


