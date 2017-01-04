// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/job_distributors/MPIWorkPoolJobDistributor.cc
/// @brief  implementation of MPIWorkPoolJobDistributor
/// @author Andy Watkins (amw579@nyu.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifdef SERIALIZATION

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

// Unit headers
#include <protocols/jd3/job_distributors/MPIWorkPoolJobDistributor.hh>

// Package headers
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/JobSummary.hh>
#include <protocols/jd3/JobResult.hh>
#include <protocols/jd3/deallocation/DeallocationMessage.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/assert.hh>
#include <utility/mpi_util.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/details/helpers.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.MPIWorkPoolJobDistributor" );

namespace protocols {
namespace jd3 {
namespace job_distributors {

using core::Size;

MPIWorkPoolJobDistributor::MPIWorkPoolJobDistributor() :
	n_archives_( basic::options::option[ basic::options::OptionKeys::jd3::n_archive_nodes ] ),
	store_on_node0_( n_archives_ == 0 && ! basic::options::option[ basic::options::OptionKeys::jd3::do_not_archive_on_node0 ]() ),
	n_results_per_archive_( n_archives_ + 1 ),
	first_call_to_determine_job_list_( true ),
	default_retry_limit_( 1 ),
	complete_( false )
{
	mpi_rank_ = utility::mpi_rank();
	mpi_nprocs_ = utility::mpi_nprocs();

	if ( mpi_rank_ == 0 ) {
		nodes_spun_down_.resize( mpi_nprocs_-1, false );
		deallocation_messages_for_node_.resize( mpi_nprocs_ );
	}

}

MPIWorkPoolJobDistributor::~MPIWorkPoolJobDistributor() {}

void
MPIWorkPoolJobDistributor::go( JobQueenOP queen )
{
	job_queen_ = queen;
	if ( mpi_rank_ == 0 ) {
		go_master();
	} else if ( mpi_rank_ <= n_archives_ ) {
		go_archive();
	} else {
		go_worker();
	}

#ifdef USEMPI
	MPI_Barrier( MPI_COMM_WORLD );
 	MPI_Finalize();
#endif

}


void
MPIWorkPoolJobDistributor::go_master()
{
	master_setup();

	while ( not_done() ) {
		int worker_node = utility::receive_integer_from_anyone();
		int message = utility::receive_integer_from_node( worker_node );
		switch ( message ) {
		case mpi_work_pool_jd_new_job_request :
			process_job_request_from_node( worker_node );
			break;
		case mpi_work_pool_jd_job_success :
			// two part process; first the worker says "I'm done"
			// and node 0 says "go archive your result to node X"
			process_job_succeeded( worker_node );
			break;
		case mpi_work_pool_jd_job_success_and_archival_complete :
			// second part of the two part process: after archival
			// has completed, the node sends the JobStatus to node 0
			// so that the JobQueen can examine it
			process_archival_complete_message( worker_node );
			break;
		case mpi_work_pool_jd_job_failed_w_message :
			process_job_failed_w_message_from_node( worker_node );
			break;
		case mpi_work_pool_jd_job_failed_do_not_retry :
			process_job_failed_do_not_retry( worker_node );
			break;
		case mpi_work_pool_jd_job_failed_bad_input :
			process_job_failed_bad_input( worker_node );
			break;
		case mpi_work_pool_jd_job_failed_retry_limit_exceeded :
			process_job_failed_retry_limit_exceeded( worker_node );
			break;
		case mpi_work_pool_jd_retrieve_job_result :
			process_retrieve_job_result_request( worker_node );
			break;
		case mpi_work_pool_jd_failed_to_retrieve_job_result :
			// Fatal error: throw an exception and shut down execution
			job_queen_->flush();
			process_failed_to_retrieve_job_result_request( worker_node );
			break;
		case mpi_work_pool_jd_error :
			job_queen_->flush();
			throw utility::excn::EXCN_Msg_Exception( "Received error from node " +
				utility::to_string( worker_node ) + ": " + utility::receive_string_from_node( worker_node ) );
			break;
		default :
			// error -- we should not have gotten here
			throw utility::excn::EXCN_Msg_Exception( "Internal Error in the MPIWorkPoolJobDistributor: recieved inappropriate signal "
				"on node 0 from node " + utility::to_string( worker_node )
				+ ": " + utility::to_string( message ) );
		}
	}

	// spin-down-signal sending
	for ( SizeList::const_iterator iter = worker_nodes_waiting_for_jobs_.begin();
			iter != worker_nodes_waiting_for_jobs_.end(); ++iter ) {
		send_spin_down_signal_to_node( *iter );
	}

	for ( int ii = 1; ii <= n_archives_; ++ii ) {
		utility::send_integer_to_node( ii, 0 ); // archive is expecting a node ID
		send_spin_down_signal_to_node( ii );
	}

	// Other nodes that have not yet sent their final job-request message to node 0.
	while ( any_nodes_not_yet_spun_down() ) {
		int worker_node = utility::receive_integer_from_anyone();
		int message = utility::receive_integer_from_node( worker_node );
		switch ( message ) {
		case mpi_work_pool_jd_new_job_request :
			process_job_request_from_node( worker_node );
			break;
		default:
			// error -- we should not have gotten here
			throw utility::excn::EXCN_Msg_Exception( "recieved inappropriate signal "
				"on node 0 from node " + utility::to_string( worker_node )
				+ ": " + utility::to_string( message ) + " in final spin-down loop" );
		}
	}

	job_queen_->flush();

}

void
MPIWorkPoolJobDistributor::go_archive()
{
	bool keep_going = true;
	while ( keep_going ) {
		int remote_node = utility::receive_integer_from_anyone();
		int message = utility::receive_integer_from_node( remote_node );
		switch ( message ) {
		case mpi_work_pool_jd_retrieve_job_result :
			process_retrieve_job_result_request( remote_node );
			break;
		case mpi_work_pool_jd_retrieve_and_discard_job_result :
			process_retrieve_and_discard_job_result_request( remote_node );
			break;
		case mpi_work_pool_jd_discard_job_result :
			process_discard_job_result_request( remote_node );
			break;
		case mpi_work_pool_jd_archive_job_result :
			process_archive_job_result_request( remote_node );
			break;
		case mpi_work_pool_jd_spin_down :
			keep_going = false; // exit the lister loop
			break;
		default :
			send_error_message_to_node0(
				"Archival node " + utility::to_string( mpi_rank_ ) + " recieved an "
				"illegal message " + utility::to_string( message ) + " from node "
				+ utility::to_string( remote_node ) );
			break;
		}
	}
}

void
MPIWorkPoolJobDistributor::go_worker()
{
	bool done = false;
	while ( ! done ) {
		int message = request_job_from_master();
		switch ( message ) {
		case mpi_work_pool_jd_new_job_available :
			retrieve_job_and_run();
			break;
		case mpi_work_pool_jd_deallocation_message :
			retrieve_deallocation_messages_and_new_job_and_run();
			break;
		case mpi_work_pool_jd_spin_down :
			done = true;
			break;
		}
	}
}

void
MPIWorkPoolJobDistributor::master_setup()
{
	job_dag_ = job_queen_->initial_job_dag();
	queue_initial_digraph_nodes_and_jobs();

	// initialize the n_results_per_archive_ heap
	bool err( false );
	if ( store_on_node0_ ) {
		n_results_per_archive_.heap_insert( 0, 0, err );
	}
	for ( int ii = 1; ii <= n_archives_; ++ii ) {
		n_results_per_archive_.heap_insert( ii, 0, err );
	}
}

void
MPIWorkPoolJobDistributor::process_job_request_from_node( int worker_node )
{
	if ( jobs_ready_to_go() ) {
		send_next_job_to_node( worker_node );
	} else if ( jobs_remain() ) {
		note_node_wants_a_job( worker_node );
	} else {
		send_spin_down_signal_to_node( worker_node );
	}
}

/// @details Retrieve the error message from the worker node, tell the JobQueen that
/// the job failed, and then note that the failed job is no longer running.
void
MPIWorkPoolJobDistributor::process_job_failed_w_message_from_node( int worker_node )
{
	Size job_id = utility::receive_size_from_node( worker_node );
	std::string error_message = utility::receive_string_from_node( worker_node );
	debug_assert( running_jobs_.count( job_id ) );
	LarvalJobOP failed_job = running_jobs_[ job_id ];
	TR.Error << "Job " << job_id << " named " << failed_job <<
		" has exited with the message:\n" << error_message << std::endl;

	job_queen_->note_job_completed( failed_job, jd3_job_status_failed_w_exception );
	note_job_no_longer_running( job_id );
}

void
MPIWorkPoolJobDistributor::process_job_failed_do_not_retry( int worker_node )
{
	Size job_id = utility::receive_size_from_node( worker_node );
	LarvalJobOP failed_job = running_jobs_[ job_id ];
	job_queen_->note_job_completed( failed_job, jd3_job_status_failed_do_not_retry );
	note_job_no_longer_running( job_id );
}

void
MPIWorkPoolJobDistributor::process_job_failed_bad_input( int worker_node )
{
	Size job_id = utility::receive_size_from_node( worker_node );
	LarvalJobOP failed_job = running_jobs_[ job_id ];
	job_queen_->note_job_completed( failed_job, jd3_job_status_inputs_were_bad );
	// TO DO: add code that purges all other jobs with the same inner job as this one
	note_job_no_longer_running( job_id );
}

void
MPIWorkPoolJobDistributor::process_job_failed_retry_limit_exceeded( int worker_node )
{
	Size job_id = utility::receive_size_from_node( worker_node );
	LarvalJobOP failed_job = running_jobs_[ job_id ];
	job_queen_->note_job_completed( failed_job, jd3_job_status_failed_max_retries );
	// TO DO: add code that purges all other jobs with the same inner job as this one
	note_job_no_longer_running( job_id );
}

void
MPIWorkPoolJobDistributor::process_job_succeeded( int worker_node )
{

	Size job_id = utility::receive_size_from_node( worker_node );

	// Now pick a node to archive the job result on, perhaps storing the
	// job result on this node.
	int archival_node = pick_archival_node();
	// tell the worker node where to send the job result
	utility::send_integer_to_node( worker_node, archival_node );

	job_result_location_map_[ job_id ] = archival_node;
	if ( archival_node == 0 ) {
		// NOTE: when worker nodes are sending the larval job / job result pair
		// they do not need to first send an archival-request message or the job id.
		// Store serialized larval job & job result pair
		std::string serialized_job_result_string =
			utility::receive_string_from_node( worker_node );
		job_results_[ job_id ] = serialized_job_result_string;
		utility::send_integer_to_node( worker_node, mpi_work_pool_jd_archival_completed );
	}
}

int
MPIWorkPoolJobDistributor::pick_archival_node()
{
	int archive_node = n_results_per_archive_.head_item();
	float nresults = n_results_per_archive_.heap_head();
	nresults += 1;
	n_results_per_archive_.reset_coval( archive_node, nresults );
	return archive_node;
}


void
MPIWorkPoolJobDistributor::process_archival_complete_message( int worker_node )
{
	// the second half of a job-completion step.
	// at this point, the job result for the completed job
	// has been archived on the archival node (perhaps on this node)
	// so it is safe to inform the JobQueen that the job has completed.
	// It is also now safe for the JobQueen to request to discard
	// or output the result, but if we were to structure the archival
	// process in a single pass manner (i.e. in the process_job_succeeded method),
	// the remote node where the archival is taking place might not have
	// completed -- or even gotten started -- by the time we ask it for the
	// result

	LarvalJobOP larval_job;
	Size job_id;

	// The remote node will also query the JobQueen to know which of the two
	// messages it should send.
	if ( job_queen_->larval_job_needed_for_note_job_completed() ||
			job_queen_->larval_job_needed_for_completed_job_summary() ) {
		std::string serialized_larval_job_string =
			utility::receive_string_from_node( worker_node );
		// Does this call to deserialize_larval_job need a try/catch block?
		// No.
		// This larval job has been deserialized previously, with no opportunity
		// for new data to have been added to it.
		larval_job = deserialize_larval_job( serialized_larval_job_string );
		job_id = larval_job->job_index();
	} else {
		job_id = utility::receive_size_from_node( worker_node );
	}

	// In either case, the remote node must send the job summary.
	std::string serialized_job_summary_string = utility::receive_string_from_node( worker_node );
	JobSummaryOP job_summary;
	try {
		job_summary = deserialize_job_summary( serialized_job_summary_string );
	} catch ( cereal::Exception & e ) {
		throw utility::excn::EXCN_Msg_Exception( "Failed to deserialize the JobSummary for job #" +
			utility::to_string( job_id ) + "\nError message from the cereal library:\n" + e.what() + "\n" );
	}

	// inform the JobQueen of the completed job
	if ( job_queen_->larval_job_needed_for_note_job_completed() ||
			job_queen_->larval_job_needed_for_completed_job_summary() ) {
		job_queen_->note_job_completed( larval_job, job_summary->status() );
		job_queen_->completed_job_summary( larval_job, job_summary );
	} else {
		job_queen_->note_job_completed( job_id, job_summary->status() );
		job_queen_->completed_job_summary( job_id, job_summary );
	}

	// This could not occur until after the archival process completed,
	// as it assumes that the JobResult is present on the node where we
	// wrote down it would be stored.
	note_job_no_longer_running( job_id );

	potentially_output_some_job_results();
	potentially_discard_some_job_results();
}

void
MPIWorkPoolJobDistributor::process_retrieve_job_result_request( int worker_node )
{
	using namespace utility;
	Size job_id = receive_size_from_node( worker_node );

	SerializedJobResultMap::const_iterator iter = job_results_.find( job_id );
	if ( iter == job_results_.end() ) {
		// fatal error: the requested job result is not present at this node.
		// inform node 0 (if this isn't node 0 )
		if ( mpi_rank_ == 0 ) {
			throw_after_failed_to_retrieve_job_result( worker_node, job_id, mpi_rank_ );
		} else {
			send_integer_to_node( worker_node, mpi_work_pool_jd_failed_to_retrieve_job_result );
			// tell node 0 about failed retrieval
			send_integer_to_node( 0, mpi_rank_ );
			send_integer_to_node( 0, mpi_work_pool_jd_failed_to_retrieve_job_result );
			send_integer_to_node( 0, worker_node );
			send_size_to_node( 0, job_id );
			return;
		}
	}

	// ok -- we have found the (serialized) JobResult
	send_integer_to_node( worker_node, mpi_work_pool_jd_job_result_retrieved );
	send_string_to_node( worker_node, iter->second );
}

void
MPIWorkPoolJobDistributor::process_failed_to_retrieve_job_result_request( int archival_node )
{
	// the node making the job-result request, which the archival node was
	// unable to fullfil
	int worker_node = utility::receive_integer_from_node( archival_node );
	Size job_id = utility::receive_size_from_node( archival_node );
	throw_after_failed_to_retrieve_job_result( worker_node, job_id, archival_node );
}

/// @details This function unlike the process_retireve_job_result_request function will only
/// be called on an archive node -- and it will be called only by the master node.
/// therefore, we do not need to handle the case when mpi_rank_ == 0.
void
MPIWorkPoolJobDistributor::process_retrieve_and_discard_job_result_request( int worker_node )
{
	using namespace utility;
	Size job_id = receive_size_from_node( worker_node );

	SerializedJobResultMap::iterator iter = job_results_.find( job_id );
	if ( iter == job_results_.end() ) {
		// fatal error: the requested job result is not present at this node.
		// inform node 0 (this isn't node 0, as this function is only called by archive nodes)
		send_integer_to_node( worker_node, mpi_work_pool_jd_failed_to_retrieve_job_result );
		return;
	}

	// ok -- we have found the (serialized) JobResult
	send_integer_to_node( worker_node, mpi_work_pool_jd_job_result_retrieved );
	send_string_to_node( worker_node, iter->second );

	job_results_.erase( iter );
}

void
MPIWorkPoolJobDistributor::process_discard_job_result_request( int remote_node )
{
	using namespace utility;
	Size job_id = receive_size_from_node( remote_node ); // node 0
	SerializedJobResultMap::iterator iter = job_results_.find( job_id );
	if ( iter == job_results_.end() ) {
		// well, something has gone wrong, so we should probably print an error, but
		// technically, we can keep going.
		std::cerr << "Failed to discard job result " << job_id << " on node " << mpi_rank_ << std::endl;
	}
	job_results_.erase( iter );
}

void
MPIWorkPoolJobDistributor::process_archive_job_result_request( int remote_node )
{
	using namespace utility;
	Size job_id = receive_size_from_node( remote_node );
	std::string serialized_larval_job_and_result = receive_string_from_node( remote_node );
	job_results_[ job_id ] = serialized_larval_job_and_result;
	utility::send_integer_to_node( remote_node, mpi_work_pool_jd_archival_completed );
}

/// @details This function is to be called only when the master node is not already in the
/// middle of communicating with this node about something; the master node should be in the
/// listening loop of go_master when this function is called (or at least it shouldn't be
/// in some other function *and* excepecting this node to send it a message).
void
MPIWorkPoolJobDistributor::send_error_message_to_node0( std::string const & error_message )
{
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_error );
	utility::send_string_to_node( 0, error_message );
}

/// @throws If the LarvalJob fails to serialize (perhaps an object stored in the
/// const_data_map of the inner-larval-job has not implemented the save & load
/// routines), then this function will throw a EXCN_Msg_Exception
void
MPIWorkPoolJobDistributor::send_next_job_to_node( int worker_node )
{
	debug_assert( ! jobs_for_current_digraph_node_.empty() );

	if ( ! deallocation_messages_for_node_[ worker_node ].empty() ) {
		send_deallocation_messages_to_node( worker_node );
	}

	// First tell the worker node that we do indeed have a job for them to work on/
	// The other option is that we are going to tell the worker to spin down.
	utility::send_integer_to_node( worker_node, mpi_work_pool_jd_new_job_available );

	// NOTE: Popping jobs_for_current_digraph_node_: we might violate the
	// invariant that this queue is never empty so long as nodes remain in the
	// digraph_nodes_ready_to_be_run_ queue.  Address that below
	LarvalJobOP larval_job = jobs_for_current_digraph_node_.front();
	jobs_for_current_digraph_node_.pop_front();

	// serialize the LarvalJob and ship it to the worker node.
	// Note that this could leave the worker hanging where it expects
	// a job, but the JD has thrown: this exception must not be caught
	// except to print its contents before exiting.
	std::string serialized_larval_job;
	try {
		serialized_larval_job = serialize_larval_job( larval_job );
	} catch ( cereal::Exception const & e ) {
		// ok, we should throw a utility::excn::EXCN_Msg_Exception letting the user know
		// that the code cannot be used in its current state because some class stored inside
		// the LarvalJob (or perhaps the descendent of the LarvalJob itself) does not implement
		// the save and load serialization funcitons.
		throw utility::excn::EXCN_Msg_Exception( "Failed to serialize LarvalJob " + utility::to_string( larval_job->job_index() )
			+ " because it holds an unserializable object.  The following error message comes from the cereal library:\n"
			+ e.what() );
	}
	utility::send_string_to_node( worker_node, serialized_larval_job );

	// Also inform the remote node where to find the JobResults that will serve as inputs to
	// this job.
	utility::vector1< core::Size > archival_nodes_for_job( larval_job->input_job_result_indices().size() );
	for ( core::Size ii = 1; ii <= archival_nodes_for_job.size(); ++ii ) {
		archival_nodes_for_job[ ii ] = job_result_location_map_[ larval_job->input_job_result_indices()[ ii ] ];
	}
	utility::send_sizes_to_node( worker_node, archival_nodes_for_job );

	running_jobs_[ larval_job->job_index() ] = larval_job;
	digraph_node_for_job_[ larval_job->job_index() ] = current_digraph_node_;
	worker_node_for_job_[ larval_job->job_index() ] = worker_node;
	if ( jobs_running_for_digraph_nodes_.count( current_digraph_node_ ) == 0 ) {
		jobs_running_for_digraph_nodes_[ current_digraph_node_ ] = JobSetOP( new JobSet );
	}
	jobs_running_for_digraph_nodes_[ current_digraph_node_ ]->insert( larval_job->job_index() );


	// Now ensure the invariant that the jobs_for_current_digraph_node_ should
	// never be empty so long as there are job-nodes listed in the
	// digraph_nodes_ready_to_be_run_ queue.  Since we just popped a job off the
	// jobs_for_current_digraph_node_ queue, we need to check if it's empty, and
	// if so, ask the JobQueen for more jobs.
	if ( jobs_for_current_digraph_node_.empty() ) {
		query_job_queen_for_more_jobs_for_current_node();
	}
}

void
MPIWorkPoolJobDistributor::send_deallocation_messages_to_node( int worker_node )
{
	for ( core::Size which_message : deallocation_messages_for_node_[ worker_node ] ) {
		utility::send_integer_to_node( worker_node, mpi_work_pool_jd_deallocation_message );
		utility::send_string_to_node( worker_node, deallocation_messages_[ which_message ] );
		core::Size n_remaining = --n_remaining_nodes_for_deallocation_message_[ which_message ];
		if ( n_remaining == 0 ) {
			deallocation_messages_[ which_message ] = std::string();
		}
	}
	deallocation_messages_for_node_[ worker_node ].clear();
}

void
MPIWorkPoolJobDistributor::note_job_no_longer_running( Size job_id )
{
	// ok, now remove the completed/failed job from the maps keeping track of
	// outstanding jobs
	running_jobs_.erase( job_id );
	Size digraph_node = digraph_node_for_job_[ job_id ];
	digraph_node_for_job_.erase( job_id );
	worker_node_for_job_.erase( job_id );

	JobSetOP digraph_nodes_remaining_jobs =
		jobs_running_for_digraph_nodes_[ digraph_node ];
	debug_assert( digraph_nodes_remaining_jobs );
	digraph_nodes_remaining_jobs->erase( job_id );

	if ( digraph_nodes_remaining_jobs->empty() ) {
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

void
MPIWorkPoolJobDistributor::potentially_output_some_job_results()
{

	// ask the job queen if she wants to output any results
	SizeList jobs_to_output = job_queen_->jobs_that_should_be_output();
	for ( SizeList::const_iterator output_iter = jobs_to_output.begin();
				output_iter != jobs_to_output.end(); ++output_iter ) {
		Size job_id = *output_iter;
		if ( job_result_location_map_.count( job_id ) == 0 ) {
			// TO DO: better diagnostic message here
			throw utility::excn::EXCN_Msg_Exception( "Failed to find job result " +
				utility::to_string( job_id ) + " for outputting as requested by the "
				+ "JobQeen. Has this job already been output or discarded?" );
		}

		Size node_housing_result = job_result_location_map_[ job_id ];
		LarvalJobAndResult job_and_result;
		if ( node_housing_result == 0 ) {
			// ok, we have the result on this node
			try {
				job_and_result = deserialize_larval_job_and_result( job_results_[ job_id ] );
			} catch ( cereal::Exception & e ) {
				throw utility::excn::EXCN_Msg_Exception( "Failed to deserialize LarvalJob & JobResult pair for job #" +
					utility::to_string( job_id ) + "\nError message from the cereal library:\n" + e.what() + "\n" );
			}
			job_results_.erase(  job_id );
		} else {
			utility::send_integer_to_node( node_housing_result, 0 ); // say "I need to talk to you"
			utility::send_integer_to_node( node_housing_result,
				mpi_work_pool_jd_retrieve_and_discard_job_result );
			utility::send_size_to_node( node_housing_result, job_id );

			int reply = utility::receive_integer_from_node( node_housing_result );
			if ( reply == mpi_work_pool_jd_failed_to_retrieve_job_result ) {
				// oh boy -- this should never happen.
				throw utility::excn::EXCN_Msg_Exception( "Failed to retrieve job result " +
					utility::to_string( job_id ) + " for outputting as requested by the "
					+ "JobQueen; it should have been archived on node " +
					utility::to_string( node_housing_result ) + " but was not found there." );
			}
			runtime_assert( reply == mpi_work_pool_jd_job_result_retrieved );

			std::string serialized_job_result_string =
				utility::receive_string_from_node( node_housing_result );
			std::istringstream iss( serialized_job_result_string );
			cereal::BinaryInputArchive arc( iss );
			arc( job_and_result );

		}

		job_queen_->completed_job_result( job_and_result.first, job_and_result.second );

		// and now erase the record of where the job had been stored, and,
		// for the sake of diagnosing errors, also store the fact that this job
		// had been output (as opposed to discarded) -- but compactly within the
		// output_jobs_ DIET
		job_result_location_map_.erase( job_id );
		output_jobs_.insert( job_id );
	}
}

void
MPIWorkPoolJobDistributor::potentially_discard_some_job_results()
{
	// ask the job queen if she wants to discard any results
	SizeList jobs_to_discard = job_queen_->job_results_that_should_be_discarded();
	if ( jobs_to_discard.empty() ) return;

	utility::vector0< int > n_discarded_jobs( n_archives_ + 1, 0 );

	for ( SizeList::const_iterator discard_iter = jobs_to_discard.begin();
				discard_iter != jobs_to_discard.end(); ++discard_iter ) {
		Size job_id = *discard_iter;
		if ( job_result_location_map_.count( job_id ) == 0 ) {
			// TO DO: better diagnostic message here
			throw utility::excn::EXCN_Msg_Exception( "Failed to find job result " +
				utility::to_string( job_id ) + " for outputting as requested by the "
				+ "JobQeen. Has this job already been output or discarded?" );
		}

		Size node_housing_result = job_result_location_map_[ job_id ];
		n_discarded_jobs[ node_housing_result ] += 1;
		LarvalJobAndResult job_and_result;
		if ( node_housing_result == 0 ) {
			// ok, we have the result on this node
			job_results_.erase( job_id );
		} else {
			utility::send_integer_to_node( node_housing_result, 0 ); // say "I need to talk to you"
			utility::send_integer_to_node( node_housing_result,
				mpi_work_pool_jd_discard_job_result );
			utility::send_size_to_node( node_housing_result, job_id );
		}

		// and now delete storage for this job result, and
		// for the sake of diagnosing errors, also store the fact that this job
		// had been discarded (as opposed to being output) -- but compactly within the
		// discarded_jobs_ DIET
		job_result_location_map_.erase( job_id );
		discarded_jobs_.insert( job_id );

	}

	// update the n_results_per_archive_ heap
	for ( int ii = 0; ii <= n_archives_; ++ii ) {
		if ( n_discarded_jobs[ ii ] == 0 ) continue;
		float new_n_results_for_archive_ii = n_results_per_archive_.coval_for_val( ii )
			- n_discarded_jobs[ ii ];
		n_results_per_archive_.heap_replace( ii, new_n_results_for_archive_ii );
	}
}

/// @details We have run out of jobs for the digraph node indicated; ask the
/// JobQueen for more jobs, and if she doesn't give us any, then consider
/// the node in the digraph exhausted.
void
MPIWorkPoolJobDistributor::query_job_queen_for_more_jobs_for_current_node()
{
	debug_assert( current_digraph_node_ );

	if ( ! first_call_to_determine_job_list_ ) {
		std::list< deallocation::DeallocationMessageOP > messages = job_queen_->deallocation_messages();
		if ( ! messages.empty() ) {
			store_deallocation_messages( messages );
		}
	}

	first_call_to_determine_job_list_ = false;

	int const maximum_jobs_to_hold_in_memory = 1000; // TODO: this becomes an option.
	LarvalJobs jobs_for_current_node = job_queen_->determine_job_list(
		current_digraph_node_, maximum_jobs_to_hold_in_memory );

	// Make sure that we haven't encountered any previous jobs that have the same job index as this job.
	// TO DO: refactor into JobDistributor base class.
	for ( LarvalJobs::const_iterator iter = jobs_for_current_node.begin();
			iter != jobs_for_current_node.end(); ++iter ) {
		if ( ! *iter ) {
			throw utility::excn::EXCN_Msg_Exception( "determine_job_list has returned a null-pointer" );
		}
		if ( job_indices_seen_.member( (*iter)->job_index() )) {
			throw utility::excn::EXCN_Msg_Exception( "determine_job_list has returned two jobs with the same job index: " +
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

/// @details Once the JobQueen has informed the JobDistributor that no more
/// jobs remain for a particular node, then we are in a position where we
/// need to check the Job DAG to see if there were nodes waiting for this particular
/// node to complete.  So we iterate across all of the edges leaving the completed
/// node, and for each node downstream, we look at all of its upstream parents.
/// If each of the upstream parents has completed (all of its jobs have completed),
/// then the node is ready to be queued.
void
MPIWorkPoolJobDistributor::mark_node_as_complete( Size digraph_node )
{

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

	bool job_node_queue_was_empty = digraph_nodes_ready_to_be_run_.empty();

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

	// now, if we found a newly active node, and we previously had no
	// active nodes, then perhaps we need to send jobs to the idle nodes
	// waiting for their next job!
	if ( job_node_queue_was_empty && found_any_nodes_ready_to_run ) {
		queue_jobs_for_next_node_to_run();
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
void
MPIWorkPoolJobDistributor::find_jobs_for_next_node()
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
}


void
MPIWorkPoolJobDistributor::queue_jobs_for_next_node_to_run()
{
	debug_assert( jobs_for_current_digraph_node_.empty() );
	debug_assert( ! digraph_nodes_ready_to_be_run_.empty() );

	find_jobs_for_next_node();

	if ( ! worker_nodes_waiting_for_jobs_.empty() &&
			! jobs_for_current_digraph_node_.empty() ) {
		assign_jobs_to_idling_nodes();
		// at the end of this call, either we have run out of jobs to assign
		// or we have run out of idling nodes
	}

	// post condition: if there are any nodes that say they are ready
	// to have their jobs run, then the jobs_for_current_digraph_node_
	// queue is not empty.  If there are any jobs that are ready to be run,
	// then we should have queued their jobs in this function call.
	debug_assert( digraph_nodes_ready_to_be_run_.empty() ||
		! jobs_for_current_digraph_node_.empty() );
}

void
MPIWorkPoolJobDistributor::queue_initial_digraph_nodes_and_jobs()
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
	queue_jobs_for_next_node_to_run();
}

bool
MPIWorkPoolJobDistributor::not_done()
{
	// we aren't done if
	// 1. There are jobs in the jobs_for_current_digraph_node_ list
	// 2. The JobQueen adds new nodes to the JobDigraph.
	// 3. There are jobs that have been started but have not completed
	//    which would be found in the jobs_running_for_digraph_nodes_ map

	// Precondition:
	// There should never be a situation in which the jobs_for_current_digraph_node_
	// list is empty, but the digraph_nodes_ready_to_be_run_ list has elements.
	// queue_jobs_for_next_node_to_run() should enforce that.

	debug_assert( digraph_nodes_ready_to_be_run_.empty() ||
		! jobs_for_current_digraph_node_.empty() );


	if ( ! jobs_for_current_digraph_node_.empty() ) return true;


	// ok -- in here, the JD asks the JQ to update the JobDigraph by adding
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
			queue_jobs_for_next_node_to_run();
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

bool
MPIWorkPoolJobDistributor::any_nodes_not_yet_spun_down()
{
	for ( int ii = 1; ii < mpi_nprocs_; ++ii ) {
		if ( ! nodes_spun_down_[ ii ] ) return true;
	}
	return false;
}



/// @brief With the invariant that the jobs_for_current_digraph_node_ queue should
/// never have anything in it if the jobs_for_current_digraph_node_ list is
/// empty, then we can simply check whether the jobs_for_current_digraph_node_ queue
/// is empty in order to decide if there are any jobs that are ready to be launched.
bool
MPIWorkPoolJobDistributor::jobs_ready_to_go()
{
	return ! jobs_for_current_digraph_node_.empty();
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
MPIWorkPoolJobDistributor::jobs_remain()
{
	return ! complete_;
}

void
MPIWorkPoolJobDistributor::note_node_wants_a_job( int worker_node )
{
	worker_nodes_waiting_for_jobs_.push_back( worker_node );
}

void
MPIWorkPoolJobDistributor::send_spin_down_signal_to_node( int worker_node )
{
	utility::send_integer_to_node( worker_node, mpi_work_pool_jd_spin_down );
	nodes_spun_down_[ worker_node ] = true;
}



//		SizeList nodes_ready_to_be_run = find_job_nodes_w_no_unfinished_ancestors();
//		for ( SizeList::const_iterator node_iter = nodes_ready_to_be_run.begin();
//				node_iter != job_total_order.end(); ++node_iter ) {
//			Size job_node = *node_iter;
//			run_jobs_for_dag_node( job_node );
//		}
//
//		// ok -- in here, the VJD asks the JQ to update the JobDigraph by adding
//		// new nodes and edges to those new nodes.  If the JQ does not add any new nodes,
//		// then the JD will exit this loop
//		JobDigraphUpdater updater( job_dag_ );
//		job_queen_->update_job_dag( updater );
//		if ( job_dag->num_nodes() == updater.orig_num_nodes() ) {
//			// no new nodes in the graph, therefore, we are done
//			break;
//		}
//	}
//}

void
MPIWorkPoolJobDistributor::assign_jobs_to_idling_nodes()
{
	while ( ! worker_nodes_waiting_for_jobs_.empty() &&
			! jobs_for_current_digraph_node_.empty() ) {
		Size worker_node = worker_nodes_waiting_for_jobs_.front();
		worker_nodes_waiting_for_jobs_.pop_front();
		send_next_job_to_node( worker_node );
	}
}

void
MPIWorkPoolJobDistributor::store_deallocation_messages( std::list< deallocation::DeallocationMessageOP > const & messages )
{
	std::ostringstream oss;
	{ // scope
		cereal::BinaryOutputArchive arc( oss );
		arc( messages );
	}

	deallocation_messages_.push_back( oss.str() );
	n_remaining_nodes_for_deallocation_message_.push_back( mpi_nprocs_ - n_archives_ - 1 );
	for ( core::Size ii = 1+n_archives_; ii < (core::Size) mpi_nprocs_; ++ii ) {
		deallocation_messages_for_node_[ ii ].push_back( deallocation_messages_.size() );
	}

}


void
MPIWorkPoolJobDistributor::throw_after_failed_to_retrieve_job_result(
	int worker_node,	Size job_id,
	int archival_node
)
{
	std::ostringstream oss;
	oss << "Failed to retrieve job result " << job_id << " which was requested ";
	oss << "from node " << archival_node << " by node " << worker_node;
	oss << " but was not present there.\n";

	if ( worker_node_for_job_.count( job_id )) {
		oss << "Internal Error in the MPIWorkPoolJobDistributor: JobDistrutor thinks that this job is still running on node " <<
			worker_node_for_job_[ job_id ] << "\n";
	} else if ( job_result_location_map_.count( job_id ) ) {
		oss << "JobDistributor on node 0 thinks the result should have been on node ";
		oss << job_result_location_map_[ job_id ] << "\n";
	} else if ( output_jobs_.member( job_id ) ) {
		oss << "JobDistributor has already output this job (and does not keep jobs ";
		oss << "that have already been output).\n";
	} else if ( discarded_jobs_.member( job_id ) ) {
		oss << "JobDistributor has discarded this job per the JobQueen's request\n";
	} else {
		oss << "Internal Error in the MPIWorkPoolJobDistributor: This job is not listed as still running, as having its JobResult stored on any"
			" archive node, as having been output already, or having been discarded already.\n";
	}

	throw utility::excn::EXCN_Msg_Exception( oss.str() );

}

int
MPIWorkPoolJobDistributor::request_job_from_master()
{
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );
	return utility::receive_integer_from_node( 0 );
}

void
MPIWorkPoolJobDistributor::retrieve_job_and_run()
{
	std::pair< LarvalJobOP, utility::vector1< JobResultCOP > > job_maturation_data =
		retrieve_job_maturation_data();

	// if the LarvalJobOP is null, then we have failed to retrieve a job result from
	// one of the archival nodes. Quit.
	if ( ! job_maturation_data.first ) return;

	// some convenient names for the two items held in the job_maturation_data class.
	LarvalJobOP & larval_job = job_maturation_data.first;
	utility::vector1< JobResultCOP > & job_results( job_maturation_data.second );

	try {

		bool communicated_w_master = false;
		// if the job's status as it exits is "fail retry", then we'll retry it.
		Size num_retries = larval_job->defines_retry_limit() ? larval_job->retry_limit() : default_retry_limit_;
		for ( core::Size ii = 1; ii <= num_retries; ++ii ) {

			JobOP job = job_queen_->mature_larval_job( larval_job, job_results );

			TR << "Starting job" << std::endl;
			CompletedJobOutput job_output = job->run();
			TR << "Job completed" << std::endl;

			if ( job_output.first->status() == jd3_job_status_failed_retry ) {
				continue;
			} else if ( job_output.first->status() == jd3_job_status_failed_do_not_retry ) {
				worker_send_fail_do_not_retry_message_to_master( larval_job );
				communicated_w_master = true;
				break;
			} else if ( job_output.first->status() == jd3_job_status_inputs_were_bad ) {
				worker_send_fail_bad_inputs_message_to_master( larval_job );
				communicated_w_master = true;
				break;
			} else if ( job_output.first->status() == jd3_job_status_success ) {
				worker_send_job_result_to_master_and_archive( larval_job, job_output );
				communicated_w_master = true;
				break;
			}
		}
		if ( ! communicated_w_master ) {
			worker_send_fail_retry_limit_exceeded( larval_job );
		}

	} catch ( std::ios_base::failure& ex ) {
		std::cerr << "std::IO error detected... exiting..." << std::endl; // We can not longer use Tracer's at this point
		// Using pure exit instead of utility_exit_with_status to avoid infinite recursion
		std::exit( basic::options::option[ basic::options::OptionKeys::out::std_IO_exit_error_code]() );
	} catch ( utility::excn::EXCN_BadInput & excn ) {
		worker_send_fail_on_bad_inputs_exception_message_to_master( larval_job, excn.msg() );
	} catch ( utility::excn::EXCN_Base & excn ) {
		worker_send_fail_w_message_to_master( larval_job, excn.msg() );
	} catch ( ... ) {
		worker_send_fail_w_message_to_master( larval_job, "Unhandled exception!" );
	}
}

void
MPIWorkPoolJobDistributor::retrieve_deallocation_messages_and_new_job_and_run()
{
	// ok, we've already received the deallocation message
	bool more_deallocation_messages_remaining( true );
	while ( more_deallocation_messages_remaining ) {
		std::list< deallocation::DeallocationMessageOP > messages;
		{ // scope
			std::string srlz_deallocation_msg = utility::receive_string_from_node( 0 );
			std::istringstream iss( srlz_deallocation_msg );
			cereal::BinaryInputArchive arc( iss );
			arc( messages );
		}

		for ( auto msg : messages ) {
			job_queen_->process_deallocation_message( msg );
		}

		core::Size next_mpi_msg = utility::receive_integer_from_node( 0 );
		if ( next_mpi_msg == mpi_work_pool_jd_new_job_available ) {
			more_deallocation_messages_remaining = false;
		} else {
			debug_assert( next_mpi_msg == mpi_work_pool_jd_deallocation_message );
		}
	}

	retrieve_job_and_run();
}


std::pair< LarvalJobOP, utility::vector1< JobResultCOP > >
MPIWorkPoolJobDistributor::retrieve_job_maturation_data()
{
	std::pair< LarvalJobOP, utility::vector1< JobResultCOP > > return_val;

	std::string larval_job_string = utility::receive_string_from_node( 0 );
	LarvalJobOP larval_job;
	try {
		larval_job = deserialize_larval_job( larval_job_string );
	} catch ( cereal::Exception & e ) {
		std::ostringstream oss;
		oss << "Failed to deserialize larval job on worker node " << mpi_rank_ << ". Exiting.\nError message from cereal library:\n";
		oss << e.what() << "\n";
		send_error_message_to_node0( oss.str() );
		return return_val;
	}

	utility::vector1< core::Size > archival_nodes_for_job_results = utility::receive_sizes_from_node( 0 );
	utility::vector1< JobResultCOP > job_results( larval_job->input_job_result_indices().size() );
	for ( core::Size ii = 1; ii <= job_results.size(); ++ii ) {
		// request job result #(larval_job->input_job_result_indices()[ ii ] )
		// from archival node #(archival_nodes_for_job_results[ii])
		core::Size archive_node = archival_nodes_for_job_results[ ii ];

		utility::send_integer_to_node( archive_node, mpi_rank_ ); // tell the other node who is talking ..
		utility::send_integer_to_node( archive_node, mpi_work_pool_jd_retrieve_job_result ); // .. and that we need a job result
		utility::send_size_to_node( archive_node, larval_job->input_job_result_indices()[ ii ] ); // .. of a particular job

		int retrieval_success = utility::receive_integer_from_node( archive_node ); // well, does the other node have that job?
		if ( retrieval_success == mpi_work_pool_jd_job_result_retrieved ) {
			std::string job_and_result_string = utility::receive_string_from_node( archive_node ); // good. Now we'll deserialize it
			LarvalJobAndResult job_and_result;
			try {
				job_and_result = deserialize_larval_job_and_result( job_and_result_string );
			} catch ( cereal::Exception & e ) {
				std::ostringstream oss;
				oss << "Failed to deserialize LarvalJob & JobResult pair from job " << larval_job->input_job_result_indices()[ ii ] <<
					" which is required as an input to job " << larval_job->job_index() << "\nError message from cereal library:\n" <<
					e.what() << "\n";
				send_error_message_to_node0( oss.str() );
				return return_val;
			}

			job_results[ ii ] = job_and_result.second;
		} else {
			// ok -- nothing more we're going to do here. We can't run this job.
			//  Node 0 is going to be notified and the job will spin down.
			// perhaps this function should return "false" to indicate that the job
			// could not be run.
			return return_val;
		}
	}
	return_val.first = larval_job;
	return_val.second = job_results;
	return return_val;
}

void
MPIWorkPoolJobDistributor::worker_send_fail_do_not_retry_message_to_master(
	LarvalJobOP larval_job
)
{
	// tell node 0 the job did not succeed.
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_job_failed_do_not_retry );
	utility::send_size_to_node( 0, larval_job->job_index() );
}

void
MPIWorkPoolJobDistributor::worker_send_fail_bad_inputs_message_to_master(
	LarvalJobOP larval_job
)
{
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_job_failed_bad_input );
	utility::send_size_to_node( 0, larval_job->job_index() );
}

void
MPIWorkPoolJobDistributor::worker_send_fail_on_bad_inputs_exception_message_to_master(
	LarvalJobOP larval_job,
	std::string const & error_message
)
{
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_job_failed_w_message );
	utility::send_size_to_node( 0, larval_job->job_index() );
	utility::send_string_to_node( 0, error_message );
}

void
MPIWorkPoolJobDistributor::worker_send_fail_w_message_to_master(
	LarvalJobOP larval_job,
	std::string const & error_message
)
{
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_job_failed_w_message );
	utility::send_size_to_node( 0, larval_job->job_index() );
	utility::send_string_to_node( 0, error_message );
}

void
MPIWorkPoolJobDistributor::worker_send_fail_retry_limit_exceeded(
	LarvalJobOP larval_job
)
{
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_job_failed_retry_limit_exceeded );
	utility::send_size_to_node( 0, larval_job->job_index() );
}

void
MPIWorkPoolJobDistributor::worker_send_job_result_to_master_and_archive(
	LarvalJobOP larval_job,
	CompletedJobOutput job_output
)
{
	std::string serialized_larval_job_and_result;
	try {
		serialized_larval_job_and_result = serialize_larval_job_and_result( std::make_pair( larval_job, job_output.second ));
	} catch ( cereal::Exception const & e ) {
		// send an error message to node 0 and return
		std::ostringstream oss;
		oss << "Failed to serialize LarvalJob or JobResult; all objects stored in these" <<
			" classes (or derived from them) must implement save & load serialization functions\nJob #" <<
			larval_job->job_index() << " failed. Error message from the cereal library:\n" << e.what() << "\n";
		send_error_message_to_node0( oss.str() );
		return;
	}

	std::string serialized_job_summary;
	try {
		serialized_job_summary = serialize_job_summary( job_output.first );
	} catch ( cereal::Exception const & e ) {
		// send an error message to node 0 and return
		std::ostringstream oss;
		oss << "Failed to serialize JobSummary; all objects stored in this class (or that are derived from" <<
			" it) must implement save & load serialization functions\nJob #" <<	larval_job->job_index() <<
			" failed. Error message from the cereal library:\n" << e.what() << "\n";
		send_error_message_to_node0( oss.str() );
		return;
	}

	// Ok -- everything that needs to be serialized can be; proceed to report a
	// successful job completion.

	// now, notify node 0 that the job is complete -- node 0 will tell us where
	// to send the JobResult for archival.
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_job_success );
	utility::send_size_to_node( 0, larval_job->job_index() );
	int archival_node = utility::receive_integer_from_node( 0 );
	if ( archival_node != 0 ) {
		utility::send_integer_to_node( archival_node, mpi_rank_ );
		utility::send_integer_to_node( archival_node, mpi_work_pool_jd_archive_job_result );
		utility::send_size_to_node( archival_node, larval_job->job_index() );
	}
	utility::send_string_to_node( archival_node, serialized_larval_job_and_result );

	// now wait to hear from the archival node before proceeding.
	// we need to ensure that the archival node has processed the archival event before
	// we send the job status to node 0.  At the time that we send the job status to node 0,
	// the JobResult's archival must have completed.  This is to ensure that no node
	// ever asks the archive for a JobResult that has not yet been archived (but that is
	// on its way to being archived).

	// int archival_completed_message =
	utility::receive_integer_from_node( archival_node );
	// debug_assert( archival_completed_message == mpi_work_pool_jd_archival_complete );

	// now go back to talk again with node 0
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
	if ( job_queen_->larval_job_needed_for_note_job_completed() ||
			job_queen_->larval_job_needed_for_completed_job_summary() ) {
		// Do we need a try/catch block around the call to serialize_larval_job?
		// No.
		// It's the same larval job that we serialized above, so if it serialized
		// before, it will serialize again.
		utility::send_string_to_node( 0, serialize_larval_job( larval_job ) );
	} else {
		utility::send_size_to_node( 0, larval_job->job_index() );
	}

	utility::send_string_to_node( 0, serialized_job_summary );
}

/// @throws Potentially throws a cereal::Exception; it is much more likely that errors are
/// hit during object serialization, rather than objet deserialization, because an object
/// that lacks the save/load routines is embedded in a LarvalJob or JobResult. However,
/// it is possible that the save routine is OK, but the load routine fails, in which case
/// the user ought to be told what was wrong.
LarvalJobOP
MPIWorkPoolJobDistributor::deserialize_larval_job(
	std::string const & larval_job_string
) const
{
	LarvalJobOP larval_job;
	std::istringstream iss( larval_job_string );
	cereal::BinaryInputArchive arc( iss );
	arc( larval_job );
	return larval_job;
}

/// @throws Potentially throws a cereal::Exception if the LarvalJob holds an
/// unserializable object or is itself unserializable
std::string
MPIWorkPoolJobDistributor::serialize_larval_job(
	LarvalJobOP larval_job
) const
{
	std::ostringstream larval_job_oss;
	{
		cereal::BinaryOutputArchive arc( larval_job_oss );
		arc( larval_job );
	}
	return larval_job_oss.str();
}

/// @throws Potentially throws a cereal::Exception; it is much more likely that errors are
/// hit during object serialization, rather than objet deserialization, because an object
/// that lacks the save/load routines is embedded in a LarvalJob or JobResult. However,
/// it is possible that the save routine is OK, but the load routine fails, in which case
/// the user ought to be told what was wrong.
MPIWorkPoolJobDistributor::LarvalJobAndResult
MPIWorkPoolJobDistributor::deserialize_larval_job_and_result(
	std::string const & job_and_result_string
) const
{
	LarvalJobAndResult job_and_result;
	std::istringstream iss( job_and_result_string );
	cereal::BinaryInputArchive arc( iss );
	arc( job_and_result );
	return job_and_result;
}

/// @throws Potentially throws a cereal::Exception if the LarvalJob or the JobResult holds an
/// unserializable object or either one is itself unserializable
std::string
MPIWorkPoolJobDistributor::serialize_larval_job_and_result(
	LarvalJobAndResult job_and_result
) const
{
	std::ostringstream oss;
	{
		cereal::BinaryOutputArchive arc( oss );
		arc( job_and_result );
	}
	return oss.str();
}

/// @throws Potentially throws a cereal::Exception; it is much more likely that errors are
/// hit during object serialization, rather than objet deserialization, because an object
/// that lacks the save/load routines is embedded in a LarvalJob or JobResult. However,
/// it is possible that the save routine is OK, but the load routine fails, in which case
/// the user ought to be told what was wrong.
JobSummaryOP
MPIWorkPoolJobDistributor::deserialize_job_summary(
	std::string const & job_summary_string
) const
{
	JobSummaryOP job_summary;
	std::istringstream iss( job_summary_string );
	cereal::BinaryInputArchive arc( iss );
	arc( job_summary );
	return job_summary;
}

/// @throws Potentially throws a cereal::Exception if the derived Summary holds an
/// unserializable object or is itself unserializable
std::string
MPIWorkPoolJobDistributor::serialize_job_summary(
	JobSummaryOP job_summary
) const
{
	std::ostringstream job_summary_oss;
	{
		cereal::BinaryOutputArchive arc( job_summary_oss );
		arc( job_summary );
	}
	return job_summary_oss.str();
}

}//job_distributors
}//jd3
}//protocols

#endif // SERIALIZATION
