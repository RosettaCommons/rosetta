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
#include <protocols/jd3/job_distributors/JobExtractor.hh>
#include <protocols/jd3/output/OutputSpecification.hh>
#include <protocols/jd3/output/ResultOutputter.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/assert.hh>
#include <utility/mpi_util.hh>
#include <utility/vector1.srlz.hh>
#include <utility/io/zipstream.hpp>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/details/helpers.hpp>

// C++
#include <utility>

static THREAD_LOCAL basic::Tracer TR( "protocols.jd3.job_distributors..MPIWorkPoolJobDistributor" );

namespace protocols {
namespace jd3 {
namespace job_distributors {

using core::Size;

MPIWorkPoolJobDistributor::MPIWorkPoolJobDistributor() :
	n_archives_( basic::options::option[ basic::options::OptionKeys::jd3::n_archive_nodes ] ),
	store_on_node0_( n_archives_ == 0 ||
	! basic::options::option[ basic::options::OptionKeys::jd3::do_not_archive_on_node0 ]() ),
	output_on_node0_( n_archives_ == 0 || ( store_on_node0_ &&
	! basic::options::option[ basic::options::OptionKeys::jd3::do_not_output_from_node0 ]() ) ),
	compress_job_results_( basic::options::option[ basic::options::OptionKeys::jd3::compress_job_results ]() ),
	archive_on_disk_( basic::options::option[ basic::options::OptionKeys::jd3::archive_on_disk ].user() ),
	archive_dir_name_( basic::options::option[ basic::options::OptionKeys::jd3::archive_on_disk ]() ),
	n_results_per_archive_( n_archives_ + 1 ),
	default_retry_limit_( 1 )
{
	mpi_rank_ = utility::mpi_rank();
	mpi_rank_string_ = utility::to_string( mpi_rank_ );
	mpi_nprocs_ = utility::mpi_nprocs();
	job_extractor_.reset( new JobExtractor );

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
	job_extractor_->set_job_queen( queen );
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
		case mpi_work_pool_jd_archive_job_result :
			process_archive_job_result_request( worker_node );
			break;
		case mpi_work_pool_jd_job_success_and_archival_complete :
			// second part of the two part process: after archival
			// has completed, the node sends the JobStatus to node 0
			// so that the JobQueen can examine it
			process_archival_complete_message( worker_node );
			break;
		case mpi_work_pool_jd_output_completed :
			process_output_complete_message( worker_node );
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
			throw utility::excn::EXCN_Msg_Exception( "Internal Error in the MPIWorkPoolJobDistributor: received inappropriate signal "
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
		default :
			// error -- we should not have gotten here
			throw utility::excn::EXCN_Msg_Exception( "received inappropriate signal "
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
		case mpi_work_pool_jd_output_job_result_already_available :
			process_output_job_result_already_available_request( remote_node );
			break;
		case mpi_work_pool_jd_accept_and_output_job_result :
			process_accept_and_output_job_result_request( remote_node );
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
				"Archival node " + utility::to_string( mpi_rank_ ) + " received an "
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
	job_dag_ = job_extractor_->create_initial_job_dag();

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
		bool job_sent = send_next_job_to_node( worker_node );
		if ( job_sent ) {
			return;
		}
	}

	if ( jobs_remain() ) {
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
	LarvalJobOP failed_job = job_extractor_->running_job( job_id );
	TR.Error << "Job " << job_id << " named " << failed_job <<
		" has exited with the message:\n" << error_message << std::endl;

	job_queen_->note_job_completed( failed_job, jd3_job_status_failed_w_exception, 1 /*dummy*/ );
	note_job_no_longer_running( job_id );
}

void
MPIWorkPoolJobDistributor::process_job_failed_do_not_retry( int worker_node )
{
	Size job_id = utility::receive_size_from_node( worker_node );
	LarvalJobOP failed_job = job_extractor_->running_job( job_id );
	job_queen_->note_job_completed( failed_job, jd3_job_status_failed_do_not_retry, 1 /*dummy*/ );
	note_job_no_longer_running( job_id );
}

void
MPIWorkPoolJobDistributor::process_job_failed_bad_input( int worker_node )
{
	Size job_id = utility::receive_size_from_node( worker_node );
	LarvalJobOP failed_job = job_extractor_->running_job( job_id );
	job_queen_->note_job_completed( failed_job, jd3_job_status_inputs_were_bad, 1 /*dummy*/ );

	// TO DO: add code that purges all other jobs with the same inner job as this one
	note_job_no_longer_running( job_id );
}

void
MPIWorkPoolJobDistributor::process_job_failed_retry_limit_exceeded( int worker_node )
{
	Size job_id = utility::receive_size_from_node( worker_node );
	LarvalJobOP failed_job = job_extractor_->running_job( job_id );
	job_queen_->note_job_completed( failed_job, jd3_job_status_failed_max_retries, 1 /*dummy*/ );

	// TO DO: add code that purges all other jobs with the same inner job as this one
	note_job_no_longer_running( job_id );
}

void
MPIWorkPoolJobDistributor::process_job_succeeded( int worker_node )
{

	Size job_id = utility::receive_size_from_node( worker_node );
	Size nresults = utility::receive_size_from_node( worker_node );

	// Now pick a node to archive the job result on, perhaps storing the
	// job result on this node.
	utility::vector1< int > archival_nodes_for_results( nresults );
	for ( Size ii = 1; ii <= nresults; ++ii ) {
		archival_nodes_for_results[ ii ] = pick_archival_node();
	}
	// tell the worker node where to send the job result
	utility::send_integers_to_node( worker_node, archival_nodes_for_results );

	for ( Size ii = 1; ii <= nresults; ++ii ) {
		JobResultID result_id = std::make_pair( job_id, ii );
		job_result_location_map_[ result_id ] = archival_nodes_for_results[ ii ];
	}

	// changing this if ( archival_node == 0 ) {
	// changing this  // NOTE: when worker nodes are sending the larval job / job result pair
	// changing this  // they do not need to first send an archival-request message or the job id.
	// changing this  // Store serialized larval job & job result pair
	// changing this  std::string serialized_job_result_string =
	// changing this   utility::receive_string_from_node( worker_node );
	// changing this  job_results_[ job_id ] = serialized_job_result_string;
	// changing this  utility::send_integer_to_node( worker_node, mpi_work_pool_jd_archival_completed );
	// changing this }
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
	// The second half of a job-completion step.
	//
	// At this point, the job results for the completed job
	// have been archived on the archival nodes (perhaps on this node)
	// so it is now safe to inform the JobQueen that the job has completed.
	// It is also safe for the JobQueen to request to discard
	// or output the results. If we were to structure the archival
	// process in a single pass manner (i.e. in the process_job_succeeded method),
	// the remote node where the archival is taking place might not have
	// completed -- or even gotten started -- by the time we ask it for the
	// result. This would produce a race condition.

	// We keep the larval jobs that are in flight, and will stop keeping them
	// at the end of this function
	Size job_id = utility::receive_size_from_node( worker_node );

	auto larval_job_iter = in_flight_larval_jobs_.find( job_id );
	debug_assert( larval_job_iter != in_flight_larval_jobs_.end() );
	LarvalJobOP larval_job = larval_job_iter->second;

	// Older implementation: do not store the in-flight jobs on the head node.
	// Is this a useful optimization? It's hard to imagine it is!
	// // The remote node will also query the JobQueen to know which of the two
	// // messages it should send.
	// if ( job_queen_->larval_job_needed_for_note_job_completed() ||
	//   job_queen_->larval_job_needed_for_completed_job_summary() ) {
	//  std::string serialized_larval_job_string =
	//   utility::receive_string_from_node( worker_node );
	//  // Does this call to deserialize_larval_job need a try/catch block?
	//  // No.
	//  // This larval job has been deserialized previously, with no opportunity
	//  // for new data to have been added to it.
	//  larval_job = deserialize_larval_job( serialized_larval_job_string );
	//  job_id = larval_job->job_index();
	// } else {
	//  job_id = utility::receive_size_from_node( worker_node );
	// }

	JobStatus status = JobStatus( utility::receive_integer_from_node( worker_node ) );
	std::string serialized_job_summaries_string = utility::receive_string_from_node( worker_node );
	utility::vector1< JobSummaryOP > job_summaries;
	try {
		job_summaries = deserialize_job_summaries( serialized_job_summaries_string );
	} catch ( cereal::Exception & e ) {
		throw utility::excn::EXCN_Msg_Exception( "Failed to deserialize the JobSummary array for job #" +
			utility::to_string( job_id ) + "\nError message from the cereal library:\n" + e.what() + "\n" );
	}

// Inform the JobQueen of the completed job.
// This step had to be delayed until after the archival process completed,
// as it assumes that the JobResult is present on the node where we
// wrote down it would be stored. Otherwise, jobs could be launched that
// would produce a race condition by asking for JobResults as inputs that
// either had or had not yet completed the archival process.
	Size const nsummaries = job_summaries.size();
	if ( job_queen_->larval_job_needed_for_note_job_completed() ||
			job_queen_->larval_job_needed_for_completed_job_summary() ) {
		job_queen_->note_job_completed( larval_job, status, nsummaries );
		for ( Size ii = 1; ii <= nsummaries; ++ii ) {
			job_queen_->completed_job_summary( larval_job, ii, job_summaries[ ii ] );
		}
	} else {
		job_queen_->note_job_completed( job_id, status, nsummaries );
		for ( Size ii = 1; ii <= nsummaries; ++ii ) {
			job_queen_->completed_job_summary( job_id, ii, job_summaries[ ii ] );
		}
	}

	note_job_no_longer_running( job_id );

	potentially_output_some_job_results();
	potentially_discard_some_job_results();
}

void
MPIWorkPoolJobDistributor::process_output_complete_message( int archival_node )
{
	// OK -- the archival node has told us that it is done outputting a job, so we can
	// decrement the number of jobs that we are saying it as doing. There are two ways
	// that the archival node may have been instructed to output a job:
	// 1. The JobResult was already being archived on that node, and now that the output
	// has been made, the archive node has discarded the result, so we decrement
	// the number of results it is holding, or
	// 2. The JobResult had previously been held on the master node, and the archive node
	// was selected to perform the output. We incremented the number of jobs it was holding
	// as a surrogate for saying the node is slightly more busy than it previously had been
	// and now that output is complete, we will decrement that count again.

	float new_result_count = n_results_per_archive_.coval_for_val( archival_node ) - 1;
	debug_assert( new_result_count >= 0 );
	n_results_per_archive_.heap_replace( archival_node, new_result_count );
}


void
MPIWorkPoolJobDistributor::process_retrieve_job_result_request( int worker_node )
{
	using namespace utility;
	Size job_id = receive_size_from_node( worker_node );
	Size result_index = receive_size_from_node( worker_node );
	JobResultID result_id = { job_id, result_index };

	bool archive_does_not_exist = false;

	if ( archive_on_disk_ ) {
		std::string const filename = filename_for_archive( job_id, result_index );
		std::ifstream in ( filename );
		if ( ! in.good() ) {
			archive_does_not_exist = true;
		} else {
			std::string const content = utility::file_contents( filename );
			send_integer_to_node( worker_node, mpi_work_pool_jd_job_result_retrieved );
			send_string_to_node( worker_node, content );
		}
		in.close();
	} else {
		SerializedJobResultMap::const_iterator iter = job_results_.find( result_id );
		if ( iter == job_results_.end() ) {
			archive_does_not_exist = true;
		} else {
			// ok -- we have found the (serialized) JobResult
			send_integer_to_node( worker_node, mpi_work_pool_jd_job_result_retrieved );
			send_string_to_node( worker_node, * iter->second );
			return;
		}
	}

	if ( archive_does_not_exist ) {
		// fatal error: the requested job result is not present at this node.
		// inform node 0 (if this isn't node 0 )
		if ( mpi_rank_ == 0 ) {
			throw_after_failed_to_retrieve_job_result( worker_node, job_id, result_index, mpi_rank_ );
		} else {
			send_integer_to_node( worker_node, mpi_work_pool_jd_failed_to_retrieve_job_result );
			// tell node 0 about failed retrieval
			send_integer_to_node( 0, mpi_rank_ );
			send_integer_to_node( 0, mpi_work_pool_jd_failed_to_retrieve_job_result );
			send_integer_to_node( 0, worker_node );
			send_size_to_node( 0, job_id );
			send_size_to_node( 0, result_index );
			return;
		}
	}
}

/// @details Called by the archive node. The archive node should already have the indicated job
/// result stored on this node.
void
MPIWorkPoolJobDistributor::process_output_job_result_already_available_request( int remote_node )
{
	// remote node should be node 0
	std::string output_spec_string = utility::receive_string_from_node( remote_node );
	output::OutputSpecificationOP spec;
	try {
		spec = deserialize_output_specification( output_spec_string );
	} catch ( cereal::Exception & e ) {
		utility::send_integer_to_node( 0, mpi_rank_ );
		utility::send_integer_to_node( 0, mpi_work_pool_jd_error );
		std::ostringstream oss;
		oss << "Error on archive node " << mpi_rank_ << " trying to deserialize output specification\n";
		oss << "Error message from cereal library:\n";
		oss << e.what() << "\n";
		utility::send_string_to_node( 0, oss.str() );
		// do not exit immediately as this will cause the MPI process to exit before
		// the head node has a chance to flush its output
		return;
	}
	JobResultID result_id = spec->result_id();

	std::string serialized_job_and_result;
	//bool archive_does_not_exist = false;
	if ( archive_on_disk_ ) {
		std::string const filename = filename_for_archive( result_id.first, result_id.second );
		std::ifstream in( filename );
		if ( ! in.good() ) {
			utility::send_integer_to_node( 0, mpi_rank_ );
			utility::send_integer_to_node( 0, mpi_work_pool_jd_error );
			std::ostringstream oss;
			oss << "Error on archive node " << mpi_rank_ << " trying to find job result archived on disk in file \"" << filename;
			oss << "\" for result id ( " << result_id.first << ", " << result_id.second << " )\n";
			utility::send_string_to_node( 0, oss.str() );
			// do not exit immediately as this will cause the MPI process to exit before
			// the head node has a chance to flush its output
			return;
		} else {
			serialized_job_and_result = utility::file_contents( filename );
			remove( filename.c_str() );
		}
		in.close();
	} else {
		auto job_and_result_iter = job_results_.find( result_id );
		if ( job_and_result_iter == job_results_.end() ) {
			utility::send_integer_to_node( 0, mpi_rank_ );
			utility::send_integer_to_node( 0, mpi_work_pool_jd_error );
			std::ostringstream oss;
			oss << "Error on archive node " << mpi_rank_ << " trying to find job result ";
			oss << "for result id ( " << result_id.first << ", " << result_id.second << " )\n";
			oss << "Job should be present in memory, but is not there. Has is already been output? ";
			oss << "Has it already been discarded?\n";
			utility::send_string_to_node( 0, oss.str() );
			// do not exit immediately as this will cause the MPI process to exit before
			// the head node has a chance to flush its output
			return;
		}
		serialized_job_and_result = std::move( *job_and_result_iter->second );//move instead of copy
		job_results_.erase( job_and_result_iter );
	}

	LarvalJobAndResult job_and_result;
	try {
		job_and_result = deserialize_larval_job_and_result( serialized_job_and_result );
	} catch ( cereal::Exception & e ) {
		utility::send_integer_to_node( 0, mpi_rank_ );
		utility::send_integer_to_node( 0, mpi_work_pool_jd_error );
		std::ostringstream oss;
		oss << "Error on archive node " << mpi_rank_ << " trying to deserialize larval job and job result";
		oss << " for result id ( " << result_id.first << ", " << result_id.second << " )\n";
		oss << "Error message from cereal library:\n";
		oss << e.what() << "\n";
		utility::send_string_to_node( 0, oss.str() );
		// do not exit immediately as this will cause the MPI process to exit before
		// the head node has a chance to flush its output
		return;
	}

	write_output( spec, job_and_result );

	// now tell node0 we are done
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_output_completed );
}

/// @details Called by the archive node. The archive node should already have the indicated job
/// result stored on this node.
void
MPIWorkPoolJobDistributor::process_accept_and_output_job_result_request( int remote_node )
{
	// remote node should be node 0
	std::string output_spec_string = utility::receive_string_from_node( remote_node );
	output::OutputSpecificationOP spec;
	try {
		spec = deserialize_output_specification( output_spec_string );
	} catch ( cereal::Exception & e ) {
		utility::send_integer_to_node( 0, mpi_rank_ );
		utility::send_integer_to_node( 0, mpi_work_pool_jd_error );
		std::ostringstream oss;
		oss << "Error on archive node " << mpi_rank_ << " trying to deserialize output specification\n";
		oss << "Error message from cereal library:\n";
		oss << e.what() << "\n";
		utility::send_string_to_node( 0, oss.str() );
		// do not exit immediately as this will cause the MPI process to exit before
		// the head node has a chance to flush its output
		return;

	}
	LarvalJobAndResult job_and_result;
	std::string job_and_result_string = utility::receive_string_from_node( remote_node );
	try {
		job_and_result = deserialize_larval_job_and_result( job_and_result_string );
	} catch ( cereal::Exception & e ) {
		utility::send_integer_to_node( 0, mpi_rank_ );
		utility::send_integer_to_node( 0, mpi_work_pool_jd_error );
		std::ostringstream oss;
		oss << "Error on archive node " << mpi_rank_ << " trying to deserialize larval job and job result";
		oss << " for result id ( " << spec->result_id().first << ", " << spec->result_id().second << " )\n";
		oss << "Error message from cereal library:\n";
		oss << e.what() << "\n";
		utility::send_string_to_node( 0, oss.str() );
		// do not exit immediately as this will cause the MPI process to exit before
		// the head node has a chance to flush its output
		return;
	}

	write_output( spec, job_and_result );

	// now tell node0 we are done
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_output_completed );
}


void
MPIWorkPoolJobDistributor::process_failed_to_retrieve_job_result_request( int archival_node )
{
	// the node making the job-result request, which the archival node was
	// unable to fullfil
	int worker_node = utility::receive_integer_from_node( archival_node );
	Size job_id = utility::receive_size_from_node( archival_node );
	Size result_index = utility::receive_size_from_node( archival_node );
	throw_after_failed_to_retrieve_job_result( worker_node, job_id, result_index, archival_node );
}

/// @details This function unlike the process_retireve_job_result_request function will only
/// be called on an archive node -- and it will be called only by the master node.
/// therefore, we do not need to handle the case when mpi_rank_ == 0.
void
MPIWorkPoolJobDistributor::process_retrieve_and_discard_job_result_request( int worker_node )
{
	using namespace utility;
	Size job_id = receive_size_from_node( worker_node );
	Size result_index = receive_size_from_node( worker_node );

	JobResultID job_result_id = { job_id, result_index };

	if ( archive_on_disk_ ) {
		std::string const filename = filename_for_archive( job_id, result_index );
		std::ifstream in ( filename );//as of now, this only checks to see if filename is good
		if ( ! in.good() ) {
			send_integer_to_node( worker_node, mpi_work_pool_jd_failed_to_retrieve_job_result );
		} else {
			std::string const content = utility::file_contents( filename );
			remove( filename.c_str() );

			send_integer_to_node( worker_node, mpi_work_pool_jd_job_result_retrieved );
			send_string_to_node( worker_node, content );
		}
		in.close();

	} else {
		SerializedJobResultMap::iterator iter = job_results_.find( job_result_id );
		if ( iter == job_results_.end() ) {
			// fatal error: the requested job result is not present at this node.
			// inform node 0 (this isn't node 0, as this function is only called by archive nodes)
			send_integer_to_node( worker_node, mpi_work_pool_jd_failed_to_retrieve_job_result );
			return;
		}

		// ok -- we have found the (serialized) JobResult
		send_integer_to_node( worker_node, mpi_work_pool_jd_job_result_retrieved );
		send_string_to_node( worker_node, * iter->second );
		job_results_.erase( iter );
	}
}

void
MPIWorkPoolJobDistributor::process_discard_job_result_request( int remote_node )
{
	using namespace utility;
	Size job_id = receive_size_from_node( remote_node ); // node 0
	Size result_index = receive_size_from_node( remote_node );

	if ( archive_on_disk_ ) {
		std::string const filename = filename_for_archive( job_id, result_index );
		if ( remove( filename.c_str() ) ) { //remove's name is not intuitive for its return value. non-zero return value means the remove failed
			std::cerr << "Failed to discard job result (" << job_id << ", " << result_index << ") on node " << mpi_rank_ << std::endl;
		}
		return;
	}

	JobResultID job_result_id = { job_id, result_index };
	SerializedJobResultMap::iterator iter = job_results_.find( job_result_id );
	if ( iter == job_results_.end() ) {
		// well, something has gone wrong, so we should probably print an error, but
		// technically, we can keep going.
		std::cerr << "Failed to discard job result (" << job_result_id.first << ", " << job_result_id.second << ") on node " << mpi_rank_ << std::endl;
	}
	job_results_.erase( iter );
}

void
MPIWorkPoolJobDistributor::process_archive_job_result_request( int remote_node )
{
	using namespace utility;
	Size job_id = receive_size_from_node( remote_node );
	Size result_index = receive_size_from_node( remote_node );
	JobResultID result_id = { job_id, result_index };
	std::string serialized_larval_job_and_result = receive_string_from_node( remote_node );
	if ( archive_on_disk_ ) {
		std::string const filename = filename_for_archive( job_id, result_index );
		std::ofstream out ( filename );
		debug_assert( out.is_open() );
		out << serialized_larval_job_and_result;
		out.close();
	} else {
		job_results_[ result_id ] = std::make_shared< std::string > ( std::move( serialized_larval_job_and_result ) );
	}
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
bool
MPIWorkPoolJobDistributor::send_next_job_to_node( int worker_node )
{
	debug_assert( ! job_extractor_->job_queue_empty() );

	if ( ! deallocation_messages_for_node_[ worker_node ].empty() ) {
		send_deallocation_messages_to_node( worker_node );
	}

	// First tell the worker node that we do indeed have a job for them to work on/
	// The other option is that we are going to tell the worker to spin down.
	utility::send_integer_to_node( worker_node, mpi_work_pool_jd_new_job_available );

	// NOTE: Popping jobs_for_current_digraph_node_: we might violate the
	// invariant that this queue is never empty so long as nodes remain in the
	// digraph_nodes_ready_to_be_run_ queue.  Address that below
	LarvalJobOP larval_job;
	while ( ! job_extractor_->job_queue_empty() ) {

		larval_job = job_extractor_->pop_job_from_queue();
		if ( job_queen_->has_job_completed( larval_job ) ) {
			job_queen_->note_job_completed( larval_job, jd3_job_previously_executed, 0 );
			job_extractor_->note_job_no_longer_running( larval_job->job_index() );
			larval_job = nullptr;
			continue;
		}
		debug_assert( in_flight_larval_jobs_.count( larval_job->job_index() ) == 0 );
		in_flight_larval_jobs_[ larval_job->job_index() ] = larval_job;
		break;
	}

	if ( ! larval_job ) {
		// i.e. the job queue is empty
		return false;
	}

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

	worker_node_for_job_[ larval_job->job_index() ] = worker_node;
	return true;
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
	worker_node_for_job_.erase( job_id );
	job_extractor_->note_job_no_longer_running( job_id );
	in_flight_larval_jobs_.erase( job_id );

	if ( job_extractor_->retrieve_and_reset_node_recently_completed() ) {
		std::list< deallocation::DeallocationMessageOP > messages = job_queen_->deallocation_messages();
		if ( ! messages.empty() ) {
			store_deallocation_messages( messages );
		}
		if ( ! worker_nodes_waiting_for_jobs_.empty() && ! job_extractor_->job_queue_empty() ) {
			assign_jobs_to_idling_nodes();
		}
	}
}

void
MPIWorkPoolJobDistributor::potentially_output_some_job_results()
{

	// ask the job queen if she wants to output any results
	std::list< output::OutputSpecificationOP > jobs_to_output =
		job_queen_->jobs_that_should_be_output();
	for ( auto spec : jobs_to_output ) {
		JobResultID result_id = spec->result_id();
		auto location_map_iter = job_result_location_map_.find( result_id );
		if ( location_map_iter == job_result_location_map_.end() ) {
			// TO DO: better diagnostic message here
			throw utility::excn::EXCN_Msg_Exception( "Failed to find job result (" +
				utility::to_string( result_id.first ) + ", " + utility::to_string( result_id.second ) +
				") for outputting as requested by the JobQeen. Has this job already been" +
				" output or discarded?" );
		}

		// TO DO: instead, wait until after the archive node has completed a job output,
		// and then erase its location from the map; this will merely require that we
		// keep track of which node is outputting which job. For now, mark the job as having
		// been erased.
		Size node_housing_result = location_map_iter->second;
		job_result_location_map_.erase( location_map_iter );


		LarvalJobAndResult job_and_result;
		if ( node_housing_result == 0 ) {
			// ok, we have the result on this node -- we optionally could ship the result
			// to another node and output it there.

			std::string serialized_job_and_result;
			if ( archive_on_disk_ ) {
				serialized_job_and_result = get_string_from_disk( result_id.first, result_id.second );
				if ( serialized_job_and_result.empty() ) {
					// OK -- exit
					std::ostringstream oss;
					oss << "Error trying to retrieve a serialized job result from disk for job result (";
					oss << result_id.first << ", " << result_id.second << " ) which should have been in ";
					oss << "\"" << filename_for_archive( result_id.first, result_id.second ) << "\"\n";
					throw utility::excn::EXCN_Msg_Exception( oss.str() );
				}
				remove( filename_for_archive( result_id.first, result_id.second ).c_str() );
			} else {
				auto job_and_result_iter = job_results_.find( result_id );
				if ( job_and_result_iter == job_results_.end() ) {
					std::ostringstream oss;
					oss << "Error trying to find job result ( " << result_id.first << ", ";
					oss << result_id.second << " ); node 0 thinks it ought to have this result.\n";
					throw utility::excn::EXCN_Msg_Exception( oss.str() );
				}
				serialized_job_and_result = std::move( *job_and_result_iter->second );
				job_results_.erase( job_and_result_iter );
			}

			if ( output_on_node0_ ) {
				try {
					job_and_result = deserialize_larval_job_and_result( serialized_job_and_result );
				} catch ( cereal::Exception & e ) {
					throw utility::excn::EXCN_Msg_Exception( "Failed to deserialize LarvalJob & JobResult pair for job #" +
						utility::to_string( result_id.first ) + " result index #" + utility::to_string( result_id.second ) +
						"\nError message from the cereal library:\n" + e.what() + "\n" );
				}
				write_output( spec, job_and_result );
			} else {
				// Pick an archival node to perform the output -- this increments the
				// number of results that we're saying this node is holding, and when the node
				// reports that it is done outputting the result, then we will decrement that count
				// again.

				Size archive_node = pick_archival_node();
				utility::send_integer_to_node( archive_node, 0 ); // say "I need to talk to you"
				utility::send_integer_to_node( archive_node,
					mpi_work_pool_jd_accept_and_output_job_result );
				utility::send_string_to_node( archive_node, serialize_output_specification( spec ) );

				utility::send_string_to_node( archive_node, serialized_job_and_result );
			}

		} else {
			utility::send_integer_to_node( node_housing_result, 0 ); // say "I need to talk to you"
			utility::send_integer_to_node( node_housing_result,
				mpi_work_pool_jd_output_job_result_already_available );
			utility::send_string_to_node( node_housing_result,
				serialize_output_specification( spec ) );
		}

		// and now erase the record of where the job had been stored, and,
		// for the sake of diagnosing errors, also store the fact that this job
		// had been output (as opposed to discarded) -- but compactly within the
		// output_jobs_ DIET -- (Note some day, I'll figure out how to do this efficiently.
		// For now, the JD does not track which jobs have already been output.)
		job_result_location_map_.erase( result_id );

		// output_jobs_.insert( result_id.first );
	}
}

void
MPIWorkPoolJobDistributor::potentially_discard_some_job_results()
{
	// ask the job queen if she wants to discard any results
	std::list< JobResultID >  jobs_to_discard = job_queen_->job_results_that_should_be_discarded();
	if ( jobs_to_discard.empty() ) return;

	utility::vector0< int > n_discarded_jobs( n_archives_ + 1, 0 );

	for ( JobResultID result_id : jobs_to_discard ) {
		if ( job_result_location_map_.count( result_id ) == 0 ) {
			// TO DO: better diagnostic message here
			throw utility::excn::EXCN_Msg_Exception( "Failed to find job result " +
				utility::to_string( result_id.first ) + ", " + utility::to_string( result_id.second ) +
				" for discarding as requested by the JobQeen. Has this job already been" +
				" output or discarded?" );
		}

		Size node_housing_result = job_result_location_map_[ result_id ];
		n_discarded_jobs[ node_housing_result ] += 1;
		LarvalJobAndResult job_and_result;
		if ( node_housing_result == 0 ) {
			// ok, we have the result on this node
			if ( archive_on_disk_ ) {
				remove( filename_for_archive( result_id.first, result_id.second ).c_str() );
			} else {
				job_results_.erase( result_id );
			}
		} else {
			utility::send_integer_to_node( node_housing_result, 0 ); // say "I need to talk to you"
			utility::send_integer_to_node( node_housing_result,
				mpi_work_pool_jd_discard_job_result );
			utility::send_size_to_node( node_housing_result, result_id.first );
			utility::send_size_to_node( node_housing_result, result_id.second );
		}

		// and now delete storage for this job result, and
		// for the sake of diagnosing errors, also store the fact that this job
		// had been discarded (as opposed to being output) -- but compactly within the
		// discarded_jobs_ DIET
		job_result_location_map_.erase( result_id );

		// discarded_jobs_.insert( job_id );

	}

	// update the n_results_per_archive_ heap
	for ( int ii = 0; ii <= n_archives_; ++ii ) {
		if ( n_discarded_jobs[ ii ] == 0 ) continue;
		float new_n_results_for_archive_ii = n_results_per_archive_.coval_for_val( ii )
			- n_discarded_jobs[ ii ];
		n_results_per_archive_.heap_replace( ii, new_n_results_for_archive_ii );
	}
}

//void
//MPIWorkPoolJobDistributor::queue_jobs_for_next_node_to_run()
//{
// debug_assert( jobs_for_current_digraph_node_.empty() );
// debug_assert( ! digraph_nodes_ready_to_be_run_.empty() );
//
// job_extractor_->queue_jobs_for_next_node_to_run();
//
// if ( ! worker_nodes_waiting_for_jobs_.empty() &&
//   ! job_extractor_->jobs_for_current_digraph_node_empty() ) {
//  assign_jobs_to_idling_nodes();
//  // at the end of this call, either we have run out of jobs to assign
//  // or we have run out of idling nodes
// }
//
//}

//void
//MPIWorkPoolJobDistributor::queue_initial_digraph_nodes_and_jobs()
//{
// runtime_assert( digraph_is_a_DAG( *job_dag_ ) );
//
// // iterate across all nodes in the job_digraph, and note every
// // node with an indegree of 0
// for ( Size ii = 1; ii <= job_dag_->num_nodes(); ++ii ) {
//  if ( job_dag_->get_node(ii)->indegree() == 0 ) {
//   digraph_nodes_ready_to_be_run_.push_back( ii );
//  }
// }
// debug_assert( ! digraph_nodes_ready_to_be_run_.empty() );
// queue_jobs_for_next_node_to_run();
//}

bool
MPIWorkPoolJobDistributor::not_done()
{
	// The JobExtractor may have found a job that's ready to run, so, if there are idling nodes
	// launch that job.
	bool we_are_not_done = job_extractor_->not_done();
	if ( ! worker_nodes_waiting_for_jobs_.empty() && ! job_extractor_->job_queue_empty() ) {
		assign_jobs_to_idling_nodes();
	}
	return we_are_not_done;
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
	return ! job_extractor_->job_queue_empty();
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
	return ! job_extractor_->complete();
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

void
MPIWorkPoolJobDistributor::assign_jobs_to_idling_nodes()
{
	while ( ! worker_nodes_waiting_for_jobs_.empty() && ! job_extractor_->job_queue_empty() ) {
		Size worker_node = worker_nodes_waiting_for_jobs_.front();
		worker_nodes_waiting_for_jobs_.pop_front();
		bool job_sent = send_next_job_to_node( worker_node );
		if ( ! job_sent ) {
			// put this node back into the queue of nodes waiting for jobs
			worker_nodes_waiting_for_jobs_.push_back( worker_node );
		}
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
	int worker_node,
	Size job_id,
	Size result_index,
	int archival_node
)
{
	std::ostringstream oss;
	oss << "Failed to retrieve job result (" << job_id << ", " << result_index << ") which was";
	oss << " requested from node " << archival_node << " by node " << worker_node;
	oss << " but was not present there.\n";

	JobResultID result_id = { job_id, result_index };

	if ( worker_node_for_job_.count( job_id ) ) {
		oss << "Internal Error in the MPIWorkPoolJobDistributor: JobDistrutor thinks that this job is still running on node " <<
			worker_node_for_job_[ job_id ] << "\n";
	} else if ( job_result_location_map_.count( result_id ) ) {
		oss << "JobDistributor on node 0 thinks the result should have been on node ";
		oss << job_result_location_map_[ result_id ] << "\n";
		//} else if ( output_jobs_.member( job_id ) ) {
		// oss << "JobDistributor has already output this job (and does not keep jobs ";
		// oss << "that have already been output).\n";
		//} else if ( discarded_jobs_.member( job_id ) ) {
		// oss << "JobDistributor has discarded this job per the JobQueen's request\n";
	} else {
		oss << "This job is not listed as still running nor as having its JobResult stored on any"
			" archive node; it perhaps has been output or discarded already.\n";
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

			if ( job_output.status == jd3_job_status_failed_retry ) {
				continue;
			} else if ( job_output.status == jd3_job_status_failed_do_not_retry ) {
				worker_send_fail_do_not_retry_message_to_master( larval_job );
				communicated_w_master = true;
				break;
			} else if ( job_output.status == jd3_job_status_inputs_were_bad ) {
				worker_send_fail_bad_inputs_message_to_master( larval_job );
				communicated_w_master = true;
				break;
			} else if ( job_output.status == jd3_job_status_success ) {
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
		utility::send_size_to_node( archive_node, larval_job->input_job_result_indices()[ ii ].first ); // .. of a particular job
		utility::send_size_to_node( archive_node, larval_job->input_job_result_indices()[ ii ].second ); // .. and a particular result from that job

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
	Size const nresults = job_output.job_results.size();
	utility::vector1< std::string > serialized_larval_job_and_results( nresults );
	try {
		for ( Size ii = 1; ii <= job_output.job_results.size(); ++ii ) {
			serialized_larval_job_and_results[ ii ] = serialize_larval_job_and_result(
				std::make_pair( larval_job, job_output.job_results[ ii ].second ));
		}
	} catch ( cereal::Exception const & e ) {
		// send an error message to node 0 and return
		std::ostringstream oss;
		oss << "Failed to serialize LarvalJob or JobResult; all objects stored in these" <<
			" classes (or derived from them) must implement save & load serialization functions\nJob #" <<
			larval_job->job_index() << " failed. Error message from the cereal library:\n" << e.what() << "\n";
		send_error_message_to_node0( oss.str() );
		return;
	}

	std::string serialized_job_summaries;
	utility::vector1< JobSummaryOP > job_summaries( nresults );
	std::transform( job_output.job_results.begin(), job_output.job_results.end(),
		job_summaries.begin(), []( SingleJobOutputPair const & summ_and_result ) { return summ_and_result.first; } );
	try {
		serialized_job_summaries = serialize_job_summaries( job_summaries );
	} catch ( cereal::Exception const & e ) {
		// send an error message to node 0 and return
		std::ostringstream oss;
		oss << "Failed to serialize JobSummary; all objects stored in this class (or that are derived from" <<
			" it) must implement save & load serialization functions\nJob #" << larval_job->job_index() <<
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
	utility::send_size_to_node( 0, nresults );
	utility::vector1< int > archival_nodes = utility::receive_integers_from_node( 0 );
	for ( Size ii = 1; ii <= nresults; ++ii ) {
		int archival_node = archival_nodes[ ii ];
		utility::send_integer_to_node( archival_node, mpi_rank_ );
		utility::send_integer_to_node( archival_node, mpi_work_pool_jd_archive_job_result );
		utility::send_size_to_node( archival_node, larval_job->job_index() );
		utility::send_size_to_node( archival_node, ii );
		utility::send_string_to_node( archival_node, serialized_larval_job_and_results[ ii ] );

		// now wait to hear from the archival node before proceeding.
		// we need to ensure that the archival nodes have processed the archival events before
		// we send the job status to node 0. Technically, we could proceed to the next archival node,
		// if we were communicating with a different archive for the next job.
		// In any case,  by the time that we send the job status to node 0,
		// the JobResult's archival must have completed to ensure that no node
		// ever asks the archive for a JobResult that has not yet been archived (but that is
		// on its way to being archived). This avoids a race condition where either the
		// archival completes first and the next job that needs that JobResult as an input
		// runs properly, or the job requests the not-yet-archived JobResult, cannot find it,
		// and exits. That'd be a pretty bad race condition to have!
		int archival_completed_message = utility::receive_integer_from_node( archival_node );
		// debug_assert( archival_completed_message == mpi_work_pool_jd_archival_completed );
		if ( archival_completed_message != mpi_work_pool_jd_archival_completed ) {
			// well, damn
			std::ostringstream oss;
			oss << "Failed to archive job result #" << ii << " on archival node " << archival_node << " for job " << larval_job->job_index() << "; exiting" << std::endl;
			send_error_message_to_node0( oss.str() );
			return;
		}
	}

	// now go back to talk again with node 0
	utility::send_integer_to_node( 0, mpi_rank_ );
	utility::send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );

	// Old Implementation: Node0 does not keep the in-flight LarvalJobOPs, so we have to send
	// it back.
	// if ( job_queen_->larval_job_needed_for_note_job_completed() ||
	//   job_queen_->larval_job_needed_for_completed_job_summary() ) {
	//  // Do we need a try/catch block around the call to serialize_larval_job?
	//  // No.
	//  // It's the same larval job that we serialized above, so if it serialized
	//  // before, it will serialize again.
	//  utility::send_string_to_node( 0, serialize_larval_job( larval_job ) );
	// } else {
	//  utility::send_size_to_node( 0, larval_job->job_index() );
	// }

	utility::send_size_to_node( 0, larval_job->job_index() );

	utility::send_integer_to_node( 0, job_output.status );
	utility::send_string_to_node( 0, serialized_job_summaries );
}

void
MPIWorkPoolJobDistributor::write_output(
	output::OutputSpecificationOP spec,
	LarvalJobAndResult const & job_and_result
)
{
	// This is a two step process. Why? Because the MPIMultithreadedJobDistributor is going
	// to deserializize the job results and perform output in separate threads, so the
	// interaction with the outputters will need to be thread safe. (The JD still guarantees
	// that only one thread will interact with the JobQueen herself at any one point in time).
	if ( n_archives_ > 1 || ( n_archives_ == 1 && ( ! store_on_node0_ || output_on_node0_ ) ) ) {
		TR << "Setting output suffix: " << mpi_rank_string_ << " n_archives: " << n_archives_ << " store_on_node0_ " << store_on_node0_ << " output_on_node0_ " << output_on_node0_ << std::endl;

		spec->jd_output_suffix( mpi_rank_string_ );
	}
	output::ResultOutputterOP outputter = job_queen_->result_outputter( *spec );
	outputter->write_output( *spec, *job_and_result.second );
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
	if ( compress_job_results_ ) {
		std::istringstream compressed_ss( job_and_result_string );
		zlib_stream::zip_istream unzipper( compressed_ss );
		cereal::BinaryInputArchive arc( unzipper );
		arc( job_and_result );
	} else {
		std::istringstream iss( job_and_result_string );
		cereal::BinaryInputArchive arc( iss );
		arc( job_and_result );
	}
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
	std::string const uncompressed_string = oss.str();

	if ( compress_job_results_ ) {
		std::ostringstream compressed_ss;
		zlib_stream::zip_ostream zipper( compressed_ss );
		zipper << uncompressed_string;
		zipper.zflush();
		return compressed_ss.str();
	} else {
		return uncompressed_string;
	}
}

/// @throws Potentially throws a cereal::Exception; it is much more likely that errors are
/// hit during object serialization, rather than objet deserialization, because an object
/// that lacks the save/load routines is embedded in a LarvalJob or JobResult. However,
/// it is possible that the save routine is OK, but the load routine fails, in which case
/// the user ought to be told what was wrong.
utility::vector1< JobSummaryOP >
MPIWorkPoolJobDistributor::deserialize_job_summaries(
	std::string const & job_summaries_string
) const
{
	utility::vector1< JobSummaryOP > job_summaries;
	std::istringstream iss( job_summaries_string );
	cereal::BinaryInputArchive arc( iss );
	arc( job_summaries );
	return job_summaries;
}

/// @throws Potentially throws a cereal::Exception if the derived Summary holds an
/// unserializable object or is itself unserializable
std::string
MPIWorkPoolJobDistributor::serialize_job_summaries(
	utility::vector1< JobSummaryOP > const & job_summaries
) const
{
	std::ostringstream job_summaries_oss;
	{
		cereal::BinaryOutputArchive arc( job_summaries_oss );
		arc( job_summaries );
	}
	return job_summaries_oss.str();
}

output::OutputSpecificationOP
MPIWorkPoolJobDistributor::deserialize_output_specification( std::string const & spec_string ) const
{
	output::OutputSpecificationOP spec;
	std::istringstream iss( spec_string );
	cereal::BinaryInputArchive arc( iss );
	arc( spec );
	return spec;

}

std::string
MPIWorkPoolJobDistributor::serialize_output_specification( output::OutputSpecificationOP output_spec ) const
{
	std::ostringstream spec_oss;
	{
		cereal::BinaryOutputArchive arc( spec_oss );
		arc( output_spec );
	}
	return spec_oss.str();
}

std::string MPIWorkPoolJobDistributor::get_string_from_disk( Size job_id, Size result_index ) const {
	return utility::file_contents( filename_for_archive( job_id, result_index ) );
}

std::string MPIWorkPoolJobDistributor::filename_for_archive( Size job_id, Size result_index ) const {
	return archive_dir_name_ + "/archive." + std::to_string( job_id ) + "." + std::to_string( result_index );
}


}//job_distributors
}//jd3
}//protocols

#endif // SERIALIZATION
