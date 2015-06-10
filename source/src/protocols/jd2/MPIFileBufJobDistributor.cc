// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MPIFileBufJobDistributor.cc
/// @brief  implementation of MPIFileBufJobDistributor
/// @author Oliver Lange olange@u.washington.edu
/// @detail freely based on the MPIWorkPoolJobDistributor from Doug

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

// Unit headers
#include <protocols/jd2/MPIFileBufJobDistributor.hh>

// Package headers
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/moves/Mover.hh>

#include <protocols/jd2/MpiFileBuffer.hh>
#include <utility/io/ozstream.hh> //to toggle MPI rerouting

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/exit.hh>

// Option headers
#include <basic/options/keys/out.OptionKeys.gen.hh>
#ifdef USEMPI
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#endif

// C++ headers
#include <string>
#include <basic/prof.hh>
// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.jd2.MPIFileBufJobDistributor" );

namespace protocols {
namespace jd2 {

using namespace core;


using namespace basic::options;
using namespace basic::options::OptionKeys;

/// @details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
///which processor it is in MPI land.
MPIFileBufJobDistributor::MPIFileBufJobDistributor() :
	JobDistributor(),
	n_rank_( 1 ),
	rank_( 0 ),
	slave_current_job_id_( 0 ),
	slave_current_batch_id_( 0 ),
	//	next_job_to_assign_( 0 ),
	bad_job_id_( 0 ),
	repeat_job_( false ),
	master_rank_( 1 ),
	file_buf_rank_( 0 ),
	min_client_rank_( 2 ),
	cumulated_runtime_( 0.0 ),
	cumulated_jobs_( 0 )
{

  // set n_rank_ and rank based on whether we are using MPI or not
#ifdef USEMPI
	int int_rank, int_n_rank;                         //don't cast pointers - copy it over instead
	MPI_Comm_rank( MPI_COMM_WORLD, &int_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &int_n_rank );
	rank_ = int_rank;
	n_rank_ = int_n_rank;
	n_worker_ = n_rank_ - min_client_rank_;
#else
	utility_exit_with_message( "ERROR ERROR ERROR: The MPIFileBufJobDistributor will not work unless you have compiled using extras=mpi" );
#endif
}

/// @details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
///which processor it is in MPI land.
MPIFileBufJobDistributor::MPIFileBufJobDistributor(
	core::Size master_rank,
	core::Size file_buf_rank,
	core::Size min_client_rank,
	bool start_empty
) :
	JobDistributor( start_empty /*call empty c'tor*/ ),
	n_rank_( 1 ),
	rank_( 0 ),
	slave_current_job_id_( 0 ),
	slave_current_batch_id_( 0 ),
	//	next_job_to_assign_( 0 ),
	bad_job_id_( 0 ),
	repeat_job_( false ),
	master_rank_( master_rank ),
	file_buf_rank_(file_buf_rank ),
	min_client_rank_( min_client_rank ),
	cumulated_runtime_( 0.0 ),
	cumulated_jobs_( 0 )
{

	// set n_rank_ and rank based on whether we are using MPI or not
#ifdef USEMPI
	//n_rank_ = MPI::COMM_WORLD.Get_size();
	//rank_ = MPI::COMM_WORLD.Get_rank();
	int int_rank, int_n_rank;                         //don't cast pointers - copy it over instead
  MPI_Comm_rank( MPI_COMM_WORLD, &int_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &int_n_rank );
  rank_ = int_rank;
  n_rank_ = int_n_rank;
  n_worker_ = n_rank_ - min_client_rank_;
#else
	utility_exit_with_message( "ERROR ERROR ERROR: The MPIFileBufJobDistributor will not work unless you have compiled using extras=mpi" );
#endif
}

/// @brief dtor
///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
MPIFileBufJobDistributor::~MPIFileBufJobDistributor()
{}

void
MPIFileBufJobDistributor::go( protocols::moves::MoverOP mover ) {
	utility::io::ozstream::enable_MPI_reroute( min_client_rank_, file_buf_rank_ );
	protocols::jd2::WriteOut_MpiFileBuffer buffer( file_buf_rank_ );
	buffer.run(); //returns immediately if not buffer_rank
	tr.Debug << "finished call to buffer.run()" << std::endl;
	if ( rank_ == master_rank_ ) {
		tr.Debug << "Master JD starts" << std::endl;
		master_go( mover );
		tr.Debug << "send STOP to FileBuffer " << std::endl;
		buffer.stop(); //this communicates to the file_buf_rank_ that it has to stop the run() loop.
	} else if ( rank_ >= min_client_rank_ ) {
			go_main( mover );
			tr.Debug << "Slave JD finished!" << std::endl;
	}

	// ideally these would be called in the dtor but the way we have the singleton pattern set up the dtors don't get
	// called
#ifdef USEMPI
 	MPI_Barrier( MPI_COMM_WORLD );
 	MPI_Finalize();
#endif
}


/// @details This is the heart of the MPIFileBufJobDistributor. It consistits of two while loops: the job
///distribution loop (JDL) and the node spin down loop (NSDL). The JDL has three functions. The first is to recieve and
///process messages from the slave nodes requesting new job ids. The second is to recieve and process messages from the
///slave nodes indicating a bad input. The third is to recive and process job_success messages from the slave nodes and
///block while the slave node is writing its output. This is prevent Sizeerleaving of output in score files and silent
///files. The function of the NSDL is to keep the head node alive while there are still slave nodes processing. Without
///the NSDL if a slave node finished its allocated job after the head node had finished handing out all of the jobs and
///exiting (a very likely scenario), it would wait indefinitely for a response from the head node when requesting a new
///job id.
void MPIFileBufJobDistributor::send_job_to_slave( Size MPI_ONLY(slave_rank) ) {
#ifdef USEMPI
	int buf[ 2 ];
	if ( rank_ == master_rank_ ) {
		buf[ 0 ] = current_job_id();
		buf[ 1 ] = current_batch_id();
		tr.Debug << "Master: send new job: " << buf[ 0 ] << " " << buf[ 1 ] << std::endl;
		MPI_Send( &buf, 2, MPI_INT, slave_rank, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
	} else {
		runtime_assert( rank_ == slave_rank );
		int buf[ 2 ];		MPI_Status status;
		MPI_Recv( &buf, 2, MPI_INT, master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD, &status );
		slave_current_job_id_ = buf[ 0 ];
		slave_current_batch_id_ = buf[ 1 ];
		tr.Debug << "Slave: receive job: " << buf[ 0 ] << " " << buf[ 1 ] << std::endl;
	}
#endif
}

/// @details messages are received constantly by Master JobDistributor and then the virtual process_message() method
/// is used to assign some action to each message ... this allows child-classes to answer to more messages or change behaviour of already known messages
bool
MPIFileBufJobDistributor::process_message( Size msg_tag, Size slave_rank, Size slave_job_id, Size slave_batch_id, core::Real runtime ) {
	switch ( msg_tag ) {
	case NEW_JOB_ID:  //slave requested a new job id ... send new job or spin-down signal
		tr.Debug << "Master Node: Sending new job id " << current_job_id() << " " << "job: " << current_job()->input_tag() << " to node " << slave_rank << std::endl;
		send_job_to_slave( slave_rank );
		if ( current_job_id() ) {
			++jobs_assigned_;
			obtain_new_job();
		} else {
			--n_nodes_left_to_spin_down_;
		}
		break;
	case BAD_INPUT: //slave reports failed job
		tr.Debug << "Master Node: Received job failure message for job id " << slave_job_id << " from node " << slave_rank << std::endl;
		++bad_input_jobs_;
		mark_job_as_bad( slave_job_id, slave_batch_id );
		++jobs_returned_;
		break;
	case JOB_SUCCESS:
		mark_job_as_completed( slave_job_id, slave_batch_id, runtime );
		++jobs_returned_;
		break;
	case JOB_FAILED_NO_RETRY :
		++jobs_returned_;
		break;
	default:
		tr.Error << "[ERROR] from " << slave_rank << " tag: "  << msg_tag << " " << slave_job_id << std::endl;
		utility_exit_with_message(" unknown tag "+ ObjexxFCL::string_of( msg_tag ) +" in master_loop of MPIFileBufJobDistributor ");
		return false;
	}
	return true;
}

/// @brief mark job as completed
void MPIFileBufJobDistributor::mark_job_as_completed( core::Size job_id, core::Size batch_id, core::Real runtime ) {
	if ( runtime > 0 ) {
		cumulated_runtime_+=runtime;
		cumulated_jobs_ += 1;
	}
	if ( batch_id == current_batch_id() ) {
		Parent::mark_job_as_completed( job_id, runtime );
	}
}

/// @brief mark job as failed --- remove future versions of same input from list
void MPIFileBufJobDistributor::mark_job_as_bad( core::Size job_id, core::Size batch_id ) {
	if ( batch_id == current_batch_id() ) {
		bad_job_id_ = job_id;
		remove_bad_inputs_from_job_list();
	}
}

/// @brief receive message of certain type -- and ignore it ... sometimes needed in communication protocol
void
MPIFileBufJobDistributor::eat_signal( Size msg_tag_in, int MPI_ONLY( source ) ) {
#ifdef USEMPI
	Size const mpi_size( 4 );
	int mpi_buf[ mpi_size ];
	while( true )  {
		MPI_Status status;
		MPI_Recv( &mpi_buf, mpi_size, MPI_INT, source, MPI_JOB_DIST_TAG, MPI_COMM_WORLD, &status);
		Size slave_rank( status.MPI_SOURCE );
		Size const msg_tag ( mpi_buf[ 0 ] );
		Size const slave_job_id( mpi_buf[ 1 ] );
		Size const slave_batch_id( mpi_buf[ 2 ]);
		Size const runtime( mpi_buf[ 3 ] );
		if ( msg_tag_in != msg_tag ) {
			tr.Debug <<" when trying to eat signal " << msg_tag_in <<" I received " << msg_tag << " going to process this now" << std::endl;
			process_message( msg_tag, slave_rank, slave_job_id, slave_batch_id, runtime );
		} else {
			break;
		}
	}
#else
	//Size const msg_tag ( 0 );
#endif

	tr.Debug << "eating expected signal " << msg_tag_in << std::endl;
	//	runtime_assert( msg_tag == msg_tag_in );
}

/// @brief the main message loop --- master cycles thru until all slave nodes have been spun down
void
MPIFileBufJobDistributor::master_go( protocols::moves::MoverOP /*mover*/ )
{
#ifdef USEMPI
	runtime_assert( rank_ == master_rank_ );

	Size const mpi_size( 4 );
	int mpi_buf[ mpi_size ];

	MPI_Status status;
	MPI_Request request;

	// set first job to assign
	obtain_new_job();

	// initialize some statistics  -- these are member variables, since they are also used in process_messages()
	jobs_assigned_ = 0;
	jobs_returned_ = 0;
	bad_input_jobs_ = 0;
	tr.Info << "Starting JobDistribution with " << n_worker() << " worker processes " << std::endl;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	n_nodes_left_to_spin_down_ = option[ OptionKeys::jd2::mpi_nowait_for_remaining_jobs ]() ? 0 : ( n_worker() );
	double timeout_wait_factor = option[ OptionKeys::jd2::mpi_timeout_factor ]();  //wait at most X-times the average job-time for stray jobs...
	// Job Distribution Loop  --- receive message and process -- repeat
	while ( current_job_id() || jobs_returned_ < jobs_assigned_ || n_nodes_left_to_spin_down_ ) {
		tr.Debug << "current_job_id: " << current_job_id() << " jobs_returned " << jobs_returned_
						 << " jobs_assigned_ " << jobs_assigned_ << " nodes_to_spin_down " << n_nodes_left_to_spin_down_ << std::endl;
		//receive message
		tr.Debug << "Master Node: Waiting for job requests..." << std::endl;
		MPI_Irecv( &mpi_buf, mpi_size, MPI_INT, MPI_ANY_SOURCE, MPI_JOB_DIST_TAG, MPI_COMM_WORLD, &request);
		int flag=0;

		//unused variable removed to silence warnings
		//bool only_some_nodes_unfinished = ( 1.0 * n_nodes_left_to_spin_down_ / n_worker() ) < 0.2;


		double timeout=MPI_Wtime();
		double timeout_limit = ( cumulated_jobs_ < 1 || timeout_wait_factor <=0 ) ? 1e9 : cumulated_runtime_/cumulated_jobs_ * timeout_wait_factor ;
		tr.Debug << "set timeout_limit to " << timeout_limit << std::endl;
		while ( flag==0 ) {
			if ( MPI_Wtime() - timeout > timeout_limit ) break;
			MPI_Test(&request, &flag, &status);
		}
		if ( flag == 0 ) { //timeout
			utility_exit_with_message("quick exit from job-distributor due to flag jd2::mpi_nowait_for_remaining_jobs and timeout of "+
				utility::to_string( MPI_Wtime()-timeout )+" seconds\n"+"increase time-out by using -jd2:mpi_timeout_factor" );
		}
		Size slave_rank( status.MPI_SOURCE );
		Size const msg_tag ( mpi_buf[ 0 ] );
		Size const slave_job_id( mpi_buf[ 1 ] );
		Size const slave_batch_id( mpi_buf[ 2 ]);
		Real const runtime( mpi_buf[ 3 ]);
		tr.Debug << "Master Node: Recieved message from  " << slave_rank << " with tag "
						 << msg_tag << " slave_jobid " << slave_job_id << " slave batchid " << slave_batch_id << std::endl;

		//process message
		process_message( msg_tag, slave_rank, slave_job_id, slave_batch_id, runtime);

	}

	//finished
	tr.Info << "Master Node: Finished sending spin down signals to slaves" << std::endl;
	tr.Info << "Master Node stats: jobs-send out: " << jobs_assigned_ << "  returned: " << jobs_returned_ << "  bad input jobs: " << bad_input_jobs_ << std::endl;

	if ( option[ OptionKeys::jd2::mpi_nowait_for_remaining_jobs ]() ) {
		utility_exit_with_message("quick exit from job-distributor due to flag jd2::mpi_nowait_for_remaining_jobs --- this is not an error " );
	}

#endif
}

/// @brief dummy for master/slave version
core::Size
MPIFileBufJobDistributor::get_new_job_id()
{
  if ( rank_ == master_rank_ ) {
    return master_get_new_job_id();
  } else {
    return slave_get_new_job_id();
  }
	return 0;
}


/// @brief work out what next job is
core::Size
MPIFileBufJobDistributor::master_get_new_job_id()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

  JobsContainer const & jobs( get_jobs() );
  JobOutputterOP outputter = job_outputter();

	//increase job-id until a new job is found
	for( core::Size next_job_to_assign = current_job_id() + 1;
       next_job_to_assign <= jobs.size(); ++next_job_to_assign ) {

    basic::show_time( tr,  "assign job "+ObjexxFCL::string_of(next_job_to_assign)+" batch: "+ObjexxFCL::lead_zero_string_of( current_batch_id(),5 ) );

		if ( jobs[ next_job_to_assign ]->bad() ) { //don't start jobs with known bad input
      tr.Debug << "Job " << ObjexxFCL::string_of(next_job_to_assign) << " being skipped due to known bad input." << std::endl;
			continue;
		} else if ( !outputter->job_has_completed( jobs[ next_job_to_assign ] ) ) { //don't start jobs which have been completed ( in previous runs )
				tr.Debug << "Master Node: Getting next job to assign from list id " << next_job_to_assign << " of " << jobs.size() << std::endl;
				return next_job_to_assign;
		} else if ( outputter->job_has_completed( jobs[ next_job_to_assign ] ) &&
                option[ out::overwrite ].value() ) {  //ignore what I just said -- we ignore previous data
			tr.Debug << "Master Node: Getting next job to assign from list, overwriting id " << next_job_to_assign << " of " << jobs.size() << std::endl;
			return next_job_to_assign;
		}

    //arrives here only if job has already been completed on the file-system
    mark_job_as_completed( next_job_to_assign, current_batch_id(), -1.0 ); //need this for the MPIArchiveJobDistributor
	}
	tr.Debug << "Master Node: No more jobs to assign, setting next job id to zero" << std::endl;
	return 0;
}

//overloaded so that slave-nodes never automatically switch to next_batch when spinning down.
bool
MPIFileBufJobDistributor::next_batch() {
	if ( rank_ != master_rank_ ) return false; //slave answer
	cumulated_jobs_ = 0;
	cumulated_runtime_ = 0;
	return Parent::next_batch();
}

core::Size
MPIFileBufJobDistributor::slave_get_new_job_id()
{
#ifdef USEMPI
	runtime_assert( rank_ != master_rank_ );

	if ( repeat_job_ == true ) {
		tr.Debug << "Slave Node " << rank_ << ": Repeating job id " << slave_current_job_id_ <<std::endl;
		repeat_job_ = false;
	}	else {
		tr.Debug << "Slave Node " << rank_ << ": Requesting new job id from master" <<std::endl;
		slave_to_master( NEW_JOB_ID );
		send_job_to_slave( rank_ );
		if ( slave_current_job_id_ ) set_batch_id( slave_current_batch_id_ );
		tr.Debug << "Slave Node " << rank_ << ": Received job id " << slave_current_job_id_
						 <<  ( slave_current_batch_id_ ? " batch: "+ get_current_batch() : "" ) << " from master" << std::endl;

	}
#endif
	return slave_current_job_id_;
}

/// @brief dummy for master/slave version
void
MPIFileBufJobDistributor::mark_current_job_id_for_repetition()
{
  if ( rank_ == master_rank_ ) {
    master_mark_current_job_id_for_repetition();
  } else {
    slave_mark_current_job_id_for_repetition();
		clear_current_job_output();
  }
}

void
MPIFileBufJobDistributor::master_mark_current_job_id_for_repetition()
{
	runtime_assert( rank_ == master_rank_ );
	tr.Debug << "Master Node: Mark current job for repetition" << std::endl;
	utility_exit_with_message( "Master Node: master_mark_current_job_id_for_repetition() should never be called" );
}

void
MPIFileBufJobDistributor::slave_mark_current_job_id_for_repetition()
{
	runtime_assert( rank_ != master_rank_ );
	tr.Debug << "Slave Node " << rank_ << ": Mark current job for repetition, id " << current_job_id() << std::endl;
	repeat_job_ = true;
}

/// @brief dummy for master/slave version
void
MPIFileBufJobDistributor::remove_bad_inputs_from_job_list()
{
  if ( rank_ == master_rank_ ) {
    master_remove_bad_inputs_from_job_list();
  } else {
    slave_remove_bad_inputs_from_job_list();
  }
}

void
MPIFileBufJobDistributor::master_remove_bad_inputs_from_job_list()
{
	runtime_assert( rank_ == master_rank_ );

	if ( tr.Debug.visible() ) {
		JobsContainer const& jobs( get_jobs() );
		std::string const & bad_job_id_input_tag( jobs[ bad_job_id_ ]->input_tag() );
		tr.Debug << "Master Node: Job id "
						 << job_outputter()->output_name( jobs[ bad_job_id_ ] )
						 << " failed, reporting bad input; other jobs of same input will be canceled: " <<  bad_job_id_input_tag << std::endl;
	}

	Parent::mark_job_as_bad( bad_job_id_ ); // this sets all jobs with this input_tag to bad!
	obtain_new_job( true /*re_consider_current_job*/ );
}

void
MPIFileBufJobDistributor::slave_to_master( Size MPI_ONLY(tag) ) {
#ifdef USEMPI
	runtime_assert( rank_ != master_rank_ );
	Size const mpi_size( 4 );
	int mpi_buf[ mpi_size ];
	mpi_buf[ 0 ] = tag;
	mpi_buf[ 1 ] = slave_current_job_id_;
	mpi_buf[ 2 ] = slave_current_batch_id_;
	mpi_buf[ 3 ] = static_cast< int > ( slave_current_runtime_ );
	MPI_Send( &mpi_buf, mpi_size, MPI_INT, master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
#endif
}

void
MPIFileBufJobDistributor::slave_remove_bad_inputs_from_job_list()
{
	slave_to_master( BAD_INPUT );
}

/// @brief dummy for master/slave version
void
MPIFileBufJobDistributor::job_succeeded(core::pose::Pose & pose, core::Real runtime, std::string const & tag)
{
  if ( rank_ == master_rank_ ) {
    master_job_succeeded( pose, tag );
  } else {
		slave_current_runtime_ = runtime;
    slave_job_succeeded( pose, tag );
  }
}

void
MPIFileBufJobDistributor::job_failed(
	core::pose::Pose &,
	bool will_retry
)
{
	assert( rank_ >= min_client_rank_ );
	if ( ! will_retry ) {
		// tell the master node that this job has failed and will not be
		// re-attempted
		slave_to_master( JOB_FAILED_NO_RETRY );
	}
}


void
MPIFileBufJobDistributor::master_job_succeeded(core::pose::Pose & /*pose*/, std::string const & /*tag*/)
{
#ifdef USEMPI
	runtime_assert( rank_ == master_rank_ );
	tr.Debug << "Master Node: Job Succeeded" << std::endl;
	utility_exit_with_message( "Master Node: master_job_succeeded() should never be called" );
#endif
}

void
MPIFileBufJobDistributor::slave_job_succeeded(core::pose::Pose &pose, std::string const & tag )
{
	runtime_assert( !( rank_ == master_rank_ ) );
	job_outputter()->final_pose( current_job(), pose, tag );
	slave_to_master( JOB_SUCCESS );
}


}//jd2
}//protocols
