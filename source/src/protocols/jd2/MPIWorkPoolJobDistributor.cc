// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MPIWorkPoolJobDistributor.cc
/// @brief  implementation of MPIWorkPoolJobDistributor
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

// Unit headers
#include <protocols/jd2/MPIWorkPoolJobDistributor.hh>

// Package headers
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <basic/message_listening/MessageListenerFactory.hh>
#include <basic/message_listening/MessageListener.hh>
#include <basic/message_listening/util.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/LargeNstructJobInputter.hh>


#include <protocols/moves/Mover.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/exit.hh>
#include <utility/assert.hh>
#include <utility/mpi_util.hh>

// Option headers
#include <basic/options/keys/out.OptionKeys.gen.hh>
#ifdef USEMPI
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#endif

// C++ headers
#include <string>

//Auto Headers
#include <utility/vector1.hh>
static thread_local basic::Tracer TR( "protocols.jd2.MPIWorkPoolJobDistributor" );

namespace protocols {
namespace jd2 {

using namespace basic::options;
using namespace basic::options::OptionKeys;

/// @details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
///which processor it is in MPI land.
MPIWorkPoolJobDistributor::MPIWorkPoolJobDistributor() :
	JobDistributor(),
	npes_( 1 ),
	rank_( 0 ),
	current_job_id_( 0 ),
	next_job_to_assign_( 0 ),
	bad_job_id_( 0 ),
	repeat_job_( false ),
	finalize_MPI_( true ),
	sequential_distribution_(false),
	starter_for_sequential_distribution_(false)
{
	// set npes and rank based on whether we are using MPI or not
#ifdef USEMPI
  //npes_ = MPI::COMM_WORLD.Get_size();
  //rank_ = MPI::COMM_WORLD.Get_rank();
  int int_npes, int_rank;                                 //don't cast pointers - copy it over instead
  MPI_Comm_rank( MPI_COMM_WORLD, &int_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &int_npes );
  rank_ = int_rank;
  npes_ = int_npes;
  if(rank_!=0) get_jobs_nonconst().set_force_job_purging(true); //If this is a slave process, the jobs list should purge old jobs even if not marked for deletion when higher-index jobs are requested.
#else
	utility_exit_with_message( "ERROR ERROR ERROR: The MPIWorkPoolJobDistributor will not work unless you have compiled using extras=mpi" );
#endif
}

/// @brief dtor
///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
MPIWorkPoolJobDistributor::~MPIWorkPoolJobDistributor()
{ }

/// @brief dummy for master/slave version
void
MPIWorkPoolJobDistributor::go( protocols::moves::MoverOP mover )
{
#ifdef USEMPI
	set_sequential_distribution( option[ OptionKeys::jd2::sequential_mpi_job_distribution	].user() ); //Set whether we're sending jobs to slaves in sequence (1,2,3,etc.) or whether slaves request jobs as they become available. 
#endif

	if ( rank_ == 0 ) {
	master_go( mover );
	} else {
	slave_go( mover );
	}

	// ideally these would be called in the dtor but the way we have the singleton pattern set up the dtors don't get
	// called
#ifdef USEMPI
	//MPI::COMM_WORLD.Barrier();
	//MPI::Finalize();
	MPI_Barrier( MPI_COMM_WORLD );
	if(finalize_MPI_)
	{
		MPI_Finalize();
	}
#endif
}


/// @details This is the heart of the MPIWorkPoolJobDistributor. It consists of two while loops: the job
///distribution loop (JDL) and the node spin down loop (NSDL). The JDL has three functions. The first is to receive and
///process messages from the slave nodes requesting new job ids. The second is to receive and process messages from the
///slave nodes indicating a bad input. The third is to receive and process job_success messages from the slave nodes and
///block while the slave node is writing its output. This is prevent interleaving of output in score files and silent
///files. The function of the NSDL is to keep the head node alive while there are still slave nodes processing. Without
///the NSDL if a slave node finished its allocated job after the head node had finished handing out all of the jobs and
///exiting (a very likely scenario), it would wait indefinitely for a response from the head node when requesting a new
///job id.
void
MPIWorkPoolJobDistributor::master_go( protocols::moves::MoverOP /*mover*/ )
{
#ifdef USEMPI
	runtime_assert( rank_ == 0 );

	core::Size slave_data( 0 );
	MPI_Status status;
	char dummy('A'); //A dummy variable for when I need to send an acknowledgement signal.  Using a char because it's guaranteed to be 1 byte, while other data types might be wider, wasting bandwidth.
	core::Size completed_job_id(0); //A storage spot for keeping track of what job has just been named as completed by the slave.

	// set first job to assign
	master_get_new_job_id();

	// Job Distribution Loop
	while ( next_job_to_assign_ != 0 ) {
		if(TR.visible()) TR << "Master Node: Waiting for job requests..." << std::endl;
		MPI_Recv( &slave_data, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if(TR.visible()) TR << "Master Node: Received message from " << status.MPI_SOURCE << " with tag " << status.MPI_TAG << std::endl;

		// decide what to do based on message tag
		//switch ( status.MPI::Status::Get_tag() ) {
		switch ( status.MPI_TAG ) {
			case NEW_JOB_ID_TAG:  //The slave node has requested a new job id.
				if(TR.visible()) TR << "Master Node: Sending new job id " << next_job_to_assign_ << " to node " << status.MPI_SOURCE << " with tag " << NEW_JOB_ID_TAG << std::endl;
				MPI_Send( &next_job_to_assign_, 1, MPI_UNSIGNED_LONG, status.MPI_SOURCE, NEW_JOB_ID_TAG, MPI_COMM_WORLD );
				master_get_new_job_id();
				break;
			case BAD_INPUT_TAG: //The slave node has reported a bad job.
				if(TR.visible()) TR << "Master Node: Received job failure message (due to bad input) for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
				bad_job_id_ = slave_data;
				master_remove_bad_inputs_from_job_list(); //This will mark the job as deletable
				break;
			case JOB_SUCCESS_TAG: //The slave node has reported a successful job completion.
				if(TR.visible()) TR << "Master Node: Received job success message for job id " << slave_data << " from node " << status.MPI_SOURCE << " blocking till output is done " << std::endl;
				completed_job_id = slave_data; //Store the id of the job just completed by the slave
				MPI_Send( &next_job_to_assign_, 1, MPI_UNSIGNED_LONG, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD );
				MPI_Recv( &slave_data, 1, MPI_UNSIGNED_LONG, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD, &status);
				if(TR.visible()) TR << "Master Node: Received job output finish message for job id " << completed_job_id << " from node " << status.MPI_SOURCE << std::endl;
				master_mark_job_as_completed( completed_job_id ); //This will mark the job as completed/deletable.
				break;
			case JOB_FAILURE_TAG:
				if(TR.visible()) {
					TR << "Master node: Received job failure message for job id " << slave_data << " from node " << status.MPI_SOURCE << "." << std::endl;
					TR.flush();
				}
				master_mark_job_as_failed( slave_data ); //This will mark the job as deletable.
				//Master needs to acknowledge the failure before the slave requests a new job.  Otherwise, master might see the slave's NEW_JOB_ID_TAG request before seeing the JOB_FAILURE_TAG.
				MPI_Send( &dummy, 1, MPI_CHAR, status.MPI_SOURCE, JOB_FAILURE_TAG, MPI_COMM_WORLD ); //Acknowledgement.
				break;
			case REQUEST_MESSAGE_TAG:
			{

				using namespace basic::message_listening;

				listener_tags listener_tag((listener_tags)slave_data);
				MessageListenerOP listener(MessageListenerFactory::get_instance()->get_listener(listener_tag));

				std::string message_data = utility::receive_string_from_node(status.MPI_SOURCE);
				std::string return_info="";
				bool request_slave_data = listener->request(message_data, return_info);
				utility::send_string_to_node(status.MPI_SOURCE, return_info);

				if(TR.visible()) TR
					<< "Master Node: node '" << status.MPI_SOURCE << "' "
					<< "requests from the message listener '" << listener_tag_to_name(listener_tag) << "' "
					<< "data on '" << message_data << "', "
					<< "respond with '" << return_info << "' "
					<< (request_slave_data ? " and requests more data." : ".") << std::endl;

				if(request_slave_data){
					message_data = utility::receive_string_from_node(status.MPI_SOURCE);
					if(TR.visible()) TR
						<< "Master Node: Received from node '" << status.MPI_SOURCE << "' "
						<< "'" << message_data << "'" << std::endl;
					listener->receive(message_data);
				}

				break;

			}
			default:
			{
				std::stringstream err_msg;
				err_msg
					<< "Received unrecognized mpi_tag '" << status.MPI_TAG << "' " << std::endl
					<< "\tfrom node '" << status.MPI_SOURCE << "' " << std::endl
					<< "\twith data '" << slave_data << "'";
				utility_exit_with_message(err_msg.str());
			}

		}
	}
	TR << "Master Node: Finished handing out jobs" << std::endl;

	core::Size n_nodes_left_to_spin_down( npes_ - 1 ); // don't have to spin down self

	// Node Spin Down loop
	while ( n_nodes_left_to_spin_down > 0 ) {
		if(TR.visible()) TR << "Master Node: Waiting for " << n_nodes_left_to_spin_down << " slaves to finish jobs" << std::endl;
		//MPI::COMM_WORLD.Recv( &slave_data, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status );
		//TR << "Master Node: Received message from  " << status.MPI::Status::Get_source() << " with tag " << status.MPI::Status::Get_tag() << std::endl;
		MPI_Recv( &slave_data, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if(TR.visible()) TR << "Master Node: Received message from  " << status.MPI_SOURCE << " with tag " << status.MPI_TAG << std::endl;

		// decide what to do based on message tag
		//switch ( status.MPI::Status::Get_tag() ) {
		switch ( status.MPI_TAG ) {
			case NEW_JOB_ID_TAG:
				//TR << "Master Node: Sending spin down signal to node " << status.MPI::Status::Get_source() << std::endl;
				//MPI::COMM_WORLD.Send( &next_job_to_assign_, 1, MPI::INT, status.MPI::Status::Get_source(), NEW_JOB_ID_TAG );
				if(TR.visible()) TR << "Master Node: Sending spin down signal to node " << status.MPI_SOURCE << std::endl;
				MPI_Send( &next_job_to_assign_, 1, MPI_UNSIGNED_LONG, status.MPI_SOURCE, NEW_JOB_ID_TAG, MPI_COMM_WORLD );
				n_nodes_left_to_spin_down--;
				break;
			case BAD_INPUT_TAG:
				if(TR.visible()) TR << "Master Node: Received job failure message (due to bad input) for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
				bad_job_id_ = slave_data;
				master_remove_bad_inputs_from_job_list(); //This will mark the job as deletable
				break;
			case JOB_SUCCESS_TAG:
				if(TR.visible()) TR << "Master Node: Received job success message for job id " << slave_data << " from node " << status.MPI_SOURCE << " blocking till output is done " << std::endl;
				completed_job_id = slave_data; //Store the id of the job just completed by the slave
				MPI_Send( &next_job_to_assign_, 1, MPI_UNSIGNED_LONG, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD );
				MPI_Recv( &slave_data, 1, MPI_UNSIGNED_LONG, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD, &status);
				if(TR.visible()) TR << "Master Node: Received job output finish message for job id " << completed_job_id << " from node " << status.MPI_SOURCE << std::endl;
				master_mark_job_as_completed( completed_job_id ); //This will mark the job as completed/deletable.
				break;
			case JOB_FAILURE_TAG:
				if(TR.visible()) TR << "Master node: Received job failure message for job id " << slave_data << " from node " << status.MPI_SOURCE << "." << std::endl;
				master_mark_job_as_failed( slave_data ); //This will mark the job as deletable.
				MPI_Send( &dummy, 1, MPI_CHAR, status.MPI_SOURCE, JOB_FAILURE_TAG, MPI_COMM_WORLD ); //Acknowledgement.
				master_mark_job_as_failed( slave_data ); //This will mark the job as deletable.
				//Master needs to acknowledge the failure before the slave requests a new job.  Otherwise, master might see the slave's NEW_JOB_ID_TAG request before seeing the JOB_FAILURE_TAG.
				MPI_Send( &dummy, 1, MPI_CHAR, status.MPI_SOURCE, JOB_FAILURE_TAG, MPI_COMM_WORLD ); //Acknowledgement.
				break;
			case REQUEST_MESSAGE_TAG:
			{
				using namespace basic::message_listening;

				listener_tags listener_tag((listener_tags)slave_data);
				MessageListenerOP listener(MessageListenerFactory::get_instance()->get_listener(listener_tag));

				std::string message_data = utility::receive_string_from_node(status.MPI_SOURCE);
				std::string return_info="";
				bool request_slave_data = listener->request(message_data, return_info);
				utility::send_string_to_node(status.MPI_SOURCE, return_info);

				if(TR.visible()) TR
					<< "Master Node: node '" << status.MPI_SOURCE << "' "
					<< "requests from the message listener '" << listener_tag_to_name(listener_tag) << "' "
					<< "data on '" << message_data << "', "
					<< "respond with '" << return_info << "' "
					<< (request_slave_data ? " and requests more data." : ".") << std::endl;

				if(request_slave_data){
					message_data = utility::receive_string_from_node(status.MPI_SOURCE);
					if(TR.visible()) TR
						<< "Master Node: Received from node '" << status.MPI_SOURCE << "' "
						<< "'" << message_data << "'" << std::endl;
					listener->receive(message_data);
				}

				break;
			}
			default:
			{
				std::stringstream err_msg;
				err_msg
					<< "Received unrecognized mpi_tag '" << status.MPI_TAG << "' " << std::endl
					<< "\tfrom node '" << status.MPI_SOURCE << "' " << std::endl
					<< "\twith data '" << slave_data << "'";
				utility_exit_with_message(err_msg.str());
			}
		}
	}
	if(TR.visible()) TR << "Master Node: Finished sending spin down signals to slaves" << std::endl;
#endif
}

void
MPIWorkPoolJobDistributor::slave_go( protocols::moves::MoverOP mover )
{
	runtime_assert( !( rank_ == 0 ) );
	
	if( sequential_distribution() && rank_ == 1 ) set_starter_for_sequential_distribution( true );
	
	go_main( mover );
}

/// @brief dummy for master/slave version
core::Size
MPIWorkPoolJobDistributor::get_new_job_id()
{
	core::Size temp( 0 );

	if ( rank_ == 0 ) {
	temp = master_get_new_job_id();
	} else {
	temp = slave_get_new_job_id();
	}

	return temp;
}

core::Size
MPIWorkPoolJobDistributor::master_get_new_job_id()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	JobsContainer & jobs( get_jobs_nonconst() ); //Must be non-const access in case jobs list needs to be updated.
	JobOutputterOP outputter = job_outputter();

	while( next_job_to_assign_ <= jobs.size()) {
		++next_job_to_assign_;
		if ( next_job_to_assign_ > jobs.size() ) {
			if(TR.visible()) TR << "Master Node: No more jobs to assign, setting next job id to zero" << std::endl;
			next_job_to_assign_ = 0;
			return 0;
		}	else if ( !outputter->job_has_completed( jobs[ next_job_to_assign_ ] ) ) {
			if(TR.visible()) TR << "Master Node: Getting next job to assign from list id " << next_job_to_assign_ << " of " << jobs.size() << std::endl;
			return next_job_to_assign_; //not used by callers
		} else if ( outputter->job_has_completed( jobs[ next_job_to_assign_ ] ) && option[ out::overwrite ].value() ) {
			if(TR.visible()) TR << "Master Node: Getting next job to assign from list, overwriting id " << next_job_to_assign_ << " of " << jobs.size() << std::endl;
			return next_job_to_assign_; //not used by callers
		}
	}

	return 0; //we won't get here
}

core::Size
MPIWorkPoolJobDistributor::slave_get_new_job_id()
{
#ifdef USEMPI
	runtime_assert( !( rank_ == 0 ) );
	
	if( sequential_distribution() && !starter_for_sequential_distribution()) wait_for_go_signal();

	if ( repeat_job_ == true ) {
		if(TR.visible()) TR << "Slave Node " << rank_ << ": Repeating job id " << current_job_id_ <<std::endl;
		repeat_job_ = false;
	}	else {
		if(TR.visible()) TR << "Slave Node " << rank_ << ": Requesting new job id from master" <<std::endl;
		core::Size empty_data( 0 );
		MPI_Status status;
		current_job_id_ = 0;
		//MPI::COMM_WORLD.Send( &empty_data, 1, MPI::INT, 0, NEW_JOB_ID_TAG );
		//MPI::COMM_WORLD.Recv( &current_job_id_, 1, MPI::INT, 0, NEW_JOB_ID_TAG );
		MPI_Send( &empty_data, 1, MPI_UNSIGNED_LONG, 0, NEW_JOB_ID_TAG, MPI_COMM_WORLD );

		MPI_Recv( &current_job_id_, 1, MPI_UNSIGNED_LONG, 0, NEW_JOB_ID_TAG, MPI_COMM_WORLD, &status );
		
		//If we're using the LargeNstructJobInputter, which limits the number of jobs held in memory at any given time, then we need to ensure that
		//the slave knows that lower-numbered jobs can be deleted:		
		JobsContainer & jobs( get_jobs_nonconst() ); //Must be nonconst to mark jobs as deletable.
		bool greater_than=false;
		if(jobs.highest_job_index() < current_job_id_) { //Mark all jobs with lower id values as deletable if we've received a higher job id.  We will definitely not receive a lower job id, I think.
			greater_than=true;
			for(core::Size i=1, imax=jobs.highest_job_index(); i<=imax; ++i) {
				if(jobs.has_job(i)) jobs[i]->set_can_be_deleted(true);
			}
		}

		if(TR.visible()) {
			TR << "Slave Node " << rank_ << ": Received job id " << current_job_id_ << " from master.";
			if(greater_than) TR << "  This index was higher than the highest index currently loaded in memory.  Purging old jobs.";
			TR << std::endl;
		}		
	}
	
	if( sequential_distribution() ) send_go_signal();
	
#endif
	return current_job_id_;
}

/// @brief dummy for master/slave version
void
MPIWorkPoolJobDistributor::mark_current_job_id_for_repetition()
{
	if ( rank_ == 0 ) {
	master_mark_current_job_id_for_repetition();
	} else {
	slave_mark_current_job_id_for_repetition();
	}
	clear_current_job_output();
}

void
MPIWorkPoolJobDistributor::master_mark_current_job_id_for_repetition()
{
	runtime_assert( rank_ == 0 );
	if(TR.visible()) TR << "Master Node: Mark current job for repetition" << std::endl;
	utility_exit_with_message( "Master Node: master_mark_current_job_id_for_repetition() should never be called" );

}

void
MPIWorkPoolJobDistributor::slave_mark_current_job_id_for_repetition()
{
	runtime_assert( !( rank_ == 0 ) );
	if(TR.visible()) TR << "Slave Node " << rank_ << ": Mark current job for repetition, id " << current_job_id_ << std::endl;
	repeat_job_ = true;
}

/// @brief dummy for master/slave version
void
MPIWorkPoolJobDistributor::remove_bad_inputs_from_job_list()
{
	if ( rank_ == 0 ) {
	master_remove_bad_inputs_from_job_list();
	} else {
	slave_remove_bad_inputs_from_job_list();
	}
}

void
MPIWorkPoolJobDistributor::master_remove_bad_inputs_from_job_list()
{
	//#ifdef USEMPI
	runtime_assert( rank_ == 0 );

	JobsContainer & jobs( get_jobs_nonconst() ); //Must be nonconst to mark jobs as deletable.

	std::string const & bad_job_id_input_tag( jobs[ bad_job_id_ ]->input_tag() );

	if(TR.visible()) TR << "Master Node: Job id " << bad_job_id_ << " failed, reporting bad input; other jobs of same input will be canceled: " << job_outputter()->output_name( jobs[ bad_job_id_ ] ) << std::endl;
	jobs[bad_job_id_]->set_can_be_deleted(true);

	while( next_job_to_assign_ <= jobs.size() && jobs[ next_job_to_assign_ ]->input_tag() == bad_job_id_input_tag ) {
		if(TR.visible()) TR << "Master Node: Job canceled without trying due to previous bad input: " << job_outputter()->output_name( jobs[ next_job_to_assign_ ] ) << " id " << next_job_to_assign_ << std::endl;
		jobs[next_job_to_assign_]->set_can_be_deleted(true);
		++next_job_to_assign_;
	}

	//iterate through for overwrite/end of vector statuses
	--next_job_to_assign_; //master_get_new_job_id() will ++ this again first thing
	master_get_new_job_id();

	//#endif
}

void
MPIWorkPoolJobDistributor::slave_remove_bad_inputs_from_job_list()
{
#ifdef USEMPI
	runtime_assert( !( rank_ == 0 ) );

	//MPI::COMM_WORLD.Send( &current_job_id_, 1, MPI::INT, 0, BAD_INPUT_TAG );
	MPI_Send( &current_job_id_, 1, MPI_UNSIGNED_LONG, 0, BAD_INPUT_TAG, MPI_COMM_WORLD );
#endif
}

/// @brief dummy for master/slave version
void
MPIWorkPoolJobDistributor::job_succeeded(core::pose::Pose &pose, core::Real /*run_time*/, std::string const & tag)
{
	if ( rank_ == 0 ) {
	master_job_succeeded( pose, tag );
	} else {
	slave_job_succeeded( pose, tag );
	}
}

/// @brief Called if the job failed.
///
void MPIWorkPoolJobDistributor::job_failed(
		core::pose::Pose & /*pose*/,
		bool /*will_retry*/
) {
#ifdef USEMPI
	if(rank_==0) {
		utility_exit_with_message( "In protocols::jd2::MPIWorkPoolJobDistributor::job_failed(): This function should never be called on the master node!\n" );
	} else {
		MPI_Status status;
		increment_failed_jobs();
		char dummy_val('D');
		if(TR.visible()) TR << "Job " << current_job_id_ << " failed!  Slave process " << rank_ << " transmitting failure to master." << std::endl;
		MPI_Send( &current_job_id_, 1, MPI_UNSIGNED_LONG, 0, JOB_FAILURE_TAG, MPI_COMM_WORLD );
		//We need to wait for acknowledgement, lest this process request a new job id (and be given one) before Master has seen the failure message and cleaned up the job list.
		MPI_Recv( &dummy_val, 1, MPI_CHAR, 0, JOB_FAILURE_TAG, MPI_COMM_WORLD, &status);
	}
	return;
#endif
}

void MPIWorkPoolJobDistributor::mpi_finalize(bool finalize)
{
	finalize_MPI_ = finalize;
}

void
MPIWorkPoolJobDistributor::master_job_succeeded(core::pose::Pose & /*pose*/, std::string const & /*tag*/)
{
#ifdef USEMPI
	runtime_assert( rank_ == 0 );
	if(TR.visible()) TR << "Master Node: Job Succeeded" << std::endl;
	utility_exit_with_message( "Master Node: master_job_succeeded() should never be called" );
#endif
}

void
MPIWorkPoolJobDistributor::slave_job_succeeded(core::pose::Pose & MPI_ONLY( pose ), std::string const & MPI_ONLY( tag ) )
{
#ifdef USEMPI
	runtime_assert( !( rank_ == 0 ) );

	if ( option[ OptionKeys::jd2::mpi_fast_nonblocking_output	].value() == true ) {
		job_outputter()->final_pose( current_job(), pose, tag );
	} else {
		core::Size empty_data( 0 );
		MPI_Status status;

		// send job success message to master
		if(TR.visible()) TR << "Slave Node " << rank_ << ": Finished job successfully! Sending output request to master." << std::endl;
		MPI_Send( &current_job_id_, 1, MPI_UNSIGNED_LONG, 0, JOB_SUCCESS_TAG, MPI_COMM_WORLD );

		// receive message from master that says is okay to write
		MPI_Recv( &empty_data, 1, MPI_UNSIGNED_LONG, 0, JOB_SUCCESS_TAG, MPI_COMM_WORLD, &status );
		if(TR.visible()) TR << "Slave Node " << rank_ << ": Received output confirmation from master. Writing output." << std::endl;
		// time and write output (pdb, silent file, score file etc.)
		clock_t starttime = clock();
		job_outputter()->final_pose( current_job(), pose, tag );
		clock_t stoptime = clock();

			// send message to master that we are done outputing
		if(TR.visible()) TR << "Slave Node " << rank_ << ": Finished writing output in " << ((double) stoptime-starttime) / CLOCKS_PER_SEC << " seconds. Sending message to master" << std::endl;
		MPI_Send( &empty_data, 1, MPI_UNSIGNED_LONG, 0, JOB_SUCCESS_TAG, MPI_COMM_WORLD );
	}
#endif
}

/// @brief Mark the job as completed/deletable in the jobs list on the master process.
///
void
MPIWorkPoolJobDistributor::master_mark_job_as_completed(
#ifdef USEMPI
	core::Size const job_index
#else
	core::Size const
#endif
) {
#ifdef USEMPI
	JobsContainer & jobs( get_jobs_nonconst() ); //Must be non-const access in case jobs list needs to be updated.
	jobs[job_index]->set_completed(true); //The job has been completed.  (Also sets it to be deletable, if we're using the LargeNstructJobInputter to save memory.)
	if(TR.visible()) {
		TR << "Master set job " << job_index << " as completed/deletable." << std::endl; //DELETE ME or switch to debug output.
		TR.flush();
	}
#endif
}

/// @brief Mark the job as deletable in the jobs list on the master process.
///
void
MPIWorkPoolJobDistributor::master_mark_job_as_failed(
#ifdef USEMPI
	core::Size const job_index
#else
	core::Size const
#endif
) {
#ifdef USEMPI
	JobsContainer & jobs( get_jobs_nonconst() ); //Must be non-const access in case jobs list needs to be updated.
	jobs[job_index]->set_can_be_deleted(true); //The job can be deleted, since it failed.
	if(TR.visible()) {
		TR << "Master set job " << job_index << " as failed/deletable." << std::endl; //DELETE ME or switch to debug output.
		TR.flush();
	}
#endif
}

/// @brief Wait for a signal from the n-1 process saying that I can proceed.
///
void MPIWorkPoolJobDistributor::wait_for_go_signal() const {
#ifdef USEMPI
	core::Size const prev_proc ( rank_ - 1 > 0 ? rank_ - 1 : npes_ - 1 );
	char dummy_signal('B');
	MPI_Status status;
	MPI_Recv( &dummy_signal, 1, MPI_CHAR, prev_proc, JOB_GO_TAG, MPI_COMM_WORLD, &status);
#endif
	return;
}

/// @brief Send a signal to the n+1 process saying that it can proceed.
/// @details This also sets starter_for_sequential_distribution_ to false, since we no longer
/// want this process to refrain from waiting.
void MPIWorkPoolJobDistributor::send_go_signal() {
#ifdef USEMPI
	core::Size const next_proc ( rank_ + 1 >= npes_ ? 1 : rank_ + 1 );
	char dummy_signal('C');
	MPI_Send( &dummy_signal, 1, MPI_CHAR, next_proc, JOB_GO_TAG, MPI_COMM_WORLD);
	set_starter_for_sequential_distribution(false);
#endif
	return;
}

}//jd2
}//protocols
