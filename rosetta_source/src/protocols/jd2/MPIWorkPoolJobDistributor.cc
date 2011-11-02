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
/// @author P. Douglas Renfrew (renfrew@unc.edu)

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

// Unit headers
#include <protocols/jd2/MPIWorkPoolJobDistributor.hh>

// Package headers
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/moves/Mover.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/exit.hh>
#include <utility/assert.hh>

// Option headers
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

// C++ headers
#include <string>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jd2/InnerJob.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobInputter.fwd.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/jd2/MPIWorkPoolJobDistributor.fwd.hh>
#include <protocols/jd2/Parser.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverStatus.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <utility>
#include <vector>
#include <basic/Tracer.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>

// AUTO-REMOVED #include <ctime>

static basic::Tracer TR("protocols.jd2.MPIWorkPoolJobDistributor");

namespace protocols {
namespace jd2 {

using namespace basic::options;
using namespace basic::options::OptionKeys;

///@details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
///which processor it is in MPI land.
MPIWorkPoolJobDistributor::MPIWorkPoolJobDistributor() :
  JobDistributor(),
  npes_( 1 ),
  rank_( 0 ),
  current_job_id_( 0 ),
	next_job_to_assign_( 0 ),
	bad_job_id_( 0 ),
	repeat_job_( false )
{
  // set npes and rank based on whether we are using MPI or not
#ifdef USEMPI
  //npes_ = MPI::COMM_WORLD.Get_size();
  //rank_ = MPI::COMM_WORLD.Get_rank();
	MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &rank_ ) );
	MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &npes_ ) );
#else
	utility_exit_with_message( "ERROR ERROR ERROR: The MPIWorkPoolJobDistributor will not work unless you have compiled using extras=mpi" );
#endif
}

///@brief dtor
///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
MPIWorkPoolJobDistributor::~MPIWorkPoolJobDistributor()
{ }

///@brief dummy for master/slave version
void
MPIWorkPoolJobDistributor::go( protocols::moves::MoverOP mover )
{
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
	MPI_Finalize();
#endif
}


///@details This is the heart of the MPIWorkPoolJobDistributor. It consists of two while loops: the job
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

	int slave_data( 0 );
	MPI_Status status;

	// set first job to assign
	master_get_new_job_id();

	// Job Distribution Loop
	while ( next_job_to_assign_ != 0 ) {
		TR << "Master Node: Waiting for job requests..." << std::endl;
		//MPI::COMM_WORLD.Recv( &slave_data, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status );
		//TR << "Master Node: Received message from  " << status.MPI::Status::Get_source() << " with tag " << status.MPI::Status::Get_tag() << std::endl;
		MPI_Recv( &slave_data, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		TR << "Master Node: Received message from " << status.MPI_SOURCE << " with tag " << status.MPI_TAG << std::endl;

		// decide what to do based on message tag
		//switch ( status.MPI::Status::Get_tag() ) {
		switch ( status.MPI_TAG ) {
		case NEW_JOB_ID_TAG:
			//TR << "Master Node: Sending new job id " << next_job_to_assign_ << " to node " << status.MPI::Status::Get_source() << " with tag " << NEW_JOB_ID_TAG << std::endl;
			//MPI::COMM_WORLD.Send( &next_job_to_assign_, 1, MPI::INT, status.MPI::Status::Get_source(), NEW_JOB_ID_TAG );
			TR << "Master Node: Sending new job id " << next_job_to_assign_ << " to node " << status.MPI_SOURCE << " with tag " << NEW_JOB_ID_TAG << std::endl;
			MPI_Send( &next_job_to_assign_, 1, MPI_INT, status.MPI_SOURCE, NEW_JOB_ID_TAG, MPI_COMM_WORLD );
			master_get_new_job_id();
			break;
		case BAD_INPUT_TAG:
			//TR << "Master Node: Received job failure message for job id " << slave_data << " from node " << status.MPI::Status::Get_source() << std::endl;
			TR << "Master Node: Received job failure message for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
			bad_job_id_ = slave_data;
			master_remove_bad_inputs_from_job_list();
			break;
		case JOB_SUCCESS_TAG:
			TR << "Master Node: Received job success message for job id " << slave_data << " from node " << status.MPI_SOURCE << " blocking till output is done " << std::endl;
			MPI_Send( &next_job_to_assign_, 1, MPI_INT, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD );
			MPI_Recv( &slave_data, 1, MPI_INT, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD, &status);
			TR << "Master Node: Received job output finish message for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
			break;
		}
	}
	TR << "Master Node: Finished handing out jobs" << std::endl;

	core::Size n_nodes_left_to_spin_down( npes_ - 1 ); // don't have to spin down self

	// Node Spin Down loop
	while ( n_nodes_left_to_spin_down > 0 ) {
		TR << "Master Node: Waiting for " << n_nodes_left_to_spin_down << " slaves to finish jobs" << std::endl;
		//MPI::COMM_WORLD.Recv( &slave_data, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status );
		//TR << "Master Node: Received message from  " << status.MPI::Status::Get_source() << " with tag " << status.MPI::Status::Get_tag() << std::endl;
		MPI_Recv( &slave_data, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		TR << "Master Node: Received message from  " << status.MPI_SOURCE << " with tag " << status.MPI_TAG << std::endl;

		// decide what to do based on message tag
		//switch ( status.MPI::Status::Get_tag() ) {
		switch ( status.MPI_TAG ) {
		case NEW_JOB_ID_TAG:
			//TR << "Master Node: Sending spin down signal to node " << status.MPI::Status::Get_source() << std::endl;
			//MPI::COMM_WORLD.Send( &next_job_to_assign_, 1, MPI::INT, status.MPI::Status::Get_source(), NEW_JOB_ID_TAG );
			TR << "Master Node: Sending spin down signal to node " << status.MPI_SOURCE << std::endl;
			MPI_Send( &next_job_to_assign_, 1, MPI_INT, status.MPI_SOURCE, NEW_JOB_ID_TAG, MPI_COMM_WORLD );
			n_nodes_left_to_spin_down--;
			break;
		case BAD_INPUT_TAG:
			break;
		case JOB_SUCCESS_TAG:
			TR << "Master Node: Received job success message for job id " << slave_data << " from node " << status.MPI_SOURCE << " blocking till output is done " << std::endl;
			MPI_Send( &next_job_to_assign_, 1, MPI_INT, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD );
			MPI_Recv( &slave_data, 1, MPI_INT, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD, &status);
			TR << "Master Node: Received job output finish message for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
			break;
		}
	}
	TR << "Master Node: Finished sending spin down signals to slaves" << std::endl;
#endif
}

void
MPIWorkPoolJobDistributor::slave_go( protocols::moves::MoverOP mover )
{
	runtime_assert( !( rank_ == 0 ) );
	go_main( mover );
}

///@brief dummy for master/slave version
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

  Jobs const & jobs( get_jobs() );
  JobOutputterOP outputter = job_outputter();

	while( next_job_to_assign_ <= jobs.size()) {
		++next_job_to_assign_;
		if ( next_job_to_assign_ > jobs.size() ) {
			TR << "Master Node: No more jobs to assign, setting next job id to zero" << std::endl;
			next_job_to_assign_ = 0;
			return 0;
		}	else if ( !outputter->job_has_completed( jobs[ next_job_to_assign_ ] ) ) {
			TR << "Master Node: Getting next job to assign from list id " << next_job_to_assign_ << " of " << jobs.size() << std::endl;
			return next_job_to_assign_; //not used by callers
		} else if ( outputter->job_has_completed( jobs[ next_job_to_assign_ ] ) && option[ out::overwrite ].value() ) {
			TR << "Master Node: Getting next job to assign from list, overwriting id " << next_job_to_assign_ << " of " << jobs.size() << std::endl;
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

	if ( repeat_job_ == true ) {
		TR << "Slave Node " << rank_ << ": Repeating job id " << current_job_id_ <<std::endl;
		repeat_job_ = false;
	}	else {
		TR << "Slave Node " << rank_ << ": Requesting new job id from master" <<std::endl;
		int empty_data( 0 );
		MPI_Status status;
		current_job_id_ = 0;
		//MPI::COMM_WORLD.Send( &empty_data, 1, MPI::INT, 0, NEW_JOB_ID_TAG );
		//MPI::COMM_WORLD.Recv( &current_job_id_, 1, MPI::INT, 0, NEW_JOB_ID_TAG );
		MPI_Send( &empty_data, 1, MPI_INT, 0, NEW_JOB_ID_TAG, MPI_COMM_WORLD );
		MPI_Recv( &current_job_id_, 1, MPI_INT, 0, NEW_JOB_ID_TAG, MPI_COMM_WORLD, &status );
		TR << "Slave Node " << rank_ << ": Received job id " << current_job_id_ << " from master" <<std::endl;
	}
#endif
	return current_job_id_;
}

///@brief dummy for master/slave version
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
	TR << "Master Node: Mark current job for repetition" << std::endl;
	utility_exit_with_message( "Master Node: master_mark_current_job_id_for_repetition() should never be called" );

}

void
MPIWorkPoolJobDistributor::slave_mark_current_job_id_for_repetition()
{
	runtime_assert( !( rank_ == 0 ) );
	TR << "Slave Node " << rank_ << ": Mark current job for repetition, id " << current_job_id_ << std::endl;
	repeat_job_ = true;
}

///@brief dummy for master/slave version
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

	Jobs const & jobs( get_jobs() );

	std::string const & bad_job_id_input_tag( jobs[ bad_job_id_ ]->input_tag() );

	TR << "Master Node: Job id " << bad_job_id_ << " failed, reporting bad input; other jobs of same input will be canceled: " << job_outputter()->output_name( jobs[ bad_job_id_ ] ) << std::endl;

	while( next_job_to_assign_ <= jobs.size() && jobs[ next_job_to_assign_ ]->input_tag() == bad_job_id_input_tag ) {
		TR << "Master Node: Job canceled without trying due to previous bad input: " << job_outputter()->output_name( jobs[ next_job_to_assign_ ] ) << " id " << next_job_to_assign_ << std::endl;
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
	MPI_Send( &current_job_id_, 1, MPI_INT, 0, BAD_INPUT_TAG, MPI_COMM_WORLD );
#endif
}

///@brief dummy for master/slave version
void
MPIWorkPoolJobDistributor::job_succeeded(core::pose::Pose & pose, core::Real /*run_time*/)
{
  if ( rank_ == 0 ) {
    master_job_succeeded( pose );
  } else {
    slave_job_succeeded( pose );
  }
}

void
MPIWorkPoolJobDistributor::master_job_succeeded(core::pose::Pose & /*pose*/)
{
#ifdef USEMPI
	runtime_assert( rank_ == 0 );
	TR << "Master Node: Job Succeeded" << std::endl;
	utility_exit_with_message( "Master Node: master_job_succeeded() should never be called" );
#endif
}

void
MPIWorkPoolJobDistributor::slave_job_succeeded(core::pose::Pose & MPI_ONLY( pose ) )
{
#ifdef USEMPI
	runtime_assert( !( rank_ == 0 ) );

	if ( option[ OptionKeys::jd2::mpi_fast_nonblocking_output	].value() == true ) {
		job_outputter()->final_pose( current_job(), pose );
	} else {
		int empty_data( 0 );
		MPI_Status status;

		// send job success message to master
		TR << "Slave Node " << rank_ << ": Finished job successfully! Sending output request to master." << std::endl;
		MPI_Send( &current_job_id_, 1, MPI_INT, 0, JOB_SUCCESS_TAG, MPI_COMM_WORLD );

		// receive message from master that says is okay to write
		TR << "Slave Node " << rank_ << ": Received output confirmation from master. Writing output." << std::endl;
		MPI_Recv( &empty_data, 1, MPI_INT, 0, JOB_SUCCESS_TAG, MPI_COMM_WORLD, &status );

		// time and write output (pdb, silent file, score file etc.)
		clock_t starttime = clock();
		job_outputter()->final_pose( current_job(), pose );
		clock_t stoptime = clock();

		// send message to master that we are done outputing
		TR << "Slave Node " << rank_ << ": Finshed writing output in " << ((double) stoptime-starttime) / CLOCKS_PER_SEC << " seconds. Sending message to master" << std::endl;
		MPI_Send( &empty_data, 1, MPI_INT, 0, JOB_SUCCESS_TAG, MPI_COMM_WORLD );
	}
#endif
}


}//jd2
}//protocols
