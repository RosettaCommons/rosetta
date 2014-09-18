// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MPIWorkPartitionJobDistributor.cc
/// @brief  implementation of MPIWorkPartitionJobDistributor - intended for MPI jobs on small numbers of nodes where the load can be balanced equally by the user
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif


// Unit headers
#include <protocols/jd2/MPIWorkPartitionJobDistributor.hh>

//Package headers
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/moves/Mover.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/exit.hh>

// Option headers
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>


///C++ headers
#include <string>

//Auto Headers
#include <utility/vector1.hh>
#include <complex>


static thread_local basic::Tracer TR( "protocols.jd2.MPIWorkPartitionJobDistributor" );

namespace protocols {
namespace jd2 {

///@details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
///which processor it is in MPI land (later used in job determination).  Note that all processors will have the same
///internal Jobs object (set by the parent class); this class merely iterates over it differently.
MPIWorkPartitionJobDistributor::MPIWorkPartitionJobDistributor() :
  JobDistributor(),
	npes_( 1 ),
  rank_( 0 ),
  job_id_start_( 0 ),
  job_id_end_( 0 ),
  next_job_to_try_assigning_( 0 )
{
  // set npes and rank based on whether we are using MPI or not

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	npes_ = option[ OptionKeys::run::nproc ](); //make this same as in "queue"-command in condor script
	rank_ = option[ OptionKeys::run::proc_id ](); //use this as -jd2:condor_rank $(PROCESS)
#ifdef USEMPI
  //npes_ = MPI::COMM_WORLD.Get_size();
  //rank_ = MPI::COMM_WORLD.Get_rank();
	int int_npes, int_rank;																	//don't cast pointers - copy it over instead
	MPI_Comm_rank( MPI_COMM_WORLD, &int_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &int_npes );
	rank_ = int_rank;
	npes_ = int_npes;
#endif

  determine_job_ids_to_run();
  next_job_to_try_assigning_ = job_id_start_;

  Jobs const & jobs( get_jobs() );
  TR << "RANK: " << rank_ << " NUM_PROCS: " << npes_ << " NUM_JOBS: " << jobs.size()
		 << " START_ID: " << job_id_start_ << " END_ID: " << job_id_end_ << std::endl;
}

///@brief dtor
///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
MPIWorkPartitionJobDistributor::~MPIWorkPartitionJobDistributor()
{}

///@details All processors will get the same Jobs object; this function determines which slice belongs to a particular
///processor determined solely by its mpi rank and the number of processors, no communication needed
/// EXAMPLE CASE: 18 jobs, 4 processors
/// processor rank   number of jobs   assigned range in Jobs vector
/// 0                5               1-5
/// 1                5               6-10
/// 2                4               11-14
/// 3                4               15-18
void
MPIWorkPartitionJobDistributor::determine_job_ids_to_run()
{
  Jobs const & jobs( get_jobs() );

  core::Size num_jobs( 0 );
  core::Size jobs_mod_procs( jobs.size() % npes_ );
  core::Real jobs_div_procs( core::Real( jobs.size() ) / core::Real( npes_ ) );

  // calculate number of jobs to run and what id to start at

	// if jobs_mod_procs == 0 (evenly divisible), an equal number of jobs go to each processor. +1 is because rank_ is
	// 0-indexed but Jobs is 1-indexed
  if( jobs_mod_procs == 0 ) {
    num_jobs =  core::Size( jobs_div_procs );
    job_id_start_ = rank_ * num_jobs + 1;
  }
	// if the rank is less than jobs%procs, the number of jobs per processor is the ceiling of jobs/processors; take that
	// many jobs.  +1 is because rank_ is 0-indexed but Jobs is 1-indexed
  else if( rank_ < jobs_mod_procs ) {
    num_jobs = (core::Size) std::ceil( jobs_div_procs );
    job_id_start_ = rank_ * num_jobs + 1;
  }
	// if the rank is more than or equal to jobs%procs, the number of jobs per processor is the floor of jobs/processors;
	// take that many jobs. rank * num jobs accounts for bulk of earlier jobs.  jobs_mod_procs accounts for all processors
	// with rank < jobs_mod_procs getting an extra job because they use ceiling instead of floor in num_jobs.  +1 is
	// because rank_ is 0-indexed but Jobs is 1-indexed
  else if( rank_ >= jobs_mod_procs ) {
    num_jobs = (core::Size) std::floor( jobs_div_procs );
		job_id_start_ = rank_ * num_jobs + jobs_mod_procs + 1;
  }
  else {
    utility_exit_with_message("ERROR: Problem determining job ids to run");
  }

  // calculate job_id_end
  job_id_end_ = job_id_start_ + num_jobs - 1;
}

void
MPIWorkPartitionJobDistributor::go( protocols::moves::MoverOP mover )
{
	go_main( mover );
#ifdef USEMPI
	//MPI::COMM_WORLD.Barrier();
	//MPI::Finalize();
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Finalize();
#endif
}

///@details determine which job to assign next: increment until we run out of available jobs
core::Size
MPIWorkPartitionJobDistributor::get_new_job_id()
{
  Jobs const & jobs( get_jobs() );
  JobOutputterOP outputter = job_outputter();

  while ( next_job_to_try_assigning_ <= job_id_end_ ) {
    if ( outputter->job_has_completed( jobs[ next_job_to_try_assigning_ ] ) &&
      !basic::options::option[ basic::options::OptionKeys::out::overwrite ].value() ) {
      ++next_job_to_try_assigning_;
    } else {
      break;
    }
  }

  if ( next_job_to_try_assigning_ <= job_id_end_ ) {
    core::Size job_to_assign = next_job_to_try_assigning_;
    ++next_job_to_try_assigning_;
    return job_to_assign;
  }

  // indicate that no jobs remain
  return 0;
}

void
MPIWorkPartitionJobDistributor::mark_current_job_id_for_repetition()
{
  runtime_assert( current_job_id() == next_job_to_try_assigning_ - 1 );
  --next_job_to_try_assigning_;
	clear_current_job_output();
}

///@details this function handles the FAIL_BAD_INPUT mover status by removing other jobs with the same input from
///consideration.  This function DOES NOT percolate across processors - so if multiple processors have jobs starting
///with the same bad input, you will get multiple hits through this function.  This is less efficient than it
///theoretically could be (but it's good enough).
void
MPIWorkPartitionJobDistributor::remove_bad_inputs_from_job_list()
{
  std::string const & current_input_tag(current_job()->input_tag());

  TR << "job failed, reporting bad input; other jobs of same input will be canceled: "
     << job_outputter()->output_name( current_job() ) << std::endl;

  Jobs const & jobs( get_jobs() );

  while(next_job_to_try_assigning_ <= job_id_end_ && //MUST BE FIRST for c++ shortcut logical evaluation
    jobs[next_job_to_try_assigning_]->input_tag() == current_input_tag) {
    TR << "job canceled without trying due to previous bad input: "
       << job_outputter()->output_name( jobs[next_job_to_try_assigning_] ) << std::endl;
    ++next_job_to_try_assigning_;
  }
}

}//jd2
}//protocols
