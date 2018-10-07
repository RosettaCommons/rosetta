// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/MPIWorkPartitionJobDistributor.cc
/// @brief  jd3 implementation of MPIWorkPartitionJobDistributor - intended for MPI jobs on small numbers of nodes where the load can be balanced equally by the user
/// @author P. Douglas Renfrew (renfrew@nyu.edu)
/// @author Andy Watkins (amw579@nyu.edu)

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif


// Unit headers
#include <protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.hh>

//Package headers
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/LarvalJob.hh>

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


static basic::Tracer TR( "protocols.jd2.MPIWorkPartitionJobDistributor" );

namespace protocols {
namespace jd3 {
namespace job_distributors {

/// @details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
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
	int int_npes, int_rank;																	//don't cast pointers - copy it over instead
	MPI_Comm_rank( MPI_COMM_WORLD, &int_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &int_npes );
	rank_ = int_rank;
	npes_ = int_npes;

#endif

	// AMW: Have to do this per round, I think.
	//determine_job_ids_to_run();
	//next_job_to_try_assigning_ = job_id_start_;

	//Jobs const & jobs( get_jobs() );
	//TR << "RANK: " << rank_ << " NUM_PROCS: " << npes_ //<< " NUM_JOBS: " << jobs.size()
	//   << " START_ID: " << job_id_start_ << " END_ID: " << job_id_end_ << std::endl;
}

/// @brief dtor
MPIWorkPartitionJobDistributor::~MPIWorkPartitionJobDistributor()
{
	// ideally these would be called in the dtor but the way we have the singleton pattern set up the dtors don't get
	// called
	// AMW: this is no longer a singleton, so this code can go here.
#ifdef USEMPI
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Finalize();
#endif
}

/// @details All processors will get the same Jobs object; this function determines which slice belongs to a particular
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
	Size total_njobs( njobs_for_round() );

	core::Size num_jobs( 0 );
	core::Size jobs_mod_procs( total_njobs % npes_ );
	core::Real jobs_div_procs( core::Real( total_njobs ) / core::Real( npes_ ) );

	// calculate number of jobs to run and what id to start at

	// if jobs_mod_procs == 0 (evenly divisible), an equal number of jobs go to each processor. +1 is because rank_ is
	// 0-indexed but Jobs is 1-indexed
	/*next_job_index_ = */

	if ( jobs_mod_procs == 0 ) {
		num_jobs = core::Size( jobs_div_procs );
	} else if ( rank_ < jobs_mod_procs ) {
		// if the rank is less than jobs%procs, the number of jobs per processor is the ceiling of jobs/processors; take that
		// many jobs.  +1 is because rank_ is 0-indexed but Jobs is 1-indexed
		num_jobs = (core::Size) std::ceil( jobs_div_procs );
	} else if ( rank_ >= jobs_mod_procs ) {
		// if the rank is more than or equal to jobs%procs, the number of jobs per processor is the floor of jobs/processors;
		// take that many jobs. rank * num jobs accounts for bulk of earlier jobs.  jobs_mod_procs accounts for all processors
		// with rank < jobs_mod_procs getting an extra job because they use ceiling instead of floor in num_jobs.  +1 is
		// because rank_ is 0-indexed but Jobs is 1-indexed
		num_jobs = job_id_start_ + (core::Size) std::floor( jobs_div_procs );
	} else {
		utility_exit_with_message("ERROR: Problem determining job ids to run");
	}

	// calculate job_id_end
	job_id_start_ = rank_ * num_jobs + 1;
	job_id_end_ = job_id_start_ + num_jobs - 1;
}

void
MPIWorkPartitionJobDistributor::go( JobQueenOP queen )
{
	set_job_queen( queen );

	// AMW: Right now, this is taken directly from JobDistributor, but
	// I think I am going to need to change stuff to account for distributing
	// over nodes.
	do {
		store_jobs_for_current_round( determine_jobs_for_next_round() );
		determine_job_ids_to_run();
		next_job_to_try_assigning_ = job_id_start_;
		TR << "RANK: " << rank_ << " NUM_PROCS: " << npes_ //<< " NUM_JOBS: " << jobs.size()
			<< " START_ID: " << job_id_start_ << " END_ID: " << job_id_end_ << std::endl;

		while ( more_jobs_in_current_round() ) {

			// select the next job to run
			LarvalJobOP larval_job = select_next_job();
			if ( ! larval_job ) break; // or if we discover there are no jobs to run, then quit

			// ask the job queen to mature the job
			JobOP mature_job = job_queen().mature_larval_job( larval_job );
			if ( ! mature_job ) {
				// signal from job_queen_ that the inputs for the job are bad
				purge_similar_jobs_which_have_bad_inputs( larval_job );
				continue;
			}

			// run the job
			JobResultOP result;
			try {
				result = mature_job->run();
			} catch ( utility::excn::Exception const & exception ) {
				process_exception_from_job( larval_job, exception );
			}
			process_job_result( larval_job, result );

		}
		note_round_completed();

	} while ( another_round_remains() );

}

LarvalJobOP
MPIWorkPartitionJobDistributor::select_next_job() {
	while ( next_job_to_try_assigning_ <= job_id_end_ ) {
		LarvalJobOP next_job = job_queen().determine_job_list()[ next_job_to_try_assigning_ ];
		++next_job_to_try_assigning_;
		if ( next_job->bad() ) continue;

		// determine if this job has already been run
		if ( job_queen().has_job_completed( next_job ) ) {
			continue;
		}
		// ok -- run this job
		job_queen().mark_job_as_having_begun( next_job );
		return next_job;
	}

	// return 0 if there are no jobs left to run
	return LarvalJobOP( 0 );
}

bool
MPIWorkPartitionJobDistributor::more_jobs_in_current_round() {
	return next_job_to_try_assigning_ <= job_id_end_;
}

/// @details determine which job to assign next: increment until we run out of available jobs
core::Size
MPIWorkPartitionJobDistributor::get_new_job_id()
{
	LarvalJobs const & larval_jobs( jobs_for_current_round() );

	while ( next_job_to_try_assigning_ <= job_id_end_ ) {
		if ( job_queen().has_job_completed( larval_jobs[ next_job_to_try_assigning_ ] ) &&
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

/*void
MPIWorkPartitionJobDistributor::mark_current_job_id_for_repetition()
{
runtime_assert( current_job_id() == next_job_to_try_assigning_ - 1 );
--next_job_to_try_assigning_;
clear_current_job_output();
}*/

/// @details this function handles the FAIL_BAD_INPUT mover status by removing other jobs with the same input from
///consideration.  This function DOES NOT percolate across processors - so if multiple processors have jobs starting
///with the same bad input, you will get multiple hits through this function.  This is less efficient than it
///theoretically could be (but it's good enough).

// AMW I think parent class handles this
/*
void
MPIWorkPartitionJobDistributor::remove_bad_inputs_from_job_list()
{
LarvalJobs const & larval_jobs( jobs_for_current_round() );
std::string const & current_input_tag(current_job()->input_tag());

TR << "job failed, reporting bad input; other jobs of same input will be canceled: "
<< job_outputter()->output_name( current_job() ) << std::endl;

// assumes JobQueen is known
Jobs const & jobs( get_jobs() );

while(next_job_to_try_assigning_ <= job_id_end_ && //MUST BE FIRST for c++ shortcut logical evaluation
jobs[next_job_to_try_assigning_]->input_tag() == current_input_tag) {
TR << "job canceled without trying due to previous bad input: "
<< job_outputter()->output_name( jobs[next_job_to_try_assigning_] ) << std::endl;
++next_job_to_try_assigning_;
}
}
*/

}//job_distributors
}//jd3
}//protocols
