// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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


static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.MPIWorkPartitionJobDistributor" );

namespace protocols {
namespace jd2 {

/// @details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
///which processor it is in MPI land (later used in job determination).  Note that all processors will have the same
///internal Jobs object (set by the parent class); this class merely iterates over it differently.
MPIWorkPartitionJobDistributor::MPIWorkPartitionJobDistributor() :
	JobDistributor(),
	npes_( 1 ),
	rank_( 0 ),
	next_job_to_try_assigning_( 1 )
{
	// set npes and rank based on whether we are using MPI or not

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	npes_ = option[ OptionKeys::run::nproc ](); //make this same as in "queue"-command in condor script
	rank_ = option[ OptionKeys::run::proc_id ](); //use this as -jd2:condor_rank $(PROCESS)
#ifdef USEMPI
  npes_ = MPI::COMM_WORLD.Get_size();
  rank_ = MPI::COMM_WORLD.Get_rank();
#endif

	next_job_to_try_assigning_ = rank_ + 1;

	TR << "RANK: " << rank_ << " NUM_PROCS: " << npes_ << " NUM_JOBS: " << get_jobs().size() << std::endl;
}

/// @brief dtor
///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
MPIWorkPartitionJobDistributor::~MPIWorkPartitionJobDistributor() = default;


void
MPIWorkPartitionJobDistributor::go( protocols::moves::MoverOP mover )
{
	go_main( mover );
#ifdef USEMPI
	MPI::COMM_WORLD.Barrier();
	MPI::Finalize();
#endif
}

/// @details determine which job to assign next: increment until we run out of available jobs
core::Size
MPIWorkPartitionJobDistributor::get_new_job_id()
{
	JobsContainer const & jobs( get_jobs() );
	JobOutputterOP outputter = job_outputter();

	// if the overwrite flag is false, advance to the next jobs we should work on
	if ( !basic::options::option[ basic::options::OptionKeys::out::overwrite ].value() ) {
		while ( next_job_to_try_assigning_ <= jobs.size() && outputter->job_has_completed( jobs[ next_job_to_try_assigning_ ] ) ) {
			next_job_to_try_assigning_ += npes_;
		}
	}

	if ( next_job_to_try_assigning_ <= jobs.size() ) {
		core::Size job_to_assign( next_job_to_try_assigning_ );
		next_job_to_try_assigning_ += npes_;
		return job_to_assign;
	}

	// no more jobs left
	return 0;
}

void
MPIWorkPartitionJobDistributor::mark_current_job_id_for_repetition()
{
	runtime_assert( current_job_id() == next_job_to_try_assigning_ - npes_ );
	next_job_to_try_assigning_ -= npes_;
	clear_current_job_output();
}

/// @details this function handles the FAIL_BAD_INPUT mover status by removing other jobs with the same input from
///consideration.  This function DOES NOT percolate across processors - so if multiple processors have jobs starting
///with the same bad input, you will get multiple hits through this function.  This is less efficient than it
///theoretically could be (but it's good enough).
void
MPIWorkPartitionJobDistributor::remove_bad_inputs_from_job_list()
{
	std::string const & current_input_tag(current_job()->input_tag());

	TR << "job failed, reporting bad input; other jobs of same input will be canceled: "
		<< job_outputter()->output_name( current_job() ) << std::endl;

	JobsContainer const & jobs( get_jobs() );

	while ( next_job_to_try_assigning_ <= jobs.size() && //MUST BE FIRST for c++ shortcut logical evaluation
			jobs[next_job_to_try_assigning_]->input_tag() == current_input_tag ) {
		TR << "job canceled without trying due to previous bad input: "
			<< job_outputter()->output_name( jobs[next_job_to_try_assigning_] ) << std::endl;
		next_job_to_try_assigning_ += npes_;
	}
}

}//jd2
}//protocols
