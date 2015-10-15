// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MPIWorkPoolJobDistributor.hh
/// @brief  header for MPIWorkPoolJobDistributor - intended for MPI jobs on large numbers of nodes where the head node is dedicated to handing out new job ids
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_protocols_jd2_MPIWorkPoolJobDistributor_hh
#define INCLUDED_protocols_jd2_MPIWorkPoolJobDistributor_hh

// Unit headers
#include <protocols/jd2/MPIWorkPoolJobDistributor.fwd.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.fwd.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <core/types.hh>

// C++ headers
#include <string>

#include <platform/types.hh>
#include <protocols/jd2/InnerJob.fwd.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/jd2/JobInputter.fwd.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/jd2/Parser.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <basic/mpi/mpi_enums.hh>
#include <algorithm>
#include <utility/assert.hh>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <vector>


namespace protocols {
namespace jd2 {


/// @details This job distributor is meant for running jobs where the machine you are using has a large number of
///processors, the number of jobs is much greater than the number of processors, or the runtimes of the individual jobs
///could vary greatly. It dedicates the head node (whichever processor gets processor rank #0) to handling job requests
///from the slave nodes (all nonzero ranks). Unlike the MPIWorkPartitionJobDistributor, this JD will not work at all
///without MPI and the implementations of all but the interface functions have been put inside of ifdef directives.
///Generally each function has a master and slave version, and the interface functions call one or the other depending
///on processor rank.
class MPIWorkPoolJobDistributor : public JobDistributor
{
protected:
	/// @brief ctor is protected; singleton pattern
	MPIWorkPoolJobDistributor();

	virtual void handle_interrupt() {}

public:
	///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
	///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
	virtual ~MPIWorkPoolJobDistributor();

	/// @brief dummy for master/slave version
	virtual
	void
	go( protocols::moves::MoverOP mover );

	/// @brief dummy for master/slave version
	virtual
	core::Size
	get_new_job_id();

	/// @brief dummy for master/slave version
	virtual
	void
	mark_current_job_id_for_repetition();


	/// @brief dummy for master/slave version
	virtual
	void
	remove_bad_inputs_from_job_list();

	/// @brief dummy for master/slave version
	virtual
	void
	job_succeeded(core::pose::Pose & pose, core::Real run_time, std::string const & tag);

	/// @brief Called if job fails.
	///
	virtual
	void job_failed(core::pose::Pose &pose, bool will_retry);

	/// @brief should the go() function call MPI_finalize()? It probably should, this is true by default.
	virtual
	void mpi_finalize(bool finalize);

	friend class JobDistributorFactory; //ctor access

protected:

	/// @brief Handles the receiving of job requests and the sending of job ids to and from slaves
	virtual
	void
	master_go( protocols::moves::MoverOP mover );

	/// @brief Proceeds to the parent class go_main() as usual
	virtual
	void
	slave_go( protocols::moves::MoverOP mover );

	/// @brief Always returns zero, simply increments next_job_to_assign_ to the next job that should be run based
	///on what has been completeted and the overwrite flags
	virtual
	core::Size
	master_get_new_job_id();

	/// @brief requests, receives, and returns a new job id from the master node or returns the current job id if the
	///repeat_job_ flag is set to true
	virtual
	core::Size
	slave_get_new_job_id();

	/// @brief This should never be called as this is handled internally by the slave nodes, it utility_exits
	virtual
	void
	master_mark_current_job_id_for_repetition();

	/// @brief Sets the repeat_job_ flag to true
	virtual
	void
	slave_mark_current_job_id_for_repetition();

	/// @brief Simply increments next_job_to_assign_ to the next job that should be run based on what has been
	///completed and if the input job tag of the job marked as having bad input
	virtual
	void
	master_remove_bad_inputs_from_job_list();

	/// @brief Sends a message to the head node that contains the id of a job that had bad input
	virtual
	void
	slave_remove_bad_inputs_from_job_list();

	/// @brief This should never be called as this is handled internally by the slave nodes, it utility_exits
	virtual
	void
	master_job_succeeded(core::pose::Pose & pose, std::string const & tag);

	/// @brief Sends a message to the head node upon successful job completion to avoid output interleaving
	virtual
	void
	slave_job_succeeded(core::pose::Pose & pose, std::string const & tag);

	/// @brief Mark the job as completed/deletable in the jobs list on the master process.
	///
	virtual
	void
	master_mark_job_as_completed( core::Size const job_index );

	/// @brief Mark the job as deletable in the jobs list on the master process.
	///
	virtual
	void
	master_mark_job_as_failed( core::Size const job_index );

	/// @brief Set whether the JobDistributor sends jobs to each slave in sequence (1, 2, 3, etc.)
	///
	virtual
	void set_sequential_distribution( bool const val ) { sequential_distribution_ = val; return; }

	/// @brief Get whether the JobDistributor sends jobs to each slave in sequence (1, 2, 3, etc.)?
	///
	virtual
	bool sequential_distribution() const { return sequential_distribution_; }

	/// @brief Set whether this is the process that should start requesting jobs (be the first for sequential distribution).
	///
	virtual
	void set_starter_for_sequential_distribution( bool const val ) { starter_for_sequential_distribution_ = val; return; }

	/// @brief Is this the process that should start requesting jobs (be the first for sequential distribution)?
	///
	virtual
	bool starter_for_sequential_distribution() const { return starter_for_sequential_distribution_; }

	/// @brief Wait for a signal from the n-1 process saying that I can proceed.
	///
	virtual
	void wait_for_go_signal() const;

	/// @brief Send a signal to the n+1 process saying that it can proceed.
	/// @details This also sets starter_for_sequential_distribution_ to false, since we no longer
	/// want this process to refrain from waiting.
	virtual
	void send_go_signal();

protected:

	/// @brief total number of processing elements
	core::Size npes_;

	/// @brief rank of the "local" instance
	core::Size rank_;

	/// @brief where slave jobs store current job id
	core::Size current_job_id_;

	/// @brief where master stores next job to assign (in a good state after get_new_job_id up until it's used)
	core::Size next_job_to_assign_;

	/// @brief where master temporarily stores id of jobs with bad input
	core::Size bad_job_id_;

	/// @brief where slave stores whether it should repeat its current job id
	bool repeat_job_;

	/// @brief should the go() function call MPI_finalize?  There are very few cases where this should be false
	bool finalize_MPI_;

	/// @brief Should the JobDistributor send jobs to each slave in sequence (1, 2, 3, etc.)?  Default false -- slaves request jobs as they
	/// become available.
	bool sequential_distribution_;

	/// @brief Is this the process that should start requesting jobs (be the first for sequential distribution)?
	/// @details Default false
	bool starter_for_sequential_distribution_;

};

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_MPIWorkPoolJobDistributor_HH
