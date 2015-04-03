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
/// @author Oliver Lange olange@u.washington.edu

#ifndef INCLUDED_protocols_jd2_MPIFileBufJobDistributor_hh
#define INCLUDED_protocols_jd2_MPIFileBufJobDistributor_hh

// Unit headers
#include <protocols/jd2/MPIFileBufJobDistributor.fwd.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.fwd.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <core/types.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

/// @brief this tag is used for all communication with JobDistributor ( use this tag to be received in the main MSG-loop of jd2 cf. process_message() )
core::Size const MPI_JOB_DIST_TAG ( 1542 ); //keep unique

/// @brief Tags used to tag messeges sent by MPI functions used to decide whether a slave is requesting a new job id or
///flagging as job as being a bad input
//FileBufjobDistributor messages range 1-100
core::Size const NEW_JOB_ID = 1;
core::Size const BAD_INPUT = 2;
core::Size const JOB_SUCCESS = 3;
core::Size const JOB_FAILED_NO_RETRY = 4;

/// @details This JobDistributor is intended for machines where you have a large number of processors.
/// two dedicated processes are used to handle JobDistribution and File-IO.
/// all other processes (higher rank ) are used for computation.
/// the file_buf_rank_ process runs the MpiFileBuffer which is at the receiving end of all ozstream output that is rerouted via MPI from the
/// slave nodes.
/// This means all slaves write to the same file without FileSystem congestion and interlacing in the file -- IO is handled from a single dedicated process
/// The other dedicated process (master_rank) runs the actual JobDistributor
/// and is only used to distribute jobs to slaves and receive their notification of successful or failed execution
/// in case you have only a small number of processors you can put say 10 MPI processes on 8 processors to get optimal CPU usage.
class MPIFileBufJobDistributor : public JobDistributor
{
	typedef JobDistributor Parent;
protected:
  /// @brief ctor is protected; singleton pattern
  MPIFileBufJobDistributor();

	/// @brief protected ctor for child-classes
	MPIFileBufJobDistributor( core::Size master_rank, core::Size file_buf_rank, core::Size min_client_rank, bool start_empty=false );

	virtual void handle_interrupt() {}

public:
	///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
	///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
  virtual ~MPIFileBufJobDistributor();

	core::Size increment_client_rank(){
		return ++min_client_rank_;
	}

	/// @brief return rank of first worker process (there might be more dedicated processes, e.g., ArchiveManager...)
	core::Size min_client_rank() const {
		return min_client_rank_;
	}


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
	job_succeeded(core::pose::Pose & pose, core::Real runtime, std::string const & tag);

	virtual
	void
	job_failed( core::pose::Pose & pose, bool will_retry );

	friend class JobDistributorFactory; //ctor access

protected:

	//return true if message was understood
	virtual bool process_message(
     core::Size msg_tag,
		 core::Size slave_rank,
		 core::Size slave_job_id,
		 core::Size slave_batch_id,
		 core::Real runtime
	);

	//overloaded so that slave-nodes never automatically switch to next_batch when spinning down.
	virtual bool next_batch();

  /// @brief Handles the receiving of job requests and the sending of job ids to and from slaves
  void master_go( protocols::moves::MoverOP mover );

  /// @brief Always returns zero, simply increments next_job_to_assign_ to the next job that should be run based
  ///on what has been completeted and the overwrite flags
  core::Size master_get_new_job_id();

  /// @brief requests, receives, and returns a new job id from the master node or returns the current job id if the
  ///repeat_job_ flag is set to true
  core::Size slave_get_new_job_id();

  /// @brief This should never be called as this is handled internally by the slave nodes, it utility_exits
  void master_mark_current_job_id_for_repetition();

  /// @brief Sets the repeat_job_ flag to true
  void slave_mark_current_job_id_for_repetition();

  /// @brief Simply increments next_job_to_assign_ to the next job that should be run based on what has been
  ///completed and if the input job tag of the job marked as having bad input
  void master_remove_bad_inputs_from_job_list();

  /// @brief Sends a message to the head node that contains the id of a job that had bad input
  void slave_remove_bad_inputs_from_job_list();

	/// @brief This should never be called as this is handled internally by the slave nodes, it utility_exits
	void master_job_succeeded(core::pose::Pose & pose, std::string const & tag);

	/// @brief Sends a message to the head node upon successful job completion to avoid output interleaving
	void slave_job_succeeded(core::pose::Pose & pose, std::string const & tag);

	/// @brief send a message to master
	void slave_to_master( core::Size tag );

	/// @brief called by master to send and by slave to receive job
	void send_job_to_slave( core::Size slave_rank );

	/// @brief return rank of this process
	core::Size rank() const {
		return rank_;
	}

	/// @brief return rank of master process ( where JobDistributor is running )
	core::Size master_rank() const {
		return master_rank_;
	}

	/// @brief return rank of file-buffer process ( where output data (via ozstream )is handled )
	core::Size file_buf_rank() const {
		return file_buf_rank_;
	}

	/// @brief how many processes --- this includes dedicated processes
	core::Size number_of_processors() {
		return n_rank_;
	}

	/// @brief how many processes --- this includes dedicated processes
	core::Size n_rank() {
		return n_rank_;
	}

	/// @brief how many workers --- important to keep track during spin-down process
	core::Size n_worker() {
		return n_worker_;
	}

	/// @brief how many workers --- important to keep track during spin-down process
	void set_n_worker( core::Size setting ) {
	  n_worker_=setting;
	}

	/// @brief marks job as completed in joblist
	virtual void mark_job_as_completed( core::Size job_id, core::Size batch_id, core::Real runtime );

	/// @brief marks job as bad in joblist
	virtual void mark_job_as_bad( core::Size job_id, core::Size batch_id );

	/// @brief receive a certain signal and ignore it.... this is needed, for instance, when MPIArchiveJobDistributor triggers an
	/// ADD_BATCH signal by sending QUEUE_EMPTY to the ArchiveManager...
	void eat_signal( core::Size signal, int source );

private:

  /// @brief total number of processing elements
	core::Size n_rank_;

	core::Size n_worker_;

  /// @brief rank of the "local" instance
	core::Size rank_;

  /// @brief where slave jobs store current job id
  core::Size slave_current_job_id_; //this overlays current_job_id_ of base class.... BAD

	/// @brief batch_id allow to run multiple batches of jobs -
	core::Size slave_current_batch_id_; //i.e. next_job_to_assign_ is from this batch (for master)

	/// @brief runtime of last job
	core::Real slave_current_runtime_; //i.e. next_job_to_assign_ is from this batch (for master)

  /// @brief where master stores next job to assign (in a good state after get_new_job_id up until it's used)
  //core::Size next_job_to_assign_;

  /// @brief where master temporarily stores id of jobs with bad input
  core::Size bad_job_id_;

	/// @brief where slave stores whether it should repeat its current job id
	bool repeat_job_;

	/// @brief keep some statistics about the jobs
	/// this is mostly just for silly tr.Info messages...
	//// but we need some of this to properly spin-down at the end

	/// @brief jobs send to slave-nodes
	core::Size jobs_assigned_;

	/// @brief jobs that have returned (either, bad or good )
	core::Size jobs_returned_;

	/// @brief jobs that have returned bad due to status BAD_INPUT
	core::Size bad_input_jobs_;

	/// @brief how many more to spin down
	core::Size n_nodes_left_to_spin_down_;

	/// @brief keep here the ranks of different functional processes

	/// @brief the job-distributor (master)
	core::Size const master_rank_; //1

	/// @brief the File-Buffer
	core::Size const file_buf_rank_; //0

	/// @brief the first slave node
	//core::Size const min_client_rank_; //2 or 3 ...
	core::Size min_client_rank_; //2 or 3 ... //ek made non-const

	/// @brief keep track of average timings for time-outs
	core::Real cumulated_runtime_;

	core::Size cumulated_jobs_;

};

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_MPIFileBufJobDistributor_HH
