// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/JobDistributor.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Base class
/// @author Andrew Leaver-Fay
/// @author Steven Lewis smlewi@gmail.com
/// @author Modified by Sergey Lyskov

#ifndef INCLUDED_protocols_jd2_JobDistributor_hh
#define INCLUDED_protocols_jd2_JobDistributor_hh

// Unit headers
#include <protocols/jd2/JobDistributor.fwd.hh>

// Package headers
#include <protocols/jd2/JobInputter.fwd.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/jd2/Parser.fwd.hh>

#include <protocols/moves/Mover.fwd.hh>

#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <core/types.hh>

#include <protocols/jd2/Job.fwd.hh>
#include <utility/vector1.hh>
#include <string>

#ifdef WIN32
	#include <protocols/jd2/Job.hh>
#endif

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <thread>
#endif
#endif

namespace protocols {
namespace jd2 {

///Enforced JobDistributor destruction turned out to cause problems - calls to Tracers and/or the Options system
///during destruction caused crashes if those systems had been destructed first.  So this is deprecated.
// simple class to ensure that JobDistributor objects are destroyed at program exit.
// class JobDistributorDestroyer {
// public:
// 	JobDistributorDestroyer(JobDistributor* = 0);
// 	~JobDistributorDestroyer();
// 	void set_job_distributor(JobDistributor* s);

// 	JobDistributor* jd_;
// };

class JobDistributor {

protected:
	/// @brief Singleton instantiation pattern; Derived classes will call default ctor, but their ctors, too must be
	/// protected (and the JDFactory must be their friend.)
	JobDistributor();

	/// @brief MPIArchiveJobDistributor starts with an empty job-list...
	JobDistributor( bool empty );

private:
	//the actual c'tor work is done here
	void init_jd();

public:

	static
	JobDistributor *
	get_instance();

	// NO LONGER called as a result of the static JobDistributorDestroyer object declared in JobDistributor.cc.
	///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
	///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
	virtual
	~JobDistributor();

public:

	/// @brief This may be overridden by derived classes.  Default implementation invokes go_main.
	virtual
	void
	go( protocols::moves::MoverOP mover );

	/// @brief invokes go, after setting JobOutputter
	void
	go( protocols::moves::MoverOP mover, JobOutputterOP jo );

	/// @brief Movers may ask their controlling job distributor for information about the current job. They may also load
	/// information into this job for later output.
	JobOP
	current_job() const;

	/// @brief Movers may ask their controlling job distributor for the output name as defined by the Job and JobOutputter.
	std::string
	current_output_name() const;

	/// @brief Movers (or derived classes) may ask for the JobOutputter
	JobOutputterOP
	job_outputter() const;

	/// @brief Movers (or derived classes) may ask for the JobOutputter
	void set_job_outputter( const JobOutputterOP &new_job_outputter );

	/// @brief JobInputter access
	JobInputterOP
	job_inputter() const;

	///@brief should the go() function call MPI_finalize()? It probably should, this is true by default.
	virtual
	void mpi_finalize(bool finalize);

	/// @brief The input source for the current JobInputter.
	JobInputterInputSource::Enum
	job_inputter_input_source() const;

	friend class JobDistributorFactory; // calls private ctor

	virtual void restart();

	core::Size total_nr_jobs() const {
		return jobs_.size();
	}


	/// @brief integer access - which job are we on?
	core::Size current_job_id() const;

protected:
	/// @brief Non-virtual get-job, run it, & output loop.  This function is pretty generic and your subclass may be able
	/// to use it. It is NOT virtual - this implementation can be shared by (at least) the simple
	/// FileSystemJobDistributor, the MPIWorkPoolJobDistributor, and the MPIWorkPartitionJobDistributor.  Do not feel that
	/// you need to use it as-is in your class - but DO plan on implementing all its functionality!

	void
	go_main( protocols::moves::MoverOP mover );

	/// Read access to private data for derived classes.

	/// @brief Jobs is the container of Job objects
	Jobs const &
	get_jobs() const;

	/// @brief Jobs is the container of Job objects
	/// need non-const to mark Jobs as completed on Master in MPI-JobDistributor
	//	Jobs&
	//	get_jobs(); //get dedicated accessor instead
	void mark_job_as_completed( core::Size job_id, core::Real run_time );

	void mark_job_as_bad( core::Size job_id );


	/// @brief Parser access
	ParserOP
	parser() const;

	void begin_critical_section();

	void end_critical_section();


protected:
	/// @brief this function updates the current_job_id_ and current_job_ fields.  The boolean return states whether or not
	///a new job was obtained (if false, quit distributing!)
	bool obtain_new_job( bool re_consider_current_job = false ); //if true we check if current_job is still selectable (after remove_bad_input)

	/// @brief Return 0 to signal that no available jobs remain.  Otherwise return an index into the Jobs object.
	virtual
	core::Size
	get_new_job_id() = 0;

	/// @brief This function is called upon a successful job completion; it has been virtualized so BOINC and MPI can delay/protect output
	///base implementation is just a call to the job outputter
	virtual
	void
	job_succeeded( core::pose::Pose & pose, core::Real run_time );

	/// @brief This function is called upon a successful job completion if there are additional poses generated by the mover
	///base implementation is just a call to the job outputter
	virtual
	void
	job_succeeded_additional_output( core::pose::Pose & pose, std::string const & tag );

	/// @brief This function is called when we give up on the job;  it has been virtualized so BOINC and MPI can delay/protect output
	///base implementation is just a call to the job outputter
	virtual
	void
	job_failed( core::pose::Pose & /*pose*/, bool /*will_retry*/ );

	/// @brief this function is called whenever a job "soft-fails" and needs to be retried.  Generally it should ensure
	///that the subsequent call to obtain_new_job returns this job over again.
	virtual
	void
	mark_current_job_id_for_repetition() = 0;

	/// @brief this function is called inside go_main if a mover returns FAIL_BAD_INPUT.  Its purpose is to remove other
	///jobs with the same input (which will still be bad) from the available list of jobs.  The default implementation is
	///a no-op and not all distributors are expected/able to implement this functionality, only those that can guaruntee
	///no other jobs of that input are currently running.
	virtual
	void
	remove_bad_inputs_from_job_list();

	/// @brief Derived classes are allowed to clean up any temporary files or data relating to the current job after the
	/// current job has completed.  Called inside go_main loop.  Default implementation is a no-op.
	virtual
	void
	current_job_finished();

	/// @brief Derived classes are allowed to perform some kind of action when the job distributor runs out of jobs to
	/// execute.  Called inside go_main.  Default implementation is a no-op.
	virtual
	void
	note_all_jobs_finished();

	void
	clear_current_job_output();

	/// @brief This function got called when job is not yet finished and got termitated abnormaly (ctrl-c, kill etc).
	///        when implimenting it in subclasses make sure to delete all in-progress-data that your job spawn.
	virtual void handle_interrupt() = 0;


private:

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static JobDistributor * create_singleton_instance();

	JobInputterOP job_inputter_;
	JobOutputterOP job_outputter_;
	ParserOP parser_;

	Jobs jobs_;
	/// @brief pointer to current job.  Information is somewhat duplicated with current_job_id_.
	JobOP current_job_;

	/// @brief access into jobs_ vector indicating current job.  Contains more information than current_job_ in that it can be incremented...
	core::Size current_job_id_;

	/// @brief access into jobs_ bector indicating the previous job.  Used with the -jd2:delete_old_poses option for deleting unnecessary poses
	core::Size last_completed_job_;

	static JobDistributor * instance_;

	///BATCH interface:
	/// @details the BATCH interface of the JobDistributor is used to allow consecutive execution of a set of jobs with different flags
	/// different batches behave like completely independent rosetta runs --- but of course a number of processes can already work on
	/// a new batch, while others are still finishing the last jobs of the previous batch.
	/// run from command-line with -run:batches flag1 flag2 flag3
	/// the flag1 flag2... point to @flag1 files that are added to all other flags ( and removed at end of batch )
	///  you can have all output in same output file or ( by redefining e.g. -out:file:silent in each batch-flag file ) in different output files

public:

	/// @brief what is the current batch ? --- name refers to the flag-file used for this batch
	std::string get_current_batch() const;

	/// @brief add a new batch ( name will be interpreted as flag_file )
	// positive id means we want to set a particular batch id ...
	// fill-up with BOGUS_BATCH_ID
	// if current_batch_id > id this will not have any effect... --> error?
	virtual void add_batch( std::string const&, core::Size id = 0 );

	/// @brief what is the current batch number ? --- refers to position in batches_
	core::Size current_batch_id() const {
				return current_batch_id_;
	}

protected:

	/// @brief set current_batch_id  --- eg for slave nodes in MPI framework
	void set_batch_id( core::Size setting );

	/// @brief switch current_batch_id_ to next batch
	virtual bool next_batch();

	/// @brief if end of batches_ reached via next_batch or set_batch_id ...
	virtual void batch_underflow() {}; //no action for base-class

	/// @brief called by next_batch() or set_batch_id() to switch-over and restart JobDistributor on new batch
	virtual void load_new_batch();

	/// @brief how many batches are in our list ... this can change dynamically
	core::Size nr_batches() const {
		return batches_.size();
	}

	/// @brief give name of batch with given id
	std::string const& batch( core::Size batch_id ) {
		return batches_[ batch_id ];
	}


protected:
	/// @brief Setting up callback function that will be call when our process is about to terminate.
    ///         This will allow us to exit propely (clean up in_progress_files/tmp files if any).
	static void setup_system_signal_handler( void (*prev_fn)(int) = jd2_signal_handler);

	/// @brief Set signal handler back to default state.
	static void remove_system_signal_handler();

	/// @brief Default callback function for signal handling
	static void jd2_signal_handler(int Signal);

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

private:

	/// @brief read -run:batches
	void populate_batch_list_from_cmd();

	/// @brief current_batch or 0
	core::Size current_batch_id_;

	/// @brief all batches if present
	utility::vector1< std::string > batches_;
};

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_JobDistributor_HH
