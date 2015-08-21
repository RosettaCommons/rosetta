// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jobdist/JobDistributors.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_jobdist_JobDistributors_hh
#define INCLUDED_protocols_jobdist_JobDistributors_hh

#define TAG_NUM_FORMAT_LENGTH 8

#ifdef USEMPI
#include <mpi.h>
#endif

// Project headers
// This has to come before boinc.hh or we get this error on VC++
// '_read' : is not a member of 'std::basic_istream<_Elem,_Traits>'

#ifdef BOINC
#include <protocols/boinc/boinc.hh>
#endif // BOINC

#include <protocols/jobdist/JobDistributors.fwd.hh>


#include <core/types.hh>
#include <basic/Tracer.hh>


#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL headers

#include <map>
#include <set>
#include <sstream>
#include <string>


// option key includes


#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


#ifdef BOINC
#ifdef USEMPI
Throw a compiler error because MPI and BOINC cannot be used together!
If you got this message, something is wrong with your build settings.
#endif
#endif


namespace protocols {
namespace jobdist {

/// @brief Coordinates processing of jobs across multiple Rosetta processes.
///
/// @details
/// Job distributors need to be customized in three different ways:
///  - by cluster architecture (none/Condor, MPI, BOINC, etc)
///  - by local test for job completion (PDB exists, tag already in silent file, etc)
///  - by type of input (single PDB file, pair of PDB files, list of tag IDs, etc)
///
/// Cluster architecture is a GLOBAL COMPILE-TIME decision:  it's handled by scons,
/// it's the same for all executables built at that time, and it should be implemented
/// using ifdef's in this base class, by modifying the next_job() method directly.
///
/// Test for job completion is a PER-EXECUTABLE decision:  it's handled
/// by subclassing BaseJobDistributor and implementing the is_finished() method.
/// BaseJobDistributor will consult is_finished() in whatever way is appropriate
/// for the current cluster architecture to ensure that jobs are not repeated.
///
/// Type of input is handled by templating the job distributor on a Job object.
/// BasicJob has been provided already, but you can subclass it if you need
/// to carry around additional information about the inputs.
///
class BaseJobDistributor : public utility::pointer::ReferenceCount
{
protected:

	typedef utility::vector1< BasicJobOP > JobVector;

public:

	BaseJobDistributor(JobVector jobs);
	BaseJobDistributor( BaseJobDistributor const & );
	virtual ~BaseJobDistributor();

	/// @brief If true, sets the next Job and nstruct number to be processed.
	/// Deliberately not virtual:  should not be overriden.  Uses the "find_available_job"
	/// method, which is common to both MPI and standard protocols, but used in slightly
	/// different manners.
	bool next_job(BasicJobOP & job, int & struct_n);

	/// @brief Must be called by client before first call to next_job().
	/// If overriden by a subclass, it MUST call the superclass implementation.
	virtual void startup();

	/// @brief Must be called by client after last call to next_job().
	/// If overriden by a subclass, it MUST call the superclass implementation.
	virtual void shutdown();

	/// @brief Signal that if at all possible, we would like to not be killed while in the critical section.
	/// If overriden by a subclass, it MUST call the superclass implementation.
	virtual void begin_critical_section();

	/// @brief Signal that if at all possible, we would like to not be killed while in the critical section.
	/// If overriden by a subclass, it MUST call the superclass implementation.
	virtual void end_critical_section();

	/// @brief Virtual function for dump_pose that is needed for main_plain_mover
	virtual void dump_pose_and_map( std::string const &, core::pose::Pose & );

	/// @brief Virtual function for temp_file main_plain_mover
	virtual void temp_file( std::string const & );

	/// @brief Virtual function for score_map that is needed for main_plain_mover
	/// sets the score_map
	virtual void score_map( std::map < std::string, core::Real> & );

	void disable_ignorefinished(){ ignorefinished_ = false; }
	void enable_ignorefinished(){ ignorefinished_ = true; }

	void disable_output(){ nooutput_ = true; }
	void enable_output(){ nooutput_ = false; }

	void disable_inprogress(){ inprogress_ = false; }
	void enable_inprogress(){ inprogress_ = true; }

	bool ignorefinished() const { return ignorefinished_; }
	bool nooutput() const { return nooutput_; }
	bool inprogress() const { return inprogress_; }

	void set_proc_id( core::Size proc_id, core::Size nproc ) {
		proc_id_ = proc_id;
		nproc_ = nproc;
	}
	/// @brief get output_tag for current job's current nstruct
	/// @details by default return current_job's current_nstruct' output_tag.
	///Overridden in derived PlainSilentFileJobDistributor to return names with "S_" at the beginning
	virtual std::string get_current_output_tag();

	virtual std::string get_output_filename();

protected:

	/// @brief Is the given nstruct number of the given input job already finished by the local process?
	/// To be implemented by subclasses.
	virtual bool is_finished(BasicJobOP const & job, int struct_n) = 0;

	/// @brief Restore state from checkpoint file, if it exists
	virtual void checkpoint_read();
	/// @brief Save state to checkpoint file, overwriting previous
	virtual void checkpoint_write();
	/// @brief Remove checkpoint file (at end of batch)
	virtual void checkpoint_clear();

	/// @brief accessor for current_nstruct_
	int current_nstruct() { return current_nstruct_; }
	/// @brief accessor for current_job owning pointer
	BasicJobOP current_job();

#ifdef USEMPI
	/// @brief Check that a call to MPI_Init has occurred already -- for use in runtime_assert statements
	bool MPI_has_been_initialized() const;

	/// @brief read access to derived classes for MPI related (const) data
	/// @details must not be called until startup has been called
	int mpi_rank() const;

	/// @brief read access to derived classes for MPI related (const) data
	/// @details must not be called until startup has been called.
	int mpi_nprocs() const;

	/// @brief Node 0 does no work, rather, it tells the other n-1 proccessors
	/// which jobs they should work on.  When all the work has been distributed, it tells the
	/// slave nodes to spin down.  Once all other nodes have been spun-down, it quits.
	void master_node_distribute_jobs();

	//// @brief Request a job from node 0. If node 0 has no work left to be done, returns false.
	bool request_job_from_master_node();

#endif

private:
	/// @brief looks for a job that has not yet been started, and stores the index for the job, and the
	/// nstruct index in the member variables.  Called by next_job() and by master_node_distribute_jobs().
	bool find_available_job();

private:

	bool const overwrite_; // = basic::options::option[ basic::options::OptionKeys::out::overwrite ];
	JobVector jobs_;
	core::Size current_job_;
	int current_nstruct_;
	bool is_started_;

	//for simple distribution  without communication
	//this Jobdistributor will only process job-nr if ( job_nr mod nprocs_ == proc_id )
	Size nproc_;
	Size proc_id_;
	Size curr_jobid_; //running number


#ifdef USEMPI
	/// Private data, readable by the derived class -- valid only after startup() has been called.
	int mpi_rank_;
	int mpi_nprocs_;
	/// Private data not to be read by the derived class
	MPI_Status stat_;
	int tag_;
#endif

protected:
	/// ignore already done jobs - redo everything !
	bool ignorefinished_;

	/// do not write files - useful such things as statistics or rescoring !
	bool nooutput_;

	/// write .in_progress files for multiple processor runs
	bool inprogress_;

	/// starttime
	int start_time_;

	/// RandomStore - needed for


	int random_counter_;
	utility::vector1< double > random_store_;

	int get_next_random_range(int low, int high);

}; // BaseJobDistributor


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Distributor for use with atomtree_diff silent files.
/// @details This class is deliberately designed for each process to work on
/// its own silent file (and preferrably in its own directory);
/// unlike Rosetta++, we DO NOT share silent files among multiple processes.
/// First, true atomic file locking is impossible on most distributed file systems (e.g. NFS),
/// meaning that files may be corrupted by multiple simultaneous writes.
/// For long-running processes, collisions may be rare,
/// but as we scale to more processors it becomes increasingly dangerous.
/// Second, many processes writing to the same file (or even directory) can cause
/// tremendous file system bottlenecks, possibly bringing down the cluster;
/// ask Keith or Chance in the Baker lab for details.
class AtomTreeDiffJobDistributor : public BaseJobDistributor
{
public:
	typedef BaseJobDistributor parent;

protected:

	typedef utility::vector1< BasicJobOP > JobVector;

public:

	AtomTreeDiffJobDistributor(JobVector jobs, std::string outfile_name);
	virtual ~AtomTreeDiffJobDistributor();

	/// @brief Appends pose to the silent file
	virtual void dump_pose(
		std::string const & tag,
		std::map< std::string, core::Real > const & scores,
		core::pose::Pose const & ref_pose,
		core::pose::Pose const & pose
	);

	/// @brief Sets number of digits used in writing atomtree diff.
	virtual void set_precision(
		int bb_precision,
		int sc_precision,
		int bondlen_precision
	);

	virtual void shutdown();

protected:

	virtual bool is_finished(BasicJobOP const & job, int struct_n );

private:

	utility::io::ozstream out_;
	std::set< std::string > used_tags_;
	core::pose::Pose const * last_ref_pose_;
	int bb_precision_, sc_precision_, bondlen_precision_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Distributor for use with plain old PDB files.
/// Use is strongly discouraged in production environments!
/// @details This class is deliberately designed for each process to write
/// its own PDB files in its own directory; it checks for pre-existing files only
/// for use by stopped and re-started jobs, NOT for coordinating between processes.
/// (To coordinate, it would have to "touch" the non-existant file in next_job()
/// before starting to process it, but see AtomtreeDiffJobDistributor for an
/// explanation of why coordinating processes via the filesystem is a bad idea.)
class PlainPdbJobDistributor : public BaseJobDistributor
{
public:
	typedef BaseJobDistributor parent;

protected:

	typedef utility::vector1< BasicJobOP > JobVector;

public:

	PlainPdbJobDistributor(JobVector jobs, std::string outfile_name="none");
	virtual ~PlainPdbJobDistributor();

	/// @brief Allows setting of inprogress.
	virtual void startup();

	using parent::get_output_filename;

	/// @brief Translates an output tag name to an output PDB file name.
	virtual std::string get_output_filename(std::string const & tag);

	/// @brief Writes pose and basic score data to a standard PDB file.
	virtual void dump_pose_and_map(
		std::string const & tag,
		core::pose::Pose & pose
	);

	/// @brief Opens a temp file (.in_progress)
	virtual void temp_file(std::string const & tag);

	void score_map( std::map < std::string, core::Real > & score_map_in ) { score_map_ = score_map_in; }
protected:

	virtual bool is_finished(BasicJobOP const & job, int struct_n );

	/// @brief Writes score data to PDB file in YAML format.
	virtual void dump_scores(
		utility::io::ozstream & out,
		std::string const & tag,
		core::pose::Pose & pose
	);

private:

	bool scorefile_;
	utility::file::FileName scorefile_name_;
	std::map< std::string, core::Real > score_map_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Distributor for use with raw files.
/// @details This class is deliberately designed for each process to work on
/// its own silent file (and preferrably in its own directory);
/// unlike Rosetta++, we DO NOT share silent files among multiple processes.
/// First, true atomic file locking is impossible on most distributed file systems (e.g. NFS),
/// meaning that files may be corrupted by multiple simultaneous writes.
/// For long-running processes, collisions may be rare,
/// but as we scale to more processors it becomes increasingly dangerous.
/// Second, many processes writing to the same file (or even directory) can cause
/// tremendous file system bottlenecks, possibly bringing down the cluster;
/// ask Keith or Chance in the Baker lab for details.
class PlainRawJobDistributor : public BaseJobDistributor
{
public:
	typedef BaseJobDistributor parent;

protected:

	typedef utility::vector1< BasicJobOP > JobVector;

public:

	PlainRawJobDistributor(JobVector jobs, std::string outfile_name);

	virtual ~PlainRawJobDistributor();

	/// @brief Writes pose and basic score data to a standard silent file.
	virtual void dump_pose_and_map(
		std::string const & tag,
		core::pose::Pose & pose
	);

	using parent::get_output_filename;

	/// @brief Translates an output tag name to an output PDB file name.
	virtual std::string get_output_filename(std::string const & tag);
	virtual std::string get_output_tag( int const & struct_n );

	void score_map( std::map < std::string, core::Real > & score_map_in ) { score_map_ = score_map_in; }
protected:

	virtual bool is_finished(BasicJobOP const & job, int struct_n );

private:

	std::string rawfile_name_;
	utility::vector1< std::string > used_tags_;
	std::map < std::string, core::Real > score_map_;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Distributor for use with silent files.
/// @details This class is deliberately designed for each process to work on
/// its own silent file (and preferrably in its own directory);
/// unlike Rosetta++, we DO NOT share silent files among multiple processes.
/// First, true atomic file locking is impossible on most distributed file systems (e.g. NFS),
/// meaning that files may be corrupted by multiple simultaneous writes.
/// For long-running processes, collisions may be rare,
/// but as we scale to more processors it becomes increasingly dangerous.
/// Second, many processes writing to the same file (or even directory) can cause
/// tremendous file system bottlenecks, possibly bringing down the cluster;
/// ask Keith or Chance in the Baker lab for details.
class PlainSilentFileJobDistributor : public BaseJobDistributor
{
public:
	typedef BaseJobDistributor parent;

protected:
	typedef utility::vector1< BasicJobOP > JobVector;

public:

	PlainSilentFileJobDistributor(JobVector jobs);

	virtual ~PlainSilentFileJobDistributor();

	/// @brief Writes pose and basic score data to a standard silent file.
	virtual void dump_pose(
		BasicJobOP const & job,
		int const & nstruct,
		bool const & fullatom,
		core::pose::Pose & pose
	);

	/// @brief Writes the silent_struct to a silen file
	void dump_silent(
		int const & struct_n,
		core::io::silent::SilentStruct & silent_struct
	);

	void dump_silent(
		core::io::silent::SilentFileData const& silent_file
	);

	virtual std::string get_output_tag( BasicJobOP const & job, int const & struct_n ) const;
	virtual std::string get_current_output_tag();
	virtual std::string get_output_filename();

	virtual void startup();

	virtual void shutdown();

protected:

	virtual bool is_finished(BasicJobOP const & job, int struct_n );

private:

	utility::vector1< std::string > used_tags_;

};


} // namespace jobdist
} // namespace protocols

#endif // INCLUDED_protocols_jobdist_JobDistributors_HH
