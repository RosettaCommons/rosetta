// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobQueen.hh
/// @brief  class declaration for JobQueen
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_JobQueen_hh
#define INCLUDED_protocols_jd3_JobQueen_hh

//unit headers
#include <protocols/jd3/JobQueen.fwd.hh>

//project headers
#include <core/types.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>

#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace jd3 {

/// @brief The %JobQueen class (think of a queen bee) has quite a few responsibilities:
/// - determining what jobs should be run and representing them as larval job objects
/// - taking a larval job object and maturing it into a fully functional Job that can be run
/// - processing the JobResult produced by a Job when it runs
/// - determining the number of rounds of larval-job creation and job execution that should
///   be performed.
///
/// @details The %JobQueen is responsible for creating all of the LarvalJobs, grouped
/// into "batches" which are dispatched by the JobDistributor; the queen then
/// matures the larval jobs into the mature form (a Job, e.g. a Mover/Pose pair) and
/// gives them up to the JobDistributor to actually execute the work.  The
/// JobDistributor then actually calls the Job's "run" method, which produces a JobResult.
/// It then hands this JobResult (e.g. a Pose or several Poses) back to the %JobQueen for
/// output.  In fact, in MPI contexts, the JobDistributor guarantees that it will hand
/// all of the completed JobResults for a single batch to a single %JobQueen, so that
/// she can compute aggregate statistics (e.g. averaging energies) or can look at the
///  results from that round of jobs and decide which to carry forward into the next
/// round of execution.
class JobQueen : public utility::pointer::ReferenceCount
{
public:

	JobQueen();
	virtual ~JobQueen();

	/// @brief Mature the input larval job into a full fledged job that will be run
	/// within this process (i.e. on this CPU), and so can hold pointers to
	/// data that may persist inside the ResourceManager or may even be used by
	/// another process running in this thread.  To pull that off, the Job object
	/// must share no non-bitwise-const data with any other Job object.  That means
	/// if the Job were to contain a Pose and a Mover (as the StandardJob does), then
	/// so each Mover should be freshly constructed, and if the Pose had been copied
	/// from an existing Pose, must use the Pose's "deep_copy" method (since Pose's
	/// sometimes share non-constant data between them, e.g. the AtomTree observer
	/// system, and the *sigh* constraints ).
	virtual JobOP mature_larval_job( LarvalJobCOP job ) = 0;

	/// @biref The JobQueen must be able to determine if a particular job has already
	/// completed (or alternatively, has already been started by another process), and
	/// if so, the JobDistributor will not re-run the job.  The JobQueen returns "true"
	/// if the job has completed (or started elsewhere), and "false" otherwise.
	virtual bool has_job_completed( LarvalJobCOP job ) = 0;

	/// @brief Some (but not all) JobDistributors mark which jobs are in the process of
	/// being run by asking the JobQueen to create a temporary file marking that
	/// fact on the file system.
	virtual void mark_job_as_having_begun( LarvalJobCOP job ) = 0;

	/// @brief The JobDistributor will call this function to inform the JobQueen that
	/// a job has "completed" -- in the sense that it will not be run in the future.
	/// It does not guarantee that the job was run on this CPU or that it successfully
	/// completed anywhere.  The JobStatus indicates whether the job completed or
	/// whether it failed.
	virtual void note_job_completed( LarvalJobCOP job, JobStatus status ) = 0;

	/// @brief The JobDistributor guarnatees that exactly one JobQueen will see every
	/// JobResult generated within a Job batch. This guarantee allows the JobQueen to
	/// aggregate data across all of the Jobs so that Rosetta is able to compute data
	/// from structures instead of forcing that computation into accessory scripts.
	virtual void completed_job_result( LarvalJobCOP job, JobResultOP result ) = 0;

	/// @brief This function determines what jobs exist.  This function neither knows nor
	/// cares what jobs are already complete on disk/memory - it just figures out what
	/// ones should exist given the input.
	virtual LarvalJobs determine_job_list() = 0;

	/// @brief The %JobQueen may indicate to the JobDistributor that multiple rounds of structure
	/// generation are desired by returning "true" to this function call (this function can constitutively
	/// return "false" and the first round will still always be executed).  After returning "true",
	/// the JobDistributor will ask the %JobQueen for another list of jobs through its determine_job_list
	/// method.  After the first round, only the %JobQueens which were given the JobResult data
	/// (one %JobQueen per Job batch) will be asked for a second job list.
	virtual bool more_jobs_remain() = 0;

	/// @brief All JobQueens should, for the sake of documentation, describe their job input
	/// XML format in the form of an XSD (XML Schema Definition). If the JobDistributor is awakened
	/// with the flag "jd3::output_xsd <output file>" on the command line, then if this function
	/// returns a non-empty string, the JobDistributor will write out the XSD to the command line.
	virtual std::string job_definition_xsd() const = 0;

}; // JobQueen



} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_JobQueen_HH
