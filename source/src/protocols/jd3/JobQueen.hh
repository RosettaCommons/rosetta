// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/JobQueen.hh
/// @brief  class declaration for JobQueen
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_JobQueen_hh
#define INCLUDED_protocols_jd3_JobQueen_hh

//unit headers
#include <protocols/jd3/JobQueen.fwd.hh>

//project headers
#include <core/types.hh>
#include <protocols/jd3/CompletedJobOutput.fwd.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/JobDigraph.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>
#include <protocols/jd3/JobSummary.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>
#include <protocols/jd3/deallocation/DeallocationMessage.fwd.hh>
#include <protocols/jd3/output/OutputSpecification.fwd.hh>
#include <protocols/jd3/output/ResultOutputter.fwd.hh>

#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <list>

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
	~JobQueen() override;

	/// @brief All JobQueens must describe their job input XML format in the form of an XSD
	/// (XML Schema Definition), and they must validate their job input files against their XSDs.
	/// If the JobDistributor is awakened with the flag "jd3::output_job_xsd <output file>" on the
	/// command line, then the JobDistributor will write out the XSD to the output file.
	virtual std::string job_definition_xsd() const = 0;

	/// @brief JobQueens may optionally define an XSD for their resource definition file
	/// which is fed to their resource manager (if they control one). If the JobDistributor
	/// is awakened with the flag "jd3::output_resource_xsd <output file>" on the command
	/// line, and if the derived JobQueen defines an XSD for resources she uses, then the
	/// JobDistributor will write out the resource definition XSD to the output file.
	/// The derived queen may return an empty string to indicate that no resources are
	/// definable for the resource manager.
	virtual std::string resource_definition_xsd() const = 0;

	/// @brief The JobQueen creates a directed acyclic graph (DAG) describing the sets of all
	/// jobs and the interdependencies between them to the JobDistributor.  The JobDistributor
	/// will allow the JobQueen to update the DAG by adding new nodes and adding edges to new
	/// nodes if the JobQueen does not know up front how many nodes will be in the graph.
	virtual JobDigraphOP initial_job_dag() = 0;

	/// @brief The JobQueen is allowed to update the JobDigraph over the course of execution
	/// but only in a particular way: by adding new nodes, and then by adding edges that land
	/// on those new nodes.  The JobQueen is not allowed to add edges that land on an old node.
	virtual void update_job_dag( JobDigraphUpdater & updater ) = 0;

	/// @brief The JobDistributor asks the JobQueen for a list of jobs that should be performed.
	/// This function neither knows nor cares what jobs are already complete on disk/memory
	/// - it just figures out what ones should exist given the input.  The JobDistributor
	/// tells the %JobQueen which node in the JobDAG it is requesting jobs for, and puts
	/// a limit on the number of LarvalJobs that the JobQueen should return.
	/// The JobDistributor will call this function repeatedly for a single node until the
	/// %JobQueen returns an empty list, at which point, the JobDistributor will consider
	/// the node's jobs exhausted.
	virtual LarvalJobs determine_job_list( Size job_dag_node_index, Size max_njobs ) = 0;

	/// @biref The JobQueen must be able to determine if a particular job has already
	/// completed (or alternatively, has already been started by another process), and
	/// if so, the JobDistributor will not re-run the job.  The JobQueen returns "true"
	/// if the job has completed (or started elsewhere), and "false" otherwise.
	virtual bool has_job_completed( LarvalJobCOP job ) = 0;

	/// @brief Some (but not all) JobDistributors mark which jobs are in the process of
	/// being run by asking the JobQueen to create a temporary file marking that
	/// fact on the file system.
	virtual void mark_job_as_having_begun( LarvalJobCOP job ) = 0;

	/// @brief Mature the input larval job into a full fledged job that will be run
	/// within this process (i.e. on this CPU), and thus can hold pointers to
	/// data that may persist inside the ResourceManager or may even be used by
	/// another process running in this thread. To pull that off, the Job object
	/// must share no non-bitwise-const data with any other Job object.  That means
	/// if the Job were to contain a Pose and a Mover (as the MoverAndPoseJob does), then
	/// so each Mover should be freshly constructed, and if the Pose had been copied
	/// from an existing Pose, must use the Pose's "deep_copy" method (since Pose's
	/// sometimes share non-constant data between them, e.g. the AtomTree observer
	/// system, and sometimes Constraints ). The JobResults vector, which is supplied
	/// by the Jobdistributor, has entries corresponding to the JobResults from
	/// JobIDs in its input_job_result_indices vector.
	virtual JobOP mature_larval_job( LarvalJobCOP job, utility::vector1< JobResultCOP > const & input_job_results ) = 0;

	/// @brief There are two interfaces to note_job_completed: one in which the
	/// job index alone is passed in, a second in which the entire LarvalJob is
	/// provided.
	virtual bool larval_job_needed_for_note_job_completed() const = 0;

	/// @brief The JobDistributor will call this function to inform the JobQueen that
	/// a job has "completed" -- in the sense that it will not be run in the future.
	/// It does not guarantee that the job was run on this CPU or that it successfully
	/// completed anywhere.  The JobStatus indicates whether the job completed or
	/// whether it failed. The nresults count lists how many JobSummary/JobResult pairs
	/// were produced by the completed job.
	virtual void note_job_completed( LarvalJobCOP job, JobStatus status, Size nresults ) = 0;

	/// @brief The JobDistributor will call this function to inform the JobQueen that
	/// a job has "completed" -- in the sense that it will not be run in the future.
	/// It does not guarantee that the job was run on this CPU or that it successfully
	/// completed anywhere.  The JobStatus indicates whether the job completed or
	/// whether it failed. The nresults count lists how many JobSummary/JobResult pairs
	/// were produced by the completed job.
	virtual void note_job_completed( core::Size job_id, JobStatus status, Size nresults ) = 0;

	/// @brief There are two interfaces to completed_job_summary: one in which the
	/// job index alone is passed in, a second in which the entire LarvalJob is
	/// provided.
	virtual bool larval_job_needed_for_completed_job_summary() const = 0;

	/// @brief The JobDistributor guarnatees that exactly one JobQueen will see every
	/// JobSummary generated within a Job batch. This guarantee allows the JobQueen to
	/// aggregate data across all of the Jobs so that Rosetta is able to compute data
	/// from structures instead of forcing that computation into accessory scripts.
	virtual void completed_job_summary( LarvalJobCOP job, Size result_index, JobSummaryOP summary ) = 0;

	/// @brief The JobDistributor guarnatees that exactly one JobQueen will see every
	/// JobSummary generated within a Job batch. This guarantee allows the JobQueen to
	/// aggregate data across all of the Jobs so that Rosetta is able to compute data
	/// from structures instead of forcing that computation into accessory scripts.
	virtual void completed_job_summary( core::Size job_id, Size result_index, JobSummaryOP summary ) = 0;

	/// @brief The JobDistributor asks the JobQueen which JobResults should be queued for output.
	/// The JobQueen should reply by returning a list of OutputSpecification objects. Each
	/// OutputSpecification indicates a particular JobResult using the JobResultID object that the
	/// Specification holds. It also indicates how the JobResult should be output by whatever object
	/// will be doing to job outputting. The JobQueen should ask for a given JobResult to be output
	/// only once. After a job result is output, the JobDistributor will discard
	/// the JobResult, so the JobQueen should not tell the JobDistributor to output a job if
	/// it will be used as an input to another job in the future. Note that this method will only
	/// be called on the master-node JobQueen (aka JQ0), but that she will not necessarily be
	/// the JobQueen to actually return the ResultOutputter for the jobs that will be output.
	/// The JobDistributor holds the option to distribute output to other nodes
	virtual std::list< output::OutputSpecificationOP >
	jobs_that_should_be_output() = 0;

	/// @brief The JobDistributor, to manage memory use, asks the JobQueen which JobResults may be
	/// discarded because they will not be used in the future.  The JobDistributor will exit with
	/// an error message if the %JobQueen gives it a LarvalJob that lists one of these discarded
	/// JobResults as a required input for that LarvalJob.
	virtual std::list< JobResultID > job_results_that_should_be_discarded() = 0;

	/// @brief The JobDistributor will ask the JobQueen to provide an outputter for a result
	/// giving it an annotated-with-job-distributor-specific-prefix OutputSpecification. The
	/// JobQueen must guarantee that a distinct ResultOutputter is given for each distinct
	/// JD prefix and that no two ResultOutputters share non-(bitwise-)constant data as the
	/// JobDistributor may use multiple threads to perform output. Note that the JobQueen
	/// that creates the outputter may not have been the JobQueen to request that output be performed.
	/// This function is thus similar to mature_larval_job: it requires that the JobQueen
	/// make few assumptions about the state of her data in order for this behavior to translate
	/// from the VanillaJobDistributor to the MPIWorkPoolJobDistributor.
	virtual output::ResultOutputterOP
	result_outputter( output::OutputSpecification const & output_specification ) = 0;

	/// @brief Send all buffered output to disk -- called by the JobDistributor right before it shuts down
	/// if it hits an error or catches an exception that it cannot ignore.
	virtual void flush() = 0;

	/// @brief Are there any deallocation messages that the JobQueen would like to have delivered
	/// to the JobQueens living on remote hosts? If there are none, the derived JobQueen may return
	/// an empty list. The JobDistributor guarantees to send the deallocation message to all remote
	/// hosts. The JobDistributor will ask the JobQueen on the head node for a set of deallocation
	/// messages after each set of jobs from determine_job_list have been distributed. (i.e., before
	/// the next call to determine_job_list is made).
	virtual
	std::list< deallocation::DeallocationMessageOP >
	deallocation_messages() = 0;

	/// @brief A deallocation message first sent to the JobDistributor on this host originating from
	/// a remote JobQueen
	virtual
	void
	process_deallocation_message( deallocation::DeallocationMessageOP message ) = 0;

}; // JobQueen



} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_JobQueen_HH
