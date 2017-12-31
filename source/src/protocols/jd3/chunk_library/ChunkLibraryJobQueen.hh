// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/chunk_library/ChunkLibraryJobQueen.hh
/// @brief  class declaration for ChunkLibraryJobQueen
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_chunk_library_ChunkLibraryJobQueen_hh
#define INCLUDED_protocols_jd3_chunk_library_ChunkLibraryJobQueen_hh

// unit headers
#include <protocols/jd3/chunk_library/ChunkLibraryJobQueen.fwd.hh>

// package headers
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/JobDigraph.fwd.hh>
#include <protocols/jd3/JobOutputIndex.fwd.hh>
#include <protocols/jd3/chunk_library_inputters/ChunkLibraryInputSource.fwd.hh>
#include <protocols/jd3/chunk_library_inputters/ChunkLibraryInputter.fwd.hh>
#include <protocols/jd3/chunk_library_inputters/ChunkLibraryInputterCreator.fwd.hh>
#include <protocols/jd3/pose_outputters/PoseOutputSpecification.fwd.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.fwd.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterCreator.fwd.hh>
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.fwd.hh>
#include <protocols/jd3/chunk_library/ChunkLibraryInnerLarvalJob.fwd.hh>

// project headers
#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/import_pose_options.fwd.hh>
#include <basic/resource_manager/JobOptions.fwd.hh>

//utility headers
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>
#include <utility/options/keys/all.fwd.hh>

// numeric headers
#include <numeric/DiscreteIntervalEncodingTree.hh>

//c++ headers
#include <string>
#include <map>

namespace protocols {
namespace jd3 {
namespace chunk_library {

class PreliminaryLarvalJob
{
public:
	PreliminaryLarvalJob();
	~PreliminaryLarvalJob();
	PreliminaryLarvalJob( PreliminaryLarvalJob const & src );
	PreliminaryLarvalJob & operator = ( PreliminaryLarvalJob const & rhs );

	ChunkLibraryInnerLarvalJobOP inner_job;
	utility::tag::TagCOP job_tag;
	utility::options::OptionCollectionCOP job_options;
	chunk_library_inputters::ChunkLibraryInputterOP chunk_library_inputter;
};

typedef std::list< PreliminaryLarvalJob > PreliminaryLarvalJobs;

/// @brief The %ChunkLibraryJobQueen is meant to handle the most common form of Rosetta jobs where
/// a protocol is applied to a single input structure to generate a single output structure.
/// Most JobQueens in Rosetta will to derive from this JobQueen.  To make the process of deriving
/// new JobQueens easy, this class's API has been carefully designed to be easy to work with and
/// the class itself has been designed to do as much heavy lifting as possible and to leave only
/// the responsibilities that the %ChunkLibraryJobQueen could not perform to the derived classes.
class ChunkLibraryJobQueen : public JobQueen
{
public:
	/// @brief The ChunkLibraryJobQueen constructor asks the ChunkLibraryInputterFactory for a ChunkLibraryInputter
	/// and creates a ResourceManager
	ChunkLibraryJobQueen();

	~ChunkLibraryJobQueen() override;


	/// @brief The %ChunkLibraryJobQueen assembles the XSD from virtual functions she invokes on
	/// the derived %JobQueen: append_job_tag_subelements, append_common_tag_subelements, and
	/// add_option/add_options.
	std::string job_definition_xsd() const override;
	std::string resource_definition_xsd() const override;

	/// @brief The %ChunkLibraryJobQueen provides an implementation of this function which returns
	/// the most straight-forward DAG representing a set of jobs that have no interdependencies:
	/// a DAG with a single node and no edges.
	JobDigraphOP initial_job_dag() override;

	/// @brief The %ChunkLibraryJobQueen's implementation is to not update the JobDAG at all: the
	/// most basic protocol defines a job DAG with only a single node.
	void update_job_dag( JobDigraphUpdater & updater ) override;

	/// @brief The %ChunkLibraryJobQueen manages the process of creating the list of LarvalJobs that
	/// will be later matured into actual jobs and run by the JobDistributor.  Classes derived
	/// from the %ChunkLibraryJobQueen need answer a few virtual functions that the %ChunkLibraryJobQueen
	/// will invoke during this process.
	LarvalJobs determine_job_list( Size job_dag_node_index, Size max_njobs ) override;

	bool has_job_completed( protocols::jd3::LarvalJobCOP job ) override;

	void mark_job_as_having_begun( protocols::jd3::LarvalJobCOP job ) override;

	protocols::jd3::JobOP
	mature_larval_job(
		protocols::jd3::LarvalJobCOP job,
		utility::vector1< JobResultCOP > const & input_job_results
	) override;

	bool larval_job_needed_for_note_job_completed() const override;
	void note_job_completed( LarvalJobCOP job, JobStatus status, Size nresults ) override;
	void note_job_completed( core::Size job_id, JobStatus status, Size nresults ) override;

	/// @brief The SJQ needs the larval job in order to create the OutputSpecification for jobs
	/// that should be output. Derived JobQueens should only override this method if they are not
	/// using the PoseOutputters that the SJQ is using or if are themselves tracking the LarvalJobOPs
	/// for currently-running jobs (a potentially memory-hogging operation, depending on the number
	/// of jobs run at once?)
	bool larval_job_needed_for_completed_job_summary() const override;

	/// @brief If the job result that is summarized here should be output, then the SJQ needs to
	/// have access to the LarvalJob in order to do that.
	void completed_job_summary( LarvalJobCOP job, Size result_index, JobSummaryOP summary ) override;

	/// @brief This method should not be called, as the SJQ needs the LarvalJob pointer when creating
	/// OutputSpecifications.
	void completed_job_summary( core::Size job_id, Size result_index, JobSummaryOP summary ) override;

	std::list< output::OutputSpecificationOP > jobs_that_should_be_output() override;
	std::list< JobResultID > job_results_that_should_be_discarded() override;

	/// @brief Return the bag of of PoseOutputters (in the form of a MultipleOutputter) for the
	/// Pose that has been requested and specified by a particular OutputSpecification; this
	/// function guarantees that for each individual outputter-name (respecting the JD-provided
	/// filename suffix) that a separate set of PoseOutputters are returned so that multiple threads
	/// can potentially output at the same time.
	output::ResultOutputterOP result_outputter( output::OutputSpecification const & spec ) override;

	/// @brief Read from an input string representing the contents of the job-definiton XML file
	/// and construct a set of LarvalJobs; this function is primarily useful for testing,
	/// but could be used to organize jobs by an enterprising job distributor or by another JobQueen.
	void determine_preliminary_job_list_from_xml_file( std::string const & job_def_string );

	void flush() override;

	std::list< deallocation::DeallocationMessageOP >
	deallocation_messages() override;

	/// @brief A deallocation message first sent to the JobDistributor on this host originating from
	/// a remote JobQueen
	void
	process_deallocation_message( deallocation::DeallocationMessageOP message ) override;

	/// @brief Naming function for the complexTypes in the job-definition XSD
	static
	std::string
	job_def_complex_type_name( std::string const & type );


protected:

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////// FUNCTIONS FOR OVERRIDING THE SET OF INPUTTERS AND OUTPUTTERS //////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// @brief The derived job queen (DJQ) may tell the %ChunkLibraryJobQueen to only accept a subset of
	/// ChunkLibraryInputters. First the DJQ says "do_not_accept_all_chunk_library_inputters_from_factory" and then she
	/// can call "allow_chunk_library_inputter." This function must only be called the DJQ's constructor.
	void
	do_not_accept_all_chunk_library_inputters_from_factory();

	/// @brief If the derived job queen (DJQ) wants to allow a ChunkLibraryInputter to be used either because
	/// it is not registered with the ChunkLibraryInputterFactory, or because she wants to allow only a
	/// subset of the inputters that are registered with the Factory, then she may indicate that
	/// to the %ChunkLibraryJobQueen with this function. If the derived job queen does not want to
	/// use all of the pose inputters provided by the factory, then she must call
	/// "do_not_accept_all_chunk_library_inputters_from_factory()" before she calls this function.
	/// This function must only be called in the DJQ's constructor. All DJQs must invoke this
	/// method in their constructors and not simply the DJQ that lives on the "head node"
	/// in order for distributed output to workk correctly.
	void
	allow_chunk_library_inputter( chunk_library_inputters::ChunkLibraryInputterCreatorOP creator );

	/// @brief The derived job queen (DJQ) may tell the %ChunkLibraryJobQueen to only accept a subset of
	/// PoseOutputters. First the DJQ says "do_not_accept_all_pose_outputters_from_factory" and then
	/// she can call "allow_pose_outputter." This function must only be called the DJQ's constructor.
	void
	do_not_accept_all_pose_outputters_from_factory();

	/// @brief If the derived job queen (DJQ) wants to allow a PoseOutputter to be used either because
	/// it is not registered with the PoseOutputterFactory, or because she wants to allow only a
	/// subset of the outputters that are registered with the Factory, then she may indicate that
	/// to the %ChunkLibraryJobQueen with this function. If the derived job queen does not want to
	/// use all of the pose outputters provided by the factory, then she must call
	/// "do_not_accept_all_pose_outputters_from_factory()" before she calls this function.
	/// If the DJQ has called the "do_not_accept_all_pose_outputters_from_factory" function, then
	/// she specifies the default outputter on the first time she calls "allow_pose_outputter."
	/// This function must only be called in the DJQ's constructor.
	void
	allow_pose_outputter( pose_outputters::PoseOutputterCreatorOP creator );

	/// @brief If the derived job queen wants the user to get a particular PoseOutputter by default
	/// (i.e. when the user hasn't specified which outputter to use in a job definition file and
	/// in the absence of a command-line flag that says which outputter to use) instead of the
	/// default PoseOutputter chosen by the PoseOutputterFactory, then the DJQ should call this function.
	/// If the DJQ has previously called "do_not_accept_all_pose_outputters_from_factory", then
	/// the first call to "allow_pose_outputter" specified the default outputter, and a separate call
	/// to this function is not absolutely necessary. This function must only be called in the
	/// DJQ's constructor.
	void
	set_default_outputter( pose_outputters::PoseOutputterCreatorOP creator );

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

protected:
	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////// FUNCTIONS FOR CONTROLLING THE STRUCTURE OF THE XML SCHEMA /////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// @brief The derived JobQueen must inform the %ChunkLibraryJobQueen of any additional tags that
	/// belong as sub elements of the <Job> tag. The %ChunkLibraryJobQueen will have already appended
	/// the <Input>, <Output>, and <Options> subtags before this function executes, so the only
	/// functions that the derived class should invoke of the XMLSchemaComplexTypeGenerator are
	/// its add_ordered_subelement_set_* functions, if you have additional subelements to add.
	/// Adding additional subelements to the complexType generator may require writing additional
	/// complexTypes to the XML Schema for the job definition file, and so that job definition xsd
	/// is passed in as well.
	virtual
	void append_job_tag_subelements(
		utility::tag::XMLSchemaDefinition & job_definition_xsd,
		utility::tag::XMLSchemaComplexTypeGenerator & ct_gen
	) const;

	/// @brief The derived JobQueen must inform the %ChunkLibraryJobQueen of any additional tags that
	/// belong as sub elements of the <Common> tag. The %ChunkLibraryJobQueen will have already appended
	/// the <Options> subtag before this function executes, so the only functions that the derived
	/// class should invoke of the XMLSchemaComplexTypeGenerator are its
	/// add_ordered_subelement_set_* functions, if you have additional subelements to add.
	/// Adding additional subelements to the complexType generator may require writing additional
	/// complexTypes to the XML Schema for the job definition file, and so that job definition xsd
	/// is passed in as well. The Tags in the <Common> block should be parsed for each Job in the
	/// mature_larval_job step; these tags are available to the derived class through the
	/// common_block_tags method. This is where data that is not constant over the course of the job
	/// can be declared (e.g. ScoreFunctions, the weights of which are often modified during
	/// the course of protocol execution), but data that is constant and that can be shared between
	/// multiple jobs ought to be initialized and held by the ResourceManager. Options provided in
	/// the <Common> block will supercede the command line, but will be overriden by any options
	/// provided for a particular job. (This is logic the %ChunkLibraryJobQueen manages).
	virtual
	void
	append_common_tag_subelements(
		utility::tag::XMLSchemaDefinition & xsd,
		utility::tag::XMLSchemaComplexTypeGenerator & ct_gen
	) const;

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

protected:

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////// FUNCTIONS FOR CONTROLLING THE STRUCTURE OF THE JOB DAG //////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Allow the derived JobQueen to tell the %ChunkLibraryJobQueen to initialize the preliminary
	/// job list; this is perhaps necessary in the context of multi-round protocols when job-definition
	/// file specifies the JobDAG.
	virtual
	void
	determine_preliminary_job_list();

	/// @brief Allow the derived job queen the opportunity to update the ChunkLibraryJobQueen's
	/// job_graph_ data member by adding nodes as well as edges that land on the new nodes.
	JobDigraphUpdater
	updater_for_sjq_job_graph();

	/// @brief Read access for the subset of nodes in the job DAG which the %ChunkLibraryJobQueen
	/// is responsible for producing the larval_jobs. They are called "preliminary" jobs because
	/// they do not depend on outputs from any previous node in the graph. (The set of job nodes
	/// that contain no incoming edges, though, could perhaps be different from the set of
	/// preliminary job nodes, so the %ChunkLibraryJobQueen requires that the deried job queen
	/// inform her of which nodes are the preliminary job nodes.)
	utility::vector1< core::Size > const &
	preliminary_job_nodes() const;

	/// @brief Returns true if all of the jobs for the given input preliminary-job node have
	/// been created and returned in a call to determine_job_list. Even if all jobs have been
	/// created, not all jobs have necessarily been completed.
	bool
	all_jobs_assigned_for_preliminary_job_node( core::Size node_id ) const;


	/// @brief The index of the first job for a particular preliminary-job node; this function
	/// returns zero if no jobs for this node have yet been created
	core::Size preliminary_job_node_begin_job_index( core::Size node_id ) const;

	/// @brief The index of the last job for a particular preliminary-job node; this function
	/// returns zero if there are some jobs for this node that have not yet been created.
	core::Size preliminary_job_node_end_job_index( core::Size node_id ) const;

	/// @brief Read access to derived JobQueens to the preliminary job list.
	/// This will return an empty list if  determine_preliminary_jobs has not yet
	/// been called.
	utility::vector1< PreliminaryLarvalJob > const &
	preliminary_larval_jobs() const;

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Read access for which jobs have completed and how; if a job-id is a member
	/// of this DIET, then it has completed (either in success or failure).
	numeric::DiscreteIntervalEncodingTree< core::Size > const & completed_jobs() const;

	/// @brief Read access for which jobs have completed and how; if a job-id is a member
	/// of this DIET, then it completed successfully.
	numeric::DiscreteIntervalEncodingTree< core::Size > const & successful_jobs() const;

	/// @brief Read access for which jobs have completed and how; if a job-id is a member
	/// of this DIET, then it completed unsuccessfully.
	numeric::DiscreteIntervalEncodingTree< core::Size > const & failed_jobs() const;

	/// @brief Read access for which jobs have completed and how; if a job-id is a member
	/// of this DIET, then the job has had all of its results output or discarded.
	numeric::DiscreteIntervalEncodingTree< core::Size > const & processed_jobs() const;

	/// @brief Ask the derived JobQueen to expand / refine a preliminary larval job, by
	/// possibly reading per-job data out of the Tag associated with each job. If there is
	/// nothing that needs to be done by the derived class, it may elect to use the base-class
	/// implementation of this function which simply returns a list of the inner-job pointers
	/// in the input list.
	/// Why all the back and forth?  So that the derived job queens can take one input structure
	/// and multiply it into several sepearate jobs; e.g. in a ddG protocol, you need to run a
	/// separate set of jobs for the WT sequence and also for the mutant sequence. The preliminary
	/// job will reflect the structure of the WT, and then that will get expanded out.
	/// This funciton is non-const so that the derived class can even track which WT sequences it
	/// has seen so that after the first WT sequence for a particular input structure, it only expands
	/// out the mutant sequences if many mutants on that input.
	virtual ChunkLibraryInnerLarvalJobs refine_preliminary_job( PreliminaryLarvalJob const & prelim_job );

	/// @brief Expand a ChunkLibraryInnerLarvalJob into a full set of LarvalJobs, creating nstruct LarvalJob objects
	/// The base class implementation of this function invokes the create_larval_job factory method so
	/// that dervied JobQueens can ensure the instantiation of derived LarvalJobs.
	virtual LarvalJobs expand_job_list( ChunkLibraryInnerLarvalJobOP inner_job, core::Size max_larval_jobs_to_create );

	/// @brief Factory method for derived classes if they wish to rely on classes derived
	/// from ChunkLibraryInnerLarvalJob.  This is invoked by the ChunkLibraryJobQueen in her determine_job_list
	/// method just as jobs are prepared. If the base ChunkLibraryInnerLarvalJob class is desired, then do
	/// not override this method.
	virtual ChunkLibraryInnerLarvalJobOP create_inner_larval_job( core::Size nstruct, core::Size prelim_job_node ) const;

	/// @brief Factory method for derived classes if they wish to rely on classes derived
	/// from LarvalJob.  This is invoked by the ChunkLibraryJobQueen in the expand_job_list method.
	/// If the base LarvalJob class is desired, then do not override this method.
	virtual LarvalJobOP create_larval_job( ChunkLibraryInnerLarvalJobOP job, core::Size nstruct_index, core::Size larval_job_index );

	/// @brief Factory method for derived classes if they wish to rely on a class besides the
	/// MoverAndPoseJob, which is returned by the %ChunkLibraryJobQueen's implementation of this
	/// function.
	virtual JobOP create_job( LarvalJobCOP job ) const;

	/// @brief The derived JobQueen must define the method that takes a larval job and the job-specific options
	/// and matures the larval job into a full job.
	virtual
	JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options,
		utility::vector1< JobResultCOP > const & input_job_results
	) = 0;

	/// @brief The ChunkLibraryJobQueen cannot readily manage the complexity of organizing
	/// larval jobs for JobDAG nodes beyond the preliminary nodes. All job-creation logic beyond
	/// the first set of nodes is the responsibility of the derived JobQueen. The %ChunkLibraryJobQueen
	/// provides a noop implementation of this function.
	///
	/// @details For the sake of outputting structures with the properties you want at the
	/// completion of a multi-stage trajectory, derived job queens with a tree-like job-progression
	/// (i.e. only one parent job for each child job) should copy the jobdef_tag into
	/// new InnerLarvalJobs from the parent job. The ChunkLibraryJobQueen does not track parentage
	/// so it cannot know that job 873 originated from job 4 when trying to figure out how to output
	/// the structure for job 873. The derived job queen should also ensure that the "outputter"
	/// data from the parent job is copied into child job.
	virtual
	LarvalJobs
	next_batch_of_larval_jobs_for_job_node( core::Size job_dag_node_index, core::Size max_njobs );

	/// @brief Create an output specification for a given job + result index and store it in the
	/// recent_successes_ list, which will be given to the JobDistributor upon request. This
	/// function invokes the factory method create_output_specification_for_job_result which
	/// derive job queens may override. This version of the function will create an OptionCollection
	/// and then invoke the version of the function that expects one.
	void
	create_and_store_output_specification_for_job_result(
		LarvalJobCOP job,
		core::Size result_index,
		core::Size nresults
	);

	/// @brief Create an output specification for a given job + result index and store it in the
	/// recent_successes_ list, which will be given to the JobDistributor upon request. This
	/// function invokes the factory method create_output_specification_for_job_result which
	/// derive job queens may override. This version of the function accepts an OptionCollection
	/// if there are multiple results for a particular job, and you can reuse the OptionCollection
	/// between calls.
	void
	create_and_store_output_specification_for_job_result(
		LarvalJobCOP job,
		utility::options::OptionCollection const & job_options,
		core::Size result_index,
		core::Size nresults
	);

	/// @brief Construct the OutputSpecification for a completed successful job. Factory method that
	/// may be overriden by the derived job queen. Invoked by
	/// create_and_store_output_specification_for_job_result which is invoked via note_job_completed
	/// by default, and may be invoked by the derived job queen as needed.
	virtual
	output::OutputSpecificationOP
	create_output_specification_for_job_result(
		LarvalJobCOP job,
		utility::options::OptionCollection const & job_options,
		core::Size result_index,
		core::Size nresults
	);

	/// @brief Construct a JobOutputIndex for a given job based on "the obvious", but giving derived
	/// classes the chance to assign their own indices via the assign_output_index function
	JobOutputIndex
	build_output_index(
		protocols::jd3::LarvalJobCOP larval_job,
		Size result_index_for_job,
		Size n_results_for_job
	);

	/// @brief The derived job queen may assign her own numbering to output Poses if she chooses.
	///
	/// @details Just before a job result is written to disk, the %ChunkLibraryJobQueen asks the derived
	/// job queen what indices should be used to identify the soon-to-be-output Pose. The default
	/// behavior is to use the nstruct index for the primary output index, and if there are multiple
	/// result Poses from a job, to use the result index (as is) for the secondary output index.
	/// These values will already have been loaded into the output_index before this function is
	/// called.
	/// This function will possibly be invoked more than a single time for a single job result, so
	/// it is important that the derived JobQueen not assume that it only happens once.
	/// This function will only be called on JQ0; it will not be called on any other JQs, so it
	/// is welcome to use information that would only be known on the "head node"
	virtual
	void
	assign_output_index(
		protocols::jd3::LarvalJobCOP larval_job,
		Size result_index_for_job,
		Size n_results_for_job,
		JobOutputIndex & output_index
	);

	/// @brief If the derived class has deallocation messages that it needs to broadcast to
	/// remote nodes and then to process on those remote nodes, then the %StanardJobQueen will
	/// pass them along to the derived job queen using this function.
	virtual
	void
	derived_process_deallocation_message(
		deallocation::DeallocationMessageOP message
	);

	/////////////////////////////////////////////////////////////////////////////////
	// The following functions are to be used by derived JobQueens to signal to the
	// ChunkLibraryJobQueen that individual options can be specified on a per-job basis
	// in the XML jobs file or that can be specified on the command line. Derived
	// JobQueens are discouraged from reading from the command line directly, since
	// that circumvents the auto-documenting features of the XSD in terms of
	// communicating to the user how Jobs are to be specified. These functions should
	// be invoked during the derived class's constructor.

	void add_options( utility::options::OptionKeyList const & opts );
	void add_option( utility::options::OptionKey const & key );

	///////////////////////////////////////////////////////////////////////////////////
	// The following functions should be used by the derived JobQueens to communicate
	// to the ChunkLibraryJobQueen the structure of the XML document that they support.
	// The ChunkLibraryJobQueen will prepare the output XSD document based on the function
	// calls here, and will also validate the input job-definition file, ensuring that
	// no options are provided to through that file that are not read.

	/// @brief Override the basic specification of input elements as either PDBs, Silent Structures,
	/// or ResourceManager defined structures. This should be called during the derived class's
	/// constructor if it is to be used.
	void remove_default_input_element();

	/// @brief Return the set of subtags that are common to the whole set of jobs in the
	/// Job definition file, if any are given.  This set of tags is read from disk at most
	/// once per execution.
	utility::tag::TagCOP common_block_tags() const;

	/// @brief Return a copy of the Pose to be used with the given job
	core::pose::PoseOP pose_for_job( LarvalJobCOP job, utility::options::OptionCollection const & options );

	// ResourceManagerOP resource_manager();

	// Of course the job inputter might vary from job to job!
	chunk_library_inputters::ChunkLibraryInputterOP
	chunk_library_inputter_for_job( ChunkLibraryInnerLarvalJob const & inner_job ) const;

	// Of course the job outputter might vary from job to job!
	pose_outputters::PoseOutputterOP
	pose_outputter_for_job( ChunkLibraryInnerLarvalJob const & innerJob );

	pose_outputters::PoseOutputterOP
	pose_outputter_for_job( ChunkLibraryInnerLarvalJob const & innerJob, utility::options::OptionCollection const & job_options );

	/// @brief Find the PoseOutputter for a job specification, respecting the Job-distributor
	/// provided output file suffix. The %ChunkLibraryJobQueen guarantees that a different Outputter
	/// is used for each unique name as long as that name is not the empty string. This allows
	/// multiple threads to write to separate files concurrently.
	pose_outputters::PoseOutputterOP
	pose_outputter_for_job( pose_outputters::PoseOutputSpecification const & spec );

	utility::vector1< pose_outputters::SecondaryPoseOutputterOP >
	secondary_outputters_for_job( ChunkLibraryInnerLarvalJob const & innerJob, utility::options::OptionCollection const & job_options );

	core::Size
	nstruct_for_job( utility::tag::TagCOP job_tag ) const;

	utility::options::OptionCollectionOP
	options_for_job( ChunkLibraryInnerLarvalJob const & inner_job ) const;

	utility::options::OptionCollectionOP
	options_from_tag( utility::tag::TagCOP job_options_tags ) const;

private:

	/// @brief After generating the job-definition XSD, construct the preliminary job
	/// list. This is invoked both from determine_job_list_from_xml_file and
	/// determine_job_list -- the latter always constructs an XSD to ensure
	/// that the derived JobQueen has properly constructed an XSD, even if
	/// a job definition file has not been provided on the command line.
	void
	determine_preliminary_job_list_from_xml_file(
		std::string const & job_def_string,
		std::string const & job_def_schema
	);

	void
	load_job_definition_file(
		std::string const & job_def_string,
		std::string const & job_def_schema
	);

	/// @brief Instead of reading a JobDefinition file, construct the set of PreliminaryLarvalJobs
	/// reading from the command line. Invoked by determine_preliminary_job_list.
	void
	determine_preliminary_job_list_from_command_line();

	LarvalJobs next_batch_of_larval_jobs_from_prelim( core::Size job_node_index, core::Size max_njobs );

	pose_outputters::SecondaryPoseOutputterOP
	secondary_outputter_for_job(
		ChunkLibraryInnerLarvalJob const & inner_job,
		utility::options::OptionCollection const & job_options,
		std::string const & secondary_outputter_type
	);

	pose_outputters::SecondaryPoseOutputterOP
	secondary_outputter_for_job(
		pose_outputters::PoseOutputSpecification const &
	);

	/// @brief The SJQ will keep track of all jobs in the processed_jobs_ diet, but requires
	/// that derived classes who are not calling the SJQ's version of completed_job_result call
	/// this function.
	void note_job_result_output_or_discarded( LarvalJobCOP job, Size result_index );

	/// @brief The SJQ will keep track of all discarded jobs in the discarded_jobs_ diet
	void note_job_result_discarded( LarvalJobCOP job, Size result_index );

	utility::options::OptionKeyList concatenate_all_options() const;

	pose_outputters::PoseOutputterOP
	get_outputter_from_job_tag( utility::tag::TagCOP tag ) const;

private:

	// ResourceManagerOP resource_manager_;

	// The set of options that the %ChunkLibraryJobQueen reads from the command line
	// and/or from the job definition XML file.
	utility::options::OptionKeyList options_;
	utility::options::OptionKeyList inputter_options_;
	utility::options::OptionKeyList outputter_options_;
	utility::options::OptionKeyList secondary_outputter_options_;

	bool use_factory_provided_chunk_library_inputters_;
	std::map< std::string, chunk_library_inputters::ChunkLibraryInputterCreatorCOP > inputter_creators_;
	std::list< chunk_library_inputters::ChunkLibraryInputterCreatorCOP > inputter_creator_list_;

	bool use_factory_provided_pose_outputters_;
	std::map< std::string, pose_outputters::PoseOutputterCreatorCOP > outputter_creators_;
	std::list< pose_outputters::PoseOutputterCreatorCOP > outputter_creator_list_;
	bool override_default_outputter_;
	pose_outputters::PoseOutputterCreatorOP default_outputter_creator_;

	// Often, you want to use the same pose outputter for multiple jobs.
	std::map< std::string, pose_outputters::PoseOutputterOP > pose_outputters_;

	bool required_initialization_performed_;
	utility::tag::TagCOP job_definition_file_tags_;
	utility::tag::TagCOP common_block_tags_;

	JobDigraphOP job_graph_;

	// The index of the last larval job that we have created. Incremented within
	// expand_job_list()
	Size larval_job_counter_;

	// For the first node in the JobDAG, the %ChunkLibraryJobQueen will spool out LarvalJobs
	// slowly to the JobDistributor (in increments of the max_njobs parameter in the call
	// to job_dag_node_index). Since max_njobs may be smaller than the nstruct parameter,
	// the %ChunkLibraryJobQueen will need to be able to interrupt the spooling of jobs until
	// the JobDistributor is ready for them. For this reason, it keeps what is effectively
	// a set of indices into a while loop for LarvalJob construction.
	bool preliminary_larval_jobs_determined_;
	utility::vector1< PreliminaryLarvalJob > preliminary_larval_jobs_;
	ChunkLibraryInnerLarvalJobs inner_larval_jobs_for_curr_prelim_job_;
	Size curr_inner_larval_job_index_;
	Size njobs_made_for_curr_inner_larval_job_;
	utility::vector1< core::Size > preliminary_job_node_inds_;
	utility::vector1< core::Size > pjn_job_ind_begin_;
	utility::vector1< core::Size > pjn_job_ind_end_;
	utility::vector1< char > preliminary_job_nodes_complete_; // complete == all jobs assigned

	typedef numeric::DiscreteIntervalEncodingTree< core::Size > SizeDIET;

	SizeDIET completed_jobs_;
	SizeDIET successful_jobs_;
	SizeDIET failed_jobs_;
	SizeDIET processed_jobs_; // jobs where all nresults have been output or discarded

	struct PartialOutputStatus {
		core::Size n_results;
		core::Size n_output_or_discarded;
		SizeDIET results_output_or_discarded;
	};

	std::map< core::Size, PartialOutputStatus > results_processed_for_job_;

	// The jobs that have completed, but have not yet been output, and the instructions for how
	// to output them.
	std::list< output::OutputSpecificationOP > recent_successes_;

	// A mapping from the outputter-type to a representative PoseOutputter/SecondaryPoseOutputter.
	typedef std::map< std::string, pose_outputters::PoseOutputterOP > RepresentativeOutputterMap;
	typedef std::map< std::string, pose_outputters::SecondaryPoseOutputterOP > SecondaryRepresentativeOutputterMap;
	RepresentativeOutputterMap representative_pose_outputter_map_;
	SecondaryRepresentativeOutputterMap representative_secondary_outputter_map_;

	// A mapping from the identifier from the (outputter-type, outputter_for_job message) pair for
	// a particular PoseOutputter/SecondaryPoseOutputter.
	typedef std::map< std::pair< std::string, std::string >, pose_outputters::PoseOutputterOP > PoseOutputterMap;
	typedef std::map< std::pair< std::string, std::string >, pose_outputters::SecondaryPoseOutputterOP > SecondaryOutputterMap;
	PoseOutputterMap pose_outputter_map_;
	SecondaryOutputterMap secondary_outputter_map_;

	// the secondary outputters that are requested given the command line flag
	std::list< pose_outputters::SecondaryPoseOutputterOP > cl_outputters_;

	// a temporary solution to the problem of not wanting to load the same pose repeatedly.
	// only works for PDBs currently -- a more general solution is required.
	// std::map< InputSourceAndImportPoseOptions, core::pose::PoseOP > previously_read_in_poses_;

	// Count how many input poses we've processed.
	core::Size input_pose_counter_;

	// The set of Poses that have been read in and are currently being held in memory.
	std::map< core::Size, core::pose::PoseOP > input_pose_index_map_;

	// For each pose that's still in memory, what is the index of the last job
	// that will use it as the starting structure?
	std::map< core::Size, core::Size > last_job_for_input_pose_;

};


} // namespace chunk_library
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_ChunkLibraryJobQueen_HH
