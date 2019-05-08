// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/standard/StandardJobQueen.hh
/// @brief  class declaration for StandardJobQueen
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) Simplification/Refactor

#ifndef INCLUDED_protocols_jd3_standard_StandardJobQueen_hh
#define INCLUDED_protocols_jd3_standard_StandardJobQueen_hh

// unit headers
#include <protocols/jd3/standard/StandardJobQueen.fwd.hh>

// package headers
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/JobDigraph.fwd.hh>
#include <protocols/jd3/JobOutputIndex.fwd.hh>
#include <protocols/jd3/JobTracker.fwd.hh>
#include <protocols/jd3/InnerLarvalJob.fwd.hh>
#include <protocols/jd3/pose_inputters/PoseInputSource.fwd.hh>
#include <protocols/jd3/pose_inputters/PoseInputter.fwd.hh>
#include <protocols/jd3/pose_inputters/PoseInputterCreator.fwd.hh>
#include <protocols/jd3/pose_outputters/PoseOutputSpecification.fwd.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.fwd.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterCreator.fwd.hh>
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.fwd.hh>
#include <protocols/jd3/standard/PreliminaryLarvalJob.fwd.hh>
#include <protocols/jd3/standard/PreliminaryLarvalJobTracker.fwd.hh>

// project headers
#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/import_pose_options.fwd.hh>

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

#ifdef PYROSETTA
#include <utility/tag/Tag.hh>
#include <utility/options/OptionCollection.hh>
#include <protocols/jd3/pose_inputters/PoseInputter.hh>
#include <protocols/jd3/pose_inputters/PoseInputterCreator.hh>
#include <protocols/jd3/pose_outputters/PoseOutputSpecification.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterCreator.hh>
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.hh>
#include <protocols/jd3/standard/PreliminaryLarvalJob.hh>
#endif

//c++ headers
#include <string>
#include <map>
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace standard {

typedef std::list< PreliminaryLarvalJob > PreliminaryLarvalJobs;

// A mapping from the outputter-type to a representative PoseOutputter/SecondaryPoseOutputter.
typedef std::map< std::string, pose_outputters::PoseOutputterOP > RepresentativeOutputterMap;
typedef std::map< std::string, pose_outputters::SecondaryPoseOutputterOP > SecondaryRepresentativeOutputterMap;

// A mapping from the identifier from the (outputter-type, outputter_for_job message) pair for
// a particular PoseOutputter/SecondaryPoseOutputter.
typedef std::map< std::pair< std::string, std::string >, pose_outputters::PoseOutputterOP > PoseOutputterMap;
typedef std::map< std::pair< std::string, std::string >, pose_outputters::SecondaryPoseOutputterOP > SecondaryOutputterMap;

typedef numeric::DiscreteIntervalEncodingTree< core::Size > SizeDIET;


/// @brief The %StandardJobQueen is meant to handle the most common form of Rosetta jobs where
/// a protocol is applied to a single input structure to generate a single output structure.
/// Most JobQueens in Rosetta will to derive from this JobQueen.  To make the process of deriving
/// new JobQueens easy, this class's API has been carefully designed to be easy to work with and
/// the class itself has been designed to do as much heavy lifting as possible and to leave only
/// the responsibilities that the %StandardJobQueen could not perform to the derived classes.
class StandardJobQueen : public JobQueen
{
public:
	/// @brief The StandardJobQueen constructor asks the PoseInputterFactory for a PoseInputter
	/// and creates a ResourceManager
	StandardJobQueen();

	~StandardJobQueen() override;


public:

	////////////////////////
	/// Starting the Job ///
	////////////////////////

	/// @brief The %StandardJobQueen's implementation is to not update the JobDAG at all: the
	/// most basic protocol defines a job DAG with only a single node.
	///
	/// Override this method if you have a need for updating the job dag during your protocol.
	void
	update_job_dag( JobDigraphUpdater & updater ) override;

	/// @brief Create a list of LarvalJobs for each Node Index.  Only can create max_njobs to save memory.
	///
	/// @brief
	/// The %StandardJobQueen manages the process of creating the list of LarvalJobs that
	///  will later be matured into actual jobs and run by the JobDistributor.
	///  It is not recommended that derived job queens override
	///  this method; doing so will mean that some of the data the SJQ relies on will not
	///  be initialized -- see comments on the SJQs data members below to understand the consequences
	///  of overriding this method.
	///
	/// It is recommended to override next_batch_of_larval_jobs_for_job_node instead of this method
	///  However, if you do override this - especially for more complicated protocols - be sure to call similar functions here
	///
	LarvalJobs
	determine_job_list( Size job_dag_node_index, Size max_njobs ) override;

	///@brief Mature the LarvalJob into a full Job that will be run on a processor.
	/// Calls complete_larval_job_maturation, which you will need to override.
	protocols::jd3::JobOP
	mature_larval_job(
		protocols::jd3::LarvalJobCOP job,
		utility::vector1< JobResultCOP > const & input_job_results
	) override;

public:

	///////////////////////////
	//// Completing the Job ///
	///////////////////////////

	///@brief The JD calls this function on completion of a LarvalJob, after updating the JobTracker
	///
	///@details
	///  If you override this method, call the SJQs version first for PJN tracking and output.
	void
	note_job_completed( LarvalJobCOP job, JobStatus status, Size nresults ) override;

	///@brief As each job completes, this function is called for each result coming from a LarvalJob.
	///
	///@details
	/// Override this method to store any results you need.
	void
	completed_job_summary( LarvalJobCOP job, Size result_index, JobSummaryOP summary ) override;

	///@brief The IDs of jobs that should be discarded, IE not kept in memory for the next set of job nodes.
	///
	///@details
	/// Override this method to note which jobs we don't need for the next set of Job Nodes.
	std::list< JobResultID >
	job_results_that_should_be_discarded() override;

	///@brief By default outputs all recently finished jobs for JD output (stored in recent_successes_).
	/// Override this method if you want to cull these lists
	std::list< output::OutputSpecificationOP >
	jobs_that_should_be_output() override;

public:


	/// @brief Return the bag of of PoseOutputters (in the form of a MultipleOutputter) for the
	/// Pose that has been requested and specified by a particular OutputSpecification;
	///
	/// @details
	///  This function guarantees that for each individual outputter-name (respecting the JD-provided
	///  filename suffix) that a separate set of PoseOutputters are returned so that multiple threads
	///  can potentially output at the same time.
	///
	/// You should not need to override this method.
	///
	output::ResultOutputterOP
	result_outputter( output::OutputSpecification const & spec ) override;

	///@brief Checks the outputter to see if the job already has been output.  Used for JD override behavior.
	///
	///@details
	/// Gets the outputter from the Job, which is cached in the SJQ.
	///
	/// You should not need to override this method
	///
	bool
	has_job_previously_been_output( protocols::jd3::LarvalJobCOP job ) override;

public:

	///////////////////////////
	/// Resource Management ///
	///////////////////////////

	std::list< deallocation::DeallocationMessageOP >
	deallocation_messages() override;

	/// @brief A deallocation message first sent to the JobDistributor on this host originating from
	/// a remote JobQueen. If a derived JobQueen has deallocation messages she needs to recieve,
	/// she should override derived_process_deallocation_message.
	void
	process_deallocation_message( deallocation::DeallocationMessageOP message ) override;

	void
	flush() override;

public:

	/////////////////////
	/// XSD Functions ///
	/////////////////////

	/// @brief The %StandardJobQueen assembles the XSD from virtual functions she invokes on
	/// the derived %JobQueen: append_job_tag_subelements, append_common_tag_subelements, and
	/// add_option/add_options.
	std::string
	job_definition_xsd() const override;

	std::string
	resource_definition_xsd() const override;

	/// @brief Read from an input string representing the contents of the job-definiton XML file
	/// and construct a set of PreliminaryLarvalJobs; this function is primarily useful for testing,
	/// but could be used to organize jobs by an enterprising job distributor or by another JobQueen.
	void
	determine_preliminary_job_list_from_xml_file( std::string const & job_def_string );

	/// @brief Naming function for the complexTypes in the job-definition XSD
	static
	std::string
	job_def_complex_type_name( std::string const & type );


	///////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////// Protected Functions ///////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

protected:

	/// @brief The job dag encodes which job nodes are independent or dependant on each other.
	///  If you have a applicationwhere all inputs are independant and only a single run of a protocol is needed,
	///  Then do not override this method.
	///
	/// @details
	///
	/// The %StandardJobQueen provides an implementation of this function which returns
	/// the most straight-forward DAG representing a set of jobs that have no interdependencies:
	/// a DAG with a single node and no edges.
	///
	/// Override this method if you have a need for any more job nodes other than what is created from input
	///  structures and the Job Definition file.
	///
	JobDigraphOP
	create_initial_job_dag() override;

	/// @brief The derived JobQueen must define the method that takes a larval job and the job-specific options
	/// and matures the larval job into a full job.
	virtual
	JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options,
		utility::vector1< JobResultCOP > const & input_job_results
	) = 0;

	///@brief Called during the SJQs determine_job_list function.
	/// Override this function to create LarvalJobs for nodes other than preliminary job nodes.
	virtual
	LarvalJobs
	next_batch_of_larval_jobs_for_job_node( Size job_dag_node_index, Size max_njobs );

	///@brief Get the current larval job index from the JobTracker.
	core::Size
	current_job_index() const;

protected:

	//////////////////
	/// XML Schema ///
	//////////////////

	/// @brief Here you can specify Schema's to use when creating your XSD for a specific JOB.
	/// AKA - you can add elements such as ResidueSelectors, task operations, etc. to parse.
	///
	/// @details
	///
	/// The derived JobQueen must inform the %StandardJobQueen of any additional tags that
	/// belong as sub elements of the <Job> tag. The %StandardJobQueen will have already appended
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

	/// @brief Here, you specify which tags that are defined in the <common> element.
	///
	/// @details
	///
	/// The derived JobQueen must inform the %StandardJobQueen of any additional tags that
	/// belong as sub elements of the <Common> tag. The %StandardJobQueen will have already appended
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
	/// provided for a particular job. (This is logic the %StandardJobQueen manages).
	virtual
	void
	append_common_tag_subelements(
		utility::tag::XMLSchemaDefinition & xsd,
		utility::tag::XMLSchemaComplexTypeGenerator & ct_gen
	) const;

	/// @author Jack Maguire, jackmaguire1444@gmail.com
	/// @brief This gives the derived JobQueen a chance to read any relevant information from the
	/// job-definition file. The second argument holds a PreliminaryLarvalJob for every job block
	/// and the PreliminaryLarvalJob itself has all of the information given in its corresponding
	/// block.
	virtual
	void
	parse_job_definition_tags(
		utility::tag::TagCOP, //common_block_tags,
		utility::vector1< PreliminaryLarvalJob > const &
	){}

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////


protected:

	///////////////////////////////
	/// Preliminary Larval Jobs ///
	///////////////////////////////

	/// @brief Initialize the PrelimininaryLarvalJob list from input poses and the Job Definition file.
	///
	/// @details
	///  Allow the derived JobQueen to tell the %StandardJobQueen to initialize the preliminary
	///  job list; this is perhaps necessary in the context of multi-round protocols when the
	///  job-definition file specifies the JobDAG. parse_job_definition_tags() is called at the
	///  end of this method.
	virtual
	void
	determine_preliminary_job_list();

	/// @brief Has determine_preliminary_job_list been called?
	bool
	get_preliminary_larval_jobs_determined() const;

	/// @brief Get write access to the PreliminaryLarvalJobTracker
	///
	PreliminaryLarvalJobTracker &
	get_prelim_larval_job_tracker();

	/// @brief Get read access to the PreliminaryLarvalJobTracker
	///
	PreliminaryLarvalJobTracker const &
	get_prelim_larval_job_tracker() const;

	/// @brief Read access to derived JobQueens to the preliminary job list.
	///
	/// @details
	///  This will return an empty list if  determine_preliminary_jobs has not yet
	///  been called.
	utility::vector1< PreliminaryLarvalJob > const &
	get_preliminary_larval_jobs() const;

	///@brief Respond to the completion of a preliminary job node in derived SJQs.
	///
	///@details
	/// Called during note_job_completed
	///
	/// The SJQ tracks which preliminary-job-nodes have had all of their jobs complete
	///  and lets the derived job queens know. DerivedJobQueens wishing to track when, e.g.,
	///  all of the jobs that might depend on a particular Resource held by the ResourceManager
	///  have completed can override this method and then be updated by the SJQ at the completion
	///  of every PJN.
	virtual void
	note_preliminary_job_node_is_complete( core::Size pjn_index );


protected:

	////////////////////////////////////////////////////////////////////
	/// Optional Larval Job Creation ///
	////////////////////////////////////////////////////////////////////

	/// @brief Ask the derived JobQueen to expand / refine a preliminary larval job, by
	/// possibly reading per-job data out of the Tag associated with each job.
	///
	/// @details
	///  If there is
	///  nothing that needs to be done by the derived class, it may elect to use the base-class
	///  implementation of this function which simply returns a list of the inner-job pointers
	///  in the input list.
	///  Why all the back and forth?  So that the derived job queens can take one input structure
	///  and multiply it into several sepearate jobs; e.g. in a ddG protocol, you need to run a
	///  separate set of jobs for the WT sequence and also for the mutant sequence. The preliminary
	///  job will reflect the structure of the WT, and then that will get expanded out.
	///  This funciton is non-const so that the derived class can even track which WT sequences it
	///  has seen so that after the first WT sequence for a particular input structure, it only expands
	///  out the mutant sequences if many mutants on that input.
	///
	/// If you derive this function, you will want to call the SJQs determine_job_list
	/// or next_batch_of_prelim_jobs within your own determine_job_list.
	///
	virtual InnerLarvalJobs
	refine_preliminary_job( PreliminaryLarvalJob const & prelim_job );

	/// @brief Expand an InnerLarvalJob into a full set of LarvalJobs, creating nstruct LarvalJob objects
	/// Increments larval_job_counter_
	virtual LarvalJobs
	expand_job_list( InnerLarvalJobOP inner_job, core::Size max_larval_jobs_to_create );

	/// @brief Adds information regarding the job definition
	/// (job tag, input source, jobdef tag, and outputter).
	/// Uses the PreliminaryInnerLarvalJob to do so.
	InnerLarvalJobOP
	create_and_init_inner_larval_job_from_preliminary( core::Size nstruct, core::Size prelim_job_node ) const;

protected:

	///////////////////////////////////////////
	/// Output Specifications and Ouputters ///
	///////////////////////////////////////////

	/// @brief Create an output specification for a given job + result index and store it in the
	/// recent_successes_ list, which will be given to the JobDistributor upon request.
	///
	/// @details
	///  This function invokes the factory method create_output_specification_for_job_result which
	///  derive job queens may override. This version of the function will create an OptionCollection
	///  and then invoke the version of the function that expects one.
	void
	create_and_store_output_specification_for_job_result(
		LarvalJobCOP job,
		core::Size result_index,
		core::Size nresults
	);

	/// @brief Create an output specification for a given job + result index and store it in the
	/// recent_successes_ list, which will be given to the JobDistributor upon request.
	///
	/// @details
	///  This function invokes the factory method create_output_specification_for_job_result which
	///  derive job queens may override. This version of the function accepts an OptionCollection
	///  if there are multiple results for a particular job, and you can reuse the OptionCollection
	///  between calls.
	void
	create_and_store_output_specification_for_job_result(
		LarvalJobCOP job,
		utility::options::OptionCollection const & job_options,
		core::Size result_index,
		core::Size nresults
	);

	/// @brief Construct the OutputSpecification for a completed successful job.
	///
	/// @details
	///  Factory method that may be overriden by the derived job queen. Invoked by
	///  create_and_store_output_specification_for_job_result which is invoked via note_job_completed
	///  by default, and may be invoked by the derived job queen as needed.
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
	/// @details Just before a job result is written to disk, the %StandardJobQueen asks the derived
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

	///@brief Create or get the PoseOutputter from the job.
	/// Cache the Outputter in the SJQ
	pose_outputters::PoseOutputterOP
	pose_outputter_for_job( InnerLarvalJob const & innerJob );

	///@brief Create or get the PoseOutputter from the job.
	/// Cache the Outputter in the SJQ
	pose_outputters::PoseOutputterOP
	pose_outputter_for_job( InnerLarvalJob const & innerJob, utility::options::OptionCollection const & job_options );

	/// @brief Find the PoseOutputter for a job specification, respecting the Job-distributor
	/// provided output file suffix.
	///
	/// @details
	///  The %StandardJobQueen guarantees that a different Outputter
	///  is used for each unique name as long as that name is not the empty string. This allows
	///  multiple threads to write to separate files concurrently.
	pose_outputters::PoseOutputterOP
	pose_outputter_for_job( pose_outputters::PoseOutputSpecification const & spec );

	utility::vector1< pose_outputters::SecondaryPoseOutputterOP >
	secondary_outputters_for_job( InnerLarvalJob const & innerJob, utility::options::OptionCollection const & job_options );

protected:

	////////////////////
	//// Pose Access ///
	////////////////////

	core::pose::PoseOP
	pose_for_inner_job( InnerLarvalJobCOP inner_job );

	/// @brief Return a copy of the Pose to be used with the given inner job
	core::pose::PoseOP
	pose_for_inner_job( InnerLarvalJobCOP inner_job, utility::options::OptionCollection const & options );

	/// @brief Return a copy of the Pose to be used with the given job. Wrapper for pose_for_inner_job()
	core::pose::PoseOP
	pose_for_job( LarvalJobCOP job, utility::options::OptionCollection const & options );


protected:

	//////////////////////////
	/// Options Collection ///
	//////////////////////////

	///@brief Create an OptionCollection from the inner_job and OptionKeys stored in the SJQ.
	utility::options::OptionCollectionOP
	options_for_job( InnerLarvalJob const & inner_job ) const;

	///@brief Create an OptionCollection from tags and OptionKeys stored in the SJQ.
	utility::options::OptionCollectionOP
	options_from_tag( utility::tag::TagCOP job_options_tags  ) const;

	// Of course the job inputter might vary from job to job!
	pose_inputters::PoseInputterOP
	pose_inputter_for_job( InnerLarvalJob const & inner_job ) const;


	/////////////////////////////////////////////////////////////////////////////////
	// The following functions are to be used by derived JobQueens to signal to the
	// StandardJobQueen that individual options can be specified on a per-job basis
	// in the XML jobs file or that can be specified on the command line. Derived
	// JobQueens are discouraged from reading from the command line directly, since
	// that circumvents the auto-documenting features of the XSD in terms of
	// communicating to the user how Jobs are to be specified. These functions should
	// be invoked during the derived class's constructor.

	void
	add_options( utility::options::OptionKeyList const & opts );

	void
	add_option( utility::options::OptionKey const & key );






	///////////////////////////////////////////////////////////////////////////////////
	// The following functions should be used by the derived JobQueens to communicate
	// to the StandardJobQueen the structure of the XML document that they support.
	// The StandardJobQueen will prepare the output XSD document based on the function
	// calls here, and will also validate the input job-definition file, ensuring that
	// no options are provided to through that file that are not read.

	/// @brief Override the basic specification of input elements as either PDBs, Silent Structures,
	/// or ResourceManager defined structures. This should be called during the derived class's
	/// constructor if it is to be used.
	void
	remove_default_input_element();

	/// @brief Return the set of subtags that are common to the whole set of jobs in the
	/// Job definition file, if any are given.  This set of tags is read from disk at most
	/// once per execution.
	utility::tag::TagCOP
	common_block_tags() const;

	///@brief By default, the common block is given before the job blocks in the job definition file.
	///Setting this value to false reverses the order and the common block goes at the end of the file.
	void set_common_block_precedes_job_blocks( bool setting ){
		common_block_precedes_job_blocks_ = setting;
	}

	//
	//
	///////////////////////////////////////////////////////////////////////////////////

protected:

	/////////////////////////////////////////////
	/// Expert Override Input/Output behavior ///
	/////////////////////////////////////////////

	////////////////////////////////////////////////////////////
	///                                                      ///
	/// Note: (Used for Complex I/O tsuch as Membrane poses) ///
	///                                                      ///
	////////////////////////////////////////////////////////////


	/// @brief The derived job queen (DJQ) may tell the %StandardJobQueen to only accept a subset of
	/// PoseInputters. First the DJQ says "do_not_accept_all_pose_inputters_from_factory" and then she
	/// can call "allow_pose_inputter." This function must only be called the DJQ's constructor.
	void
	do_not_accept_all_pose_inputters_from_factory();

	/// @details If the derived job queen (DJQ) wants to allow a PoseInputter to be used either because
	/// it is not registered with the PoseInputterFactory, or because she wants to allow only a
	/// subset of the inputters that are registered with the Factory, then she may indicate that
	/// to the %StandardJobQueen with this function. If the derived job queen does not want to
	/// use all of the pose inputters provided by the factory, then she must call
	/// "do_not_accept_all_pose_inputters_from_factory()" before she calls this function.
	/// This function must only be called in the DJQ's constructor. All DJQs must invoke this
	/// method in their constructors and not simply the DJQ that lives on the "head node"
	/// in order for distributed output to workk correctly.
	void
	allow_pose_inputter( pose_inputters::PoseInputterCreatorOP creator );

	/// @details The derived job queen (DJQ) may tell the %StandardJobQueen to only accept a subset of
	/// PoseOutputters. First the DJQ says "do_not_accept_all_pose_outputters_from_factory" and then
	/// she can call "allow_pose_outputter." This function must only be called the DJQ's constructor.
	void
	do_not_accept_all_pose_outputters_from_factory();

	/// @details If the derived job queen (DJQ) wants to allow a PoseOutputter to be used either because
	/// it is not registered with the PoseOutputterFactory, or because she wants to allow only a
	/// subset of the outputters that are registered with the Factory, then she may indicate that
	/// to the %StandardJobQueen with this function. If the derived job queen does not want to
	/// use all of the pose outputters provided by the factory, then she must call
	/// "do_not_accept_all_pose_outputters_from_factory()" before she calls this function.
	/// If the DJQ has called the "do_not_accept_all_pose_outputters_from_factory" function, then
	/// she specifies the default outputter on the first time she calls "allow_pose_outputter."
	/// This function must only be called in the DJQ's constructor.
	void
	allow_pose_outputter( pose_outputters::PoseOutputterCreatorOP creator );

	/// @details If the derived job queen wants the user to get a particular PoseOutputter by default
	/// (i.e. when the user hasn't specified which outputter to use in a job definition file and
	/// in the absence of a command-line flag that says which outputter to use) instead of the
	/// default PoseOutputter chosen by the PoseOutputterFactory, then the DJQ should call this function.
	/// If the DJQ has previously called "do_not_accept_all_pose_outputters_from_factory", then
	/// the first call to "allow_pose_outputter" specified the default outputter, and a separate call
	/// to this function is not absolutely necessary. This function must only be called in the
	/// DJQ's constructor.
	void
	set_default_outputter( pose_outputters::PoseOutputterCreatorOP creator );


	//////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// Private Functions ///////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////


private:


	/// @brief After generating the job-definition XSD, construct the preliminary job
	/// list.
	///
	/// @details
	///  This is invoked both from determine_job_list_from_xml_file and
	///  determine_job_list -- the latter always constructs an XSD to ensure
	///  that the derived JobQueen has properly constructed an XSD, even if
	///  a job definition file has not been provided on the command line.
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

	pose_inputters::PoseInputSourcesAndInputters
	read_input_poses_from_command_line();

	///@brief This function is called in determine_job_list.  It is responsible for creating the set of LarvalJobs
	/// for ALL PreliminaryLarvalJobs.  Call this function if you derive determine_job_list.
	///
	///@details
	///  Updates PreliminaryLarvalJobTracker
	///
	LarvalJobs
	next_batch_of_larval_jobs_from_prelim( core::Size job_node_index, core::Size max_njobs );

	///@brief Get the SecondardyPoseOutputter for a job.  An example of an SPO is a scorefile.
	///
	pose_outputters::SecondaryPoseOutputterOP
	secondary_outputter_for_job(
		InnerLarvalJob const & inner_job,
		utility::options::OptionCollection const & job_options,
		std::string const & secondary_outputter_type,
		utility::tag::TagCOP outputter_tag
	);

	///@brief Get the SecondardyPoseOutputter for a job.  An example of an SPO is a scorefile.
	///
	pose_outputters::SecondaryPoseOutputterOP
	secondary_outputter_for_job(
		pose_outputters::PoseOutputSpecification const &
	);

	/// @brief The SJQ will keep track of all discarded jobs in the discarded_jobs_ diet
	void
	note_job_result_discarded( LarvalJobCOP job, Size result_index );

	utility::options::OptionKeyList
	concatenate_all_options() const;

	pose_outputters::PoseOutputterOP
	get_outputter_from_job_tag( utility::tag::TagCOP tag ) const;

private:

	// ResourceManagerOP resource_manager_;

	// The set of options that the %StandardJobQueen reads from the command line
	// and/or from the job definition XML file.
	utility::options::OptionKeyList options_;
	utility::options::OptionKeyList inputter_options_;
	utility::options::OptionKeyList outputter_options_;
	utility::options::OptionKeyList secondary_outputter_options_;

	bool use_factory_provided_pose_inputters_;
	std::map< std::string, pose_inputters::PoseInputterCreatorCOP > inputter_creators_;
	std::list< pose_inputters::PoseInputterCreatorCOP > inputter_creator_list_;

	bool use_factory_provided_pose_outputters_;
	std::map< std::string, pose_outputters::PoseOutputterCreatorCOP > outputter_creators_;
	std::list< pose_outputters::PoseOutputterCreatorCOP > outputter_creator_list_;
	bool override_default_outputter_;
	pose_outputters::PoseOutputterCreatorOP default_outputter_creator_;

	// Often, you want to use the same pose outputter for multiple jobs, E.g.
	// if you are writing multiple Poses to a single silent file and are buffering
	// the output in the Outputter to reduce the number of times you push data to disk
	std::map< std::string, pose_outputters::PoseOutputterOP > pose_outputters_;

	// By default, the <Common> block goes first in the job definition file, but derived
	// job queens may request that the jobs are listed first and the <Common> block loast.
	bool common_block_precedes_job_blocks_;

	// All job queens, even those that spawn on the remote nodes, must perform some degree of
	// initialization. In particular, they need to read the job definition file to obtain
	// the <Common> block data
	bool required_initialization_performed_;
	utility::tag::TagCOP job_definition_file_tags_;
	utility::tag::TagCOP common_block_tags_;

	// For the first node in the JobDAG, the %StandardJobQueen will spool out LarvalJobs
	// slowly to the JobDistributor (in increments of the max_njobs parameter in the call
	// to job_dag_node_index). Since max_njobs may be smaller than the nstruct parameter,
	// the %StandardJobQueen will need to be able to interrupt the spooling of jobs until
	// the JobDistributor is ready for them. For this reason, it keeps what is effectively
	// a set of indices into a while loop for LarvalJob construction.
	bool preliminary_larval_jobs_determined_;
	utility::vector1< PreliminaryLarvalJob > preliminary_larval_jobs_;

	PreliminaryLarvalJobTrackerOP prelim_job_tracker_;

	InnerLarvalJobs inner_larval_jobs_for_curr_prelim_job_;
	Size curr_inner_larval_job_index_;
	Size njobs_made_for_curr_inner_larval_job_;

	// The jobs that have completed, but have not yet been output, and the instructions for how
	// to output them.
	std::list< output::OutputSpecificationOP > recent_successes_;

	RepresentativeOutputterMap representative_pose_outputter_map_;
	SecondaryRepresentativeOutputterMap representative_secondary_outputter_map_;

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

	bool load_starting_poses_only_once_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace standard
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_standard_StandardJobQueen )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_StandardJobQueen_HH
