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

#ifndef INCLUDED_protocols_jd3_standard_StandardJobQueen_hh
#define INCLUDED_protocols_jd3_standard_StandardJobQueen_hh

// unit headers
#include <protocols/jd3/standard/StandardJobQueen.fwd.hh>

// package headers
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/JobDigraph.fwd.hh>
#include <protocols/jd3/PoseInputSource.fwd.hh>
#include <protocols/jd3/pose_inputters/PoseInputter.fwd.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.fwd.hh>
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.fwd.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.fwd.hh>

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
namespace standard {

class PreliminaryLarvalJob
{
public:
	PreliminaryLarvalJob();
	~PreliminaryLarvalJob();
	PreliminaryLarvalJob( PreliminaryLarvalJob const & src );
	PreliminaryLarvalJob & operator = ( PreliminaryLarvalJob const & rhs );

	StandardInnerLarvalJobOP inner_job;
	utility::tag::TagCOP job_tag;
	utility::options::OptionCollectionCOP job_options;
	pose_inputters::PoseInputterOP pose_inputter;
};

typedef std::list< PreliminaryLarvalJob > PreliminaryLarvalJobs;

///// @brief A small class to figure out which Poses do not need to be loaded
///// a second time because they have already been read in. This class mainly
///// stores the PoseInputSource and the set of options that are used to turn
///// a file into a Pose. I don't imagine this class sticking around that long;
///// it is only truely compatible with PDB/mmCIF file reading.
//class InputSourceAndImportPoseOptions
//{
//public:
// InputSourceAndImportPoseOptions();
// InputSourceAndImportPoseOptions(
//  PoseInputSource const & input_source,
//  core::import_pose::ImportPoseOptions const & options
// );
// InputSourceAndImportPoseOptions( InputSourceAndImportPoseOptions const & );
//
// virtual ~InputSourceAndImportPoseOptions();
//
// InputSourceAndImportPoseOptions & operator = ( InputSourceAndImportPoseOptions const & rhs );
//
// bool operator == ( InputSourceAndImportPoseOptions const & ) const;
// bool operator <  ( InputSourceAndImportPoseOptions const & ) const;
//
// PoseInputSource const & input_source() const;
// void input_source( PoseInputSource const & setting );
//
// core::import_pose::ImportPoseOptions const & import_pose_options() const;
// void import_pose_options( core::import_pose::ImportPoseOptions const & setting );
//
//private:
// PoseInputSourceOP input_source_;
// core::import_pose::ImportPoseOptionsOP import_pose_options_;
//
//};

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


	/// @brief The %StandardJobQueen assembles the XSD from virtual functions she invokes on
	/// the derived %JobQueen: append_job_tag_subelements, append_common_tag_subelements, and
	/// add_option/add_options.
	std::string job_definition_xsd() const override;
	std::string resource_definition_xsd() const override;

	/// @brief The %StandardJobQueen provides an implementation of this function which returns
	/// the most straight-forward DAG representing a set of jobs that have no interdependencies:
	/// a DAG with a single node and no edges.
	JobDigraphOP initial_job_dag() override;

	/// @brief The %StandardJobQueen's implementation is to not update the JobDAG at all: the
	/// most basic protocol defines a job DAG with only a single node.
	void update_job_dag( JobDigraphUpdater & updater ) override;

	/// @brief The %StandardJobQueen manages the process of creating the list of LarvalJobs that
	/// will be later matured into actual jobs and run by the JobDistributor.  Classes derived
	/// from the %StandardJobQueen need answer a few virtual functions that the %StandardJobQueen
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
	void note_job_completed( LarvalJobCOP job, JobStatus status ) override;
	void note_job_completed( core::Size job_id, JobStatus status ) override;

	bool larval_job_needed_for_completed_job_summary() const override;
	void completed_job_summary( LarvalJobCOP job, JobSummaryOP summary ) override;
	void completed_job_summary( core::Size job_id, JobSummaryOP summary ) override;

	std::list< core::Size > jobs_that_should_be_output() override;
	std::list< core::Size > job_results_that_should_be_discarded() override;
	void completed_job_result( LarvalJobCOP job, JobResultOP job_result ) override;


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

protected:

	/// @brief The derived JobQueen must inform the %StandardJobQueen of any additional tags that
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

	/// @brief The derived JobQueen must inform the %StandardJobQueen of any additional tags that
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

	/// @brief Allow the derived JobQueen to tell the %StandardJobQueen to initialize the preliminary
	/// job list; this is perhaps necessary in the context of multi-round protocols when job-definition
	/// file specifies the JobDAG.
	virtual
	void
	determine_preliminary_job_list();

	/// @brief Allow the derived job queen the opportunity to update the StandardJobQueen's
	/// job_graph_ data member by adding nodes as well as edges that land on the new nodes.
	JobDigraphUpdater
	updater_for_sjq_job_graph();

	/// @brief Read access for the subset of nodes in the job DAG which the %StandardJobQueen
	/// is responsible for producing the larval_jobs. They are called "preliminary" jobs because
	/// they do not depend on outputs from any previous node in the graph. (The set of job nodes
	/// that contain no incoming edges, though, could perhaps be different from the set of
	/// preliminary job nodes, so the %StandardJobQueen requires that the deried job queen
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

	/////// @brief Allow the derived job queen to specify a node in the JobDAG as being
	/////// "preliminary" in the sense a) that the %StandardJobQueen is responsible for creating the
	/////// list of larval jobs for this node, and b) there are no nodes that this node depends
	/////// on having completed before it can run.
	////virtual
	////void
	////declare_job_node_to_be_preliminary( core::Size job_node_index );

	/// @brief Read access to derived JobQueens to the preliminary job list.
	/// This will return an empty list if  determine_preliminary_jobs has not yet
	/// been called.
	utility::vector1< PreliminaryLarvalJob > const &
	preliminary_larval_jobs() const;

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
	/// of this DIET, then the job has already been output.
	numeric::DiscreteIntervalEncodingTree< core::Size > const & output_jobs() const;

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
	virtual StandardInnerLarvalJobs refine_preliminary_job( PreliminaryLarvalJob const & prelim_job );

	/// @brief Expand a StandardInnerLarvalJob into a full set of LarvalJobs, creating nstruct LarvalJob objects
	/// The base class implementation of this function invokes the create_larval_job factory method so
	/// that dervied JobQueens can ensure the instantiation of derived LarvalJobs.
	virtual LarvalJobs expand_job_list( StandardInnerLarvalJobOP inner_job, core::Size max_larval_jobs_to_create );

	/// @brief Factory method for derived classes if they wish to rely on classes derived
	/// from StandardInnerLarvalJob.  This is invoked by the StandardJobQueen in her determine_job_list
	/// method just as jobs are prepared. If the base StandardInnerLarvalJob class is desired, then do
	/// not override this method.
	virtual StandardInnerLarvalJobOP create_inner_larval_job( core::Size nstruct, core::Size prelim_job_node ) const;

	/// @brief Factory method for derived classes if they wish to rely on classes derived
	/// from LarvalJob.  This is invoked by the StandardJobQueen in the expand_job_list method.
	/// If the base LarvalJob class is desired, then do not override this method.
	virtual LarvalJobOP create_larval_job( StandardInnerLarvalJobOP job, core::Size nstruct_index, core::Size larval_job_index );

	/// @brief Factory method for derived classes if they wish to rely on a class besides the
	/// MoverAndPoseJob, which is returned by the %StandardJobQueen's implementation of this
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

	/// @brief The StandardJobQueen cannot readily manage the complexity of organizing
	/// larval jobs for JobDAG nodes beyond the preliminary nodes. All job-creation logic beyond
	/// the first set of nodes is the responsibility of the derived JobQueen. The %StandardJobQueen
	/// provides a noop implementation of this function.
	///
	/// @details For the sake of outputting structures with the properties you want at the
	/// completion of a multi-stage trajectory, derived job queens with a tree-like job-progression
	/// (i.e. only one parent job for each child job) should copy the jobdef_tag into
	/// new InnerLarvalJobs from the parent job. The StandardJobQueen does not track parentage
	/// so it cannot know that job 873 originated from job 4 when trying to figure out how to output
	/// the structure for job 873. The derived job queen should also ensure that the "outputter"
	/// data from the parent job is copied into child job.
	virtual
	LarvalJobs
	next_batch_of_larval_jobs_for_job_node( core::Size job_dag_node_index, core::Size max_njobs );

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
	// StandardJobQueen that individual options can be specified on a per-job basis
	// in the XML jobs file or that can be specified on the command line. Derived
	// JobQueens are discouraged from reading from the command line directly, since
	// that circumvents the auto-documenting features of the XSD in terms of
	// communicating to the user how Jobs are to be specified. These functions should
	// be invoked during the derived class's constructor.

	void add_options( utility::options::OptionKeyList const & opts );
	void add_option( utility::options::OptionKey const & key );

	///////////////////////////////////////////////////////////////////////////////////
	// The following functions should be used by the derived JobQueens to communicate
	// to the StandardJobQueen the structure of the XML document that they support.
	// The StandardJobQueen will prepare the output XSD document based on the function
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
	pose_inputters::PoseInputterOP
	pose_inputter_for_job( StandardInnerLarvalJob const & inner_job ) const;

	// Of course the job outputter might vary from job to job!
	pose_outputters::PoseOutputterOP
	pose_outputter_for_job( StandardInnerLarvalJob const & innerJob );

	pose_outputters::PoseOutputterOP
	pose_outputter_for_job( StandardInnerLarvalJob const & innerJob, utility::options::OptionCollection const & job_options );

	std::list< pose_outputters::SecondaryPoseOutputterOP >
	secondary_outputters_for_job( StandardInnerLarvalJob const & innerJob, utility::options::OptionCollection const & job_options );

	core::Size
	nstruct_for_job( utility::tag::TagCOP job_tag ) const;

	utility::options::OptionCollectionOP
	options_for_job( StandardInnerLarvalJob const & inner_job ) const;

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
		StandardInnerLarvalJob const & inner_job,
		utility::options::OptionCollection const & job_options,
		std::string const & secondary_outputter_type
	);

	/// @brief The SJQ will keep track of all output jobs in the output_jobs_ diet, but requires
	/// that derived classes who are not calling the SJQ's version of completed_job_result call
	/// this function.
	void note_job_result_output( LarvalJobCOP job );

private:

	// ResourceManagerOP resource_manager_;

	// The set of options that the %StandardJobQueen reads from the command line
	// and/or from the job definition XML file.
	utility::options::OptionKeyList options_;

	// Often, you want to use the same pose outputter for multiple jobs.
	std::map< std::string, pose_outputters::PoseOutputterOP > pose_outputters_;

	bool required_initialization_performed_;
	utility::tag::TagCOP job_definition_file_tags_;
	utility::tag::TagCOP common_block_tags_;

	JobDigraphOP job_graph_;

	// The index of the last larval job that we have created. Incremented within
	// expand_job_list()
	Size larval_job_counter_;

	// For the first node in the JobDAG, the %StandardJobQueen will spool out LarvalJobs
	// slowly to the JobDistributor (in increments of the max_njobs parameter in the call
	// to job_dag_node_index). Since max_njobs may be smaller than the nstruct parameter,
	// the %StandardJobQueen will need to be able to interrupt the spooling of jobs until
	// the JobDistributor is ready for them. For this reason, it keeps what is effectively
	// a set of indices into a while loop for LarvalJob construction.
	bool preliminary_larval_jobs_determined_;
	utility::vector1< PreliminaryLarvalJob > preliminary_larval_jobs_;
	StandardInnerLarvalJobs inner_larval_jobs_for_curr_prelim_job_;
	Size curr_inner_larval_job_index_;
	Size njobs_made_for_curr_inner_larval_job_;
	utility::vector1< core::Size > preliminary_job_node_inds_;
	utility::vector1< core::Size > pjn_job_ind_begin_;
	utility::vector1< core::Size > pjn_job_ind_end_;
	utility::vector1< char > preliminary_job_nodes_complete_; // complete == all jobs assigned

	numeric::DiscreteIntervalEncodingTree< core::Size > completed_jobs_;
	numeric::DiscreteIntervalEncodingTree< core::Size > successful_jobs_;
	numeric::DiscreteIntervalEncodingTree< core::Size > failed_jobs_;
	numeric::DiscreteIntervalEncodingTree< core::Size > output_jobs_;
	std::list< core::Size > recent_successes_;

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


} // namespace standard
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_StandardJobQueen_HH
