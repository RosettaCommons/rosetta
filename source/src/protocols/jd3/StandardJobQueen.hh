// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/StandardJobQueen.hh
/// @brief  class declaration for StandardJobQueen
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_StandardJobQueen_hh
#define INCLUDED_protocols_jd3_StandardJobQueen_hh

// unit headers
#include <protocols/jd3/StandardJobQueen.fwd.hh>

// package headers
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/InnerLarvalJob.fwd.hh>
#include <protocols/jd3/pose_inputters/PoseInputter.fwd.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.fwd.hh>

// project headers
#include <core/pose/Pose.fwd.hh>
#include <basic/resource_manager/JobOptions.fwd.hh>


//utility headers
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>
#include <utility/options/keys/all.fwd.hh>

//c++ headers
#include <string>
#include <map>

namespace protocols {
namespace jd3 {

class PreliminaryLarvalJob
{
public:
	PreliminaryLarvalJob();
	~PreliminaryLarvalJob();
	PreliminaryLarvalJob( PreliminaryLarvalJob const & src );
	PreliminaryLarvalJob & operator = ( PreliminaryLarvalJob const & rhs );

	InnerLarvalJobOP inner_job;
	utility::tag::TagCOP job_tag;
	utility::options::OptionCollectionCOP job_options;
};

typedef std::list< PreliminaryLarvalJob > PreliminaryLarvalJobs;

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

	std::string job_definition_xsd() const override;
	std::string resource_definition_xsd() const override;

	/// @brief The %StandardJobQueen manages the process of creating the list of LarvalJobs that
	/// will be later matured into actual jobs and run by the JobDistributor.  Classes derived
	/// from the %StandardJobQueen need answer a few virtual functions that the %StandardJobQueen
	/// will invoke during this process.
	///
	/// @details The process begins by first constructing the job definition and resource definition
	/// XSDs.  With these schemas, the %StandardJobQueen validates the input XML files (if present).
	/// The %StandardJobQueen then populates preliminary versions of LarvalJob objects./ If the XSD
	/// includes "command line options" (which may be specified either from the command line or in
	/// the <options> section of the Job XML file), the %StandardJobQueen loads the preliminary
	/// LarvalJob objects with the options. These preliminary LarvalJob objects will not have been
	/// nstruct expanded (i.e. if there are 100 nstruct for each of 5 different jobs, then there will
	/// only be 5 preliminary larval jobs created). It then passes the preliminary LarvalJob list and
	/// the TagOP objects for each preliminary LarvalJob to the derived class through the
	/// refine_job_list method.
	LarvalJobs determine_job_list() override;

	bool has_job_completed( protocols::jd3::LarvalJobCOP job ) override;


	protocols::jd3::JobOP
	mature_larval_job( protocols::jd3::LarvalJobCOP job ) override;

	/// @brief Read from an input string representing the contents of the job-definiton XML file
	/// and construct a set of LarvalJobs; this function is primarily useful for testing,
	/// but could be used to organize jobs by an enterprising job distributor or by another JobQueen.
	LarvalJobs determine_job_list_from_xml_file( std::string const & job_def_string );

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
	virtual InnerLarvalJobs refine_preliminary_job( PreliminaryLarvalJob const & prelim_job );

	/// @brief Expand an InnerLarvalJob into a full set of LarvalJobs, creating nstruct LarvalJob objects
	/// The base class implementation of this function invokes the create_larval_job factory method so
	/// that dervied JobQueens can ensure the instantiation of derived LarvalJobs.
	virtual LarvalJobs expand_job_list( InnerLarvalJobOP inner_job ) const;

	/// @brief Factory method for derived classes if they wish to rely on classes derived
	/// from InnerLarvalJob.  This is invoked by the StandardJobQueen in her determine_job_list
	/// method just as jobs are prepared. If the base InnerLarvalJob class is desired, then do
	/// not override this method.
	virtual InnerLarvalJobOP create_inner_larval_job() const;

	/// @brief Factory method for derived classes if they wish to rely on classes derived
	/// from LarvalJob.  This is invoked by the StandardJobQueen in the expand_job_list method.
	/// If the base LarvalJob class is desired, then do not override this method.
	virtual LarvalJobOP create_larval_job( InnerLarvalJobOP job, core::Size nstruct_index ) const;

	/// @brief Factory method for derived classes if they wish to rely on a class besides the
	/// MoverAndPoseJob, which is returned by the %StandardJobQueen's implementation of this
	/// function.
	virtual JobOP create_job( LarvalJobCOP job ) const;

	virtual
	JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options
	) const = 0;

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
	core::pose::PoseOP pose_for_job( LarvalJobCOP job, utility::options::OptionCollection const & options ) const;

	// ResourceManagerOP resource_manager();

	// Of course the job inputter might vary from job to job!
	pose_inputters::PoseInputterOP
	pose_inputter_for_job( InnerLarvalJob const & inner_job ) const;

	// Of course the job outputter might vary from job to job!
	pose_outputters::PoseOutputterOP
	pose_outputter_for_job( InnerLarvalJob const & innerJob ) const;

	core::Size
	nstruct_for_job( InnerLarvalJob const & inner_job ) const;

	utility::options::OptionCollectionOP
	options_for_job( InnerLarvalJob const & inner_job ) const;

	utility::options::OptionCollectionOP
	options_from_tag( utility::tag::TagCOP job_options_tags ) const;

private:

	/// @brief After generating the job-definition XSD, construct the job list.
	/// This is invoked both from determine_job_list_from_xml_file and
	/// determine_job_list -- the latter always constructs an XSD to ensure
	/// that the derived JobQueen has properly constructed an XSD, even if
	/// a job definition file has not been provided on the command line.
	LarvalJobs
	determine_job_list_from_xml_file(
		std::string const & job_def_string,
		std::string const & job_def_schema
	);

	/// @brief Instead of reading a JobDefinition file, construct the set of LarvalJobs
	/// reading from the command line. Invoked by determine_job_list.
	LarvalJobs
	determine_job_list_from_command_line();


	/// @brief After constructing a PreliminaryLarvalJob, ask the derived JobQueen
	/// to refine that job into (potentially more) InnerLarvalJobs, and to then
	/// construct a full complement of LarvalJob objects and splice them into the input
	/// LarvalJobs list.
	void
	expand_preliminary_larval_job(
		PreliminaryLarvalJob const & prelim_job,
		pose_outputters::PoseOutputterOP outputter,
		utility::options::OptionCollectionCOP job_options,
		utility::tag::TagCOP output_tag, // the <Output> subtag of the <Job> tag, if present
		LarvalJobs & jobs
	);


private:

	// ResourceManagerOP resource_manager_;

	// The set of options that the %StandardJobQueen reads from the command line
	// and/or from the job definition XML file.
	utility::options::OptionKeyList options_;

	// Often, you want to use the same pose outputter for multiple jobs.
	std::map< std::string, pose_outputters::PoseOutputterOP > pose_outputters_;

	utility::tag::TagCOP common_block_tags_;
};



} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_StandardJobQueen_HH
