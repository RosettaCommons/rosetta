// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <protocols/jd3/PoseInputter.fwd.hh>
#include <protocols/jd3/PoseOutputter.fwd.hh>

// project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/options/keys/all.fwd.hh>

//c++ headers
#include <string>

namespace protocols {
namespace jd3 {

class StandardJobQueen : public JobQueen
{
public:
	/// @brief The StandardJobQueen constructor asks the PoseInputterFactory for a PoseInputter
	/// and creates a ResourceManager
	StandardJobQueen();

	virtual ~StandardJobQueen();

	virtual std::string job_definition_xsd() const;

protected:

	/// @brief Read from the command line and from the job definition file and prepare a preliminary
	/// list of jobs.  This function will verify that the input job definition file meets the requirement
	/// of the standard-form job definition XSD.  This code does not perform the formal verification
	/// performed by a true XSD validator, but it will make sure that the input job-definition XML file does
	/// not contain any attributes or options that were not given to the base class.  After the InnerLarvalJob
	/// objects are created, the derived JobQueen may create more InnerJobs or may leave them as-is, before
	/// then calling the convenience function expand_job_list.
	InnerLarvalJobs prepare_preliminary_job_list();

	/// @brief Expand the list of InnerJobs into a full set of Jobs, creating nstruct Job objects for each
	/// inner job.  This is a convenience function.
	LarvalJobs expand_job_list( InnerLarvalJobs const & inner_jobs ) const;

	/// @brief Factory method for derived classes if they wish to rely on classes derived
	/// from LarvalJob.  This is invoked by the StandardJobQueen in the expand_job_list method.
	virtual LarvalJobOP create_job( InnerLarvalJobOP job, core::Size nstruct_index ) const;

	/////////////////////////////////////////////////////////////////////////////////
	// The following functions are to be used by derived JobQueens to signal to the
	// StandardJobQueen that individual options can be specified on a per-job basis
	// in the XML jobs file or that can be specified on the command line. Derived
	// JobQueens are discouraged from reading from the command line directly, since
	// that circumvents the auto-documenting features of the XSD in terms of
	// communicating to the user how Jobs are to be specified.

	/// @brief Indicate to the StandardJobInputter that a Boolean option can be specified on a per-job basis.
	void add_option( utility::options::BooleanOptionKey const & key );
	/// @brief Indicate to the StandardJobInputter that a BooleanVector option can be specified on a per-job basis.
	void add_option( utility::options::BooleanVectorOptionKey const & key );
	/// @brief Indicate to the StandardJobInputter that a File option can be specified on a per-job basis.
	void add_option( utility::options::FileOptionKey const & key );
	/// @brief Indicate to the StandardJobInputter that a FileVector option can be specified on a per-job basis.
	void add_option( utility::options::FileVectorOptionKey const & key );
	/// @brief Indicate to the StandardJobInputter that an Integer option can be specified on a per-job basis.
	void add_option( utility::options::IntegerOptionKey const & key );
	/// @brief Indicate to the StandardJobInputter that an IntegerVector option can be specified on a per-job basis.
	void add_option( utility::options::IntegerVectorOptionKey const & key );
	/// @brief Indicate to the StandardJobInputter that a Path option can be specified on a per-job basis.
	void add_option( utility::options::PathOptionKey const & key );
	/// @brief Indicate to the StandardJobInputter that a PathVector option can be specified on a per-job basis.
	void add_option( utility::options::PathVectorOptionKey const & key );
	/// @brief Indicate to the StandardJobInputter that a Real option can be specified on a per-job basis.
	void add_option( utility::options::RealOptionKey const & key );
	/// @brief Indicate to the StandardJobInputter that a RealVector option can be specified on a per-job basis.
	void add_option( utility::options::RealVectorOptionKey const & key );
	/// @brief Indicate to the StandardJobInputter that a String option can be specified on a per-job basis.
	void add_option( utility::options::StringOptionKey const & key );
	/// @brief Indicate to the StandardJobInputter that a StringVector option can be specified on a per-job basis.
	void add_option( utility::options::StringVectorOptionKey const & key );

	///////////////////////////////////////////////////////////////////////////////////
	// The following functions should be used by the derived JobQueens to communicate
	// to the StandardJobQueen the structure of the XML document that they support.
	// The StandardJobQueen will prepare the output XSD document based on the function
	// calls here, and will also validate the input job-definition file, ensuring that
	// no options are provided to through that file that are not read.

	/// @brief Override the basic specification of input elements as either PDBS, Silent Structures,
	/// or ResourceManager defined structures
	void remove_default_input_element();

	///// @brief Add a sub-element to the element given by the string-separated "nesting level" element (e.g. "" if top level,
	///// or "FastRelax MoveMap" to append beneath the already-declared MoveMap element that itself lives beneath the already-declared
	///// FastRelax element".  This sub-element is an "all" sub-element, which means it is required, but can appear in any order among
	///// the block of "all" sub-elements (which by convention will appear after the "choice" sub-elements).
	//void add_all_sub_element(
	// std::string const & nesting_level,
	// std::string const & element_name );
	//
	///// @brief Add a "choice" sub-element -- one among the set of choice sub-elements that are in the same block.
	///// Multiple "choice blocks" can be put into the same element.  Choice blocks should start indexing at 1
	///// and be incremented in order (e.g. an element for choice block i-1 should have been inserted before an
	///// element for choice block i for all i > 1).
	//void add_choice_sub_element(
	// std::string const & nesting_level,
	// core::Size choice_block_index,
	// std::string const & element_name );
	//
	///// @brief Add an attribute to a given element.
	//void add_attribute(
	// std::string const & nested_element_name,
	// std::string const & attribute_nanme,
	// XMLAttributeType attribute_type );
	//
	///// @brief Add an attribute to a given element.
	//void add_attribute(
	// std::string const & nested_element_name,
	// std::string const & attribute_nanme,
	// XMLAttributeType attribute_type,
	// std::string const & default_value
	// );
	//
	///// @brief Add an attribute to a given element and indicate that it represents
	///// the name of a managed resource and thus that the StandardJobQueen, when it
	///// encounters this attribute, should check with the ResourceManager that a
	///// resource with that name exists and to throw an exception if it does not.
	///// This allows bad inputs to be caught very early.
	//void add_managed_resource_attribute(
	// std::string const & nested_element_name,
	// std::string const & attribute_name );

	// ResourceManagerOP resource_manager();

	/// @brief Access the pose inputter
	PoseInputter & pose_inputter();

	/// @brief Access the pose outputter
	PoseOutputter & pose_outputter();

	/// @brief Return a copy of the Pose to be used with the given job.
	core::pose::PoseOP pose_for_job( LarvalJobCOP job );


private:

	// std::map< std::string, XMLElementOP > job_subtags_by_name_;
	// std::list< XMLElementOP > job_subtags_;

	// ResourceManagerOP resource_manager_;

	PoseInputterOP  pose_inputter_;
	PoseOutputterOP pose_outputter_;
};



} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_StandardJobQueen_HH
