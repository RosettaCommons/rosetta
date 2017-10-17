// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/SilentFilePoseOutputter.hh
/// @brief  Definition of the %SilentFilePoseOutputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Andy Watkins (amw579@stanford.edu)


#ifndef INCLUDED_protocols_jd3_pose_outputters_SilentFilePoseOutputter_HH
#define INCLUDED_protocols_jd3_pose_outputters_SilentFilePoseOutputter_HH

//unit headers
#include <protocols/jd3/pose_outputters/SilentFilePoseOutputter.fwd.hh>

//package headers
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentFileOptions.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/file/FileName.fwd.hh>

//project headers

namespace protocols {
namespace jd3 {
namespace pose_outputters {

/// @brief The %SilentFilePoseOutputter
class SilentFilePoseOutputter : public PoseOutputter
{
public:

	SilentFilePoseOutputter();
	virtual ~SilentFilePoseOutputter();

	static
	bool
	outputter_specified_by_command_line();

	void
	determine_job_tag(
		utility::tag::TagCOP output_tag,
		utility::options::OptionCollection const & job_options,
		InnerLarvalJob & job
	) const override;

	std::string
	outputter_for_job(
		utility::tag::TagCOP output_tag,
		utility::options::OptionCollection const & job_options,
		InnerLarvalJob const & job
	) const override;

	std::string
	outputter_for_job(
		PoseOutputSpecification const & spec
	) const override;

	bool job_has_already_completed( LarvalJob const & job, utility::options::OptionCollection const & options ) const override;

	void mark_job_as_having_started( LarvalJob const & job, utility::options::OptionCollection const & options ) const override;

	/// @brief Create the PoseOutputSpecification for a particular job
	PoseOutputSpecificationOP
	create_output_specification(
		LarvalJob const & job,
		JobOutputIndex const & output_index,
		utility::options::OptionCollection const & options,
		utility::tag::TagCOP tag // possibly null-pointing tag pointer
	) override;

	/// @brief Write a pose out to permanent storage (whatever that may be).
	void write_output(
		output::OutputSpecification const & spec,
		JobResult const & result
	) override;

	/// @brief Write out all cached silent files -- may happen after an exception is caught by
	/// the JobDistributor and the process is about to be spun down.
	void flush() override;

	/// @brief Return the string used by the SilentFilePoseOutputterCreator to identify the
	/// class in an XML file (which is the string returned by the static keyname method)
	std::string
	class_key() const override;

	std::string
	output_silent_name( LarvalJob const & job ) const;

	static
	std::string
	keyname();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static
	void
	list_options_read( utility::options::OptionKeyList & read_options );

private:
	std::string fname_out_;
	core::Size buffer_limit_;
	core::io::silent::SilentFileOptionsOP opts_;
	utility::vector1< core::io::silent::SilentStructOP > buffered_structs_;

};

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_SilentFilePoseOutputter_HH
