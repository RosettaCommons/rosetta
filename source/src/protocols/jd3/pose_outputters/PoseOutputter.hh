// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_outputters/PoseOutputter.hh
/// @brief  Definition of the %PoseOutputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_outputters_PoseOutputter_HH
#define INCLUDED_protocols_jd3_pose_outputters_PoseOutputter_HH

//unit headers
#include <protocols/jd3/pose_outputters/PoseOutputter.fwd.hh>

//package headers
#include <protocols/jd3/LarvalJob.fwd.hh>
#include <protocols/jd3/InnerLarvalJob.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

// utility headers
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace jd3 {
namespace pose_outputters {

/// @brief The %PoseOutputter
class PoseOutputter : utility::pointer::ReferenceCount
{
public:

	PoseOutputter();
	virtual ~PoseOutputter();

	/// @brief Determine the inner-larval job's "job_tag" from the <Output> tag / per-job options
	virtual
	void determine_job_tag(
		utility::tag::TagCOP output_tag,
		utility::options::OptionCollection const & job_options,
		InnerLarvalJob & job
	) const = 0;

	/// @brief Return an identifier string for the specific instance of the %PoseOutputter that ought to be used
	/// for a particular job so that the %PoseOutputter can e.g. aggregate all of the outputs for a group of jobs
	/// and output them all at once when flush is called.  The outputter may return the empty string if all
	/// outputters (of the same type) are interchangable (e.g. the PDBPoseOutputter).
	/// e.g., the SilentFilePoseOutputter returns the name of the file that it sends its outputs to.
	virtual
	std::string
	outputter_for_job(
		utility::tag::TagCOP output_tag,
		utility::options::OptionCollection const & job_options,
		InnerLarvalJob const & job
	) const = 0;

	virtual
	bool job_has_already_completed( LarvalJob const & job, utility::options::OptionCollection const & options ) const = 0;

	virtual
	void mark_job_as_having_started( LarvalJob const & job, utility::options::OptionCollection const & options ) const = 0;

	virtual
	void write_output_pose(
		LarvalJob const & job,
		utility::options::OptionCollection const & options,
		utility::tag::TagCOP tag, // possibly null-pointing tag pointer
		core::pose::Pose const & pose
	) = 0;

	/// @brief Output from a pose outputter may be held back and only flushed when requested
	/// by the JobQueen; I/O can be expensive, so it's a good idea to gather up the
	/// results of many outputs before flushing them to disk.
	virtual
	void flush() = 0;

	/// @brief Return the stiring used by the PoseOutputterCreator for this class
	virtual
	std::string
	class_key() const = 0;

};

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_pose_outputters_PoseOutputter_HH
