// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/pose_inputters/PoseInputter.hh
/// @brief  Declaration of the %PoseInputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_inputters_PoseInputter_HH
#define INCLUDED_protocols_jd3_pose_inputters_PoseInputter_HH

//unit headers
#include <protocols/jd3/pose_inputters/PoseInputter.fwd.hh>

// Package headers
#include <protocols/jd3/PoseInputSource.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace jd3 {
namespace pose_inputters {

/// @brief The %PoseInputter is responsible for reading from the command line a set of structures
/// that are each to be run through a protocol, where each input struture will be the starting
/// point for some number of jobs (where that number is at the JobQueen's discretion).  The
/// %PoseInputter is responsible for two things:
/// - for creating a list of PoseInputSource objects
/// - for turning a PoseInputSource object into a Pose on demand.
///
/// @details In addition to the virtual functions, each derived PoseInputter should define two
/// static methods that will be invoked by its corresponding PoseInputterCreator:
/// - static std::string keyname();
/// - static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
class PoseInputter : public utility::pointer::ReferenceCount
{
public:

	PoseInputter();
	virtual ~PoseInputter();

	virtual bool job_available_on_command_line() const = 0;

	/// @brief Construct a list of PoseInputSource objects from the command line that will
	/// be used by the JobQueen to construct a list of LarvalJobs.
	virtual PoseInputSources pose_input_sources_from_command_line() const = 0;

	/// @brief Construct a list of PoseInputSource objects from a Tag object that will
	/// be used by the JobQueen to construct a list of LarvalJobs.
	virtual PoseInputSources pose_input_sources_from_tag( utility::tag::TagCOP ) const = 0;

	/// @brief Convert a single PoseInputSource into a Pose that will be used to
	/// initialize a Job.  The PoseInputSource object must have originated from
	/// this %PoseInputter in the prior call to initialize_pose_input_sources.
	virtual
	core::pose::PoseOP
	pose_from_input_source(
		PoseInputSource const & input_source,
		utility::options::OptionCollection const & options
	) const = 0;


}; // PoseInputter

} // namespace pose_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_pose_inputters_PoseInputter_HH
