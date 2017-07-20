// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/full_model_inputters/FullModelInputter.hh
/// @brief  Declaration of the %FullModelInputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_full_model_inputters_FullModelInputter_HH
#define INCLUDED_protocols_jd3_full_model_inputters_FullModelInputter_HH

//unit headers
#include <protocols/jd3/full_model_inputters/FullModelInputter.fwd.hh>

// Package headers
#include <protocols/jd3/full_model_inputters/FullModelInputSource.fwd.hh>

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
namespace full_model_inputters {

/// @brief The %FullModelInputter is responsible for reading from the command line a set of structures
/// that are each to be run through a protocol, where each input struture will be the starting
/// point for some number of jobs (where that number is at the JobQueen's discretion).  The
/// %FullModelInputter is responsible for two things:
/// - for creating a list of FullModelInputSource objects
/// - for turning a FullModelInputSource object into a Pose on demand.
///
/// @details In addition to the virtual functions, each derived FullModelInputter should define two
/// static methods that will be invoked by its corresponding FullModelInputterCreator:
/// - static std::string keyname();
/// - static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
class FullModelInputter : public utility::pointer::ReferenceCount
{
public:

	FullModelInputter();
	virtual ~FullModelInputter();

	virtual bool job_available_on_command_line() const = 0;

	/// @brief Construct a list of FullModelInputSource objects from the command line that will
	/// be used by the JobQueen to construct a list of LarvalJobs.
	virtual FullModelInputSources full_model_input_sources_from_command_line() = 0;

	/// @brief Construct a list of FullModelInputSource objects from a Tag object that will
	/// be used by the JobQueen to construct a list of LarvalJobs.
	virtual
	FullModelInputSources
	full_model_input_sources_from_tag(
		utility::options::OptionCollection const & opts,
		utility::tag::TagCOP
	) = 0;

	/// @brief Convert a single FullModelInputSource into a Pose that will be used to
	/// initialize a Job.  The FullModelInputSource object must have originated from
	/// this %FullModelInputter in the prior call to initialize_full_model_input_sources.
	virtual
	core::pose::PoseOP
	full_model_from_input_source(
		FullModelInputSource const & input_source,
		utility::options::OptionCollection const & options,
		utility::tag::TagCOP tag // possibly null-pointing tag pointer
	) = 0;


}; // FullModelInputter

} // namespace full_model_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_full_model_inputters_FullModelInputter_HH
