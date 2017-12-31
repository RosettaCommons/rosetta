// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/chunk_library_inputters/ChunkLibraryInputter.hh
/// @brief  Declaration of the %ChunkLibraryInputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_chunk_library_inputters_ChunkLibraryInputter_HH
#define INCLUDED_protocols_jd3_chunk_library_inputters_ChunkLibraryInputter_HH

//unit headers
#include <protocols/jd3/chunk_library_inputters/ChunkLibraryInputter.fwd.hh>

// Package headers
#include <protocols/jd3/chunk_library_inputters/ChunkLibraryInputSource.fwd.hh>

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
namespace chunk_library_inputters {

/// @brief The %ChunkLibraryInputter is responsible for reading from the command line a set of structures
/// that are each to be run through a protocol, where each input struture will be the starting
/// point for some number of jobs (where that number is at the JobQueen's discretion).  The
/// %ChunkLibraryInputter is responsible for two things:
/// - for creating a list of ChunkLibraryInputSource objects
/// - for turning a ChunkLibraryInputSource object into a Pose on demand.
///
/// @details In addition to the virtual functions, each derived ChunkLibraryInputter should define two
/// static methods that will be invoked by its corresponding ChunkLibraryInputterCreator:
/// - static std::string keyname();
/// - static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
class ChunkLibraryInputter : public utility::pointer::ReferenceCount
{
public:

	ChunkLibraryInputter();
	virtual ~ChunkLibraryInputter();

	virtual bool job_available_on_command_line() const = 0;

	/// @brief Construct a list of ChunkLibraryInputSource objects from the command line that will
	/// be used by the JobQueen to construct a list of LarvalJobs.
	virtual ChunkLibraryInputSources chunk_library_input_sources_from_command_line() = 0;

	/// @brief Construct a list of ChunkLibraryInputSource objects from a Tag object that will
	/// be used by the JobQueen to construct a list of LarvalJobs.
	virtual
	ChunkLibraryInputSources
	chunk_library_input_sources_from_tag(
		utility::options::OptionCollection const & opts,
		utility::tag::TagCOP
	) = 0;

	/// @brief Convert a single ChunkLibraryInputSource into a Pose that will be used to
	/// initialize a Job.  The ChunkLibraryInputSource object must have originated from
	/// this %ChunkLibraryInputter in the prior call to initialize_chunk_library_input_sources.
	virtual
	core::pose::PoseOP
	chunk_library_from_input_source(
		ChunkLibraryInputSource const & input_source,
		utility::options::OptionCollection const & options,
		utility::tag::TagCOP tag // possibly null-pointing tag pointer
	) = 0;


}; // ChunkLibraryInputter

} // namespace chunk_library_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_chunk_library_inputters_ChunkLibraryInputter_HH
