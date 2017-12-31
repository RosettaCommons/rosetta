// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/chunk_library_inputters/PDBChunkLibraryInputter.hh
/// @brief  Declaration of the %PDBChunkLibraryInputter class for initializing Poses from .pdb or .pdb.gz files
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_protocols_jd3_chunk_library_inputters_PDBChunkLibraryInputter_HH
#define INCLUDED_protocols_jd3_chunk_library_inputters_PDBChunkLibraryInputter_HH

// Unit headers
#include <protocols/jd3/chunk_library_inputters/PDBChunkLibraryInputter.fwd.hh>

// Package headers
#include <protocols/jd3/chunk_library_inputters/ChunkLibraryInputter.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace jd3 {
namespace chunk_library_inputters {

/// @brief This is the simplest implementation of ChunkLibraryInputter, which reads from -s/-l and PDB files.
class PDBChunkLibraryInputter : public ChunkLibraryInputter
{
public:

	PDBChunkLibraryInputter();
	virtual ~PDBChunkLibraryInputter();

	bool job_available_on_command_line() const override;

	/// @brief Constructs a list of ChunkLibraryInputSource objects reading from the
	/// -s or -l command line flags. This stores the names of the PDBs that
	/// are to be read in, and it initializes the input tags based on the pdb
	/// names, stripping the path and the extension from the file name.
	ChunkLibraryInputSources chunk_library_input_sources_from_command_line() override;

	ChunkLibraryInputSources chunk_library_input_sources_from_tag(
		utility::options::OptionCollection const & opts,
		utility::tag::TagCOP tag ) override;

	/// @brief Takes a ChunkLibraryInputSource object previously initialized in the
	/// call to initialize_chunk_library_input_sources()
	core::pose::PoseOP
	chunk_library_from_input_source(
		ChunkLibraryInputSource const &,
		utility::options::OptionCollection const &,
		utility::tag::TagCOP tag // possibly null-pointing tag pointer
	) override;

	/// @brief returns the name for the element that will be used in a job-definition
	/// file for a structure originating from a .pdb file: "PDB"
	static std::string keyname();

	/// @brief returns the schema for the PDB element used in a job-definition file
	/// including all options that govern how a PDB is loaded.
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static void list_options_read( utility::options::OptionKeyList & read_options );


}; // PDBChunkLibraryInputter

} // namespace chunk_library_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_PDBChunkLibraryInputter_HH
