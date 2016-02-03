// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/pose_inputters/PDBPoseInputter.cc
/// @brief  %PDBPoseInputter class definition
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)

// Unit headers
#include <protocols/jd3/pose_inputters/PDBPoseInputter.hh>

// Package headers
#include <protocols/jd3/PoseInputSource.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/jd2/util.hh>

//utility headers
#include <utility/assert.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

namespace protocols {
namespace jd3 {
namespace pose_inputters {


PDBPoseInputter::PDBPoseInputter() {}
PDBPoseInputter::~PDBPoseInputter() {}

PoseInputSources
PDBPoseInputter::initialize_pose_input_sources()
{
	utility::vector1< utility::file::FileName > filenames_from_command_line = jd2::input_pdb_files_from_command_line();
	PoseInputSources input_sources;
	input_sources.reserve( filenames_from_command_line.size() );
	for ( core::Size ii = 1; ii <= filenames_from_command_line.size(); ++ii ) {
		utility::file::FileName ii_filename = filenames_from_command_line[ ii ];
		PoseInputSourceOP ii_source( new PoseInputSource( piso_command_line ));
		ii_source->input_tag( ii_filename.base() );
		ii_source->store_string_pair( "filename", ii_filename.name() );
		ii_source->input_kind( pik_pdb_file );
		input_sources.push_back( ii_source );
	}
	return input_sources;
}

/// @details This is only a stub implementation -- the logic here ought to be robust enough
/// to handle non-fullatom poses, but it currently isn't.  This will improve.
core::pose::PoseOP
PDBPoseInputter::pose_from_input_source( PoseInputSource const & input_source )
{
	assert( input_source.string_string_map().find( "filename" ) != input_source.string_string_map().end() );
	return core::import_pose::pose_from_file( input_source.string_string_map().find( "filename" )->second , core::import_pose::PDB_file);
}

// PDBPoseInputterCreator::PDBPoseInputterCreator() {}
// PDBPoseInputterCreator::~PDBPoseInputterCreator() {}

} // namespace pose_inputters
} // namespace jd3
} // namespace protocols

