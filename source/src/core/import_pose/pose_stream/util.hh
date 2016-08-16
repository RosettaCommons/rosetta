// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/pose_stream/util.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_core_import_pose_pose_stream_util_hh
#define INCLUDED_core_import_pose_pose_stream_util_hh

#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <utility/vector1.hh>


namespace core {
namespace import_pose {
namespace pose_stream {

/// @brief Get all input streams based on command-line input.
/// @details If renumber_decoys is true, silent file decoys are sorted in alphabetical order of tags.
MetaPoseInputStream streams_from_cmd_line( bool const do_renumber_decoys);

/// @brief Get all input streams based on command-line input, sorting silent file decoys in alphabetical order of tags.
///
MetaPoseInputStream streams_from_cmd_line();

} // pose_stream
} // import_pose
} // core

#endif
