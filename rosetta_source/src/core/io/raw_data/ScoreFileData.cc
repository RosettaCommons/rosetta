// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/raw_data/ScoreFileData.cc
///
/// @brief a way to write out just the score using the SilentFile stuff
/// @author Monica Berrondo

// C++ Headers
#include <iostream>
#include <string>
#include <map>

// mini headers
#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>
#include <core/io/raw_data/ScoreStruct.hh>
#include <core/io/raw_data/ScoreFileData.hh>

#include <core/pose/Pose.hh>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace core {
namespace io {
namespace raw_data {

	/// @brief write the given Score_Silent_Struct to the supplied filename.
	bool ScoreFileData::write_struct(
		ScoreStruct s,
		std::map < std::string, core::Real > const & score_map,
        std::map < std::string, std::string > const & string_map
	) {
		bool success = false;

		utility::io::ozstream output;
		if ( !utility::file::file_exists( filename_ ) ) {
			output.open( filename_ );
			s.print_header( output, score_map, string_map );
		} else {
			output.open_append( filename_ );
		}

		s.print_scores      ( output, score_map, string_map );

		return success;
	}


	/// @brief write the given Score_Silent_Struct to the supplied filename.
	bool ScoreFileData::write_pose(
		const core::pose::Pose & pose,
		std::map < std::string, core::Real > const & score_map,
		std::string tag = "empty_tag",
        std::map < std::string, std::string > const & string_map
	) {
		ScoreStruct s( pose, tag );
		bool success = write_struct( s, score_map, string_map );
		return success;
	}

} // namespace silent
} // namespace io
} // namespace core
