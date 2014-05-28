// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/raw_data/ScoreFileData.hh
///
/// @brief a way to write out just the score using the SilentFile stuff
/// @author Monica Berrondo

#ifndef INCLUDED_core_io_raw_data_ScoreFileData_hh
#define INCLUDED_core_io_raw_data_ScoreFileData_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/raw_data/Raw.fwd.hh>
#include <core/io/raw_data/RawFileData.hh>

// C++ Headers
// AUTO-REMOVED #include <string>
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace raw_data {
	class ScoreFileData : public RawFileData {
		public:
			///////////////////////////////////////////////////////////////////////////
			// constructor
			ScoreFileData(std::string filename) : filename_(filename) {}

			bool write_struct(
				const RawStructOP s,
				std::map < std::string, core::Real > const & score_map,
                std::map < std::string, std::string > const & string_map = ( std::map < std::string, std::string > () )
			);

			bool write_pose(
				const core::pose::Pose & pose,
				std::map < std::string, core::Real > const & score_map,
                std::string tag,
                std::map < std::string, std::string > const & string_map = ( std::map < std::string, std::string > () )
			);

		private:
			std::string filename_;
	};

} // namespace silent
} // namespace io
} // namespace core

#endif
