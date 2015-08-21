// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/DecoyFileData.hh
///
/// @brief Output a decoy in "silent" format
/// @author James Thompson, Monica Berrondo

#ifndef INCLUDED_core_io_raw_data_DecoyFileData_hh
#define INCLUDED_core_io_raw_data_DecoyFileData_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/io/raw_data/RawFileData.hh>

// C++ Headers
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace raw_data {
class DecoyFileData : public RawFileData {
public:
	///////////////////////////////////////////////////////////////////////////
	// constructor
	DecoyFileData(std::string filename) : filename_(filename) {}

	bool write_struct(
		const DecoyStruct s,
		std::map < std::string, core::Real > const & score_map
	);

	bool write_pose(
		const core::pose::Pose & pose,
		std::map < std::string, core::Real > const & score_map,
		std::string tag,
		bool fullatom
	);

private:
	std::string filename_;
};

} // namespace silent
} // namespace io
} // namespace core

#endif
