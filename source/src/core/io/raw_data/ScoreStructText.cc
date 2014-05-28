// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/raw_data/ScoreStructText.cc
///
/// @brief Write out only score information
/// @author Monica Berrondo

// C++ Headers
// AUTO-REMOVED #include <string>
#include <map>

// mini headers
#include <core/pose/Pose.hh>
#include <core/io/raw_data/ScoreStructText.hh>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace raw_data {

	ScoreStructText::ScoreStructText() {
		decoy_tag_   = "empty";
	}

	ScoreStructText::ScoreStructText(
		core::pose::Pose, // pose,
		std::string tag
	) {
		// tag information
		decoy_tag_ = tag;
	}
	// @brief Fill a Pose with the data in this ScoreStructText.
	void ScoreStructText::fill_pose(
		core::pose::Pose &
	) {
//		basic::T("core.io.silent.ScoreStructText") << "Error: don't have a conformation! (called fill_pose())";
	} // fill_pose

	void ScoreStructText::fill_pose(
		core::pose::Pose &,
		core::chemical::ResidueTypeSet const&
	) {
//		basic::T("core.io.silent.ScoreStructText") << "Error: don't have a conformation! (called fill_pose())";
	} // fill_pose

	/// @brief Print the conformation information contained in this object to the given ozstream.
	void ScoreStructText::print_conformation( std::ostream& ) const {
//		basic::T("core.io.silent.ScoreStructText") << "Error: don't have a conformation (called print_conformation())!";
	} // print_conformation

	Real ScoreStructText::get_debug_rmsd() {
//		basic::T("core.io.silent.ScoreStructText") << "Error: don't have a conformation (called get_debug_rmsd())!";
		return 0;
	}

} // namespace silent
} // namespace io
} // namespace core
