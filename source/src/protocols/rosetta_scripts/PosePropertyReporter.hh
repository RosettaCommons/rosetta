// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/rosetta_scripts/PosePropertyReporter.hh
/// @brief  Base class for pose selectors used by MultiplePoseMover
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_rosetta_scripts_PosePropertyReporter_hh
#define INCLUDED_protocols_rosetta_scripts_PosePropertyReporter_hh

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <string>


namespace protocols {
namespace rosetta_scripts {

class PosePropertyReporter : public utility::pointer::ReferenceCount {

protected:
	PosePropertyReporter();
	virtual ~PosePropertyReporter();

public:
	static std::string name() {
		return "UNDEFINED NAME";
	}

	/// @brief Each derived selector must specify its name
	virtual std::string get_name() const = 0;

	/// @brief Reports a specific property for the provided pose
	virtual core::Real report_property( core::pose::Pose & p ) const;

	/// @brief Reports a specific property for the combination of two provided poses
	virtual core::Real report_property( core::pose::Pose & p1, core::pose::Pose & p2 ) const;

	/// @brief Called by PosePropertyReporterFactory when constructing new PosePropertyReporters. Takes care of the specific selector's parsing.
	virtual
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

protected:

private:

}; // PosePropertyReporter

} // rosetta_scripts
} // protocols

#endif //INCLUDED_protocols_rosetta_scripts_PosePropertyReporter_HH
