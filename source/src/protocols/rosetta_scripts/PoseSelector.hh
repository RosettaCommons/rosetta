// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/rosetta_scripts/PoseSelector.hh
/// @brief  Base class for pose selectors used by MultiplePoseMover
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_rosetta_scripts_PoseSelector_hh
#define INCLUDED_protocols_rosetta_scripts_PoseSelector_hh

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.fwd.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>


namespace protocols {
namespace rosetta_scripts {

enum PoseSelectorFlags {
	PSF_NONE = 0,
	PSF_NEED_FULL_POSE_SET = (1 << 0)
};

class PoseSelector : public utility::pointer::ReferenceCount {

protected:
	PoseSelector();
	virtual ~PoseSelector();

public:
	static std::string name() {
		return "UNDEFINED NAME";
	}
	
	/// @brief Each derived selector must specify its name
	// Why do we have get_name() and name() ?
	virtual std::string get_name() const = 0;	

	/// @brief Report selector flags
	virtual PoseSelectorFlags get_flags() const = 0;

	/// @brief Called by PoseSelectorFactory when constructing new PoseSelectors. Takes care of the specific selector's parsing.
	virtual
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	/// @brief Select poses from provided vector by applying the selection criteria parsed from Tag
	virtual utility::vector1<bool> select_poses( utility::vector1< core::pose::PoseOP > poses ) = 0;

private:
	PosePropertyReporterCOP reporter_;

}; // PoseSelector

} // rosetta_scripts
} // protocols

#endif //INCLUDED_protocols_rosetta_scripts_PoseSelector_HH
