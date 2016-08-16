// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ContactMapEvaluator.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_protocols_simple_filters_ContactMapEvaluator_hh
#define INCLUDED_protocols_simple_filters_ContactMapEvaluator_hh

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <algorithm>

#include <utility/vector1.hh>

namespace protocols  {
namespace simple_filters {

class ContactMapEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {

public:
	ContactMapEvaluator(
		core::pose::Pose const & native_pose,
		core::Real const max_dist,
		core::Size const min_seqsep
	);

	~ContactMapEvaluator() {}

	virtual void apply(
		core::pose::Pose & pose,
		std::string tag,
		core::io::silent::SilentStruct & ss
	) const;

	virtual core::Real apply(
		core::pose::Pose & /*pose*/
	) const {
		return 0;
	}

private:
	core::Real const max_dist_;
	core::Size const min_seqsep_;
	core::pose::Pose native_;
};

} // simple_filter
} // protocols

#endif
