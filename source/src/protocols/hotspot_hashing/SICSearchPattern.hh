// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_protocols_hotspot_hashing_SICSearchPattern_hh
#define INCLUDED_protocols_hotspot_hashing_SICSearchPattern_hh

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <core/kinematics/Stub.hh>

#include <protocols/sic_dock/SICFast.hh>
#include <protocols/hotspot_hashing/SearchPattern.hh>

namespace protocols {
namespace hotspot_hashing {

class SICPatternAtTransform : public SearchPattern
{
public:
	const static core::Size default_displacement = 100;

	SICPatternAtTransform(
		core::pose::Pose const & source_pose,
		core::pose::Pose const & placed_pose,
		SearchPatternOP slide_pattern,
		SearchPatternOP source_pattern,
		core::Real starting_displacement = default_displacement);

	SICPatternAtTransform(
		core::pose::Pose const & source_pose,
		core::pose::Pose const & placed_pose,
		SearchPatternOP slide_pattern,
		core::Real starting_displacement = default_displacement);

	utility::vector1<core::kinematics::Stub> Searchpoints() override;

private:
	protocols::sic_dock::SICFast sic_fast_;
	core::Real starting_displacement_;
	SearchPatternOP slide_pattern_;
	SearchPatternOP source_pattern_;
};

}
}

#endif
