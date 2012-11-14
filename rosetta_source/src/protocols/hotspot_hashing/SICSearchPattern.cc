// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author 

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>

#include <core/pack/task/TaskFactory.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>

#include <protocols/sic_dock/SICFast.hh>

#include <protocols/hotspot_hashing/SearchPattern.hh>
#include <protocols/hotspot_hashing/SICSearchPattern.hh>

namespace protocols {
namespace hotspot_hashing {

static basic::Tracer TR( "protocols.hotspot_hashing.SICSearchPattern" ); 

SICPatternAtTransform::SICPatternAtTransform(
    core::pose::Pose const & source_pose,
    core::pose::Pose const & placed_pose,
    SearchPatternOP slide_pattern,
    SearchPatternOP source_pattern,
		core::Real starting_displacement) :
  sic_fast_(),
	starting_displacement_(starting_displacement),
	slide_pattern_(slide_pattern),
	source_pattern_(source_pattern)
{
	TR.Debug << "Initializing SICPatternAtTransform."<< std::endl;
  sic_fast_.init(placed_pose, source_pose);
}

SICPatternAtTransform::SICPatternAtTransform(
    core::pose::Pose const & source_pose,
    core::pose::Pose const & placed_pose,
    SearchPatternOP slide_pattern,
		core::Real starting_displacement) :
  sic_fast_(),
	starting_displacement_(starting_displacement),
	slide_pattern_(slide_pattern),
	source_pattern_(new ConstPattern())
{
	TR.Debug << "Initializing SICPatternAtTransform with no source pattern."<< std::endl;
  sic_fast_.init(placed_pose, source_pose);
}

utility::vector1<core::kinematics::Stub> SICPatternAtTransform::Searchpoints()
{
	utility::vector1<core::kinematics::Stub> result_transforms;

	utility::vector1<core::kinematics::Stub> slide_locations = slide_pattern_->Searchpoints();
	utility::vector1<core::kinematics::Stub> source_transforms = source_pattern_->Searchpoints();

	result_transforms.reserve(slide_locations.size() * source_transforms.size());

	for (core::Size i = 1; i <= slide_locations.size(); i++)
	{
		for (core::Size j = 1; j <= source_transforms.size(); j++)
		{
			result_transforms.push_back(
					slideSourceAtLocation(slide_locations[i], source_transforms[j]));
		}
	}

	return result_transforms;
}

core::kinematics::Stub SICPatternAtTransform::slideSourceAtLocation(core::kinematics::Stub slide_location, core::kinematics::Stub source_transform)
{
  using core::kinematics::Stub;
  using core::kinematics::RT;

	// Get transform relating the source transform to the global coord frame,
	// then apply that rt at the slide location frame to produce new stub
	Stub source_at_slide_location;
	RT(core::kinematics::default_stub, source_transform).make_jump(slide_location, source_at_slide_location);

	Vector slide_vector = slide_location.M * Stub::Vector(1,0,0);

	source_at_slide_location.v -= slide_vector * starting_displacement_;

	core::Real sic_distance = sic_fast_.slide_into_contact(
      source_at_slide_location,
      core::kinematics::default_stub,
			-slide_vector);

	source_at_slide_location.v += -slide_vector * sic_distance;

	return source_at_slide_location;
}

}
}
