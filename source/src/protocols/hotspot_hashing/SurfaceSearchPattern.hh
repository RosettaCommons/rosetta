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


#ifndef INCLUDED_protocols_hotspot_hashing_SurfaceSearchPattern_hh
#define INCLUDED_protocols_hotspot_hashing_SurfaceSearchPattern_hh

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <core/kinematics/RT.hh>
#include <protocols/hotspot_hashing/SearchPattern.hh>

namespace protocols {
namespace hotspot_hashing {

class SurfaceSearchPattern : public SearchPattern
{
public:
	SurfaceSearchPattern(
		core::pose::Pose const & source_pose,
		core::pack::task::TaskFactoryOP surface_selection,
		core::Real surface_density);

	utility::vector1<core::kinematics::Stub> Searchpoints() override;
	utility::vector1<VectorPair> surface_vectors() { return surface_vectors_; }

private:
	utility::vector1<VectorPair> surface_vectors_;

	core::Real surface_density_;
};
}
}

#endif
