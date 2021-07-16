// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file source/src/protocols/electron_density/DockIntoDensityUtils.hh
/// @brief Utils for folding into density
/// @details
/// @author Frank DiMaio
/// @author Danny Farrell


#ifndef INCLUDED_protocols_electron_density_DockIntoDensityUtils_hh
#define INCLUDED_protocols_electron_density_DockIntoDensityUtils_hh

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray3D.hh>


namespace protocols {
namespace electron_density {


struct PoseSphericalSamplesOptions {
	core::Size B_, nRsteps_;
	core::Real delRsteps_, laplacian_offset_;
	bool center_on_middle_ca_;
};

/// @brief  step 0: map pose to spherically sampled density + mask
void
pose_spherical_samples(
	core::pose::Pose const &pose,
	ObjexxFCL::FArray3D< core::Real > &sigR,
	ObjexxFCL::FArray3D< core::Real > &epsR,
	PoseSphericalSamplesOptions const & params);

}  // electron_density
}  // protocols
#endif
