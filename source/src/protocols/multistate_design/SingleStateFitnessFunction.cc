// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SingleStateFitnessFunction.cc
/// @brief
/// @author Colin A. Smith

#include <protocols/multistate_design/SingleStateFitnessFunction.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace multistate_design {

core::Real
SingleStateFitnessFunction::evaluate(
	core::pose::Pose const & pose
) const
{
	// the base class fitness is just the pose score
	runtime_assert(pose.energies().energies_updated());
	return pose.energies().total_energy();
}

} // namespace multistate_design
} // namespace protocols
