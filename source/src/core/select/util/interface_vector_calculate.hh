// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/util/interface_vector_calculate.hh
/// @brief calculates an interface definition based on a short distance
/// @author Ben Stranges (stranges@unc.edu)

#ifndef INCLUDED_core_select_util_interface_vector_calculate_hh
#define INCLUDED_core_select_util_interface_vector_calculate_hh


#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <set>

#include <utility/vector1.hh>

namespace core {
namespace select {
namespace util {

/// @details Calculates the residues at an interface between two protein chains or jump.
/// The calculation is done in the following manner.  First the point graph
/// is used to find all residues within some big cutoff of residues on the other chain.
/// For these residues near the interface, two metrics are used to decide if they are actually
/// possible interface residues.  The first metric is to iterate through all the side chain
/// atoms in the residue of interest and check to see if their distance is less than the nearby
/// atom cutoff, if so then they are an interface residue.  If a residue does not pass that
/// check, then two vectors are drawn, a CA-CB vector and a vector from CB to a CB atom on the
/// neighboring chain.  The dot product between these two vectors is then found and if the angle
/// between them is less than some cutoff then they are classified as interface.

/// @details minimal chain number definition
utility::vector1_bool
calc_interface_vector( core::pose::Pose const & pose, core::Size const chain1_number, core::Size const chain2_number );

/// @details full runner that takes all of the inputs for chains
utility::vector1_bool
calc_interface_vector(
	core::pose::Pose const & pose,
	core::Size const chain1_number, core::Size const chain2_number,
	core::Real const CB_dist_cutoff, core::Real const nearby_atom_cutoff,
	core::Real const vector_angle_cutoff, core::Real const vector_dist_cutoff );

/// @details full runner that takes the jump
utility::vector1_bool
calc_interface_vector(
	core::pose::Pose const & pose,
	int const interface_jump,
	core::Real const CB_dist_cutoff,
	core::Real const nearby_atom_cutoff,
	core::Real const vector_angle_cutoff,
	core::Real const vector_dist_cutoff );

/// @details minimal jump runner
utility::vector1_bool
calc_interface_vector( core::pose::Pose const & pose, int const interface_jump );

/// @brief calculate the same thing but don't require an interface
utility::vector1_bool
calc_interacting_vector(
	core::pose::Pose const & pose,
	std::set< core::Size > const & part1res,
	std::set< core::Size > const & part2res,
	core::Real const CB_dist_cutoff,
	core::Real const nearby_atom_cutoff,
	core::Real const vector_angle_cutoff,
	core::Real const vector_dist_cutoff
);

}//end namespace util
}//end namespace select
}//end namespace core

#endif //INCLUDED_core_select_util_vector_calculate_HH
