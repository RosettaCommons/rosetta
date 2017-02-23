// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/util.hh
/// @brief Utilities to help in selecting residues.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_select_util_hh
#define INCLUDED_core_select_util_hh

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <core/types.hh>


namespace core {
namespace select {

/// @brief Get a vector1 of the indexes corresponding to True or False in the ResidueSubset.
utility::vector1< core::Size >
get_residues_from_subset( utility::vector1< bool > const & subset, bool select = true);


/// @brief Get a boolean vector of neighbor residues given some distance < 10A
///
utility::vector1< bool >
get_neighbor_residues(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & residue_positions,
	core::Real neighbor_dis
);

/// @brief
///  Fill a boolean vector of residues with neighboring residues.
/// @details
///  Uses 10A neighbor graph for neighbors,
///  THEN takes distance from selection to neighbors to trim neighbors.
///  residue_positions includes positions and neighbors.
///  Use a NOT selector to get only the neighbors.
void
fill_neighbor_residues(
	core::pose::Pose const & pose,
	utility::vector1<bool> & residue_positions,
	core::Real neighbor_dis = 10.0
);



/// @brief filter neighbors to a certain CB distance.
///  Takes distance from selection to neighbors to trim neighbors.
utility::vector1< bool >
trim_neighbors_by_distance(
	core::pose::Pose const & pose,
	utility::vector1<bool> const & selection,
	utility::vector1<bool> const & selection_and_neighbors,
	core::Real & dist_cutoff
);

/// @brief filter neighbors to a certain CB distance.
///  Takes distance from selection to neighbors to trim neighbors.
void
filter_neighbors_by_distance(
	core::pose::Pose const & pose,
	utility::vector1<bool> const & selection,
	utility::vector1<bool> & selection_and_neighbors,
	core::Real & dist_cutoff
);


/// @brief Get neighbor residues within 10 A CB distance cutoff using the 10 A neighbor graph in energies.
///  Returns ONLY the neighbor residues.
utility::vector1< bool >
get_tenA_neighbor_residues(
	core::pose::Pose const & pose,
	utility::vector1<bool> const & residue_positions
);

/// @brief Fill neighbor residues within 10 A CB distance cutoff using the 10 A neighbor graph in energies into the boolean vector.
void
fill_tenA_neighbor_residues(
	pose::Pose const & pose,
	utility::vector1<bool> & residue_positions
);


} //core
} //select


#endif //core/select_util_hh

