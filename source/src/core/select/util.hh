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

#include <core/select/residue_selector/ResidueIndexSelector.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/types.hh>
#include <set>

namespace core {
namespace select {

/// @brief Get a vector1 of the indexes which are the same value as 'select' in the ResidueSubset.
utility::vector1< core::Size >
get_residues_from_subset( utility::vector1< bool > const & subset, bool select = true);

/// @brief Get a std::set of the indexes which are the same value as 'select' in the ResidueSubset.
std::set< core::Size >
get_residue_set_from_subset( utility::vector1< bool > const & subset, bool select = true);

/// @brief Get a vector1 of true/false (a ResidueSubset) from a list of residue numbers.
/// @details If invert is true, residues *not* in the selection are true in the returned vector.
utility::vector1< bool >
get_subset_from_residues( utility::vector1< core::Size > const & selection, core::Size total_nres, bool invert = false);

/// @brief Convert a residue subset back into a residue selector.
/// @detail Not only does this allow one to go from a residue subset back to a selector,
///          but by applying a selector to a pose and using this to go back to a selector,
///          one can turn context-dependent residue selectors into context-independent residue selectors
core::select::residue_selector::ResidueIndexSelectorOP
get_residue_selector_from_subset(
	core::select::residue_selector::ResidueSubset subset );

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

///@brief Return a pymol selection of a set of atoms.
std::string
get_pymol_selection_for_atoms(pose::Pose const & pose, utility::vector1< id::AtomID > const & atoms,std::string const & sele_name, bool skip_virts=true );

///@brief Turns off all residues that are not part of the master subunit.
utility::vector1< bool >
get_master_subunit_selection(pose::Pose const & pose, utility::vector1<bool> const & full_subset);

} //core
} //select


#endif //core/select_util_hh

