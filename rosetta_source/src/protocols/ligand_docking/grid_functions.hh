// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/grid_functions.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_ligand_docking_grid_functions_hh
#define INCLUDED_protocols_ligand_docking_grid_functions_hh

// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/grid/CartGrid.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

#include <iostream>

//Auto Headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/grid/CartGrid.hh>

namespace protocols {
namespace ligand_docking {


/// @brief Sum the grid values for all heavy atoms in the residue
int grid_score(
	core::grid::CartGrid<int> const & grid,
	core::conformation::Residue const & rsd,
	int max_score = 9999
);


/// @brief Sum the grid values for all heavy atoms in the residue
void grid_score_atr_rep(
	core::grid::CartGrid<int> const & grid,
	core::conformation::Residue const & rsd,
	int & atr_out, //< sum of negative grid values
	int & rep_out, //< sum of positive grid values
	int max_rep = 9999
);

/// @brief Sum the grid values for all heavy atoms in the residue

void rb_grid_score_atr_rep(
	core::grid::CartGrid<int> const & grid,
	core::pose::Pose const & pose,
	core::Size begin,
	core::Size const end,
	int & atr_out, //< sum of negative grid values
	int & rep_out, //< sum of positive grid values
	int max_rep = 9999
);

///@brief a cleaner implementation of rb_grid_score_atr_rep
std::pair<int, int> get_rb_atr_and_rep_scores(
		core::grid::CartGrid<int> const & grid,
		core::pose::Pose const & pose,
		core::Size begin,
		core::Size end
);

/// @brief Try all rotamers for the specified residue and install the first one
/// that minimizes the grid score.  Only tested with ligand residues w/ a conformer library.
void grid_rotamer_trials(
	core::grid::CartGrid<int> const & grid,
	core::pose::Pose & pose,
	core::Size rsd_no,
	int const min_score = 0
);

void rb_grid_rotamer_trials_atr_rep(
		core::grid::CartGrid<int> const & grid,
		core::pose::Pose & pose,
		core::Size begin,
		core::Size end
);

/// @brief Try all rotamers for the specified residue and install the first one
/// that minimizes the repulsive score, breaking ties by the attractive score.
/// Only tested with ligand residues w/ a conformer library.
void grid_rotamer_trials_atr_rep(
	core::grid::CartGrid<int> const & grid,
	core::pose::Pose & pose,
	core::Size rsd_no
);

/// @brief Internal helper function for rotamer trials;  fills conformers_out.
void rotamers_for_trials(
	core::pose::Pose & pose,
	core::Size rsd_no,
	utility::vector1< core::conformation::ResidueOP > & conformers_out
);


/// @brief Set the value for all grid boxes whose centers fall inside the sphere.
/* void set_sphere(
	core::grid::CartGrid<int> const & grid,
	core::Vector const & center,
	core::Real radius,
	int value
); */


/// @brief Make a grid around the specified point with attractive (negative)
/// and repulsive (positive) values for the protein backbone.
utility::pointer::owning_ptr<core::grid::CartGrid<int> > make_atr_rep_grid(
	core::pose::Pose const & pose,
	core::Vector const & center
);

/// @brief Make a grid around the specified point with attractive (negative)
/// and repulsive (positive) values for all heavy atoms not in ligand_chain_id_to_exclude
/* utility::pointer::owning_ptr<core::grid::CartGrid<int> > make_atr_rep_grid_with_ligands(
	core::pose::Pose const & pose,
	core::Vector const & center,
	core::Size const & ligand_chain_id_to_exclude
); */

/// @brief Make a grid around the specified point with attractive (negative)
/// and repulsive (positive) values for all heavy atoms not in ligand_chain_id_to_exclude
utility::pointer::owning_ptr<core::grid::CartGrid<int> > make_atr_rep_grid_without_ligands(
	core::pose::Pose const & pose,
	core::Vector const & center,
	utility::vector1<core::Size> ligand_chain_ids_to_exclude
);

void grid_to_kin(
	std::ostream & out,
	core::grid::CartGrid<int> const & grid,
	int min_val,
	int max_val,
	int stride = 1
);


} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_grid_functions_HH
