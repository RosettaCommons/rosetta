// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/grid_functions.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_ligand_docking_grid_functions_hh
#define INCLUDED_protocols_ligand_docking_grid_functions_hh


#include <iostream>

#include <core/conformation/Residue.fwd.hh>
#include <core/grid/CartGrid.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

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

/// @brief a cleaner implementation of rb_grid_score_atr_rep
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
utility::pointer::shared_ptr<core::grid::CartGrid<int> > make_atr_rep_grid(
	core::pose::Pose const & pose,
	core::Vector const & center
);

/// @brief Make a grid around the specified point with attractive (negative)
/// and repulsive (positive) values for all heavy atoms not in ligand_chain_id_to_exclude
utility::pointer::shared_ptr<core::grid::CartGrid<int> > make_atr_rep_grid_without_ligand(
	core::pose::Pose const & pose,
	core::Vector const & center,
	core::Size const & ligand_chain_id_to_exclude
);

/// @brief Make a grid around the specified point with attractive (negative)
/// and repulsive (positive) values for all heavy atoms not in ligand_chain_ids_to_exclude
utility::pointer::shared_ptr<core::grid::CartGrid<int> > make_atr_rep_grid_without_ligands(
	core::pose::Pose const & pose,
	core::Vector const & center,
	utility::vector1<core::Size> ligand_chain_ids_to_exclude
);

/// @details Just writes the points -- you have to write @ dotlist, etc.
template<typename T>
void grid_to_kin(
	std::ostream & out,
	core::grid::CartGrid<T> const & grid,
	T min_val,
	T max_val,
	int stride /*= 1*/
)
{
	typedef core::grid::CartGrid<int>::GridPt GridPt;
	using core::Vector;

	core::Size nx(0), ny(0), nz(0); // grid points in each dimension
	grid.getNumberOfPoints(nx, ny, nz);
	for ( core::Size i = 0, i_end = nx; i < i_end; i += stride ) {
		for ( core::Size j = 0, j_end = ny; j < j_end; j += stride ) {
			for ( core::Size k = 0, k_end = nz; k < k_end; k += stride ) {
				GridPt grid_pt(i, j, k);
				Vector box_ctr = grid.coords( grid_pt );
				T value = grid.getValue(grid_pt);
				if ( min_val <= value && value <= max_val ) {
					out << '{' << i << ' ' << j << ' ' << k << "}U "
						<< box_ctr.x() << ' ' << box_ctr.y() << ' ' << box_ctr.z() << '\n';
				}
			}
		}
	}
}


} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_grid_functions_HH
