// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/sasa.hh
/// @brief  routines which calculate solvent accessible surface area
/// @author Jeff Gray

#ifndef INCLUDED_core_scoring_sasa_hh
#define INCLUDED_core_scoring_sasa_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>

// Utility headers

// ObjexxFCL header

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>

// C++

namespace core {
namespace scoring {


// DO not use these functions/this file.  They are deprecated and will be removed soon.
// Please use scoring/sasa/util for functions; and scoring/sasa/SasaCalc class for sasa calculation.
// JAB 2/18/14


void input_sasa_dats();

void get_overlap( Real const radius_a, Real const radius_b, Real const distance_ijxyz, int & degree_of_overlap );
void get_orientation( Vector const & a_xyz, Vector const & b_xyz, int & phi_index, int & theta_index, Real distance_ijxyz );

/// @brief Return total SASA
Real calc_per_atom_sasa( pose::Pose const & pose, id::AtomID_Map< Real > & atom_sasa, utility::vector1< Real > & rsd_sasa,
	Real const probe_radius, bool const use_big_polar_H = false );

//@brief Return total SASA for side chains
Real calc_per_atom_sasa_sc( pose::Pose const & pose, utility::vector1< Real > & rsd_sasa, bool normalize);

/// @brief Get the area of the sidechain.
/// @details Threadsafe now, but these values are suspect.
Real normalizing_area(char const res);

/// @brief Given a one-letter code for a canonical amino acid, return
/// its total surface area.
/// @details Threadsafe now, but these values are suspect.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
Real normalizing_area_total (char const res);

/// @brief Given a one-letter code for a canonical amino acid, return
/// its total surface area, computed only using hydrophobic atoms.
/// @details Threadsafe now.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
Real normalizing_area_total_hydrophobic_atoms_only ( char const res );

/// @brief Given a one-letter code for a canonical amino acid, return
/// its total surface area, computed only using polar atoms.
/// @details Threadsafe.  Based on Gabe Rocklin's values (grocklin@gmail.com).
/// @author Vikram K. Mulligan (vmullig@uw.edu).
Real normalizing_area_total_polar_atoms_only ( char const res );

/// returns total sasa
//Real
//calc_per_atom_sasa( pose::Pose const & pose, id::AtomID_Map< Real > & atom_sasa, utility::vector1< Real > & rsd_sasa,
// Real const probe_radius, bool const use_big_polar_H, id::AtomID_Map< bool > & atom_subset );

Real
calc_per_atom_sasa(
	pose::Pose const & pose,
	id::AtomID_Map< Real > & atom_sasa,
	utility::vector1< Real > & rsd_sasa,
	Real const probe_radius,
	bool const use_big_polar_H,
	id::AtomID_Map< bool > & atom_subset,
	bool const use_naccess_sasa_radii = false,
	bool const expand_polar_radii = false,
	Real const polar_expansion_radius = 1.0,
	bool const include_probe_radius_in_atom_radii = true,
	bool const use_lj_radii = false
);

void
calc_atom_masks(
	core::conformation::Residue const & irsd,
	core::conformation::Residue const & jrsd,
	Real const probe_radius,
	Real const cutoff_distance,
	utility::vector1< Real > const & radii,
	id::AtomID_Map< bool > const & atom_subset,
	core::id::AtomID_Map< utility::vector1< ObjexxFCL::ubyte > > & atom_mask
);


/// returns total sasa
Real calc_total_sasa( pose::Pose const & pose, Real const probe_radius );

int get_num_bytes();
ObjexxFCL::FArray2D_int const & get_angles();
ObjexxFCL::FArray2D_ubyte const & get_masks();

Real calc_per_res_hydrophobic_sasa( pose::Pose const & pose,
	utility::vector1< Real > & rsd_sasa,
	utility::vector1< Real > & rsd_hydrophobic_sasa,
	Real const probe_radius,
	bool use_naccess_sasa_radii = false);

// Undefined, commenting out to fix PyRosetta build  void print_dot_bit_string( utility::vector1< ObjexxFCL::ubyte > & values );

} // namespace scoring
} // namespace core


#endif //
