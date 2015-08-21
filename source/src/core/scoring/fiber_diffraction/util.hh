// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionEnergyDens.hh
/// @brief  FiberDiffraction utility functions
/// @author Wojciech Potrzebowski Ingemar Andre


#ifndef INCLUDED_core_scoring_fiber_diffraction_util_hh
#define INCLUDED_core_scoring_fiber_diffraction_util_hh

#include <core/scoring/fiber_diffraction/xray_scattering.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/chemical/ChemicalManager.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <string>

namespace core {
namespace scoring {
namespace fiber_diffraction {

void setup_cylindrical_coords(
	pose::Pose const & pose,
	core::Size & natoms,
	utility::vector1< Size > & atom_type_number,
	std::map<  core::id::AtomID, core::Size > & AtomID_to_atomnbr,
	utility::vector1< Real > & phi,
	utility::vector1< Real > & z,
	utility::vector1< Real > & r,
	utility::vector1< Real > & bfactor
);

utility::vector0< utility::vector1< utility::vector1< core::Real > > > setup_form_factors(
	pose::Pose & pose,
	core::Size const & lmax,
	utility::vector0 < utility::vector1< core::Real > >::iterator const & layer_lines_R,
	core::Real const & c,
	core::Real const & B_factor,
	core::Real const & B_factor_solv,
	core::Real const & Ksolv
);

void find_pitch(
	pose::Pose const & pose,
	core::Real & pitch
);

void
find_min_xyz(
	pose::Pose const & pose,
	core::Real &minX,
	core::Real &minY,
	core::Real &minZ,
	core::Real &maxX,
	core::Real &maxY,
	core::Real &maxZ
);

void
find_num_scattering_atoms(
	pose::Pose & pose,
	core::Size & nscatterers
);

void
centroid_scatter(
	std::string const & res_name,
	OneGaussianScattering & sig_centroid
);

utility::vector1< OneGaussianScattering > setup_centroid_scatter(
	pose::Pose & pose
);

bool isPowerOfTwo(int n);


void find_max_r(
	pose::Pose const & pose,
	core::Real & maxR
);

void generate_shannon_points(
	core::Size const & lmax,
	core::Real const & dmax,
	utility::vector0< utility::vector1< core::Real > >::iterator const & layer_lines_R,
	utility::vector0< core::Size > ::iterator & sampling_points_l,
	utility::vector0< utility::vector1< core::Real > >::iterator & shanon_points_lS,
	utility::vector0< core::Size > ::iterator & lowest_bessel_orders_l,
	utility::vector0< utility::vector0 < int > >::iterator const & nvals
);

void bessel_roots(
	core::Size const & lmax,
	core::Real const & c,
	core::Real const & res_cutoff_high,
	core::Real const & res_cutoff_low,
	core::Real & structure_cutoff,
	utility::vector0 < utility::vector1< core::Real > > & bessel_roots_lE,
	utility::vector0 < core::Size >  & sampling_points_l,
	utility::vector0 < core::Size >  & lowest_bessel_orders_l,
	utility::vector0 < core::Size >  & highest_resolution_l,
	utility::vector0 < core::Size >  & lowest_resolution_l,
	utility::vector0 < utility::vector1< core::Real > >::iterator const & layer_lines_R,
	utility::vector0 < utility::vector0 < int > >::iterator const & nvals
);

void interpolate_sampled_to_grid(
	core::Size const & lmax,
	utility::vector0 < utility::vector1< core::Real > > const & bessel_roots_lE,
	utility::vector0 < core::Size > const & sampling_points_l,
	utility::vector0 < core::Size > const & highest_resolution_l,
	utility::vector0 < core::Size > const & lowest_resolution_l,
	utility::vector0 < utility::vector1< core::Real > >::iterator const & layer_lines_R,
	utility::vector0 < utility::vector1< core::Size > > & selected_R_l,
	utility::vector0 < utility::vector1< core::Real > > & selected_Rinv_l
);

void calculate_I_of_E(
	core::Size const & lmax,
	core::Size const & k_iteration,
	utility::vector0 < utility::vector1< core::Real > > const & sampling_points_lE,
	core::Size const & natoms,
	core::Size const & c_,
	utility::vector0< utility::vector0 < int > >::iterator const & nvals,
	utility::vector1< Size > const & atom_type_number,
	utility::vector1< Real > const & phi,
	utility::vector1< Real > const & z,
	utility::vector1< Real > const & r,
	utility::vector0< utility::vector1< utility::vector1< core::Real > > >::iterator const & form_factors,
	utility::vector0< utility::vector1 < Real > > & I_E
);

core::Real calculate_chi_k(
	core::Size const & lmax,
	core::Size const & k_iteration,
	utility::vector0< utility::vector1< core::Size > > const & selected_R_l,
	utility::vector0< utility::vector1< core::Real > > const & selected_Rinv_l,
	utility::vector0< utility::vector1< core::Real > >::iterator const & layer_lines_I,
	utility::vector0< utility::vector1< core::Real > > const & I_interpolated
);

core::Real calculate_chi2_free(
	pose::Pose & pose,
	core::Size const &chi_free_iterations_,
	core::Size const & lmax,
	//utility::vector0< utility::vector1< core::Real > >::iterator const & sampling_points_lE,
	utility::vector0< utility::vector1< core::Real > >::iterator const & layer_lines_I,
	utility::vector0< utility::vector1< core::Real > >::iterator const & layer_lines_R,
	core::Size const & natoms,
	core::Size const & c_,
	utility::vector0< utility::vector0 < int > >::iterator const & nvals,
	utility::vector1< Size > const & atom_type_number,
	utility::vector1< Real > const & phi,
	utility::vector1< Real > const & z,
	utility::vector1< Real > const & r,
	core::Real const b_factor_,
	core::Real const b_factor_solv,
	core::Real const b_factor_solv_K
);

void sample_layer_lines_from_bins(
	core::Size const & lmax,
	core::Size const & k_iteration,
	utility::vector0< utility::vector1< core::Real > >::iterator const & layer_lines_R,
	utility::vector0< utility::vector1< core::Size > > & selected_R_l,
	utility::vector0< utility::vector1< core::Real > > & selected_Rinv_l
);

void rootj(
	int const & N,
	core::Real const & CUTOFF,
	utility::vector1< core::Real >  & zeroj,
	core::Size & npoints
);

void secant(
	int const & N,
	int const & NITMX,
	core::Real TOL,
	core::Real *ZEROJ,
	int *IER
);

}
}
}

#endif
