// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionKernelGpu.hh
/// @brief  FiberDiffraction kernel module for GPU computations
/// @author Wojciech Potrzebowski and Ingemar Andre

#ifndef INCLUDED_core_scoring_fiber_diffraction_kernel_gpu_hh
#define INCLUDED_core_scoring_fiber_diffraction_kernel_gpu_hh

#include <core/types.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace fiber_diffraction {

void calculate_intensity_gpu(
	Size const l_max,
	Size const natoms,
	utility::vector0< utility::vector0 < int > >::iterator & nvals,
	utility::vector0< utility::vector1< core::Real > >::iterator & layer_lines_R,
	utility::vector0< utility::vector1< core::Real > > & I,
	utility::vector0< utility::vector1< utility::vector1< core::Real > > >::iterator & form_factors,
	utility::vector1< Real > & phi,
	utility::vector1< Real > & z,
	utility::vector1< Real > & r,
	utility::vector1< Size > & atom_type_number,
	Real const c_,
	Real const res_cutoff_low_,
	Real const res_cutoff_high_,
	int const gpu_processor_);

void  calculate_derivatives_gpu(
	Size const l_max,
	Size const natoms,
	utility::vector0< utility::vector0 < int > >::iterator & nvals,
	utility::vector0< utility::vector1< Real > >::iterator & layer_lines_R,
	utility::vector0< utility::vector1< Real > >::iterator & layer_lines_I,
	utility::vector0< utility::vector1< Real > > & I,
	utility::vector0< utility::vector1< utility::vector1< Real > > >::iterator & form_factors,
	utility::vector1< Real > & phi,
	utility::vector1< Real > & z,
	utility::vector1< Real > & r,
	utility::vector1< Size > & atom_type_number,
	utility::vector1< numeric::xyzVector< core::Real > > & dchi2_d,
	utility::vector1< numeric::xyzVector< core::Real > > & dchi2_d_cross_R,
	Real const c_,
	Real const res_cutoff_low_,
	Real const res_cutoff_high_,
	Real const scale_factor_,
	Real const square_obs_,
	int const gpu_processor_,
	bool rfactor_refinement);
}
}
}
#endif
