// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/util.cc
/// @brief  Utility functions for beta-barrel construction.
/// @author Andy Watkins

#include <protocols/beta_barrel/util.hh>

#include <numeric/constants.hh>
#include <utility/exit.hh>

namespace protocols {
namespace beta_barrel {

/// @brief Compute per-strand azimuthal position, axial stagger, and inversion flag
/// for each strand in a beta-barrel, given global barrel parameters.
void compute_strand_positions(
	core::Size n_strands,
	core::Size shear_number,
	core::Real z1,
	bool antiparallel,
	utility::vector1< core::Real > & delta_omega0_values,
	utility::vector1< core::Real > & delta_z0_values,
	utility::vector1< bool > & invert_values
) {
	runtime_assert_string_msg( n_strands >= 2,
		"Error in protocols::beta_barrel::compute_strand_positions(): At least two strands are required." );

	core::Real const coiling_angle( numeric::constants::d::pi_2 / static_cast< core::Real >( n_strands ) );

	delta_omega0_values.resize( n_strands );
	delta_z0_values.resize( n_strands );
	invert_values.resize( n_strands );

	for ( core::Size j = 1; j <= n_strands; ++j ) {
		delta_omega0_values[ j ] = static_cast< core::Real >( j - 1 ) * coiling_angle;
		delta_z0_values[ j ] = static_cast< core::Real >( j - 1 ) * static_cast< core::Real >( shear_number ) * z1 / static_cast< core::Real >( n_strands );
		invert_values[ j ] = antiparallel && ( j % 2 == 0 );
	}
}

} //namespace beta_barrel
} //namespace protocols
