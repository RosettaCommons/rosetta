// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/magnesium/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_magnesium_util_HH
#define INCLUDED_core_scoring_magnesium_util_HH

#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh> // silly, for GaussianParameter
#include <core/types.hh>

namespace core {
namespace scoring {
namespace magnesium {

Real
get_cos_theta( core::conformation::Residue const & rsd1,
	Size const i, Vector const & j_xyz,
	Size const i_base = 0 );

Real
get_cos_theta( core::conformation::Residue const & rsd1,
	Size const i,  Vector const & j_xyz,
	Size & i_base,
	Vector & xyz_base );

Real
get_gaussian_potential_score(
	core::chemical::rna::GaussianParameter const & mg_potential_gaussian_parameter,
	Vector const & pos1,
	Vector const & pos2 );

Real
get_gaussian_score(
	core::chemical::rna::GaussianParameter const & mg_potential_gaussian_parameter,
	Real const & d );

Real
get_gaussian_deriv(
	core::chemical::rna::GaussianParameter const & mg_potential_gaussian_parameter,
	Real const & d );

Real
get_cos_angle_to_closest_orbital_axis( core::conformation::Residue const & rsd1, Vector const & j_xyz );

Size
get_closest_orbital_axis( core::conformation::Residue const & mg_rsd, Vector const & j_xyz );

} //magnesium
} //scoring
} //core

#endif
