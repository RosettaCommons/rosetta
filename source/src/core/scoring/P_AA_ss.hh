// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/P_AA.hh
/// @brief  Amino acid probability arrays and functions
/// @author FD

#ifndef INCLUDED_core_scoring_P_AA_ss_hh
#define INCLUDED_core_scoring_P_AA_ss_hh

// Unit Headers
#include <core/scoring/P_AA.fwd.hh>

// Package headers
#include <core/scoring/types.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/TorsionID.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/interpolation/spline/Bicubic_spline.hh>

namespace core {
namespace scoring {

class P_AA_ss : public utility::pointer::ReferenceCount
{
private:
	utility::vector1< core::Real > p_L_, p_H_, p_E_;
	core::Real p0_L_, p0_H_, p0_E_;

public:
	P_AA_ss();

	/// @brief Read the ss-dep amino acid probability file into P_AA_ss
	void
	read_P_AA_ss();

public:
	/// @brief Probability energies from P(aa|ss)
	core::Real
	P_AA_ss_energy( chemical::AA aa, char ss ) const;
};

} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_P_AA_HH
