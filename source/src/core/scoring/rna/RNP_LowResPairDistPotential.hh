// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RNP_LowResPairDistPotential.hh
/// @brief  Statistically derived RNP low resolution protential, simple distance-based
/// @author Kalli Kappel

#ifndef INCLUDED_core_scoring_RNP_LowResPairDistPotential_hh
#define INCLUDED_core_scoring_RNP_LowResPairDistPotential_hh

#include <core/types.hh>

// Unit headers
#include <core/scoring/rna/RNP_LowResPairDistPotential.fwd.hh>

// Package headers
#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers

#include <ObjexxFCL/FArray3D.hh>

#include <utility/vector1.hh>


// C++


namespace core {
namespace scoring {
namespace rna {

////////////////////////////////////////////////////////////////////////////////////////////////////
class RNP_LowResPairDistPotential : public utility::pointer::ReferenceCount {

public:
	RNP_LowResPairDistPotential();

	void
	initialize_rnp_pair_dist();

	void
	evaluate_rnp_pair_dist_score(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Real & rnp_pair_dist_score
	) const;


	Real
	calc_rnp_pair_dist_score(
		conformation::Residue const & protein_rsd,
		conformation::Residue const & RNA_rsd ) const;
	Size
	convert_aa_to_index( char const c ) const;

private: // data

	Size const max_aa_ = 25;
	Size const max_base_ = 8;
	Size const num_dbins_ = 10;
	bool use_actual_centroid_;
	utility::vector1< std::string > RNA_atoms_;
	utility::vector1< std::string > protein_atoms_;

	ObjexxFCL::FArray3D< Real > rnp_pair_dist_potential_;

};

} // rna
} // scoring
} // core

#endif
