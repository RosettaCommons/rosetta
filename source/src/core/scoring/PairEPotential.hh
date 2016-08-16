// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/PairEPotential.fwd.hh
/// @brief  pairE knowledge-based potential class delcaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_PairEPotential_hh
#define INCLUDED_core_scoring_PairEPotential_hh

#include <core/types.hh>
#include <core/scoring/types.hh>

//ObjexxFCL
#include <ObjexxFCL/FArray5D.hh>

#include <core/conformation/Residue.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace scoring {

class PairEPotential : public utility::pointer::ReferenceCount
{
public:
	PairEPotential();

	Energy
	pair_term_energy(
		conformation::Residue const & res1,
		int res1_num_10A_neighbors,
		conformation::Residue const & res2,
		int res2_num_10A_neighbors
	) const;

	Energy
	pair_term_energy_and_deriv(
		conformation::Residue const & res1,
		int res1_num_10A_neighbors,
		conformation::Residue const & res2,
		int res2_num_10A_neighbors,
		EnergyDerivative & dpairE_dr
	) const;

	bool
	pair_term_energy_exists( conformation::Residue const & rsd ) const;

	Real
	range() const {
		return 5.5; // TEMP LIE
		return (max_bin_+1)*pair_score_bin_range_ + pair_score_bin_base_; // 9.0 A if max_bin_ == 3; 7.5 if max_bin == 2
	}

private:

	Energy
	pair_term_energy(
		conformation::Residue const & res1,
		int res1_num_10A_neighbors,
		conformation::Residue const & res2,
		int res2_num_10A_neighbors,
		Probability & pair_lhood_ratio,
		Probability & pair_lhood_ratio_high,
		Probability & pair_lhood_ratio_low ) const;

	FArray5D_TableProbability  pair_corr_;
	Size pair_score_min_sep_;
	int pair_score_cb_thresh_;
	Distance pair_score_bin_range_;
	Distance pair_score_bin_base_;
	int max_bin_;

};

}
}

#endif
