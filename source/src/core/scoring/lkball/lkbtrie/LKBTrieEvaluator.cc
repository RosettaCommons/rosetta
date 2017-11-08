// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/etrie/EtableAtom.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Package Headers
#include <core/scoring/lkball/lkbtrie/LKBTrieEvaluator.hh>
#include <core/scoring/lkball/LK_BallEnergy.hh>
#include <core/scoring/etable/Etable.hh>

// Project Headers
#include <core/types.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

namespace core {
namespace scoring {
namespace lkball {
namespace lkbtrie {

LKBTrieEvaluator::LKBTrieEvaluator(
	core::Real wt_lk_ball, core::Real wt_lk_ball_iso, core::Real wt_lk_ball_wtd,
	core::Real wt_lk_ball_bridge, core::Real wt_lk_ball_bridge_uncpl,
	core::scoring::lkball::LK_BallEnergy const &lkb,
	core::scoring::etable::EtableCOP etable
) :
	wt_lk_ball_(wt_lk_ball),
	wt_lk_ball_iso_(wt_lk_ball_iso),
	wt_lk_ball_wtd_(wt_lk_ball_wtd),
	wt_lk_ball_bridge_(wt_lk_ball_bridge),
	wt_lk_ball_bridge_uncpl_(wt_lk_ball_bridge_uncpl),
	lkb_(lkb),
	etable_(etable)
{}

LKBTrieEvaluator::~LKBTrieEvaluator() {}

Energy
LKBTrieEvaluator::heavyatom_heavyatom_energy(
	LKBAtom const & at1,
	LKBAtom const & at2,
	DistanceSquared & d2,
	Size & /*path_dist*/
) const
{
	d2 = at1.atom( ).xyz().distance_squared( at2.atom( ).xyz() );

	core::Real lk_desolvation_of_atom1_by_atom2, lk_desolvation_of_atom2_by_atom1;
	etable_->analytic_lk_energy(
		at1.atom( ), at2.atom( ), lk_desolvation_of_atom1_by_atom2, lk_desolvation_of_atom2_by_atom1 );

	core::Real lk_desolvation_of_atom1_by_atom2_lkb = lk_desolvation_of_atom1_by_atom2 *
		lkb_.get_lk_fractional_contribution( at2.atom().xyz(), at2.atom().type(),
		at1.n_attached_waters(), at1.waters() );
	core::Real lk_desolvation_of_atom2_by_atom1_lkb = lk_desolvation_of_atom2_by_atom1 *
		lkb_.get_lk_fractional_contribution( at1.atom().xyz(), at1.atom().type(),
		at2.n_attached_waters(), at2.waters() );

	core::Real lk_ij = 0.0;
	if ( at1.n_attached_waters() != 0 ) lk_ij += lk_desolvation_of_atom1_by_atom2;
	if ( at2.n_attached_waters() != 0 ) lk_ij += lk_desolvation_of_atom2_by_atom1;
	core::Real lkbr_ij = lk_desolvation_of_atom1_by_atom2_lkb+lk_desolvation_of_atom2_by_atom1_lkb;

	core::Real lkbridge_frac = 0.0, lkbr_wt = 0.0, lkbr_uncpl_wt = 0.0;
	if ( wt_lk_ball_bridge_ != 0.0 || wt_lk_ball_bridge_uncpl_ != 0.0 ) {
		lkbridge_frac = lkb_.get_lkbr_fractional_contribution(
			at1.atom().xyz(), at2.atom().xyz(),
			at1.n_attached_waters(), at2.n_attached_waters(),
			at1.waters(), at2.waters(),
			lkbr_ij,
			lkbr_wt, lkbr_uncpl_wt );
	}

	core::Real score =
		wt_lk_ball_iso_ * lk_ij +
		wt_lk_ball_  * lkbr_ij +
		wt_lk_ball_bridge_  * lk_ij * lkbridge_frac +
		wt_lk_ball_bridge_uncpl_  * lkbridge_frac +
		wt_lk_ball_wtd_ * (
		at1.atom_weights()[1] * lk_desolvation_of_atom1_by_atom2 +
		at1.atom_weights()[2] * lk_desolvation_of_atom1_by_atom2_lkb +
		at2.atom_weights()[1] * lk_desolvation_of_atom2_by_atom1 +
		at2.atom_weights()[2] * lk_desolvation_of_atom2_by_atom1_lkb
	);

	return score;
}

core::Real
LKBTrieEvaluator::hydrogen_interaction_cutoff2() const {
	//return etable_->hydrogen_interaction_cutoff2();

	//fd uncpl makes this a bit different
	core::Real maxHdist = (lkb_.atomic_interaction_cutoff());
	return (maxHdist*maxHdist);
}

} // namespace lkbtrie
} // namespace lkball
} // namespace scoring
} // namespace core

