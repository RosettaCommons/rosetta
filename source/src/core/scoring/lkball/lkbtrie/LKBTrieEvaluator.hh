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

#ifndef INCLUDED_core_scoring_lkball_lkbtrie_LKBTrieEvaluator_hh
#define INCLUDED_core_scoring_lkball_lkbtrie_LKBTrieEvaluator_hh

// Package Headers
#include <core/scoring/lkball/LK_BallEnergy.fwd.hh>
#include <core/scoring/lkball/lkbtrie/LKBAtom.hh>
#include <core/scoring/etable/Etable.fwd.hh>

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

class LKBTrieEvaluator
{
public:
	LKBTrieEvaluator(
		core::Real wt_lk_ball, core::Real wt_lk_ball_iso, core::Real wt_lk_ball_wtd,
		core::Real wt_lk_ball_bridge, core::Real wt_lk_ball_bridge_uncpl,
		core::scoring::lkball::LK_BallEnergy const &lkb,
		core::scoring::etable::EtableCOP etable
	);

	~LKBTrieEvaluator();


	Energy
	heavyatom_heavyatom_energy(
		LKBAtom const & at1,
		LKBAtom const & at2,
		DistanceSquared & d2,
		Size & /*path_dist*/
	) const;

	inline
	Energy heavyatom_hydrogenatom_energy(
		LKBAtom const & /*at1*/,
		LKBAtom const & /*at2*/,
		Size & /*path_dist*/
	) const
	{
		return 0.0;
	}

	inline
	Energy hydrogenatom_heavyatom_energy(
		LKBAtom const & /*at1*/,
		LKBAtom const & /*at2*/,
		Size & /*path_dist*/
	) const
	{
		return 0.0;
	}

	inline
	Energy hydrogenatom_hydrogenatom_energy(
		LKBAtom const & /*at1*/,
		LKBAtom const & /*at2*/,
		Size & /*path_dist*/
	) const
	{
		return 0.0;
	}

	core::Real
	hydrogen_interaction_cutoff2() const;

private:
	core::Real wt_lk_ball_, wt_lk_ball_iso_, wt_lk_ball_wtd_, wt_lk_ball_bridge_, wt_lk_ball_bridge_uncpl_;
	core::scoring::lkball::LK_BallEnergy const & lkb_; // store reference to energy method (which does the heavy lifting)
	core::scoring::etable::EtableCOP etable_; // pointer to etable
};


} // namespace lkbtrie
} // namespace lkball
} // namespace scoring
} // namespace core

#endif
