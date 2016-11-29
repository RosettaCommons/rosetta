// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/EtableAtom.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_elec_electrie_ElecTrieEvaluator_hh
#define INCLUDED_core_scoring_elec_electrie_ElecTrieEvaluator_hh

// Package Headers
#include <core/scoring/elec/FA_ElecEnergy.fwd.hh>
#include <core/scoring/elec/electrie/ElecAtom.hh>

// Project Headers
#include <core/types.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

namespace core {
namespace scoring {
namespace elec {
namespace electrie {

class ElecTrieEvaluator
{
public:
	ElecTrieEvaluator(
		core::Real wt_bb_bb, core::Real wt_bb_sc, core::Real wt_sc_bb,
		core::Real wt_sc_sc,
		core::scoring::elec::FA_ElecEnergy const &elec
	);

	~ElecTrieEvaluator();

	Real
	elec_weight( bool at1isbb, bool at2isbb ) const;

	Energy
	heavyatom_heavyatom_energy(
		ElecAtom const & at1,
		ElecAtom const & at2,
		DistanceSquared & d2,
		Size & /*path_dist*/
	) const;

	Energy
	heavyatom_hydrogenatom_energy(
		ElecAtom const & at1,
		ElecAtom const & at2,
		Size & /*path_dist*/
	) const;

	Energy
	hydrogenatom_heavyatom_energy(
		ElecAtom const & at1,
		ElecAtom const & at2,
		Size & /*path_dist*/
	) const;

	Energy
	hydrogenatom_hydrogenatom_energy(
		ElecAtom const & at1,
		ElecAtom const & at2,
		Size & /*path_dist*/
	) const;

	core::Real
	hydrogen_interaction_cutoff2() const;

private:
	core::Real wbb_bb_, wbb_sc_, wsc_bb_, wsc_sc_;
	core::scoring::elec::FA_ElecEnergy const & elec_; // store reference to energy method (which does the heavy lifting)
};


} // namespace lkbtrie
} // namespace lkball
} // namespace scoring
} // namespace core

#endif
