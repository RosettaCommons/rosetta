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

// Package Headers
#include <core/scoring/elec/electrie/ElecTrieEvaluator.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>

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

Real
ElecTrieEvaluator::elec_weight( bool at1isbb, bool at2isbb ) const {
	return ( at1isbb ? ( at2isbb ? wbb_bb_ : wbb_sc_ ) : ( at2isbb ? wsc_bb_ : wsc_sc_ ) ); // fpd add asymmetry
}


ElecTrieEvaluator::ElecTrieEvaluator(
		core::Real wt_bb_bb,
		core::Real wt_bb_sc,
		core::Real wt_sc_bb,
		core::Real wt_sc_sc,
		core::scoring::elec::FA_ElecEnergy const &elec
) :
	wbb_bb_(wt_bb_bb),
	wbb_sc_(wt_bb_sc),
	wsc_bb_(wt_sc_bb),
	wsc_sc_(wt_sc_sc),
	elec_(elec)
{}

ElecTrieEvaluator::~ElecTrieEvaluator() {}

Energy
ElecTrieEvaluator::heavyatom_heavyatom_energy(
	ElecAtom const & at1,
	ElecAtom const & at2,
	DistanceSquared & d2,
	Size & /*path_dist*/
) const
{
	Real bbscwt = elec_weight(at1.isbb(),at2.isbb());
	Real score = bbscwt * elec_.coulomb().eval_atom_atom_fa_elecE( at1.xyz(), at1.charge(), at2.xyz(), at2.charge(), d2  );
	return score;
}

Energy
ElecTrieEvaluator::heavyatom_hydrogenatom_energy(
	ElecAtom const & at1,
	ElecAtom const & at2,
	Size & /*path_dist*/
) const
{
	Real d2;
	Real bbscwt = elec_weight(at1.isbb(),at2.isbb());
	Real score = bbscwt * elec_.coulomb().eval_atom_atom_fa_elecE( at1.xyz(), at1.charge(), at2.xyz(), at2.charge(), d2  );
	return score;
}

Energy
ElecTrieEvaluator::hydrogenatom_heavyatom_energy(
	ElecAtom const & at1,
	ElecAtom const & at2,
	Size & /*path_dist*/
) const
{
	Real d2;
	Real bbscwt = elec_weight(at1.isbb(),at2.isbb());
	Real score = bbscwt * elec_.coulomb().eval_atom_atom_fa_elecE( at1.xyz(), at1.charge(), at2.xyz(), at2.charge(), d2  );

	return score;
}

Energy
ElecTrieEvaluator::hydrogenatom_hydrogenatom_energy(
	ElecAtom const & at1,
	ElecAtom const & at2,
	Size & /*path_dist*/
) const
{
	Real d2;
	Real bbscwt = elec_weight(at1.isbb(),at2.isbb());
	Real score = bbscwt * elec_.coulomb().eval_atom_atom_fa_elecE( at1.xyz(), at1.charge(), at2.xyz(), at2.charge(), d2  );
	return score;
}

core::Real
ElecTrieEvaluator::hydrogen_interaction_cutoff2() const {
	return elec_.hydrogen_interaction_cutoff2();
}

} // namespace lkbtrie
} // namespace lkball
} // namespace scoring
} // namespace core

