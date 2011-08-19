// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPair1BC3.cc
/// @brief  Count pair for residue pairs connected with one bond, where the
/// crossover from excluding to counting atom pair interactions is at 3 bonds.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#include <core/scoring/etable/count_pair/CountPair1BC3.hh>
#include <core/scoring/etable/atom_pair_energy_inline.hh>

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {


/// @brief take a row from the path distances table
/// to retrieve the lower and upper path distances
/// for all atoms in each residue
CountPair1BC3::CountPair1BC3(
	conformation::Residue const & res1,
	Size const res1_connect_atom,
	conformation::Residue const & res2,
	Size const res2_connect_atom
) :
	parent(),
	res1_conn_dist_( res1.path_distance( res1_connect_atom )),
	res2_conn_dist_( res2.path_distance( res2_connect_atom ))
{
}

CountPair1BC3::~CountPair1BC3() {}

bool
CountPair1BC3::count(
	int const at1,
	int const at2,
	Real & w,
	Size & path_dist
) const
{
	return operator() ( at1, at2, w, path_dist );
}


void
CountPair1BC3::residue_atom_pair_energy(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::EtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy( res1, res2, etable_energy, *this, emap );
}


void
CountPair1BC3::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::EtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_backbone( res1, res2, etable_energy, *this, emap );
}


void
CountPair1BC3::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::EtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_whole( res1, res2, etable_energy, *this, emap );
}

//XRW_B_T1
/*

void
CountPair1BC3::residue_atom_pair_energy(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::CoarseEtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy( res1, res2, etable_energy, *this, emap );
}


void
CountPair1BC3::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::CoarseEtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_backbone( res1, res2, etable_energy, *this, emap );
}


void
CountPair1BC3::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::CoarseEtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_whole( res1, res2, etable_energy, *this, emap );
}

*/
//XRW_E_T1

} // namespace count_pair
} // namespace etable {
} // namespace scoring
} // namespace core
