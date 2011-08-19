// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPairIntraResC4.cc
/// @brief  Count pair for intra-residue LJ/LK, where the
/// crossover from excluding to counting atom pair interactions is at 4 bonds.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#include <core/scoring/etable/count_pair/CountPairIntraResC4.hh>
#include <core/scoring/etable/atom_pair_energy_inline.hh>

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {


/// @brief take a row from the path distances table
/// to retrieve the lower and upper path distances
/// for all atoms in each residue
CountPairIntraResC4::CountPairIntraResC4(
	conformation::Residue const & res
) :
	parent(),
	path_dists_( res.path_distances())
{
}

CountPairIntraResC4::~CountPairIntraResC4() {}

bool
CountPairIntraResC4::count(
	int const at1,
	int const at2,
	Real & w,
	Size & path_dist
) const
{
	return operator() ( at1, at2, w, path_dist );
}


void
CountPairIntraResC4::residue_atom_pair_energy(
	conformation::Residue const & res,
	conformation::Residue const & /*res2*/,
	etable::EtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_intraresidue_atom_pair_energy( res, etable_energy, *this, emap );
}


void
CountPairIntraResC4::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & /*res1*/,
	conformation::Residue const & /*res2*/,
	etable::EtableEnergy const & /*etable_energy*/,
	EnergyMap & /*emap*/
) const
{
	std::cerr << "Error: illegal call to CountPairIntraResC4::residue_atom_pair_energy_sidechain_backbone" << std::endl;
	utility_exit();
}


void
CountPairIntraResC4::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & ,
	conformation::Residue const & ,
	etable::EtableEnergy const & ,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraResC4::residue_atom_pair_energy_sidechain_whole" << std::endl;
	utility_exit();
}

//XRW_B_T1
/*

void
CountPairIntraResC4::residue_atom_pair_energy(
	conformation::Residue const & res,
	conformation::Residue const & res2,
	etable::CoarseEtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_intraresidue_atom_pair_energy( res, etable_energy, *this, emap );
}


void
CountPairIntraResC4::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::CoarseEtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	std::cerr << "Error: illegal call to CountPairIntraResC4::residue_atom_pair_energy_sidechain_backbone" << std::endl;
	utility_exit();
}


void
CountPairIntraResC4::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & ,
	conformation::Residue const & ,
	etable::CoarseEtableEnergy const & ,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraResC4::residue_atom_pair_energy_sidechain_whole" << std::endl;
	utility_exit();
}

*/
//XRW_E_T1

} // namespace count_pair
} // namespace etable {
} // namespace scoring
} // namespace core
