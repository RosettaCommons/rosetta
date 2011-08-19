// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPairAll.cc
/// @brief  Count pair for residues where all atom pairs should be counted.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/atom_pair_energy_inline.hh>

//Auto Headers
//XRW_B_T1
//#include <core/scoring/etable/CoarseEtableEnergy.hh>
//XRW_E_T1
#include <core/scoring/etable/EtableEnergy.hh>


namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

bool
CountPairAll::count(
	int const at1,
	int const at2,
	Real & w,
	Size & path_dist
) const
{
	return operator() ( at1, at2, w, path_dist );
}

void
CountPairAll::residue_atom_pair_energy(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::EtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy( res1, res2, etable_energy, *this, emap );
}


void
CountPairAll::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::EtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_backbone( res1, res2, etable_energy, *this, emap );
}


void
CountPairAll::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::EtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_whole( res1, res2, etable_energy, *this, emap );
}

void
CountPairAll::residue_atom_pair_energy_backbone_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::EtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_backbone_backbone( res1, res2, etable_energy, *this, emap );
}



void
CountPairAll::residue_atom_pair_energy_sidechain_sidechain(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::EtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_sidechain( res1, res2, etable_energy, *this, emap );
}


//XRW_B_T1
/*

void
CountPairAll::residue_atom_pair_energy(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::CoarseEtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy( res1, res2, etable_energy, *this, emap );
}


void
CountPairAll::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::CoarseEtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_backbone( res1, res2, etable_energy, *this, emap );
}


void
CountPairAll::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::CoarseEtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_whole( res1, res2, etable_energy, *this, emap );
}

void
CountPairAll::residue_atom_pair_energy_backbone_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::CoarseEtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_backbone_backbone( res1, res2, etable_energy, *this, emap );
}



void
CountPairAll::residue_atom_pair_energy_sidechain_sidechain(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::CoarseEtableEnergy const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_sidechain( res1, res2, etable_energy, *this, emap );
}

*/
//XRW_E_T1

} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core
