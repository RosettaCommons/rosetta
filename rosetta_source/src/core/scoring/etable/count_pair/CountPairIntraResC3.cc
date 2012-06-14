// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPairIntraResC3.cc
/// @brief  Count pair for residue pairs connected with one bond, where the
/// crossover from excluding to counting atom pair interactions is at 3 bonds.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#include <core/scoring/etable/count_pair/CountPairIntraResC3.hh>
#include <core/scoring/etable/atom_pair_energy_inline.hh>

#include <core/scoring/etable/EtableEnergy.hh>

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {


/// @brief take a reference to the path distances table
CountPairIntraResC3::CountPairIntraResC3(
	conformation::Residue const & res
) :
	parent(),
	path_dists_( res.path_distances() )
{
}

CountPairIntraResC3::~CountPairIntraResC3() {}

bool
CountPairIntraResC3::count(
	int const at1,
	int const at2,
	Real & w,
	Size & path_dist
) const
{
	return operator() ( at1, at2, w, path_dist );
}


void
CountPairIntraResC3::residue_atom_pair_energy(
	conformation::Residue const & res,
	conformation::Residue const & ,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_intraresidue_atom_pair_energy( res, etable_energy, *this, emap );
}


void
CountPairIntraResC3::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & ,
	conformation::Residue const & ,
	etable::TableLookupEvaluator const & ,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraResC3::residue_atom_pair_energy_sidechain_backbone" << std::endl;
	utility_exit();
}


void
CountPairIntraResC3::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & ,
	conformation::Residue const & ,
	etable::TableLookupEvaluator const & ,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraResC3::residue_atom_pair_energy_sidechain_whole" << std::endl;
	utility_exit();
}

void
CountPairIntraResC3::residue_atom_pair_energy(
	conformation::Residue const & res,
	conformation::Residue const & ,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_intraresidue_atom_pair_energy( res, etable_energy, *this, emap );
}


void
CountPairIntraResC3::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & ,
	conformation::Residue const & ,
	etable::AnalyticEtableEvaluator const & ,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraResC3::residue_atom_pair_energy_sidechain_backbone" << std::endl;
	utility_exit();
}


void
CountPairIntraResC3::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & ,
	conformation::Residue const & ,
	etable::AnalyticEtableEvaluator const & ,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraResC3::residue_atom_pair_energy_sidechain_whole" << std::endl;
	utility_exit();
}

} // namespace count_pair
} // namespace etable {
} // namespace scoring
} // namespace core
