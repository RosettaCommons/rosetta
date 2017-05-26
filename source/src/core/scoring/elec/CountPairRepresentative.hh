// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/elec/CountPairRepresentative.hh
/// @brief  Count pair function for use in initializing the NBList that uses the
/// FA_ElecEnergy::get_countpair_representative_atom function to determine the representative
/// for an atom, and then to ask the inner CountPairFunction whether the two representatives
/// should have their interaction counted. Note: this count pair function is not set up to be
/// used efficiently in the atom_pair_energy_inline functions.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_elec_CountPairRepresentative_hh
#define INCLUDED_core_scoring_elec_CountPairRepresentative_hh

#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <core/scoring/elec/FA_ElecEnergy.hh>

#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>


namespace core {
namespace scoring {
namespace elec {

class CountPairRepresentative : public etable::count_pair::CountPairFunction
{
public:
	typedef etable::count_pair::CountPairFunction parent;

public:

	CountPairRepresentative(
		FA_ElecEnergy const & energy,
		conformation::Residue const & r1,
		conformation::Residue const & r2,
		etable::count_pair::CountPairFunctionCOP inner_cpfxn
	) :
		energy_( energy ),
		r1_type_( r1.type() ),
		r2_type_( r2.type() ),
		inner_cpfxn_( inner_cpfxn )
	{}

	virtual ~CountPairRepresentative() {}

	bool
	count(
		int const at1,
		int const at2,
		Real & weight,
		Size & path_dist
	) const override
	{
		Size rep1 = energy_.get_countpair_representative_atom( r1_type_, at1 );
		Size rep2 = energy_.get_countpair_representative_atom( r2_type_, at2 );
		return inner_cpfxn_->count( rep1, rep2, weight, path_dist );
	}

	void
	residue_atom_pair_energy(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const override {}

	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const override {}

	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const override {}

	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const override {}

	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const override {}

	void
	residue_atom_pair_energy(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const override {}

	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const override {}

	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const override {}

	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const override {}

	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const override {}

private:
	FA_ElecEnergy const & energy_;
	chemical::ResidueType const & r1_type_;
	chemical::ResidueType const & r2_type_;
	etable::count_pair::CountPairFunctionCOP inner_cpfxn_;
};


} // namespace ele
} // namespace scoring
} // namespace core

#endif
