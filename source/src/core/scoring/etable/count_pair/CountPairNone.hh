// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPairNone.hh
/// @brief  Count pair for residues where all atom pairs should be counted.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_etable_count_pair_CountPairNone_hh
#define INCLUDED_core_scoring_etable_count_pair_CountPairNone_hh

// Unit Headers
#include <core/scoring/etable/count_pair/CountPairNone.fwd.hh>

// Package Headers
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

class CountPairNone : public CountPairFunction
{
public:
	public:
	typedef CountPairFunction parent;

public:
	CountPairNone();
	virtual ~CountPairNone();

	/// @brief function required by templated functions in atom_pair_energy_inline
	// -- those functions are in fact never reached since residue_atom_pair_energy simply returns
	// without ever doing any work.  This function is merely to mimic the style in the
	// other count pair classes where the nblist needs to ask which pairs count: the nblist
	// uses the count() method, which in tern invokes operator().
	inline
	bool
	operator () (
		int const /*at1*/,
		int const /*at2*/,
		Real & /*weight*/,
		Size & /*path_dist*/
	) const
	{
		return false;
	}

	virtual
	bool
	count(
		int const at1,
		int const at2,
		Real & w,
		Size & path_dist
	) const;

	/// Type resolution functions
	virtual
	void
	residue_atom_pair_energy(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const;

	virtual
	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const;

	virtual
	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const;

	virtual
	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const;

	/// Type resolution functions
	virtual
	void
	residue_atom_pair_energy(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const;

	virtual
	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const;

	virtual
	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const;

	virtual
	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const;

};

} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

#endif
