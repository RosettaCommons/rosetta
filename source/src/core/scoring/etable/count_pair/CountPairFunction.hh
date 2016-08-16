// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/count_pair/CountPairFunction.hh
/// @brief  Count pair base class interface
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_etable_count_pair_CountPairFunction_hh
#define INCLUDED_core_scoring_etable_count_pair_CountPairFunction_hh

// Unit headers
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>


// C++

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

class CountPairFunction : public utility::pointer::ReferenceCount
{

public:
	CountPairFunction() {} // inlined when declared on the stack
	virtual ~CountPairFunction() {} // inlined when declared on the stack

	// each derived class must override this (non virtual) function
	// in order to work with residue_atom_pair_energy< T > defined in
	// atom_pair_energy_inline.hh
	bool
	operator () (
		int const at1,
		int const at2,
		Real &,
		Size & path_dist
	) const;

	// virtual method calls operator(); for use when type resolution is non-critical
	virtual
	bool
	count(
		int const at1,
		int const at2,
		Real &,
		Size & path_dist
	) const = 0;

	virtual
	void
	residue_atom_pair_energy(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const = 0;

	virtual
	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const = 0;


	virtual
	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const = 0;

	virtual
	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const = 0;


	virtual
	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const = 0;

	virtual
	void
	residue_atom_pair_energy(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const = 0;

	virtual
	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const = 0;


	virtual
	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const = 0;

	virtual
	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const = 0;


	virtual
	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const = 0;

	/*static
	count_pair_type
	find_count_pair_type( std::string const & ); */

	static Real const cp_half;

};

} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core


#endif
