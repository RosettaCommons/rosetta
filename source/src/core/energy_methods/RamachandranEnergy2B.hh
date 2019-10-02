// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/RamachandranEnergy2B.hh
/// @brief  Ramachandran energy method class declaration
/// @author Guoli Wang

#ifndef INCLUDED_core_scoring_methods_RamachandranEnergy2B_hh
#define INCLUDED_core_scoring_methods_RamachandranEnergy2B_hh

// Unit headers
#include <core/energy_methods/RamachandranEnergy2B.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/Ramachandran2B.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class RamachandranEnergy2B : public ContextIndependentTwoBodyEnergy  {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;
public:

	/// ctor
	RamachandranEnergy2B();

	/// clone
	EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentTwoBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap & emap
	) const override;

	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const override;

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const &, // unused,
		ScoreFunction const &, // unused,
		EnergyMap & emap
	) const override;


	Distance
	atomic_interaction_cutoff() const override;


	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const &,// dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		ScoreFunction const &,// sfxn,
		EnergyMap const & weights
	) const;

	/// @brief Ramachandran Energy is context independent and thus indicates that no context graphs need to
	/// be maintained by class Energies
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const override;

	// data
private:
	Ramachandran2B const & potential_;
	core::Size version() const override;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
