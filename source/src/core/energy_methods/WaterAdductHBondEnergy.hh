// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/WaterAdductHBondEnergy.hh
/// @brief  Energy potential for water mediated hydrogen bonds
///         involving adduct-placed water molecules
/// @author Jim Havranek


#ifndef INCLUDED_core_scoring_methods_WaterAdductHBondEnergy_hh
#define INCLUDED_core_scoring_methods_WaterAdductHBondEnergy_hh

// Unit headers
#include <core/energy_methods/WaterAdductHBondEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/WaterAdductHBondPotential.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class WaterAdductHBondEnergy : public ContextIndependentTwoBodyEnergy  {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;

public:

	/// ctor
	WaterAdductHBondEnergy();

	/// clone
	EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentTwoBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const override;

	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const override;

	void
	setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const override;


	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const override {}


	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;

	void
	get_atom_h2o_hbond_derivative(
		id::AtomID const & atom,
		hbonds::HBondSet const & hbond_set,
		EnergyMap const & weights,
		Vector & f1,
		Vector & f2
	) const;

	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const override { return false; }

	Distance
	atomic_interaction_cutoff() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const override;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	core::scoring::WaterAdductHBondPotential const & potential_;
	core::Size version() const override;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_methods_WaterAdductHBondEnergy_HH
