// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/energy_methods/FaMPEnvEnergy.hh
///
/// @brief  LK-Type Membrane Environment Energy
///
/// @author  Patrick Barth (Original)
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPEnvEnergy_hh
#define INCLUDED_core_scoring_membrane_FaMPEnvEnergy_hh

// Unit headers
#include <core/energy_methods/FaMPEnvEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Project Headers

// Package headers
#include <core/scoring/memb_etable/MembEtable.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/conformation/Atom.fwd.hh>

#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <ObjexxFCL/FArray1.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh> // DO NOT AUTO-REMOVE

// C++ Headers

namespace core {
namespace energy_methods {

/// @brief Fullatom Membrane Environment Energy
class FaMPEnvEnergy : public core::scoring::methods::ContextDependentOneBodyEnergy {

public:

	typedef core::scoring::methods::ContextDependentOneBodyEnergy parent;

	/// @brief Construct Energy Method from Etable
	FaMPEnvEnergy( core::scoring::etable::MembEtableCAP memb_etable_in );

	/// @brief Clone Energy Method
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/// @brief Compute Per-Residue Energies
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const override;

	/// @brief Fianlzie Total Per-Residue Energies
	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	/// @brief Setup for Computing Derivatives
	void
	setup_for_derivatives(
		pose::Pose & pose,
		core::scoring::ScoreFunction const & scfxn
	) const override;

	/// @brief Evaluate Per-Atom Derivatives
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & emap,
		Vector & F1,
		Vector & F2
	) const override;

	/// @brief Fa_MbenvEnergy is context independent
	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const override;

	/// @brief Setup Method for initial scoring
	void
	setup_for_scoring(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &
	) const override;

private: // helper methods

	/// @brief Evaluate Per-Atom Env term
	Real
	eval_fa_mbenv(
		conformation::Atom const & atom1,
		Real const & f1
	) const;

	/// @brief Versioning
	core::Size version() const override;

	/// @brief Initialize Energy Method data for derivatives
	void
	init( pose::Pose & pose ) const;

	/// @brief Allocate memory for derivatives
	void setup_for_fullatom( pose::Pose & pose ) const;

private:

	// Store a copy of the etable in construction
	core::scoring::etable::MembEtableCAP memb_etable_;

	// Make copies from the etable
	ObjexxFCL::FArray1< Real > const & lk_dgrefce_;
	ObjexxFCL::FArray1< Real > const & memb_lk_dgrefce_;

	// Store mbenv weight when computing derivatives
	mutable Real fa_mbenv_weight_;

	// Arrays used for computing derivatives
	mutable utility::vector1 < utility::vector1 < Real > > fa_proj_;
	mutable utility::vector1 < utility::vector1 < Vector > > fa_f1_;
	mutable utility::vector1 < utility::vector1 < Vector > > fa_f2_;

};

} // scoring
} // core

#endif // INCLUDED_core_energy_methods_FaMPEnvEnergy_hh
