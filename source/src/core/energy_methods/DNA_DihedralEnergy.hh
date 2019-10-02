// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/energy_methods/DNA_DihedralEnergy.hh
/// @brief  dna scoring
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_DNA_DihedralEnergy_HH
#define INCLUDED_core_scoring_methods_DNA_DihedralEnergy_HH

// Unit headers
#include <core/energy_methods/DNA_DihedralEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/dna/DNA_DihedralPotential.fwd.hh>
//#include <core/scoring/dna/DNA_DihedralPotential.fwd.hh>


namespace core {
namespace scoring {
namespace methods {

///
class DNA_DihedralEnergy : public ContextDependentOneBodyEnergy {
public:
	typedef ContextDependentOneBodyEnergy  parent;
public:

	/// ctor
	DNA_DihedralEnergy();
	DNA_DihedralEnergy( DNA_DihedralEnergy const & src);

	/// clone
	EnergyMethodOP
	clone() const override;

	bool
	defines_score_for_residue(
		conformation::Residue const &
	) const override;

	bool
	defines_dof_derivatives( pose::Pose const & ) const override { return true; }

	void
	configure_from_options_system();

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	///
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const override;

	//  virtual
	//  void
	//  residue_energy_old(
	//   conformation::Residue const & rsd,
	//   pose::Pose const & pose,
	//   EnergyMap & emap
	//  ) const;


	///
	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const &,// dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		ScoreFunction const &,// sfxn,
		EnergyMap const & weights
	) const;


	Real
	eval_residue_dof_derivative(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const &, // min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights
	) const override;


	/// this function is used for the sugar derivs, which dont match up to a "torsion" in the current scheme
	void
	eval_residue_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const override;

	//  virtual
	//  Real
	//  eval_dof_derivative_old(
	//   id::DOF_ID const &,// dof_id,
	//   id::TorsionID const & tor_id,
	//   pose::Pose const & pose,
	//   ScoreFunction const &,// sfxn,
	//   EnergyMap const & weights
	//  ) const;

	/// @brief DNA_Dihedral Energy is context dependent, but indicates that no context graphs need to
	/// be maintained by class Energies
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const override;

	core::Size version() const override { return 1; }

	// data
private:
	dna::DNA_DihedralPotential const & potential_;
	bool score_delta_;
	bool score_chi_;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
