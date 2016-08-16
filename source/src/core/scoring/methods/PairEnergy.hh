// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/PairEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_PairEnergy_hh
#define INCLUDED_core_scoring_methods_PairEnergy_hh

// Unit Headers
#include <core/scoring/methods/PairEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/PairEPotential.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace methods {


class PairEnergy : public ContextDependentTwoBodyEnergy  {
public:
	typedef ContextDependentTwoBodyEnergy  parent;

public:


	PairEnergy();


	/// clone
	virtual
	EnergyMethodOP
	clone() const;


	virtual
	void
	setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;


	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;


	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & set
	) const;

	virtual
	void
	update_residue_for_packing(
		pose::Pose &,
		Size resid
	) const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	virtual
	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const;

	/*virtual
	void
	eval_atom_derivative_for_residue_pair(
	Size const atom_index,
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const & minsingle_data1,
	ResSingleMinimizationData const & minsingle_data2,
	ResPairMinimizationData const & minpair_data,
	pose::Pose const & pose, // provides context
	kinematics::DomainMap const & domain_map,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
	) const;*/

	virtual
	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	/// @brief Batch computation of rotamer pair energies.  Need not be overriden in
	/// derived class -- by default, iterates over all pairs of rotamers,
	/// and calls derived class's residue_pair_energy method.  Since short range rotamer pairs
	/// may not need calculation, the default method looks at blocks of residue type pairs
	/// and only calls the residue_pair_energy method if the rotamer pairs are within range
	virtual
	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const;


	/// @brief Batch computation of rotamer/background energies.  Need not be overriden
	/// in derived class -- by default, iterates over all rotamers in the set, and calls
	/// derived class's residue_pair_energy method for each one against the background rotamer
	/// Since short range rotamer pairs may not need calculation, the default method
	/// looks at blocks of residue type pairs and only calls the residue_pair_energy method
	/// if the rotamer pairs are within range
	virtual
	void
	evaluate_rotamer_background_energies(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		utility::vector1< core::PackerEnergy > & energy_vector
	) const;

	/// @brief Batch computation of rotamer/background energies.  Need not be overriden
	/// in derived class -- by default, iterates over all rotamers in the set, and calls
	/// derived class's residue_pair_energy method for each one against the background rotamer
	/// Since short range rotamer pairs may not need calculation, the default method
	/// looks at blocks of residue type pairs and only calls the residue_pair_energy method
	/// if the rotamer pairs are within range
	virtual
	void
	evaluate_rotamer_background_energy_maps(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		utility::vector1< EnergyMap > & emaps
	) const;


	/*virtual
	void
	eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
	) const;*/

	virtual
	Distance
	atomic_interaction_cutoff() const;

	/// @details non-virtual accessor for speed
	Distance
	interaction_cutoff() const;


	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	// const-ref to scoring database
	PairEPotential const & potential_;
	virtual
	core::Size version() const;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
