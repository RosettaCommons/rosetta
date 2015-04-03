// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ProClosureEnergy.hh
/// @brief  Proline ring closure energy method class declaration
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_methods_ProClosureEnergy_hh
#define INCLUDED_core_scoring_methods_ProClosureEnergy_hh

// Unit headers
#include <core/scoring/methods/ProClosureEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>

#include <core/kinematics/DomainMap.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class ProClosureEnergy : public ContextIndependentTwoBodyEnergy  {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;

public:

	/// ctor
	ProClosureEnergy();

	// dstor
	~ProClosureEnergy();

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentTwoBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	/// @brief Pro-closure terms only apply between bonded residues where i+1 is
	/// proline -- skip residue pairs that don't apply during minimization.
	virtual
	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const;

	/// @brief Evaluate the interaction between a given residue pair
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	/// @brief Evaluate the derivative for the input residue pair.
	/// This will only be called if rsd1 and rsd2 are bonded and one of them is a proline
	/// because ProClosureEnergy defines the "defines_score_for_residue_pair" method.
	/*virtual
	void
	eval_atom_derivative_for_residue_pair(
		Size const atom_index,
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const & minsingle_data1,
		ResSingleMinimizationData const & minsingle_data2,
		ResPairMinimizationData const & min_data,
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
		ResPairMinimizationData const &,
		pose::Pose const &,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	/// @brief Non-virtual interface; takes only the needed parameters.
	/*void
	eval_atom_derivative_for_residue_pair2(
		Size const atom_index,
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;*/

	/// @brief Penalize the pucker-up residue type if its chi1 is positive;
	/// penalize the pucker-down residue type if its chi1 is negative.  Only
	/// applies this penalty when the other_residue is the next polymeric residue
	/// after pro_residue (i+1), unless pro_residue is an upper_term,
	/// in which case it applies the penalty for pro_residue's previous polymeric
	/// residue.
	virtual
	void
	bump_energy_full(
		conformation::Residue const & pro_residue,
		conformation::Residue const & other_residue,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const;

	/// @brief Penalize the pucker-up residue type if its chi1 is positive;
	/// penalize the pucker-down residue type if its chi1 is negative.  Only
	/// applies this penalty when the other_residue is the next polymeric residue
	/// after pro_residue (i+1), unless pro_residue is an upper_term,
	/// in which case it applies the penalty for pro_residue's previous polymeric
	/// residue.
	virtual
	void
	bump_energy_backbone(
		conformation::Residue const & pro_residue,
		conformation::Residue const & other_residue,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const;


	virtual
	bool
	defines_intrares_energy( EnergyMap const & weights ) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;


	/// @brief Returns false if res is not a proline.
	bool
	defines_intrares_energy_for_residue(
		conformation::Residue const & res
	) const;

	/// @brief This should only be handed a proline.
	virtual
	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;

	/// @brief ProClosure Energy is context independent and thus
	/// indicates that no context graphs need to
	/// be maintained by class Energies
	virtual
	void indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const;

	virtual
	Distance
	atomic_interaction_cutoff() const
	{ return 4.0; }

private:

	/// @brief measure in radians the chi4 between two residues;
	/// upper_residue must be a proline chi4 is wrapped to be in
	/// the range [-pi_over_2, 3/2*pi )
	Real
	measure_chi4(
		conformation::Residue const & lower_residue,
		conformation::Residue const & upper_residue
	) const;

	Real
	chi4E(
		Real chi4
	) const;

	Real
	dchi4E_dchi4(
		Real chi4
	) const;

	// data
private:
	Real const n_nv_dist_sd_; // coorinate variation stdev

	Real const trans_chi4_mean_;
	Real const trans_chi4_sd_;
	Real const cis_chi4_mean_;
	Real const cis_chi4_sd_;

	std::string const bbN_;
	std::string const scNV_;
	std::string const scCD_;
	std::string const bbC_;
	std::string const bbO_;
virtual
core::Size version() const;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
