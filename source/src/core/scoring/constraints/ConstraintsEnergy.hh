// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/ConstraintsEnergy.hh
/// @brief  Constraints Energy Method declaration
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_constraints_ConstraintsEnergy_hh
#define INCLUDED_core_scoring_constraints_ConstraintsEnergy_hh

// Unit headers
#include <core/scoring/constraints/ConstraintsEnergy.fwd.hh>

// Package headers

#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/Methods.hh>

#include <core/kinematics/MinimizerMapBase.fwd.hh>


// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace constraints {


class ConstraintsEnergy : public methods::ContextIndependentLRTwoBodyEnergy {
public:
	typedef methods::ContextIndependentLRTwoBodyEnergy parent;
	typedef methods::EnergyMethodOP EnergyMethodOP;

public:
	ConstraintsEnergy();
	virtual ~ConstraintsEnergy();

	virtual EnergyMethodOP clone() const;

	/// @brief Identification for this LR2B energy that links it with the
	/// long-range energy container that it stores in the Energies object
	virtual methods::LongRangeEnergyType long_range_type() const;

	/// @brief Interface from the LongRangeTwoBodyEnery base class; returns "true" if there's any non-zero
	/// or potentially non-zero interaction between a pair of residues in a pose.
	virtual
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const;


	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentTwoBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const;

	/// @brief Returns true as the ConstraintsEnergy caches its residue pair constraints in
	/// the ResPairMinimicationData objects stored on edges in the MinimizationGraph
	virtual
	bool
	use_extended_residue_pair_energy_interface() const;

	/// @brief The ConstraintsEnergy caches all residue-pair constraints operating between a particular
	/// pair of residues in the min_data object so that it does not need to search for those constraints
	/// when scoring during minimization.
	virtual
	void
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/// @brief Returns false if residue1 and residue2 have no inter-residue pair constraints
	/// or if the two residues are not moving wrt each other.
	virtual
	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const;


	/// @brief Cache the intra-residue constraints in the ResSingleMinimizationData object
	/// for rapid retrieval during minimization and allow the constraints to store data
	/// in the res_data_cache if they need to.
	virtual
	void
	setup_for_minimizing_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		ResSingleMinimizationData & res_data_cache
	) const;

	/// @brief Cache the residue-pair constraints in the ResPairMinimizationData object
	/// for rapid retrieval during minimization and allow constraints to cache useful
	/// data in the respair_data_cache if they need to
	virtual
	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		ResSingleMinimizationData const & res1_data_cache,
		ResSingleMinimizationData const & res2_data_cache,
		ResPairMinimizationData & respair_data_cache
	) const;


	/// @brief Request the opportunity to setup for scoring
	virtual
	bool
	requires_a_setup_for_scoring_for_residue_pair_opportunity( pose::Pose const & pose ) const;

	/// @brief Do any setup work should the coordinates of a pair of residues, who are still guaranteed to be
	/// of the same residue type as when setup_for_minimizing_for_residue was called, have changed so dramatically
	/// as to possibly require some amount of setup work before scoring should proceed
	virtual
	void
	setup_for_scoring_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const & minsingle_data1,
		ResSingleMinimizationData const & minsingle_data2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		ResPairMinimizationData & data_cache
	) const;

	/// @brief Ask for the opportunity to setup for derivative evaluation
	virtual
	bool
	requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & pose ) const;

	/// @brief Do any setup work necessary before evaluating the derivatives for this residue pair
	virtual
	void
	setup_for_derivatives_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const & minsingle_data1,
		ResSingleMinimizationData const & minsingle_data2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		ResPairMinimizationData & data_cache
	) const;

	/// @brief Evaluate the derivative for an atom in rsd1 with respect to rsd2 in the context
	/// of a particular pose, and increment the F1 and F2 vectors.
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


	void
	prepare_constraints_energy_container( pose::Pose & pose ) const;

	// check compatibility with atomtypeset and store a constraint graph in
	// the pose.energies object
 	virtual
 	void
 	setup_for_scoring( pose::Pose &pose, ScoreFunction const &scfxn ) const;

	/// @brief Make sure that the ConstraintsEnergyContainer is ready for packing.
	virtual
	void
	setup_for_packing( 
		pose::Pose & pose,
		utility::vector1< bool > const & residues_repacking,
		utility::vector1< bool > const & residues_designing
	) const;

	// call the cst setup_for_derivatives wrapper
 	virtual
 	void
 	setup_for_derivatives( pose::Pose &pose, ScoreFunction const &scfxn ) const;

	/// @brief Returns true if there are three-body or higher-body constraints that the user
	/// has defined.  High-order terms will not be correctly evaluated in RTMin or other
	/// extra-pose techniques for minimization (but will work correctly when minimizing an entire Pose).
	bool
	defines_high_order_terms( pose::Pose const & ) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & weights ) const;

	/// @brief Evaluate the intra-residue constraint energy for a given residue
	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/// @brief request of minimization routines that they use the extended intraresidue energy
	/// interface
	virtual
	bool
	use_extended_intrares_energy_interface() const;

	/// @brief Evaluate the intra-residue energies using ConstraintCOPs cached in the data_cache object
	virtual
	void
	eval_intrares_energy_ext(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & data_cache,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/// @brief Ask for the opportunity to setup for scoring for evaluating 1-body constraints
	virtual
	bool
	requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const & pose ) const;

	/// @brief Allow the component constraints to setup for scoring for a particular residue
	virtual
	void
	setup_for_scoring_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		ResSingleMinimizationData & min_data
	) const;

	/// @brief Does this EnergyMethod require the opportunity to examine each residue before derivative evaluation begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residue pairs that are uninterested
	/// in doing so.
	virtual
	bool
	requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const & pose ) const;

	/// @brief Do any setup work necessary before evaluating the derivatives for this residue
	virtual
	void
	setup_for_derivatives_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		ResSingleMinimizationData & min_data
	) const;

	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;

	/// @brief The DunbrackConstraint, which adds a "constant term" to the energy for a particular
	/// rotamer so that the particular rotamer matches the score of the best rotamer at a particular
	/// phi/psi, defines derivatives for phi and psi. or rather -- it doesn't, but it should.
	virtual
	bool
	defines_intrares_dof_derivatives( pose::Pose const & p ) const;

	/// @brief Evaluate the DOF derivative for a particular residue.  The Pose merely serves as context,
	/// and the input residue is not required to be a member of the Pose.
	virtual
	Real
	eval_intraresidue_dof_derivative(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights
	) const;

	/// called at the end of energy evaluation
	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;


	/// called during gradient-based minimization inside dfunc
	/**
		 F1 and F2 are not zeroed -- contributions from this atom are
		 just summed in
	**/
	virtual
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;


	/// @brief Evaluate dof derivatives when DOF_constraints are in use
	Real
	eval_dof_derivative(
		id::DOF_ID const & id,
		id::TorsionID const & tor,
		pose::Pose const & pose,
		ScoreFunction const & scorefxn,
		EnergyMap const & weights
	) const;


	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const;
virtual
core::Size version() const;

};


} // constraints
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
