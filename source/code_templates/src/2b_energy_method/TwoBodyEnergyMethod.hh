// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.hh
/// @brief --brief--
/// @author --name-- (--email--)

#ifndef INCLUDED_--path_underscore--_--class--_hh
#define INCLUDED_--path_underscore--_--class--_hh

// Unit headers
#include <--path--/--class--.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.hh> 
#include <core/scoring/methods/EnergyMethod.hh> 

// Core headers
#include <core/pose/Pose.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>

#include <core/scoring/DerivVectorPair.fwd.hh>
#include <utility/vector1.hh>

--namespace--

///@brief --brief--
/// The Two Body Energy Method specifies an interface for all two body methods: both
/// long and short range, both context dependent and independent.  Any two body method
/// must implement this interface as well as the EnergyMethod interface.
class --class-- : public core::scoring::methods::TwoBodyEnergy {

public:

	--class--();

	// copy constructor (not needed unless you need deep copies)
	//--class--( --class-- const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~--class--();


	/// @brief Evaluate the interaction between a given residue pair
	/// accumulating the unweighted energies in an EnergyMap
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;


	/// @brief During minimization, energy methods are allowed to decide that they say nothing
	/// about a particular residue pair (e.g. no non-zero energy) and as a result they will not be queried for
	/// a derivative or an energy.  The default implementation returns "true" for all residue pairs.
	/// Context-dependent two-body energies have the option of behaving as if they are context-independent
	/// by returning "false" for residue pairs that do no move wrt each other.
	virtual
	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const;


	/// @brief Rely on the extended version of the residue_pair_energy function during score-function
	/// evaluation in minimization? The extended version (below) takes a ResPairMinimizationData in which
	/// the derived base class has (or should have) cached a piece of data that will make residue-pair
	/// energy evaluation faster than its absense (e.g. a neighbor list). Derived energy methods should
	/// return 'true' from this function to use the extended interface. The default method implemented
	/// in this class returns 'false'
	virtual
	bool
	use_extended_residue_pair_energy_interface() const;

	/// @brief Evaluate the two-body energies for a particular residue, in the context of a
	/// given Pose, and with the help of a piece of cached data for minimization, increment those
	/// two body energies into the input EnergyMap.  The calling function must guarantee that this
	/// EnergyMethod has had the opportunity to update the input ResPairMinimizationData object
	/// for the given residues in a call to setup_for_minimizing_for_residue_pair before this function is
	/// invoked. This function should not be called unless the use_extended_residue_pair_energy_interface()
	/// method returns "true".  Default implementation provided by this base class calls
	/// utility::exit().
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

	/// @brief Called at the beginning of minimization, allowing this energy method to cache data
	/// pertinent for a single residue in the the ResPairMinimizationData that is used for a
	/// particular residue in the context of a particular Pose.  This base class provides a noop
	/// implementation for this function if there is nothing that the derived class needs to perform
	/// in this setup phase.
	virtual
	void
	setup_for_minimizing_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		ResSingleMinimizationData & res_data_cache
	) const;

	/// @brief Called at the beginning of minimization, allowing this energy method to cache data
	/// pertinent for a single residue in the the ResPairMinimizationData that is used for a
	/// particular residue in the context of a particular Pose.  This base class provides a noop
	/// implementation for this function if there is nothing that the derived class needs to perform
	/// in this setup phase.
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
		ResPairMinimizationData & data_cache
	) const;

	/// @brief Does this EnergyMethod require the opportunity to examine the residue before scoring begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residues that are uninterested
	/// in doing so.
	virtual
	bool
	requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const & pose ) const;

	/// @brief Do any setup work should the coordinates of this residue (who is still guaranteed to be
	/// of the same residue type as when setup_for_minimizing_for_residue was called) have changed so dramatically
	/// as to possibly require some amount of setup work before scoring should proceed.
	/// This function is used for both intra-residue setup and pre-inter-residue setup
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

	/// @brief Does this EnergyMethod require the opportunity to examine each residue pair before scoring begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residue pairs that are uninterested
	/// in doing so.
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

	/// @brief Does this EnergyMethod require the opportunity to examine each residue pair before derivative evaluation begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residue pairs that are uninterested
	/// in doing so.
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

	/// @brief Evaluate the derivatives for all atoms on rsd1 and rsd2 with respect
	/// to each other and increment the derivatives in atom-derivatives vector1s.
	/// The calling function must guarantee that the r1_atom_derivs vector1 holds at
	/// least as many entries as there are atoms in rsd1, and that the r2_atom_derivs
	/// vector1 holds at least as many entries as there are atoms in rsd2.
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

	/// @brief Evaluate the interaction between the backbone of rsd1 and the
	/// backbone of rsd2 and accumulate the unweighted energies.  The sum
	/// bb_bb(r1,r2) + bb_sc(r1,r2) + bb_sc(r2,r1) + sc_sc( r1,r2) must
	/// equal the weighted result of a call to residue_pair_energy.
	/// By default, bb_bb & bb_sc return 0 and sc_sc returns
	/// residue pair energy.
	virtual
	void
	backbone_backbone_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;


	/// @brief Evaluate the interaction between the backbone of rsd1 and the
	/// sidechain of rsd2 and accumulate the unweighted energies.  The sum
	/// bb_bb(r1,r2) + bb_sc(r1,r2) + bb_sc(r2,r1) + sc_sc( r1,r2) must
	/// equal the unweighted result of a call to residue_pair_energy.
	/// By default, bb_bb & bb_sc return 0 and sc_sc returns
	/// residue pair energy.
	virtual
	void
	backbone_sidechain_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/// @brief Evaluate the interaction between the sidechain of rsd1 and the
	/// sidechain of rsd2 and accumulate the unweighted energies.  The sum
	/// bb_bb(r1,r2) + bb_sc(r1,r2) + bb_sc(r2,r1) + sc_sc( r1,r2) must
	/// equal the unweighted result of a call to residue_pair_energy.
	/// By default, bb_bb & bb_sc return 0 and sc_sc returns
	/// residue pair energy.
	virtual
	void
	sidechain_sidechain_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/// @brief Two body energies are able to define intra-residue energies, and to do so
	/// only in the presence of certain non-zero weights.  The ScoreFunction will hand over its
	/// weight set as it asks whether the energy method defines an intraresidue energy or not.
	///
	/// For example, the Etable method defines intra-residue energies only when one or more
	/// of the fa_intra_{atr,rep,sol} weights are non-zero.
	virtual
	bool
	defines_intrares_energy( EnergyMap const & weights ) const = 0;

	/// @brief Evaluate the intra-residue energy for a given residue
	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const = 0;

	/// @brief If a score function defines no intra-residue scores for a particular
	/// residue, then it may opt-out of being asked during minimization to evaluate
	/// the score for this residue.
	virtual
	bool
	defines_intrares_energy_for_residue(
		conformation::Residue const & res
	) const;

	/// @brief Derived classes wishing to invoke the alternate, extended interface for eval_intrares_energy
	/// during minimization routines should return "true" when this function is invoked on them.  This
	/// class provides a default "return false" implementation so that classes not desiring to take advantage
	/// of this alternate interface need to do nothing.
	virtual
	bool
	use_extended_intrares_energy_interface() const;

	/// @brief Evaluate the intra-residue energy for a given residue using the data held within the
	/// ResSingleMinimizationData object.  This function should be invoked only on derived instances
	/// of this class if they return "true" in a call to their use_extended_intrares_energy_interface
	/// method.  This base class provides a noop implementation for classes that do not implement this
	/// interface, or that do not define intrares energies.
	virtual
	void
	eval_intrares_energy_ext(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & data_cache,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;


	/// @brief Evaluate the derivative for the intra-residue component of this energy method
	/// for all the atoms in a residue in the context of a particular pose,
	/// and increment the F1 and F2 vectors held in the atom_derivs vector1.
	/// This base class provides a default noop implementation
	/// of this function. The calling function must guarantee that this EnergyMethod has had the
	/// opportunity to update the input ResSingleMinimizationData object for the given residue
	/// in a call to prepare_for_minimization before this function is invoked.
	/// The calling function must also guarantee that there are at least as many entries
	/// in the atom_derivs vector1 as there are atoms in the input rsd.
	virtual
	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;

	/// @brief Use the dof_derivative interface for this energy method when
	/// calculating derivatives?  It is possible to define both dof_derivatives and
	/// atom-derivatives; they are not mutually exclusive.
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


	// @brief Low fidelity evaluation of sc/bb + sc/sc energy.
	// Do not define in derived class if that class should not be
	// used in the packer's bump-check phase.
	virtual
	void
	bump_energy_full(
		conformation::Residue const &,
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const;

	// @brief Low fidelity evaluation of sc/bb energy.
	// Do not define in derived class if that class should not be
	// used in the packer's bump-check phase.
	virtual
	void
	bump_energy_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const;

	/// @brief Batch computation of rotamer intrares energies.  Need not be overriden in
	/// derived class -- by default, iterates over all rotamers,
	/// and calls derived class's intrares _energy method.
	virtual
	void
	evaluate_rotamer_intrares_energies(
		conformation::RotamerSetBase const & set,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		utility::vector1< core::PackerEnergy > & energies
	) const;

	/// @brief Batch computation of rotamer intrares energy map.  Need not be overriden in
	/// derived class -- by default, iterates over all rotamers,
	/// and calls derived class's intrares _energy method.
	virtual
	void
	evaluate_rotamer_intrares_energy_maps(
		conformation::RotamerSetBase const & set,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		utility::vector1< EnergyMap > & emaps
	) const;

	/// @brief Batch computation of rotamer pair energies.  Need not be overriden in
	/// derived class -- by default, iterates over all pairs of rotamers,
	/// and calls the derived class's residue_pair_energy method.
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
	/// derived class's residue_pair_energy method for each one against the background rotamr
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
	/// derived class's residue_pair_energy method for each one against the background rotamr
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

private:

};

std::ostream &
operator<<( std::ostream & os, --class-- const & mover );

--end_namespace--

#endif //--path--_--class--_hh
