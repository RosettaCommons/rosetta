// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/energy_methods/RNA_CoarseDistEnergy.hh
/// @brief Score two-body energies in coarse RNA poses between P, S, and CEN using a statistical potential.
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_core_energy_methods_RNA_CoarseDistEnergy_hh
#define INCLUDED_core_energy_methods_RNA_CoarseDistEnergy_hh

// Unit headers
#include <core/energy_methods/RNA_CoarseDistEnergy.fwd.hh>
#include <core/energy_methods/RNA_CoarseDistEnergyCreator.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/OneDDistPotential.fwd.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/elec/CPRepMapType.fwd.hh>

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

namespace core {
namespace energy_methods {

///@brief Score two-body energies in coarse RNA poses between P, S, and CEN using a statistical potential.
/// The Two Body Energy Method specifies an interface for all two body methods: both
/// long and short range, both context dependent and independent.  Any two body method
/// must implement this interface as well as the EnergyMethod interface.
class RNA_CoarseDistEnergy : public core::scoring::methods::ContextIndependentTwoBodyEnergy {

public:

	//Simply to keep the interface clean

	typedef core::pose::Pose Pose;
	typedef core::conformation::Residue Residue;
	typedef conformation::RotamerSetBase RotamerSetBase;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::ResPairMinimizationData ResPairMinimizationData;

public:

	RNA_CoarseDistEnergy() = delete;

	RNA_CoarseDistEnergy( core::scoring::methods::EnergyMethodOptions const & options );

	// copy constructor (not needed unless you need deep copies)
	//RNA_CoarseDistEnergy( RNA_CoarseDistEnergy const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	~RNA_CoarseDistEnergy() override;

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const override {}

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const override;

	//inline
	Real
	score_atom_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Size const at1,
		Size const at2
	) const;//inline

	Real
	deriv_atom_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Size const at1,
		Size const at2
	) const;

	/// @brief Evaluate the interaction between a given residue pair
	/// accumulating the unweighted energies in an EnergyMap
	void
	residue_pair_energy(
		Residue const & rsd1,
		Residue const & rsd2,
		Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;


	/// @brief During minimization, energy methods are allowed to decide that they say nothing
	/// about a particular residue pair (e.g. no non-zero energy) and as a result they will not be queried for
	/// a derivative or an energy.  The default implementation returns "true" for all residue pairs.
	/// Context-dependent two-body energies have the option of behaving as if they are context-independent
	/// by returning "false" for residue pairs that do no move wrt each other.
	bool
	defines_score_for_residue_pair(
		Residue const & res1,
		Residue const & res2,
		bool res_moving_wrt_eachother
	) const override;


	/// @brief Rely on the extended version of the residue_pair_energy function during score-function
	/// evaluation in minimization? The extended version (below) takes a ResPairMinimizationData in which
	/// the derived base class has (or should have) cached a piece of data that will make residue-pair
	/// energy evaluation faster than its absense (e.g. a neighbor list). Derived energy methods should
	/// return 'true' from this function to use the extended interface. The default method implemented
	/// in this class returns 'false'
	bool
	use_extended_residue_pair_energy_interface() const override;

	/// @brief Returns a regular count-pair function as opposed to a CountPairRepresentative function
	core::scoring::etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2
	) const;


	/// @brief Evaluate the two-body energies for a particular residue, in the context of a
	/// given Pose, and with the help of a piece of cached data for minimization, increment those
	/// two body energies into the input EnergyMap.  The calling function must guarantee that this
	/// EnergyMethod has had the opportunity to update the input ResPairMinimizationData object
	/// for the given residues in a call to setup_for_minimizing_for_residue_pair before this function is
	/// invoked. This function should not be called unless the use_extended_residue_pair_energy_interface()
	/// method returns "true".  Default implementation provided by this base class calls
	/// utility::exit().
	void
	residue_pair_energy_ext(
		Residue const & rsd1,
		Residue const & rsd2,
		core::scoring::ResPairMinimizationData const & min_data,
		Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;

	/// @brief Called at the beginning of minimization, allowing this energy method to cache data
	/// pertinent for a single residue in the the ResPairMinimizationData that is used for a
	/// particular residue in the context of a particular Pose.  This base class provides a noop
	/// implementation for this function if there is nothing that the derived class needs to perform
	/// in this setup phase.
	void
	setup_for_minimizing_for_residue(
		Residue const & rsd,
		Pose const & pose,
		ScoreFunction const & sfxn,
		core::kinematics::MinimizerMapBase const & minmap,
		basic::datacache::BasicDataCache & residue_data_cache,
		core::scoring::ResSingleMinimizationData & res_data_cache
	) const override;

	/// @brief Called at the beginning of minimization, allowing this energy method to cache data
	/// pertinent for a single residue in the the ResPairMinimizationData that is used for a
	/// particular residue in the context of a particular Pose.  This base class provides a noop
	/// implementation for this function if there is nothing that the derived class needs to perform
	/// in this setup phase.
	void
	setup_for_minimizing_for_residue_pair(
		Residue const & rsd1,
		Residue const & rsd2,
		Pose const & pose,
		ScoreFunction const & sfxn,
		core::kinematics::MinimizerMapBase const & minmap,
		core::scoring::ResSingleMinimizationData const & res1_data_cache,
		core::scoring::ResSingleMinimizationData const & res2_data_cache,
		core::scoring::ResPairMinimizationData & data_cache
	) const override;

	/// @brief Does this EnergyMethod require the opportunity to examine the residue before scoring begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residues that are uninterested
	/// in doing so.
	bool
	requires_a_setup_for_scoring_for_residue_opportunity_during_minimization( Pose const & pose ) const override;

	/// @brief Do any setup work should the coordinates of this residue (who is still guaranteed to be
	/// of the same residue type as when setup_for_minimizing_for_residue was called) have changed so dramatically
	/// as to possibly require some amount of setup work before scoring should proceed.
	/// This function is used for both intra-residue setup and pre-inter-residue setup
	void
	setup_for_scoring_for_residue(
		Residue const & rsd,
		Pose const & pose,
		ScoreFunction const & sfxn,
		core::scoring::ResSingleMinimizationData & min_data
	) const override;

	/// @brief Does this EnergyMethod require the opportunity to examine each residue before derivative evaluation begins?  Not
	/// all energy methods would.  The ScoreFunction will not ask energy methods to examine residue pairs that are uninterested
	/// in doing so.
	bool
	requires_a_setup_for_derivatives_for_residue_opportunity( Pose const & pose ) const override;

	/// @brief Do any setup work necessary before evaluating the derivatives for this residue
	void
	setup_for_derivatives_for_residue(
		Residue const & rsd,
		Pose const & pose,
		ScoreFunction const & sfxn,
		core::scoring::ResSingleMinimizationData & min_data,
		basic::datacache::BasicDataCache & res_data_cache
	) const override;

	/// @brief Evaluate the derivatives for all atoms on rsd1 and rsd2 with respect
	/// to each other and increment the derivatives in atom-derivatives vector1s.
	/// The calling function must guarantee that the r1_atom_derivs vector1 holds at
	/// least as many entries as there are atoms in rsd1, and that the r2_atom_derivs
	/// vector1 holds at least as many entries as there are atoms in rsd2.
	void
	eval_residue_pair_derivatives(
		Residue const & rsd1,
		Residue const & rsd2,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData const & min_data,
		Pose const & pose, // provides context
		EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const override;

	/// @brief Evaluate the interaction between the backbone of rsd1 and the
	/// backbone of rsd2 and accumulate the unweighted energies.  The sum
	/// bb_bb(r1,r2) + bb_sc(r1,r2) + bb_sc(r2,r1) + sc_sc( r1,r2) must
	/// equal the weighted result of a call to residue_pair_energy.
	/// By default, bb_bb & bb_sc return 0 and sc_sc returns
	/// residue pair energy.
	void
	backbone_backbone_energy(
		Residue const & rsd1,
		Residue const & rsd2,
		Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;


	/// @brief Evaluate the interaction between the backbone of rsd1 and the
	/// sidechain of rsd2 and accumulate the unweighted energies.  The sum
	/// bb_bb(r1,r2) + bb_sc(r1,r2) + bb_sc(r2,r1) + sc_sc( r1,r2) must
	/// equal the unweighted result of a call to residue_pair_energy.
	/// By default, bb_bb & bb_sc return 0 and sc_sc returns
	/// residue pair energy.
	void
	backbone_sidechain_energy(
		Residue const & rsd1,
		Residue const & rsd2,
		Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;

	/// @brief Evaluate the interaction between the sidechain of rsd1 and the
	/// sidechain of rsd2 and accumulate the unweighted energies.  The sum
	/// bb_bb(r1,r2) + bb_sc(r1,r2) + bb_sc(r2,r1) + sc_sc( r1,r2) must
	/// equal the unweighted result of a call to residue_pair_energy.
	/// By default, bb_bb & bb_sc return 0 and sc_sc returns
	/// residue pair energy.
	void
	sidechain_sidechain_energy(
		Residue const & rsd1,
		Residue const & rsd2,
		Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;

	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const override { return false; }

	void
	eval_intrares_energy( const core::conformation::Residue&, const core::pose::Pose&, const core::scoring::ScoreFunction&, core::scoring::EnergyMap&) const override { }

	core::Size
	version() const override { return 1; }

	core::Real
	atomic_interaction_cutoff() const override { return 20; }

private:

	core::scoring::OneDDistPotential const & P_P_potential_;
	core::scoring::OneDDistPotential const & P_S_potential_;
	core::scoring::OneDDistPotential const & P_CEN_potential_;
	core::scoring::OneDDistPotential const & S_S_potential_;
	core::scoring::OneDDistPotential const & S_CEN_potential_;
	core::scoring::OneDDistPotential const & CEN_CEN_potential_;

};


} //energy_methods
} //core

#endif //core_energy_methods_RNA_CoarseDistEnergy_hh
