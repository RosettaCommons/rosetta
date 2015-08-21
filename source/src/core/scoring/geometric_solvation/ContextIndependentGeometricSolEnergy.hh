// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ContextIndependentGeometricSolEnergy.hh
/// @brief  Geometric solvation energy.
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_core_scoring_geometric_solvation_ContextIndependentGeometricSolEnergy_hh
#define INCLUDED_core_scoring_geometric_solvation_ContextIndependentGeometricSolEnergy_hh

// Unit Headers
#include <core/scoring/geometric_solvation/ContextIndependentGeometricSolEnergy.fwd.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.fwd.hh>
#include <core/types.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/kinematics/MinimizerMapBase.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>

//Auto Headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>


namespace core {
namespace scoring {
namespace geometric_solvation {


class ContextIndependentGeometricSolEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;
public:


	ContextIndependentGeometricSolEnergy( methods::EnergyMethodOptions const & options );

	/// @brief copy c-tor
	ContextIndependentGeometricSolEnergy( ContextIndependentGeometricSolEnergy const & src );

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	/// attempt to precalculate backbone/backbone energies in advance
	virtual
	void
	setup_for_packing(
		pose::Pose & pose,
		utility::vector1< bool > const &,
		utility::vector1< bool > const & ) const;

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn ) const;

	virtual
	bool
	defines_score_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		bool res_moving_wrt_eachother
	) const;

	virtual
	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		Size const res1,
		Size const res2,
		pose::Pose const & pose,
		ScoreFunction const &
	) const;

	virtual
	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2
	) const;

	virtual
	etable::count_pair::CountPairFunctionCOP
	get_intrares_countpair(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &
	) const;

	virtual
	bool
	use_extended_residue_pair_energy_interface() const;

	virtual
	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData & pair_data
	) const;

	virtual
	void
	eval_residue_pair_derivatives(
		conformation::Residue const & ires,
		conformation::Residue const & jres,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	virtual
	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const &,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;

	virtual
	bool
	requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const;

	virtual
	void
	setup_for_minimizing_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		ResSingleMinimizationData & res_data_cache
	) const;

	virtual
	void
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	/// This evaluates everything for now,
	/// but eventually may want to split this
	/// based on backbone/backbone vs. others,
	/// as is carried out in HBondEnergy.cc
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & pose ) const;

	// Undefined, commenting out to fix PyRosetta build
	/* void
	eval_atom_derivative_intra_RNA(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
	) const;
	*/


	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	/// f1 and f2 are zeroed
	// virtual
	// void
	// eval_atom_derivative(
	//  id::AtomID const & atom_id,
	//  pose::Pose const & pose,
	//  kinematics::DomainMap const &,
	//  ScoreFunction const &,
	//  EnergyMap const & weights,
	//  Vector & F1,
	//  Vector & F2
	// ) const;


	virtual
	bool
	defines_intrares_energy( EnergyMap const & weights ) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & ,
		EnergyMap & emap
	) const;

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	core::Size version() const;

	/// @brief GeometricSolEnergy is context sensitive
	virtual
	void indicate_required_context_graphs(
		utility::vector1< bool > & context_graphs_required ) const;

	void
	precalculate_bb_bb_energy_for_design(
		pose::Pose const &
	) const;


private:

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
	methods::EnergyMethodOptions const & options_;

	GeometricSolEnergyEvaluatorOP evaluator_;

	mutable Real precalculated_bb_bb_energy_;

	mutable bool using_extended_method_;


};

} // hbonds
} // scoring
} // core

#endif
