// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/FullatomDisulfideEnergy.fwd.hh
/// @brief  Disulfide Energy class forward declaration
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_disulfides_FullatomDisulfideEnergy_hh
#define INCLUDED_core_scoring_disulfides_FullatomDisulfideEnergy_hh

// Unit headers
#include <core/scoring/disulfides/FullatomDisulfideEnergy.fwd.hh>

// Package headers
#include <core/scoring/disulfides/FullatomDisulfidePotential.fwd.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace disulfides {

class FullatomDisulfideEnergy : public methods::ContextIndependentLRTwoBodyEnergy {
public:
	typedef methods::ContextIndependentLRTwoBodyEnergy parent;

public:

	FullatomDisulfideEnergy( FullatomDisulfidePotential const & potential );
	virtual ~FullatomDisulfideEnergy();

	// EnergyMethod Methods:
	virtual
	methods::EnergyMethodOP
	clone() const;

	/// @brief check that the fullatom disulfid energy container is the right size, and the
	/// set of disulfides it holds corresponds correctly to the set of disulfides in the Pose.
	void
	ensure_lrenergy_container_is_up_to_date(
		pose::Pose & pose
	) const;

	virtual
	void
	setup_for_scoring( pose::Pose &, ScoreFunction const & ) const;

	/// @brief Make sure that the FullatomDisulfideEnergyContainer is ready for packing.
	virtual
	void
	setup_for_packing(
		pose::Pose & pose,
		utility::vector1< bool > const & residues_repacking,
		utility::vector1< bool > const & residues_designing
	) const;

	/// @brief Returns true only for disulfide-bonded residue pairs
	virtual
	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const;

	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const;

	/// @brief During minimization, access atom-index information from the ResPairMinimizationData
	virtual
	bool
	use_extended_residue_pair_energy_interface() const;

	/// @brief Pull the atom-index information out of min_data object, and use those indices to
	/// score the disulfide bond
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

	/// @brief Initialize the atom-index information for a particular residue pair and store those
	/// indices in the ResPairMinimizationData data_cache
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

	/// @brief Retrieve the atom-index information for this residue pair from the minpair_data object
	/// and evaluate the derivatives for a particular atom.
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


	virtual
	void
	old_eval_atom_derivative(
		id::AtomID const &,
		pose::Pose const &,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const &,
		Vector &,// F1,
		Vector & // F2
	) const;


	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const &,
		id::TorsionID const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap const &
	) const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

	// TwoBodyEnergy Methods:
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
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

	// LongRangeTwoBodyEnergy methods
	virtual methods::LongRangeEnergyType long_range_type() const;

	virtual
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const;

private:
	FullatomDisulfidePotential const & potential_;
	virtual
	core::Size version() const;

};


} // namespace disulfides
} // namespace scoring
} // namespace core

#endif
