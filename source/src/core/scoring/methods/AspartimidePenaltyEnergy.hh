// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/AspartimidePenaltyEnergy.hh
/// @brief  This is a score term that penalizes sequences that are likely to result in aspartimide formation during peptide synthesis.
/// @details This is intended for peptide design applications only.  Sequences penalized are LASP-D*, LASP-LSER, LASP-LTHR, LASP-LGLN,
/// and the mirror-image equivalents (DASP-L*, DASP-DSER, DASP-DTHR, DASP-DGLN), plus D/LASP-GLY.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_scoring_methods_AspartimidePenaltyEnergy_hh
#define INCLUDED_core_scoring_methods_AspartimidePenaltyEnergy_hh

// Unit headers
#include <core/scoring/methods/AspartimidePenaltyEnergy.fwd.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/util.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/types.hh>

// C++ headers
#include <iostream>
#include <map>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class AspartimidePenaltyEnergy : public ContextIndependentLRTwoBodyEnergy {
public:
	typedef ContextIndependentLRTwoBodyEnergy  parent;

public:

	/// @brief Constructor.
	///
	AspartimidePenaltyEnergy( );

	/// @brief Constructor that sets penalty value.
	/// @details The penalty value is the energetic hit for each aspartimide in the
	/// sequence (which will be multiplied by the score function's weight, of course).
	AspartimidePenaltyEnergy( core::Real const &penalty_value );

	/// @brief Destructor.
	///
	~AspartimidePenaltyEnergy( );

	/// @brief Copy this energy object and return an owning pointer to the copy.
	///
	virtual EnergyMethodOP clone() const;

	/// @brief Method called before scoring a pose.
	///
	virtual void setup_for_scoring( pose::Pose & pose, ScoreFunction const &sfxn ) const;

	/// @brief Are the two residues (rsd1, rsd2) two residues that should be scored by this scorefunction?
	/// @details Returns true only if rsd2 is connected to the C-terminus of rsd1 by its N-terminal connection, or
	/// vice versa.
	virtual
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size rsd1,
		Size rsd2
	) const;

	/// @brief Score the residues (rsd1, rsd2) and put the energy in the EnergyMap.
	///
	virtual void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	/// @brief Does nothing, since there is no one-body energy associated with this term.
	///
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const { }

	/// @brief Does nothing, since this term is spatially invariant (i.e. has no DoF derivatives, because
	/// the value depends only on residue identities).
	virtual
	Real
	eval_intraresidue_dof_derivative(
		conformation::Residue const & /*rsd*/,
		ResSingleMinimizationData const & /*min_data*/,
		id::DOF_ID const & /*dof_id*/,
		id::TorsionID const & /*tor_id*/,
		pose::Pose const & /*pose*/,
		ScoreFunction const & /*sfxn*/,
		EnergyMap const & /*weights*/
	) const;

	/// @brief Returns false -- there's no one-body energy defined here.
	///
	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	/// @brief Returns false -- there are no derivatives defined here.
	/// @details This score term's value depends only on residue identities, not on geometric DoFs.
	virtual
	bool
	defines_intrares_dof_derivatives( pose::Pose const & ) const { return false; }

	virtual
	Distance
	atomic_interaction_cutoff() const { return 0.0; }

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const { }

	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	methods::LongRangeEnergyType
	long_range_type() const;

private: // Private functions:

	/// @brief Is the given AA for the first residue one of the possible types that this score term scores?
	/// @details Returns true for aa_asp/aa_das, false otherwise.
	bool first_res_types( core::chemical::AA const aa) const;

	/// @brief Is the given AA for the second residue one of the possible types that this score term scores?
	/// @details Returns true for aa_gly/aa_asn/aa_ser/aa_thr/aa_dan/aa_dse/aa_dth, false otherwise.
	bool second_res_types( core::chemical::AA const aa) const;

	virtual
	core::Size version() const;

private: // Private variables:

	/// @brief The point value of the penalty applied for each aspartimide-forming sequence encountered.
	/// @details Set to 25.0 in default constructor; read from options system in constructor that's actually called.
	core::Real aspartimide_penalty_value_;

};


} // namespace methods
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_methods_AspartimidePenaltyEnergy_HH
