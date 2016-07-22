// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RamaPreProEnergy.hh
/// @brief  A variation on the Ramachandran scorefunction that has separate probability tables for residues that precede prolines.
/// @author Frank DiMaio
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- modified this to work with D-amino acids, BACKBONE_AA amino acids, and cyclic geometry.

#ifndef INCLUDED_core_scoring_methods_RamaPreProEnergy_hh
#define INCLUDED_core_scoring_methods_RamaPreProEnergy_hh

// Unit headers
#include <core/scoring/methods/RamaPreProEnergy.fwd.hh>
#include <core/scoring/RamaPrePro.fwd.hh>

// Package headers
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


class RamaPreProEnergy : public ContextIndependentLRTwoBodyEnergy {
public:
	typedef ContextIndependentLRTwoBodyEnergy  parent;

public:
	RamaPreProEnergy( );

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	virtual void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size rsd1,
		Size rsd2
	) const;


	virtual void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const { }

	virtual
	Real
	eval_intraresidue_dof_derivative(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & /*min_data*/,
		id::DOF_ID const & /*dof_id*/,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		ScoreFunction const & /*sfxn*/,
		EnergyMap const & weights
	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	bool
	defines_intrares_dof_derivatives( pose::Pose const & ) const { return true; }

	virtual
	Distance
	atomic_interaction_cutoff() const { return 0.0; }

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const { }

	//fpd  use the new minimizer interface
	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	methods::LongRangeEnergyType
	long_range_type() const;

private:
	RamaPrePro const & potential_;

	virtual
	core::Size version() const;
};


} // namespace methods
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_methods_RamaPreProEnergy_HH
