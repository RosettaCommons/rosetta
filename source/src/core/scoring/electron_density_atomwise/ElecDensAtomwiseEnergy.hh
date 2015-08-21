// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/electron_density_atomwise/ElecDensAtomwiseEnergy.hh
/// @brief  Declaration for elec_dens_atomwise scoring method
/// @author Fang-Chieh Chou

#ifndef INCLUDED_core_scoring_electron_density_atomwise_ElecDensAtomwiseEnergy_HH
#define INCLUDED_core_scoring_electron_density_atomwise_ElecDensAtomwiseEnergy_HH

// Package headers
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>


namespace core {
namespace scoring {
namespace electron_density_atomwise {

class ElecDensAtomwiseEnergy : public methods::ContextIndependentLRTwoBodyEnergy {
public:
	typedef methods::ContextIndependentLRTwoBodyEnergy parent;

public:


	ElecDensAtomwiseEnergy();


	virtual ~ElecDensAtomwiseEnergy();

	virtual
	methods::EnergyMethodOP
	clone() const;


	virtual
	void indicate_required_context_graphs ( utility::vector1< bool > & ) const {};

	/////////////////////////////////////////////////////////////////////////////

	methods::LongRangeEnergyType
	long_range_type() const;


	virtual
	bool
	defines_residue_pair_energy (
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const;


	virtual
	bool
	defines_intrares_energy ( EnergyMap const & ) const {
		return true;
	}


	virtual
	void
	residue_pair_energy (
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;


	virtual
	void
	eval_intrares_energy (
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {
		return;
	}


	virtual
	void
	eval_atom_derivative (
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;


	virtual
	void
	setup_for_scoring (
		pose::Pose & pose,
		ScoreFunction const &
	) const;

	using parent::finalize_total_energy;


	virtual
	void
	finalize_total_energy (
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {
		return;
	}

private:
	mutable bool pose_is_proper;
	virtual
	core::Size version() const;
};

} // electron_density_atomwise
} // scoring
} // core


#endif
