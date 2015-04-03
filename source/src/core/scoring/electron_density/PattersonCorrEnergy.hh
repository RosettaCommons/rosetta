// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  Scoring a structure's fit to a patterson map
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_electron_density_PattersonCorrEnergy_hh
#define INCLUDED_core_scoring_electron_density_PattersonCorrEnergy_hh

#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers

namespace core {
namespace scoring {
namespace electron_density {


class PattersonCorrEnergy : public methods::ContextDependentLRTwoBodyEnergy  {
public:
	typedef methods::ContextDependentLRTwoBodyEnergy  parent;

public:


	PattersonCorrEnergy();


	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	methods::LongRangeEnergyType
	long_range_type() const;

	virtual	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const;

	virtual	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sfxn ) const;

	virtual	void
	finalize_after_derivatives( pose::Pose &, ScoreFunction const &  ) const;

	virtual	bool
	defines_intrares_energy( EnergyMap const & ) const { return false; }

	virtual	void
	eval_intrares_energy(
		conformation::Residue const & ,
		pose::Pose const & ,
		ScoreFunction const & ,
		EnergyMap &
	) const { }

	virtual	void
	setup_for_packing( pose::Pose &, utility::vector1< bool > const &, utility::vector1< bool > const & ) const { isRepacking=true; };

	virtual	void
	update_residue_for_packing(
		pose::Pose &,
		Size resid
	) const;

	virtual	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	using methods::ContextDependentLRTwoBodyEnergy::finalize_total_energy;

	/// called at the end of energy evaluation
	virtual	void
	finalize_total_energy(
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;


	virtual	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;


	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const {};

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////


private:
	bool map_loaded;
	bool scoreRepacks;
	mutable bool isRepacking;
	mutable bool pose_is_proper;
	mutable double pcc_structure;
	mutable int nreses;
virtual
core::Size version() const;
};


}
}
}

#endif

