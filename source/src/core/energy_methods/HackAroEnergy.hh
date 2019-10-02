// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/HackAroEnergy.hh
/// @brief  Electrostatics for RNA
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_hackaro_HackAroEnergy_hh
#define INCLUDED_core_scoring_hackaro_HackAroEnergy_hh

/// Unit Headers
#include <core/energy_methods/HackAroEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/kinematics/Stub.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace hackaro {


class HackAroEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;

public:


	HackAroEnergy();

	/// clone
	methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const override;

	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const override;

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const override;


	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const override {}


	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const & scorefxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;

	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const override { return false; }

	Distance
	atomic_interaction_cutoff() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const override;

private:

	Vector
	get_centroid( conformation::Residue const & rsd ) const;

	kinematics::Stub
	get_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid ) const;

	core::Real
	get_aro_axis_score_ANGLE(
		Real const cos_theta,
		Real & deriv ) const;

	core::Real
	get_aro_axis_score_DIST(
		Real const dist,
		Real & deriv ) const;

	void
	residue_pair_energy_aro_aro(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		EnergyMap & emap) const;


	void
	eval_atom_derivative_aro_aro(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;
	core::Size version() const override;

private:

};


}
}
}

#endif
