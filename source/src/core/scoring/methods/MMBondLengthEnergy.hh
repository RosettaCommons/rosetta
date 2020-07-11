// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MMBondLengthEnergy.hh
/// @brief  Molecular mechanics bond length score class
/// @author Frank DiMaio (based on Colin Smith's MMBondAngle potential)

#ifndef INCLUDED_core_scoring_methods_MMBondLengthEnergy_hh
#define INCLUDED_core_scoring_methods_MMBondLengthEnergy_hh

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/mm/MMBondLengthScore.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/types.hh>

// C++ headers
#include <iosfwd>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class MMBondLengthEnergy : public ContextIndependentTwoBodyEnergy  {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;
public:


	MMBondLengthEnergy( methods::EnergyMethodOptions const & options );


	MMBondLengthEnergy( MMBondLengthEnergy const & src );

	~MMBondLengthEnergy() override;

	/// clone
	EnergyMethodOP
	clone() const override;


	void
	setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const override;


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

	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const override ;

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;

	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & emap,
		Vector & F1,
		Vector & F2
	) const override;


	/// @brief MMBondLengthEnergy does not have an atomic interation threshold
	Distance
	atomic_interaction_cutoff() const override;

	/// @brief MMBondLengthEnergy is context independent; indicates that no
	/// context graphs are required
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;

private:
	core::scoring::mm::MMBondLengthScore potential_;

	core::Size version() const override;

};

} // namespace methods
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_methods_MMBondLengthEnergy_HH
