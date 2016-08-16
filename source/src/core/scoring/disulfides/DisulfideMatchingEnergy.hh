// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/disulfides/DisulfideMatchingEnergy.hh
/// @brief  Centroid Disulfide Energy class declaration
/// @author rvernon@u.washington.edu
/// @date   02/09/10

#ifndef INCLUDED_core_scoring_disulfides_DisulfideMatchingEnergy_hh
#define INCLUDED_core_scoring_disulfides_DisulfideMatchingEnergy_hh

// Unit headers
#include <core/scoring/disulfides/DisulfideMatchingEnergy.fwd.hh>

// Package headers
#include <core/scoring/disulfides/DisulfideMatchingPotential.fwd.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace disulfides {

class DisulfideMatchingEnergy : public methods::ContextIndependentLRTwoBodyEnergy {
public:
	typedef methods::ContextIndependentLRTwoBodyEnergy parent;

public:
	DisulfideMatchingEnergy( DisulfideMatchingPotential const & potential );
	virtual ~DisulfideMatchingEnergy();

	// EnergyMethod Methods:
	virtual
	methods::EnergyMethodOP
	clone() const;

	virtual
	void
	setup_for_scoring( pose::Pose &, ScoreFunction const & ) const;

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
	DisulfideMatchingPotential const & potential_;
	virtual
	core::Size version() const;
};


} // namespace disulfides
} // namespace scoring
} // namespace core

#endif //INCLUDED_core_scoring_disulfides_DisulfideMatchingEnergy_HH
