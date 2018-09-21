// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/ApproximateBuriedUnsatPenalty.hh
/// @brief  Guidance term that gives a quadratic approximation to no buried unsats
/// @author Brian Coventry (bcov@uw.edu)


#ifndef INCLUDED_core_pack_guidance_scoreterms_approximate_buried_unsat_penalty_ApproximateBuriedUnsatPenalty_hh
#define INCLUDED_core_pack_guidance_scoreterms_approximate_buried_unsat_penalty_ApproximateBuriedUnsatPenalty_hh

// Unit headers
#include <core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/ApproximateBuriedUnsatPenalty.fwd.hh>

// Package headers
#include <core/pack_basic/RotamerSetsBase.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <basic/datacache/CacheableUint64MathMatrixFloatMap.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace approximate_buried_unsat_penalty {



class ApproximateBuriedUnsatPenalty : public scoring::methods::ContextDependentTwoBodyEnergy {

public:
	ApproximateBuriedUnsatPenalty( core::scoring::methods::EnergyMethodOptions const & options );

	scoring::methods::EnergyMethodOP
	clone() const override;

	// This is so we know that packing has started
	void
	setup_for_packing(
		pose::Pose &,
		utility::vector1< bool > const &,
		utility::vector1< bool > const &
	) const override;

	// This is where we evaluate the 3-body terms
	void
	setup_for_packing_with_rotsets(
		pose::Pose & pose,
		pack_basic::RotamerSetsBaseOP const & rotsets,
		scoring::ScoreFunction const & sfxn
	) const override;

	// This is where we add the one-body terms in
	void
	evaluate_rotamer_intrares_energies(
		conformation::RotamerSetBase const & set,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		utility::vector1< core::PackerEnergy > & energies
	) const override;

	// This is where we add the one-body terms in
	void
	evaluate_rotamer_intrares_energy_maps(
		conformation::RotamerSetBase const & set,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		utility::vector1< scoring::EnergyMap > & emaps
	) const override;

	// This is where we add the two-body terms in
	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		scoring::EnergyMap const & weights,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const override;

	// This is where we evaluate 3-body terms during scoring
	void
	setup_for_scoring(
		pose::Pose &,
		scoring::ScoreFunction const &
	) const override;

	// So we know when we're done scoring
	void
	finalize_total_energy(
		pose::Pose & pose,
		scoring::ScoreFunction const & sfxn,
		scoring::EnergyMap & total_energy
	) const override;

	// So we can disable ourselves during minimization
	void
	setup_for_minimizing(
		pose::Pose & ,
		scoring::ScoreFunction const & ,
		kinematics::MinimizerMapBase const &
	) const override;

	// So we can re-enable ourselves after minimization
	void
	finalize_after_minimizing(
		pose::Pose & pose
	) const override;

	// This throws an exception if we aren't disabled. Does nothing otherwise
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		scoring::EnergyMap & emap
	) const override;

	// Always returns true
	bool
	defines_intrares_energy( scoring::EnergyMap const & ) const override { return true; }

	// This throws an exception if we aren't disabled. Returns 0 otherwise
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		scoring::EnergyMap & emap
	) const override;


	// We do this in intrares energy. Not even trying to follow the rules here
	void
	evaluate_rotamer_background_energies(
		conformation::RotamerSetBase const &,
		conformation::Residue const &,
		pose::Pose const &,
		scoring::ScoreFunction const &,
		scoring::EnergyMap const &,
		utility::vector1< core::PackerEnergy > &
	) const override {}

	// We do this in intrares energy. Not even trying to follow the rules here
	void
	evaluate_rotamer_background_energy_maps(
		conformation::RotamerSetBase const &,
		conformation::Residue const &,
		pose::Pose const &,
		scoring::ScoreFunction const &,
		scoring::EnergyMap const &,
		utility::vector1< scoring::EnergyMap > &
	) const override {}

	// We don't want any of these
	void indicate_required_context_graphs(
		utility::vector1< bool > &
	) const override {}

	// Version 1
	core::Size version() const override { return 1; }

	// This is a little tricky. Our direct interactions are simply h-bonds
	//  which means that we should set this to 4.35 like the h-bond terms.
	//  But then we have those weird 3-body terms... What's the max distance
	//  two atoms can be apart and both be h-bonding to the same atom?
	//    Maybe 4.35 * 2 ???
	//  The etable energy uses 6.0. So let's just go with that so we don't
	//  make anything slower
	Distance
	atomic_interaction_cutoff() const override { return 6.0; }

private:

	basic::datacache::CacheableUint64MathMatrixFloatMapCOP
	get_energies_cache( pose::Pose const & pose ) const;


private:

	enum Mode {
		IDLE,
		SCORING,
		PACKING,    // There's no "done" call for packing so we can get stuck here
		MINIMIZING
	};

	mutable Mode mode_;

	// We need this score function to find hbonds. In an absolute sense, time
	//  could be saved by scavenging the hbonds calculated by the normal score function
	//  but those get calculated after we do our calculation, so I don't really see
	//  any way to make that happen.
	scoring::ScoreFunctionOP scorefxn_sc_;
	scoring::ScoreFunctionOP scorefxn_bb_;
	scoring::ScoreFunctionOP scorefxn_hbond_;

	Real hbond_energy_threshold_;
	Real burial_atomic_depth_;
	Real burial_probe_radius_;
	Real burial_resolution_;
	Real oversat_penalty_;

	// This scoreterm can either assume a fixed backbone, or can simulate
	//   a changing backbone (proline, n-methyl, etc.). The tradeoff is that
	//   simulating the changing backbone makes it slower.
	bool assume_const_backbone_;

};

} // approximate_buried_unsat_penalty
} // guidance_scoreterms
} // pack
} // core

#endif
