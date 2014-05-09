// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange ( oliver.lange@tum.de )

#ifndef INCLUDED_protocols_canonical_sampling_BiasedMonteCarlo_hh
#define INCLUDED_protocols_canonical_sampling_BiasedMonteCarlo_hh

// Unit Headers
#include <protocols/canonical_sampling/BiasedMonteCarlo.fwd.hh>

#include <protocols/canonical_sampling/BiasEnergy.fwd.hh>

#include <protocols/moves/MonteCarlo.hh>

// Utility Headers
#include <core/types.hh>

namespace protocols {
namespace canonical_sampling {

///@details
class BiasedMonteCarlo : public protocols::moves::MonteCarlo {
	typedef protocols::moves::MonteCarlo Parent;
public:

	/// @brief Constructs a useable MonteCarlo object
	///
	/// mc = MonteCarlo( init_pose , scorefxn , temp )
	///
	/// Pose           init_pose   /manipulated during the simulation
	/// ScoreFunction  scorefxn    /evaluates pose scores
	/// Real (float)   temp        /used in the Metropolis Criterion
	BiasedMonteCarlo(
		Pose const & init_pose, // PoseCOP init_pose,
		ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
		Real const temperature,
		BiasEnergyOP bias_energy
	);


	/// @brief Constructor without Pose -- call reset(pose) before first use
	BiasedMonteCarlo(
		ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
		Real const temperature,
		BiasEnergyOP bias_energy
	);

	BiasedMonteCarlo( BiasedMonteCarlo const& );

	virtual
	protocols::moves::MonteCarloOP clone() {
		return new BiasedMonteCarlo( *this );
	}
	//BiasedMonteCarlo& operator=( BiasedMonteCarlo const& );

	/// @brief Applies the Metropolis Criterion on pose based on
	/// the ScoreFunction, temperature, and the last accepted
	/// pose. This method evaluates the change in score, compares
	/// the trial pose to the last accepted pose, and updates the
	/// pose structure and simulation statistics appropriately
	///
	/// example(s):
	///     mc.boltzmann( pose )
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accepted_score
	///     MonteCarlo.lowest_score
	virtual bool
	boltzmann(
		Pose & pose,//PoseOP pose,
		std::string const & move_type = "unk",
		core::Real const proposal_density_ratio = 1,
		core::Real const inner_score_temperature_delta = 0
	);

	virtual
	void
	reset( Pose const & pose );

	virtual
	void
	score_function( ScoreFunction const & scorefxn ); // ScoreFunctionCOP scorefxn )

	virtual void
	set_temperature( Real const temp );

	void set_bias_energy( BiasEnergyOP );

	BiasEnergyOP bias_energy() { return bias_energy_; };

private:
	BiasEnergyOP bias_energy_;
};

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_BiasedMonteCarlo_HH
