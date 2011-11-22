// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/moves/ThermodynamicMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_moves_ThermodynamicMover_hh
#define INCLUDED_protocols_moves_ThermodynamicMover_hh

// Unit Headers
#include <protocols/moves/ThermodynamicMover.fwd.hh>

// Project Headers
#include <protocols/moves/MetropolisHastingsMover.fwd.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>

#include <core/id/DOF_ID_Range.fwd.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/id/TorsionID_Range.fwd.hh>

namespace protocols {
namespace moves {

///@details
class ThermodynamicMover : public protocols::moves::Mover {

public:

	///@brief
	ThermodynamicMover();

	virtual
	~ThermodynamicMover();

	/// @brief callback executed before any Monte Carlo trials
	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief callback executed after the Metropolis criterion is evaluated
	virtual
	void
	observe_after_metropolis(
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief callback for proposal density ratio of last apply method
	virtual
	core::Real
	last_proposal_density_ratio();

	/// @brief callback executed after all Monte Carlo trials
	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief get whether detailed balance is preserved (i.e. no branch angle optimization during moves)
	virtual
	bool
	preserve_detailed_balance() const = 0;

	/// @brief set whether detailed balance is preserved (i.e. no branch angle optimization during moves)
	virtual
	void
	set_preserve_detailed_balance(
		bool preserve_detailed_balance
	) = 0;

	/// @brief determine whether the move performs multiple trials on a single apply
	virtual
	bool
	is_multi_trial();

	/// @brief get change in internal score/temperature for last apply method of multiple trial movers
	virtual
	core::Real
	last_inner_score_temperature_delta();

	/// @brief get the MetropolisHastingsMover for multiple trial movers
	virtual
	protocols::moves::MetropolisHastingsMoverAP
	metropolis_hastings_mover();

	/// @brief set the MetropolisHastingsMover for multiple trial movers
	virtual
	void
	set_metropolis_hastings_mover(
		protocols::moves::MetropolisHastingsMoverAP metropolis_hastings_mover
	);

	/// @brief get the TorsionIDs perturbed by the mover during moves, along with their ranges
	virtual
	utility::vector1<core::id::TorsionID_Range>
	torsion_id_ranges(
		core::pose::Pose & pose
	) = 0;

	/// @brief get the DOF_IDs perturbed by the mover during moves (not including those returned by torsion_id_ranges(), along with their ranges
	virtual
	utility::vector1<core::id::DOF_ID_Range>
	dof_id_ranges(
		core::pose::Pose & pose
	);

private:

}; //end ThermodynamicMover

} //namespace moves
} //namespace protocols

#endif // INCLUDED_protocols_moves_ThermodynamicMover_HH
