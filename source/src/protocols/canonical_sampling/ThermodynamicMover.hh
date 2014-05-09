// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/canonical_sampling/ThermodynamicMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_ThermodynamicMover_hh
#define INCLUDED_protocols_canonical_sampling_ThermodynamicMover_hh

// Unit Headers
#include <protocols/canonical_sampling/ThermodynamicMover.fwd.hh>

// Project Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
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
namespace canonical_sampling {

/// @brief Base class for moves that can obey detailed balance.
///
/// @details In order to sample a thermodynamic state using a Monte Carlo 
/// simulation, the moves must obey detailed balance.  This base class provides 
/// a framework for writing moves that can obey this condition.  One 
/// interesting method is set_preserve_detailed_balance(), which indicates 
/// whether or not detailed balance needs to be obeyed.  This flag makes it 
/// possible to implement fancy (but biased) features for use in contexts where 
/// rigorous thermodynamic sampling isn't needed.  If the move requires a 
/// non-unity proposal ratio to obey detailed balance, it can reimplement 
/// last_proposal_density_ratio().  Support for movers that make multiple trial 
/// moves under the hood is provided by is_multi_trial() and its related 
/// methods.  A number of callbacks, including  initialize_simulation(), 
/// observe_after_metropolis(), and finalize_simulation(), are also defined to 
/// let the mover react to certain milestones in the simulation.

class ThermodynamicMover : public protocols::moves::Mover {

public:

	/// @brief Default constructor.
	ThermodynamicMover();

	/// @brief Default destructor.
	virtual
	~ThermodynamicMover();

	/// @brief Callback executed before any Monte Carlo trials are attempted.
	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	);

	/// @brief Callback executed after the Metropolis criterion is evaluated.
	virtual
	void
	observe_after_metropolis(
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief Return the proposal density ratio for last apply method.
	virtual
	core::Real
	last_proposal_density_ratio();

	/// @brief Callback executed after all Monte Carlo trials are completed.
	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief Return true if detailed balance is being preserved (i.e. no branch 
	/// angle optimization).
	virtual
	bool
	preserve_detailed_balance() const = 0;

	/// @brief Set to true if detailed balance should be preserved (i.e. no 
	/// branch angle optimization).  This will be set to true for all movers used 
	/// by MetropolisHastingsMover.
	virtual
	void
	set_preserve_detailed_balance(
		bool preserve_detailed_balance
	) = 0;

	/// @brief Return true if the move performs multiple trials on each apply.
	/// @see last_inner_score_temperature_delta()
	/// @see metropolis_hastings_mover()
	/// @see set_metropolis_hastings_mover()
	virtual
	bool
	is_multi_trial();

	/// @brief If this is a multi-trial move, return the change in internal 
	/// score/temperature caused by the last call to apply().
	/// @see is_multi_trial()
	virtual
	core::Real
	last_inner_score_temperature_delta();

	/// @brief If this is a multi-trial move, return the MetropolisHastingsMover 
	/// being used internally.
	/// @see is_multi_trial()
	virtual
	protocols::canonical_sampling::MetropolisHastingsMoverAP
	metropolis_hastings_mover();

	/// @brief If this is a multi-trial move, set the MetropolisHastingsMover to 
	/// be used internally.
	/// @see is_multi_trial()
	virtual
	void
	set_metropolis_hastings_mover(
		protocols::canonical_sampling::MetropolisHastingsMoverAP metropolis_hastings_mover
	);

	/// @brief Return a list specifying which torsions may be perturbed by 
	/// apply(), and the in what range each perturbation may be.
	/// @details This method should probably not be pure virtual, and in fact 
	/// should probably not even exist.  I searched most of the codebase, and 
	/// could only find it being used in one pilot app.  It is also a somewhat 
	/// difficult method to write, which means that most of the implementations 
	/// are either untested or no-ops.  It might be better to remove the method 
	/// altogether and implement it on a class-by-class basis as necessary.
	virtual
	utility::vector1<core::id::TorsionID_Range>
	torsion_id_ranges(
		core::pose::Pose & pose
	) = 0;

	/// @brief Return a list specifying which degrees of freedom may be perturbed 
	/// by apply(), and the in what range each perturbation may be.  
	virtual
	utility::vector1<core::id::DOF_ID_Range>
	dof_id_ranges(
		core::pose::Pose & pose
	);

private:

}; //end ThermodynamicMover

} //namespace canonical_sampling
} //namespace protocols

#endif // INCLUDED_protocols_canonical_sampling_ThermodynamicMover_HH
