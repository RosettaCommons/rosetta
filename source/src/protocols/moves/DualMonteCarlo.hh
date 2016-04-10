// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/DualMonteCarlo.hh
/// @brief "dual" MonteCarlo header - wraps MonteCarlo to allow for centroid MonteCarlo scoring of a fullatom pose (or other similar situations)
/// @author Steven Lewis


#ifndef INCLUDED_protocols_moves_DualMonteCarlo_hh
#define INCLUDED_protocols_moves_DualMonteCarlo_hh

// Unit headers
#include <protocols/moves/DualMonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarlo.hh> //composition needs full header (?)

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

// these are not used here but I've left them for reference
//  mc_accepted
//   3 = accepted:score beat low score and last_accepted score
//   2 = accepted:score beat last_accepted score
//   1 = thermally accepted: score worse than last_accepted score
//   0 = not accepted
// enum MCA {
//  MCA_accepted_score_beat_low=3,
//  MCA_accepted_score_beat_last=2,
//  MCA_accepted_thermally=1,
//  MCA_rejected=0
// }; // mc_accepted state

/// @details DualMonteCarlo is a wrapper class around MonteCarlo.  Its original purpose is to allow for Boltzmann scoring a pose in centroid mode while simultaneously tracking a fullatom equivalent.  It should work for any paired poses.  Generally, "DMC_pose" refers to the pose this class tracks, whereas "MC_pose" refers to the one used by the underlying MC object.  While this class contains a MonteCarlo object (by composition, not inheritance); you cannot pass in a MC object to it, nor can you get an OP to that object.  const access to the underlying MC is provided at this time; if you want nonconst access you can write it in but beware that the user can then ruin the bookkeeping!  Also note that I (SML) left the framework in place in comments for a lot of MC's functions within DMC like checkpointing and counters.  Flesh these out if you want them...
class DualMonteCarlo : public utility::pointer::ReferenceCount {

public:

	DualMonteCarlo(
		core::pose::Pose const & DMC_pose,
		core::pose::Pose const & MC_pose,
		core::scoring::ScoreFunction const & DMC_scorefunction,
		core::scoring::ScoreFunction const & MC_scorefunction,
		core::Real const temperature
	);

	~DualMonteCarlo();

	bool
	boltzmann(
		core::pose::Pose & DMC_pose,
		core::pose::Pose & MC_pose,
		std::string const & move_type = "unk"
	);

	/// @brief sets lowest and last_accepted poses to pose; calls MC's reset too
	void
	reset( core::pose::Pose const & DMC_pose, core::pose::Pose const & MC_pose );

	/// @brief this gets you only the DMC's pose; call MC directly for its.  Note that this pose is accepted based on the partner pose in the underlying MC.
	core::pose::Pose const &
	last_accepted_pose() const
	{
		return *last_accepted_pose_;
	}

	//@brief this gets you only the DMC's pose; call MC directly for its.  Note that this is not necessarily lower in score than the last_accepted; it changes not based on score but based on the partner pose.
	core::pose::Pose const &
	lowest_score_pose() const
	{
		return *lowest_score_pose_;
	}

	// for recovering from a checkpoint
	//void set_last_accepted_pose( core::pose::Pose const & pose );

	// for recovering from a checkpoint
	//void set_lowest_score_pose( core::pose::Pose const & pose );

	// void
	//  attach_observer_to_last_accepted_conformation( core::conformation::ConformationObserverOP observer );

	//  void
	//  attach_observer_to_lowest_score_conformation( core::conformation::ConformationObserverOP observer );

	//  void
	//  attach_observer_to_last_accepted_pose( core::pose::PoseObserver * observer );

	//  void
	//  attach_observer_to_lowest_score_pose( core::pose::PoseObserver * observer );

	/// @brief Copies the lowest_score_ pose into the input pose, as well as into
	/// last_accepted_pose_.  Copies action for underlying MC but does not return that pose.
	void
	recover_low( core::pose::Pose & DMC_pose );

	/// @brief Copies the lowest_score_ pose into the input pose, as well as into
	/// last_accepted_pose_.  Always copies action for underlying MC.  This version returns the MC's pose.
	void
	recover_low( core::pose::Pose & DMC_pose, core::pose::Pose & MC_pose );

	//  /// @brief change the scorefxn,  re-scores last-accepted and lowest-score pose.  Will NOT change their identities.
	//  void
	//  score_function( core::scoring::ScoreFunction const & DMC_scorefunction );


	//  /// @brief change the scorefxns,  re-scores last-accepted and lowest-score pose, calls MC's score_function, and may change the identity of the poses if MC does too.
	//  void
	//  score_function(
	//          core::scoring::ScoreFunction const & DMC_scorefunction,
	//          core::scoring::ScoreFunction const & MC_scorefunction)

	core::scoring::ScoreFunction const &
	score_function() const;

	void
	show_scores() const;

	//  void
	//  reset_counters();

	//  void
	//  show_counters() const;

	//  //@brief returns number in all trial-counters
	//  Size total_trials() const;

	core::Real
	last_accepted_score() const;

	core::Real
	lowest_score() const;

	/// @brief const access to MC object
	protocols::moves::MonteCarlo const &
	MC() const;

private: //functions
	void reset_DMC( core::pose::Pose const & DMC_pose );

private: //magic 4 functions, don't call these
	/// @brief declared only for the #inclusion drive.  If you need this function implement it and move it to public.
	DualMonteCarlo & operator= (DualMonteCarlo const & rhs);

	/// @brief declared only for the #inclusion drive.  If you need this function implement it and move it to public.
	DualMonteCarlo();

	/// @brief declared only for the #inclusion drive.  If you need this function implement it and move it to public.
	DualMonteCarlo( DualMonteCarlo const & rhs);

private:

	/// @brief accepted structure
	core::pose::PoseOP last_accepted_pose_;

	/// @brief lowest energy structure we've seen
	core::pose::PoseOP lowest_score_pose_;

	/// @brief our own internal scoring function - passed in and out by reference to force calls to score_function, enforcing bookkeeping
	core::scoring::ScoreFunctionOP score_function_;

	/// @brief our MC; not intended to be dropped into EnergyCutRotamerTrialsMover, TrialMover, etc.
	protocols::moves::MonteCarlo MC_;

};

} // moves
} // protocols

#endif //INCLUDED_protocols_moves_DualMonteCarlo_HH
