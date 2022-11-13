// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/MixedMonteCarlo.hh
/// @brief A hybrid monte carlo mover
/// @author AmeyaHarmalkar (harmalkar.ameya24@gmail.com)


#ifndef INCLUDED_protocols_moves_MixedMonteCarlo_hh
#define INCLUDED_protocols_moves_MixedMonteCarlo_hh

#include <protocols/moves/MixedMonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

namespace protocols {
namespace moves {

/// @brief A wrapper class around MonteCarlo that uses both centroid and FA scorefxns.
/// @details A Hybrid Monte Carlo mover that employs both centroid and ful-atom score functions. The purpose
/// of this object is to use a combination of centroid and full-atom scorefxn for Boltzmann scoring but only
/// use only the full-atom pose for bookkeeping. This object is particularly important for Resolution Exchange
/// MC based docking.
/// @author AmeyaHarmalkar (harmalkar.ameya24@gmail.com)

class MixedMonteCarlo : public utility::VirtualBase {

public:

	/// @brief Default constructor.
	MixedMonteCarlo(
		core::pose::Pose const & low_pose,
		core::pose::Pose const & high_pose,
		core::scoring::ScoreFunction const & low_scorefxn,
		core::scoring::ScoreFunction const & high_scorefxn,
		core::Real const tuning_param,
		core::Real const temperature
	);

	/// @brief Destructor.
	~MixedMonteCarlo() override;


	/// @brief This method should evaluate net score from CG and AA modes
	bool
	boltzmann(
		core::pose::Pose & low_pose,
		core::pose::Pose & high_pose,
		std::string const & move_type = "unk"
	);


	/// @brief Sets the lowest and last accepted poses to pose and calls MC reset
	void
	reset( core::pose::Pose const & low_pose, core::pose::Pose const & high_pose);


	/// @brief this gets you only the low-res pose
	core::pose::Pose const &
	last_accepted_pose() const{
		return *last_accepted_pose_;
	}


	/// @brief this gets you only the low-res pose
	core::pose::Pose const &
	lowest_score_pose() const{
		return *lowest_score_pose_;
	}


	/// @brief Copies the lowest_score_ pose into the input pose as well as in the
	/// last_accepted_ pose
	void
	recover_low( core::pose::Pose & high_pose );


	/// @brief Copies the lowest_score_ pose into the input pose as well as in the
	/// last_accepted_ pose. Returns the highres pose
	void
	recover_low( core::pose::Pose & low_pose, core::pose::Pose & high_pose );


	/// @brief Returns the highres scorefxn
	core::scoring::ScoreFunction const &
	highres_score_function() const;


	/// @brief Returns the lowres scorefxn
	core::scoring::ScoreFunction const &
	lowres_score_function() const;


	/// @brief Displays the last accepted and lowest scores
	void
	show_scores() const;

	/// @brief Resets the mover counters
	//  void
	//  reset_counters();


	/// @brief Displays the number of trials performed, fraction
	/// of trial moves accepted, and the average energy drop per
	/// accepted trial by mover types applied (unknown movers or
	/// perturbations are listed as "unktrials")
	//  void
	//  show_counters() const;


	/// @brief Returns the total number of trials since the last reset
	/// @note: MonteCarlo.boltzmann(pose) updates the number of trials
	//  core::Size total_trials() const;


	/// @brief Returns the score value of the last accepted pose
	core::Real
	last_accepted_score() const;


	/// @brief Returns the score value of the lowest score pose encountered
	core::Real
	lowest_score() const;


	/// @brief Const access to the MC object
	protocols::moves::MonteCarlo const &
	MC() const;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	// MixedMonteCarloOP clone() const;

private: /// Methods

	/// @brief Sets the lowest and last accepted low pose to pose
	void
	reset_mmc( core::pose::Pose const & high_pose );


	/// @brief Returns the mixed resolution score
	core::Real
	score_mixed_res(
		core::pose::Pose const & low_pose,
		core::pose::Pose const & high_pose
	);


	/// @brief Assignment operator
	MixedMonteCarlo & operator = ( MixedMonteCarlo const & rhs );

	/// @brief Declared only for the #inclusion drive.
	MixedMonteCarlo();

	/// @brief Declared only for the #inclusion drive.
	MixedMonteCarlo( MixedMonteCarlo const & rhs );


private:

	/// @brief accepted structure
	core::pose::PoseOP last_accepted_pose_;

	/// @brief Lowest structure we have seen
	core::pose::PoseOP lowest_score_pose_;

	/// @brief Internal scoring functions for force calls if not set
	core::scoring::ScoreFunctionOP lowres_scorefxn_;
	core::scoring::ScoreFunctionOP highres_scorefxn_;

	/// @brief Tuning parameter for scorefxns
	core::Real tuning_param_;

	/// @brief Getting the MC object
	protocols::moves::MonteCarlo MC_;


};

} //moves
} //protocols

#endif //INCLUDED_protocols_moves_MixedMonteCarlo_HH
