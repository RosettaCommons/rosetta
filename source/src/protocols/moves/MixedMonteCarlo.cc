// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/MixedMonteCarlo.cc
/// @brief A hybrid monte carlo mover
/// @author AmeyaHarmalkar (harmalkar.ameya24@gmail.com)

// Project headers:
#include <protocols/moves/MixedMonteCarlo.hh>
#include <protocols/moves/MonteCarlo.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>


// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR( "protocols.moves.MixedMonteCarlo" );


namespace protocols {
namespace moves {


/////////////////////
/// Constructors  ///
/////////////////////

/// @details Default Constructor. This would call the MC Constructor as well.
MixedMonteCarlo::MixedMonteCarlo(
	core::pose::Pose const & low_pose,
	core::pose::Pose const & high_pose,
	core::scoring::ScoreFunction const & low_scorefxn,
	core::scoring::ScoreFunction const & high_scorefxn,
	core::Real const tuning_param,
	core::Real const temperature
) :
	last_accepted_pose_( utility::pointer::make_shared< core::pose::Pose >() ),
	lowest_score_pose_( utility::pointer::make_shared< core::pose::Pose >() ),
	MC_( high_pose, high_scorefxn, temperature )
{
	highres_scorefxn_ = high_scorefxn.clone();
	lowres_scorefxn_ = low_scorefxn.clone();
	tuning_param_ = tuning_param ;
	reset( low_pose, high_pose );
}


/// @brief Destructor.
MixedMonteCarlo::~MixedMonteCarlo() = default;


/// @brief
void
MixedMonteCarlo::reset( core::pose::Pose const & low_pose, core::pose::Pose const & high_pose)
{
	reset_mmc( high_pose );
	core::Real score = score_mixed_res(  low_pose, high_pose);
	// Here we are passing the highres pose and the cumulative score to the MC object.
	MC_.reset( low_pose, score );
}


/// @brief Bookkeeping with highres pose
void
MixedMonteCarlo::reset_mmc( core::pose::Pose const & high_pose ){

	// I believe we are just setting the pose here for book-keeping. Is scoring necessary? Maybe not?
	*last_accepted_pose_ = high_pose;
	( *highres_scorefxn_ )( *last_accepted_pose_ );
	*lowest_score_pose_ = *last_accepted_pose_;

}


/// @brief Copies the lowest_score_ pose into the input pose as well as in the last_accepted_ pose.
/// Copies action for underlying MC but does not return that pose
void
MixedMonteCarlo::recover_low( core::pose::Pose & low_pose ){
	core::pose::Pose toss;
	recover_low( low_pose, toss );
}


/// @brief Returns simulation state to the lowest energy structure that has been observed
void
MixedMonteCarlo::recover_low( core::pose::Pose & low_pose, core::pose::Pose & high_pose )
{
	high_pose  = *lowest_score_pose_;
	*last_accepted_pose_ = *lowest_score_pose_;
	// Maybe there should be a method in MC that takes in both score and pose??
	MC_.recover_low( low_pose );
}


/// @brief Shows the current scores
void
MixedMonteCarlo::show_scores() const{
	TR << "Last accepted score : " << last_accepted_score() << std::endl;
	TR << "Lowest score : " << lowest_score() << std::endl;
}

/// @brief Boltzmann for MixedMC. Gets the score and passes to the MC boltzman
/// @details Returns true for an accept, and false otherwise
/// MC_accepted
/// 3 = accepted : score beats low_score and last_accepted score
/// 2 = accepted : score beats last_accepted score
/// 1 = thermally accepted: score worse than last_accepted score
/// 0 = not accepted
bool
MixedMonteCarlo::boltzmann(
	core::pose::Pose & low_pose,
	core::pose::Pose & high_pose,
	std::string const & move_type// = "unk"
){

	core::Real score = score_mixed_res( low_pose, high_pose );
	// TR << "Mixed Res Score : " << score << std::endl;
	// Calculate the mixed res score above and pass it on to the boltzmann of the MC
	const bool accepted( MC_.boltzmann( score, move_type ) );

	switch( MC_.mc_accepted() ){
	case MCA_accepted_score_beat_low :
		reset_mmc( high_pose );
		break;
	case MCA_accepted_score_beat_last :
	case MCA_accepted_thermally :
		// PROF_START( basic::MC_ACCEPT );
		*last_accepted_pose_ = high_pose;
		//PROF_STOP( basic::MC_ACCEPT );
		( *highres_scorefxn_ )( *last_accepted_pose_ );
		break;
	case MCA_rejected :
		high_pose = *last_accepted_pose_;
		break;
	default :
		break;
	}
	return accepted;
}


core::scoring::ScoreFunction const &
MixedMonteCarlo::highres_score_function() const
{
	return *highres_scorefxn_;
}


core::scoring::ScoreFunction const &
MixedMonteCarlo::lowres_score_function() const
{
	return *lowres_scorefxn_;
}


core::Real
MixedMonteCarlo::last_accepted_score() const
{
	// As of now, just returning the high res scores, but this should return the
	// cumulative score somehow, right?
	// return last_accepted_pose_->energies().total_energy();
	return MC_.last_accepted_score();
}


core::Real
MixedMonteCarlo::lowest_score() const
{
	// Just return the highres energies
	// return lowest_score_pose_->energies().total_energy();
	return MC_.lowest_score();
}


/// @brief Const access to the MC object
protocols::moves::MonteCarlo const &
MixedMonteCarlo::MC() const
{
	return MC_;
}


core::Real
MixedMonteCarlo::score_mixed_res(
	core::pose::Pose const & low_pose,
	core::pose::Pose const & high_pose
){
	core::pose::Pose lowres_pose = low_pose;
	core::pose::Pose highres_pose = high_pose;
	core::Real const lowres_score_( (*lowres_scorefxn_)( lowres_pose ) );
	core::Real const highres_score_( (*highres_scorefxn_)( highres_pose ) );
	core::Real score = tuning_param_ * lowres_score_ + ( 1 - tuning_param_ ) * highres_score_;
	return score;
}


/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
// MixedMonteCarloOP
// MixedMonteCarlo::clone() const {
// return utility::pointer::make_shared< MixedMonteCarlo >( *this );
//}

} //moves
} //protocols
